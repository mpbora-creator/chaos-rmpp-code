"""Optional bridge from a resolved PIC IICSI mode to the unchanged nonlinear
charging relaxator used by the theoretical SL->vdPM model.

The plug-in never changes the field solve or particle equations. It obtains
two causal IICSI quadratures from the selected PIC observable, advances only
the published charging coordinate z, and returns a bounded multiplier for the existing physical
photoelectron source.  `enabled=false` returns exactly 1.0 and consumes no RNG.
"""
mutable struct PhotoemissionSLvdPMState
    z::Float64
    x::Float64
    y::Float64
    multiplier::Float64
end

"Checkpointable causal measurement state, kept separate for old-checkpoint compatibility."
mutable struct PhotoemissionIICSIFilterState
    initialized::Bool
    raw_signal::Float64
    baseline::Float64
    filtered_signal::Float64
    filtered_velocity::Float64
end

PhotoemissionIICSIFilterState() =
    PhotoemissionIICSIFilterState(false, 0.0, 0.0, 0.0, 0.0)

PhotoemissionSLvdPMState(cfg) = PhotoemissionSLvdPMState(
    cfg.initial_z, 0.0, 0.0, 1.0)

default_photoemission_sl_vdpm_config() = (
    enabled=false,
    measurement=:local_electric_field_bandpass,
    electric_field_reference=0.01,
    band_frequency=(0.015, 0.05),
    baseline_relaxation_time_omega_pe=12.0,
    mode_number=40,
    window_fraction=(0.10, 0.90),
    mode_reference_amplitude=0.01,
    quadrature_clip=2.5,
    tau_rate_per_omega_pe=0.03,
    coefficients=(
        Lambda=0.19072799832134626,
        C=-2.376110373607035,
        D=-0.9267200175037225,
        E=-0.7870886774439758,
        G=1.183328923259103,
    ),
    initial_z=0.0,
    mean_rate_multiplier=1.0,
    modulation_depth=0.80,
    output_scale=1.0,
    minimum_rate_multiplier=0.05,
    maximum_rate_multiplier=2.0,
)

function photoemission_sl_vdpm_config(cfg)
    defaults = default_photoemission_sl_vdpm_config()
    isnothing(cfg) && return defaults
    coefficients = merge(defaults.coefficients,
        hasproperty(cfg, :coefficients) ? cfg.coefficients : (;))
    return merge(defaults, cfg, (coefficients=coefficients,))
end

"Local debris electric field, normalized before causal baseline/band filtering."
function photoemission_local_field_sample(grid, cfg, bridge)
    center = cfg.debris.center_fraction * grid.length
    scalar_at_position(grid.electric_field, center, grid, cfg) /
        bridge.electric_field_reference
end

function advance_photoemission_iicsi_filter!(state::PhotoemissionIICSIFilterState,
        raw, bridge, dt)
    if !state.initialized
        state.initialized = true
        state.raw_signal = raw
        state.baseline = raw
        return 0.0, 0.0
    end
    low, high = bridge.band_frequency
    frequency = sqrt(low*high)
    omega = 2pi*frequency
    damping = (high-low)/(2frequency)
    tau0 = bridge.baseline_relaxation_time_omega_pe
    # RK4 for baseline plus a causal second-order band-pass resonator.  The
    # returned velocity/omega is the phase quadrature of the filtered field.
    function rhs(u)
        baseline, q, velocity = u
        highpass = raw-baseline
        ((raw-baseline)/tau0,
         velocity,
         omega^2*(highpass-q)-2damping*omega*velocity)
    end
    u = (state.baseline, state.filtered_signal, state.filtered_velocity)
    k1=rhs(u)
    u2=ntuple(i->u[i]+0.5dt*k1[i],3); k2=rhs(u2)
    u3=ntuple(i->u[i]+0.5dt*k2[i],3); k3=rhs(u3)
    u4=ntuple(i->u[i]+dt*k3[i],3); k4=rhs(u4)
    updated=ntuple(i->u[i]+dt*(k1[i]+2k2[i]+2k3[i]+k4[i])/6,3)
    state.raw_signal=raw
    state.baseline=updated[1]
    state.filtered_signal=updated[2]
    state.filtered_velocity=updated[3]
    x=clamp(state.filtered_signal,-bridge.quadrature_clip,bridge.quadrature_clip)
    y=clamp(state.filtered_velocity/omega,-bridge.quadrature_clip,
            bridge.quadrature_clip)
    return x,y
end

"Measure cosine/sine quadratures of a resolved ion-density Fourier mode."
function photoemission_iicsi_quadratures(grid, cfg, bridge)
    left, right = bridge.window_fraction
    selected = findall(x -> left*grid.length <= x <= right*grid.length, grid.x)
    isempty(selected) && return (0.0, 0.0)
    signal = grid.ni[selected]
    signal = signal .- sum(signal)/length(signal)
    phase = 2pi * bridge.mode_number .* grid.x[selected] ./ grid.length
    normalization = 2 / (length(selected)*bridge.mode_reference_amplitude)
    x = clamp(normalization * sum(signal .* cos.(phase)),
              -bridge.quadrature_clip, bridge.quadrature_clip)
    y = clamp(-normalization * sum(signal .* sin.(phase)),
              -bridge.quadrature_clip, bridge.quadrature_clip)
    return x, y
end

@inline function photoemission_relaxator_rhs(z, x, y, bridge)
    p = bridge.coefficients
    p.C*x + p.D*y - p.Lambda*z + p.E*x*z - p.G*z^3
end

"Advance the exact first-order nonlinear charging equation with held PIC input."
function advance_photoemission_sl_vdpm!(state::PhotoemissionSLvdPMState,
        filter_state::PhotoemissionIICSIFilterState, grid, cfg, bridge, dt)
    bridge.enabled || return 1.0
    x,y = if bridge.measurement === :local_electric_field_bandpass
        raw=photoemission_local_field_sample(grid,cfg,bridge)
        advance_photoemission_iicsi_filter!(filter_state,raw,bridge,dt)
    elseif bridge.measurement === :ion_density_mode
        photoemission_iicsi_quadratures(grid,cfg,bridge)
    else
        error("Unsupported photoemission feedback measurement=$(bridge.measurement)")
    end
    h = bridge.tau_rate_per_omega_pe * dt
    z = state.z
    k1 = photoemission_relaxator_rhs(z, x, y, bridge)
    k2 = photoemission_relaxator_rhs(z + 0.5h*k1, x, y, bridge)
    k3 = photoemission_relaxator_rhs(z + 0.5h*k2, x, y, bridge)
    k4 = photoemission_relaxator_rhs(z + h*k3, x, y, bridge)
    state.z = z + h*(k1 + 2k2 + 2k3 + k4)/6
    state.x, state.y = x, y
    raw = bridge.mean_rate_multiplier *
        (1 + bridge.modulation_depth*tanh(state.z/bridge.output_scale))
    state.multiplier = clamp(raw, bridge.minimum_rate_multiplier,
                             bridge.maximum_rate_multiplier)
    return state.multiplier
end

photoemission_sl_vdpm_history_path(history_file) =
    joinpath(dirname(history_file), "photoemission_sl_vdpm_feedback.tsv")

function append_photoemission_sl_vdpm_history!(path, step, time,
        state::PhotoemissionSLvdPMState,
        filter_state::PhotoemissionIICSIFilterState; new_file=false)
    mkpath(dirname(path))
    open(path, new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\tiicsi_x\tiicsi_y\trelaxator_z\tphoto_rate_multiplier\tlocal_E_normalized\tlocal_E_baseline")
        println(io, join((step, time, state.x, state.y, state.z,
                          state.multiplier,filter_state.raw_signal,
                          filter_state.baseline), '\t'))
    end
    return nothing
end
