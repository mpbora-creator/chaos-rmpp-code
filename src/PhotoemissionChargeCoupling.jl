"""Optional coupling between the external RF charge drive and photoemission.

This plug-in never creates particles itself. It supplies a bounded multiplier
to the existing, charge-conserving `emit_debris_photoelectrons!` routine.  Its
state represents a finite surface excitation time, cubic saturation, and the
change in photoelectron escape probability caused by the debris potential.

Set `enabled=false` for an exact multiplier of one. To recover the original
IICSI trajectory, also set `photoemission.enabled=false` and
`external_debris_charge_driver.enabled=false`; no base PIC equation is replaced.
"""
default_photoemission_charge_coupling_config() = (
    enabled=false,
    frequency_omega_pe=0.022,
    phase_lag=pi/2,
    modulation_depth=0.80,
    ramp_start_time=50.0,
    ramp_duration=100.0,
    surface_response_time_omega_pe=20.0,
    nonlinear_saturation=0.30,
    potential_baseline_time_omega_pe=100.0,
    potential_scale_normalized=0.50,
    potential_coupling_strength=0.35,
    minimum_rate_multiplier=0.20,
    maximum_rate_multiplier=1.80,
)

photoemission_charge_coupling_config(raw=nothing)=
    isnothing(raw) ? default_photoemission_charge_coupling_config() :
    merge(default_photoemission_charge_coupling_config(),raw)

mutable struct PhotoemissionChargeCouplingState
    initialized::Bool
    surface_modulation::Float64
    potential_baseline::Float64
    escape_gate::Float64
    multiplier::Float64
end
PhotoemissionChargeCouplingState()=
    PhotoemissionChargeCouplingState(false,0.0,0.0,1.0,1.0)

function photo_charge_coupling_ramp(time,cfg)
    time<=cfg.ramp_start_time && return 0.0
    cfg.ramp_duration==0 && return 1.0
    u=clamp((time-cfg.ramp_start_time)/cfg.ramp_duration,0.0,1.0)
    return u^2*(3-2u)
end

"Advance the surface memory and potential-dependent escape gate."
function advance_photoemission_charge_coupling!(state,grid,pic_cfg,cfg,time,dt)
    cfg.enabled || return 1.0
    center=pic_cfg.debris.center_fraction*pic_cfg.domain_length
    potential=scalar_at_position(grid.phi,center,grid,pic_cfg)
    if !state.initialized
        state.initialized=true
        state.potential_baseline=potential
    end
    state.potential_baseline += dt*(potential-state.potential_baseline)/
        cfg.potential_baseline_time_omega_pe

    envelope=photo_charge_coupling_ramp(time,cfg)
    phase=2pi*cfg.frequency_omega_pe*time-cfg.phase_lag
    drive=envelope*cfg.modulation_depth*sin(phase)
    x=state.surface_modulation
    # First-order surface relaxation with a cubic, non-emissive saturation.
    rhs=(drive-x-cfg.nonlinear_saturation*x^3)/
        cfg.surface_response_time_omega_pe
    state.surface_modulation=x+dt*rhs

    # A more-positive-than-baseline grain recaptures a larger fraction of its
    # cold photoelectrons; a more-negative grain assists their escape.
    exponent=-cfg.potential_coupling_strength*
        (potential-state.potential_baseline)/cfg.potential_scale_normalized
    state.escape_gate=clamp(exp(clamp(exponent,-20.0,20.0)),
                            cfg.minimum_rate_multiplier,
                            cfg.maximum_rate_multiplier)
    raw=(1+state.surface_modulation)*state.escape_gate
    state.multiplier=clamp(raw,cfg.minimum_rate_multiplier,
                           cfg.maximum_rate_multiplier)
    return state.multiplier
end

photoemission_charge_coupling_history_path(history_file)=
    joinpath(dirname(history_file),"photoemission_charge_coupling.tsv")

function append_photoemission_charge_coupling_history!(path,step,time,state;
        new_file=false)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\tsurface_modulation\tpotential_baseline\tescape_gate\tphoto_rate_multiplier")
        println(io,join((step,time,state.surface_modulation,
            state.potential_baseline,state.escape_gate,state.multiplier),'\t'))
    end
    return nothing
end
