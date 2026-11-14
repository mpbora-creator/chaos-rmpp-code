"""Opt-in dielectric surface-polarization feedback for a collecting debris.

`polarization_charge_e` is the signed elementary-charge excess placed on the
right face; the left face receives its exact opposite.  It therefore changes
the local sheath without changing the net debris or system charge.
"""
mutable struct DielectricPolarizationState
    polarization_charge_e::Float64
    electron_source_this_step_e::Float64
    ion_source_this_step_e::Float64
    cumulative_electron_source_e::Float64
    cumulative_ion_source_e::Float64
    cumulative_relaxation_loss_e::Float64
    cumulative_nonlinear_loss_e::Float64
    cumulative_capacity_loss_e::Float64
end

DielectricPolarizationState() = DielectricPolarizationState(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

default_dielectric_polarization_config() = (
    enabled=false,
    source_model=:expected_impact_asymmetry,
    relaxation_time_omega_pe=20.0,
    polarization_capacity_electrons=2.0e7,
    deposition_fraction=0.05,
    nonlinear_saturation_rate=0.10,
    modulate_see_yield=false,
    see_yield_coupling=0.15,
    minimum_yield_multiplier=0.80,
    maximum_yield_multiplier=1.20,
)

function dielectric_polarization_config(cfg)
    defaults = default_dielectric_polarization_config()
    isnothing(cfg) && return defaults
    return merge(defaults, cfg)
end

"Face-local trapped-charge correction to the existing material SEE curve."
function polarization_see_yield_multiplier(state::DielectricPolarizationState,
        cfg, face)
    (!cfg.enabled || !cfg.modulate_see_yield) && return 1.0
    face_charge_e = sign(face)*state.polarization_charge_e
    raw = 1-cfg.see_yield_coupling*face_charge_e/
          cfg.polarization_capacity_electrons
    return clamp(raw,cfg.minimum_yield_multiplier,
                    cfg.maximum_yield_multiplier)
end

function begin_dielectric_polarization_step!(state::DielectricPolarizationState)
    state.electron_source_this_step_e = 0.0
    state.ion_source_this_step_e = 0.0
    return nothing
end

"Accumulate the expected captured charge from one crossing, without RNG use."
function deposit_expected_polarization!(state::DielectricPolarizationState,
        cfg, species_charge, capture_probability, face,
        physical_particles_per_macroparticle)
    cfg.enabled || return nothing
    cfg.source_model === :expected_impact_asymmetry ||
        error("Unsupported dielectric polarization source_model=$(cfg.source_model)")
    # P is the charge excess on the right face. A right-face impact contributes
    # its signed charge to P; a left-face impact contributes the opposite sign.
    source_e = cfg.deposition_fraction * sign(face) * species_charge *
               capture_probability * physical_particles_per_macroparticle
    if species_charge < 0
        state.electron_source_this_step_e += source_e
    else
        state.ion_source_this_step_e += source_e
    end
    return nothing
end

"Apply impact drive, linear dielectric relaxation, and odd cubic saturation."
function advance_dielectric_polarization!(state::DielectricPolarizationState,
        cfg, dt)
    cfg.enabled || return nothing
    electron_source = state.electron_source_this_step_e
    ion_source = state.ion_source_this_step_e
    state.cumulative_electron_source_e += electron_source
    state.cumulative_ion_source_e += ion_source

    capacity = cfg.polarization_capacity_electrons
    driven = state.polarization_charge_e + electron_source + ion_source
    limited = clamp(driven, -capacity, capacity)
    state.cumulative_capacity_loss_e += abs(driven-limited)

    after_linear = limited * exp(-dt/cfg.relaxation_time_omega_pe)
    state.cumulative_relaxation_loss_e += abs(limited-after_linear)

    occupancy2 = (after_linear/capacity)^2
    after_nonlinear = cfg.nonlinear_saturation_rate == 0 ? after_linear :
        after_linear/sqrt(1 + 2cfg.nonlinear_saturation_rate*dt*occupancy2)
    state.cumulative_nonlinear_loss_e += abs(after_linear-after_nonlinear)
    state.polarization_charge_e = after_nonlinear
    return nothing
end

dielectric_polarization_history_path(history_file) =
    joinpath(dirname(history_file), "dielectric_polarization.tsv")

function append_dielectric_polarization_history!(path, step, time,
        state::DielectricPolarizationState, cfg; new_file=false)
    mkpath(dirname(path))
    open(path, new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\tpolarization_charge_e\t" *
            "electron_source_this_step_e\tion_source_this_step_e\t" *
            "cumulative_electron_source_e\tcumulative_ion_source_e\t" *
            "cumulative_relaxation_loss_e\tcumulative_nonlinear_loss_e\t" *
            "cumulative_capacity_loss_e\tsee_multiplier_left\t" *
            "see_multiplier_right")
        multiplier_left=polarization_see_yield_multiplier(state,cfg,-1.0)
        multiplier_right=polarization_see_yield_multiplier(state,cfg,1.0)
        println(io, join((step,time,state.polarization_charge_e,
            state.electron_source_this_step_e,state.ion_source_this_step_e,
            state.cumulative_electron_source_e,state.cumulative_ion_source_e,
            state.cumulative_relaxation_loss_e,
            state.cumulative_nonlinear_loss_e,
            state.cumulative_capacity_loss_e,multiplier_left,
            multiplier_right),'\t'))
    end
    return nothing
end
