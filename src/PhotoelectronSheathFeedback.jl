# =============================================================================
# OPTIONAL PHOTOELECTRON-SHEATH ESCAPE FEEDBACK
# =============================================================================
# The ordinary photoemitter already launches real electron macroparticles;
# those particles enter Poisson's equation and may escape or be recollected by
# the debris.  This plug-in supplies only the unresolved, sub-cell surface
# response: recently emitted low-energy electrons form a finite-memory sheath
# which reduces the escape probability of subsequent photoelectrons.
#
# With enabled=false every function below is an exact deterministic no-op.
# The established PIC collection, reservoir, sponge and photoelectron launch
# paths are unchanged.

default_photoelectron_sheath_feedback_config() = (
    enabled=false,
    model=:space_charge_return,
    activation_time_omega_pe=300.0,
    residence_time_omega_pe=8.0,
    return_time_omega_pe=4.0,
    capacity_electrons=1.0e7,
    trapped_fraction=0.50,
    potential_baseline_time_omega_pe=100.0,
    potential_scale_normalized=0.50,
    potential_return_coupling=0.75,
    space_charge_strength=1.0,
    minimum_escape_multiplier=0.10,
)

photoelectron_sheath_feedback_config(raw=nothing) = isnothing(raw) ?
    default_photoelectron_sheath_feedback_config() :
    merge(default_photoelectron_sheath_feedback_config(),raw)

mutable struct PhotoelectronSheathFeedbackState
    initialized::Bool
    cloud_electrons::Float64
    occupation::Float64
    potential_baseline::Float64
    potential_excursion::Float64
    escape_gate::Float64
    space_charge_factor::Float64
    multiplier::Float64
    emitted_electrons_this_step::Float64
    cumulative_emitted_electrons::Float64
end

PhotoelectronSheathFeedbackState() = PhotoelectronSheathFeedbackState(
    false,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0)

"Advance the unresolved near-surface cloud and return its net escape factor."
function advance_photoelectron_sheath_feedback!(state,grid,pic_cfg,cfg,time,dt)
    cfg.enabled || return 1.0
    center=pic_cfg.debris.center_fraction*pic_cfg.domain_length
    phi=scalar_at_position(grid.phi,center,grid,pic_cfg)
    if !state.initialized
        state.potential_baseline=phi
        state.initialized=true
    end
    if time<cfg.activation_time_omega_pe
        # Learn the charging-transient baseline without feeding it back. This
        # prevents an imposed initial Qd from being interpreted as a physical
        # photoelectron-return event.
        state.potential_baseline=phi
        state.potential_excursion=0.0
        state.cloud_electrons=0.0
        state.occupation=0.0
        state.escape_gate=1.0
        state.space_charge_factor=1.0
        state.multiplier=1.0
        state.emitted_electrons_this_step=0.0
        return 1.0
    end
    alpha=1-exp(-dt/cfg.potential_baseline_time_omega_pe)
    state.potential_baseline += alpha*(phi-state.potential_baseline)
    state.potential_excursion=phi-state.potential_baseline

    # A surface excursion toward positive potential attracts emitted electrons
    # back to the debris.  The logistic gate is normalized to unity at the
    # slowly varying floating potential and has no discontinuous clipping.
    argument=clamp(cfg.potential_return_coupling*state.potential_excursion/
                   cfg.potential_scale_normalized,-50.0,50.0)
    state.escape_gate=2/(1+exp(argument))

    # Cloud electrons leave either by escape or by return.  Both channels are
    # physical losses from the unresolved sheath state; real emitted PIC
    # particles continue to evolve and recollect through the ordinary solver.
    escape_rate=state.escape_gate/cfg.residence_time_omega_pe
    return_rate=(2-state.escape_gate)/cfg.return_time_omega_pe
    state.cloud_electrons *= exp(-dt*(escape_rate+return_rate))
    state.occupation=state.cloud_electrons/cfg.capacity_electrons

    # Smooth space-charge-limited escape.  The 3/2 power is the planar
    # space-charge scaling; its local Taylor expansion supplies nonlinear
    # saturation without inserting a polynomial or an upper hard clamp.
    state.space_charge_factor=(1+cfg.space_charge_strength*state.occupation)^(-1.5)
    state.multiplier=max(cfg.minimum_escape_multiplier,
        state.escape_gate*state.space_charge_factor)
    state.emitted_electrons_this_step=0.0
    return state.multiplier
end

# Compatibility helper used by disabled-path tests and archived callers.
advance_photoelectron_sheath_feedback!(state,grid,pic_cfg,cfg,dt)=
    cfg.enabled ? error("enabled photoelectron sheath requires the current simulation time") : 1.0

"Register the real photoelectrons launched by the existing PIC emitter."
function register_photoelectron_emission!(state,number,debris_cfg,cfg)
    cfg.enabled || return nothing
    emitted=number*debris_cfg.physical_particles_per_macroparticle
    state.emitted_electrons_this_step=emitted
    state.cumulative_emitted_electrons += emitted
    state.cloud_electrons += cfg.trapped_fraction*emitted
    state.occupation=state.cloud_electrons/cfg.capacity_electrons
    return nothing
end

photoelectron_sheath_feedback_history_path(history_file)=
    joinpath(dirname(history_file),"photoelectron_sheath_feedback.tsv")

function append_photoelectron_sheath_feedback_history!(path,step,time,state;
        new_file=false)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file && println(io,"step\ttime_omega_pe\tcloud_electrons\toccupation\t" *
            "potential_baseline\tpotential_excursion\tescape_gate\t" *
            "space_charge_factor\tescape_multiplier\t" *
            "emitted_electrons_this_step\tcumulative_emitted_electrons")
        println(io,join((step,time,state.cloud_electrons,state.occupation,
            state.potential_baseline,state.potential_excursion,state.escape_gate,
            state.space_charge_factor,state.multiplier,
            state.emitted_electrons_this_step,
            state.cumulative_emitted_electrons),'\t'))
    end
    return nothing
end
