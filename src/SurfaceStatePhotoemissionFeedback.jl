# =============================================================================
# OPTIONAL IMPACT-DRIVEN SURFACE-STATE PHOTOEMISSION FEEDBACK
# =============================================================================
# Physical interpretation: incident electrons populate a finite surface-state
# reservoir.  The population relaxes and undergoes non-emissive cubic
# quenching.  Its occupation changes the yield of the already existing
# photoelectron source, including the escape probability set by the resolved
# debris potential.  It never adds charge or particles directly.
#
# With enabled=false, every entry point immediately returns a multiplier of
# one, performs no particle loop, changes no state, and consumes no RNG.

default_surface_state_photoemission_feedback_config() = (
    enabled=false,
    source_model=:expected_incident_energy,
    relaxation_time_omega_pe=5.0,
    capacity_electrons=2.0e9,
    deposition_efficiency=0.50,
    reference_impact_energy_ev=100.0,
    nonlinear_quenching_rate=1.0,
    photoemission_gain=4.0,
    potential_baseline_time_omega_pe=50.0,
    potential_scale_normalized=0.50,
    potential_escape_coupling=0.30,
    minimum_rate_multiplier=0.50,
    maximum_rate_multiplier=1.50,
)

surface_state_photoemission_feedback_config(raw=nothing) =
    isnothing(raw) ? default_surface_state_photoemission_feedback_config() :
    merge(default_surface_state_photoemission_feedback_config(),raw)

mutable struct SurfaceStatePhotoemissionFeedbackState
    initialized::Bool
    occupation::Float64
    potential_baseline::Float64
    expected_impact_energy_ev::Float64
    deposited_electrons::Float64
    nonemissive_loss_electrons::Float64
    escape_gate::Float64
    multiplier::Float64
end

SurfaceStatePhotoemissionFeedbackState() =
    SurfaceStatePhotoemissionFeedbackState(false,0.0,0.0,0.0,0.0,0.0,1.0,1.0)

"Expected captured impact energy from crossings during the current PIC step."
function expected_debris_impact_energy_ev(electrons,old_x,debris_cfg,center,
        radius,reference_temperature_ev)
    left_edge,right_edge=center-radius,center+radius
    total=0.0
    @inbounds for p in eachindex(electrons.x)
        crossed=max(old_x[p],electrons.x[p])>=left_edge &&
                min(old_x[p],electrons.x[p])<=right_edge
        crossed||continue
        kinetic=0.5*electrons.mass*electrons.v[p]^2
        total += capture_probability(kinetic,debris_cfg.electron_capture)*
                 kinetic*reference_temperature_ev
    end
    total
end

"Advance the physical surface reservoir and return its photo-rate multiplier."
function advance_surface_state_photoemission_feedback!(state,electrons,old_x,
        debris,grid,pic_cfg,cfg,dt)
    cfg.enabled||return 1.0
    center=pic_cfg.debris.center_fraction*pic_cfg.domain_length
    radius=pic_cfg.debris.size_debye
    potential=scalar_at_position(grid.phi,center,grid,pic_cfg)
    if !state.initialized
        state.initialized=true
        state.potential_baseline=potential
    end
    state.potential_baseline += dt*(potential-state.potential_baseline)/
        cfg.potential_baseline_time_omega_pe

    expected_energy=expected_debris_impact_energy_ev(electrons,old_x,
        pic_cfg.debris,center,radius,
        pic_cfg.secondary_electron_emission.reference_temperature_ev)
    proposed=cfg.deposition_efficiency*expected_energy/
        cfg.reference_impact_energy_ev*
        pic_cfg.debris.physical_particles_per_macroparticle
    accepted=proposed*clamp(1-state.occupation,0.0,1.0)
    state.expected_impact_energy_ev=expected_energy
    state.deposited_electrons=accepted

    occupied_electrons=state.occupation*cfg.capacity_electrons+accepted
    after_relaxation=occupied_electrons*
        exp(-dt/cfg.relaxation_time_omega_pe)
    occupation_after_relaxation=after_relaxation/cfg.capacity_electrons
    occupation_after_quenching=cfg.nonlinear_quenching_rate==0 ?
        occupation_after_relaxation :
        occupation_after_relaxation/sqrt(1+2cfg.nonlinear_quenching_rate*dt*
                                          occupation_after_relaxation^2)
    state.nonemissive_loss_electrons=max(0.0,
        occupied_electrons-occupation_after_quenching*cfg.capacity_electrons)
    state.occupation=clamp(occupation_after_quenching,0.0,1.0)

    exponent=-cfg.potential_escape_coupling*
        (potential-state.potential_baseline)/cfg.potential_scale_normalized
    state.escape_gate=exp(clamp(exponent,-20.0,20.0))
    raw=1+cfg.photoemission_gain*state.occupation*state.escape_gate
    state.multiplier=clamp(raw,cfg.minimum_rate_multiplier,
                           cfg.maximum_rate_multiplier)
    state.multiplier
end

surface_state_photoemission_feedback_history_path(history_file)=
    joinpath(dirname(history_file),"surface_state_photoemission_feedback.tsv")

function append_surface_state_photoemission_feedback_history!(path,step,time,
        state;new_file=false)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        new_file&&println(io,
            "step\ttime_omega_pe\toccupation\texpected_impact_energy_ev_this_step\t" *
            "deposited_electrons_this_step\tnonemissive_loss_electrons_this_step\t" *
            "potential_baseline\tescape_gate\tphoto_rate_multiplier")
        println(io,join((step,time,state.occupation,
            state.expected_impact_energy_ev,state.deposited_electrons,
            state.nonemissive_loss_electrons,state.potential_baseline,
            state.escape_gate,state.multiplier),'\t'))
    end
    nothing
end
