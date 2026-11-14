# =============================================================================
# OPTIONAL ELECTRON-IMPACT SECONDARY ELECTRON EMISSION FROM DEBRIS
# =============================================================================
# This plug-in is called only when enabled=true. The disabled branch remains
# in the established collection routine and therefore consumes no additional
# random numbers or changes the legacy IICSI trajectory.

mutable struct SecondaryEmissionState
    macroparticle_remainder::Float64
    captured_primary_macroparticles::Int
    emitted_secondary_macroparticles::Int
    incident_energy_ev::Float64
end

SecondaryEmissionState() = SecondaryEmissionState(0.0, 0, 0, 0.0)

"State of the optional trap-mediated delayed SEE channel."
mutable struct DelayedSecondaryEmissionState
    trapped_left_macroparticles::Float64
    trapped_right_macroparticles::Float64
    emission_remainder_left::Float64
    emission_remainder_right::Float64
    deposited_excitation_macroparticles::Float64
    released_excitation_macroparticles::Float64
    nonemissive_loss_macroparticles::Float64
    emitted_secondary_macroparticles::Int
end

DelayedSecondaryEmissionState() = DelayedSecondaryEmissionState(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)

"Continuous surface-excitation state, stored in physical-electron units."
mutable struct ContinuousSurfaceEmissionState
    excitation_left_electrons::Float64
    excitation_right_electrons::Float64
    emission_remainder_left_electrons::Float64
    emission_remainder_right_electrons::Float64
    deposited_excitation_electrons::Float64
    deposited_impact_energy_ev::Float64
    released_excitation_electrons::Float64
    nonemissive_loss_electrons::Float64
    emitted_secondary_macroparticles::Int
end

ContinuousSurfaceEmissionState() = ContinuousSurfaceEmissionState(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)

default_continuous_surface_config() = (
    source_model=:captured_events,
    energy_weighted_deposition=true,
    reference_impact_energy_ev=300.0,
    conservative_emission_remainder=true,
)

default_secondary_emission_delay_config() = (
    enabled=false,
    model=:discrete_trap,
    delayed_fraction=0.15,
    relaxation_time_omega_pe=100.0,
    trap_capacity_electrons=1.0e10,
    release_efficiency=1.0,
    nonlinear_saturation_rate=0.01,
    continuous_surface=default_continuous_surface_config(),
)

default_secondary_emission_config() = (
    enabled=false,
    model=:electron_impact,
    yield_curve=:sternglass,
    yield_scale=1.0,
    peak_yield=1.0,
    peak_energy_ev=300.0,
    threshold_energy_ev=0.0,
    maximum_yield=5.0,
    emitted_temperature_ev=2.0,
    reference_temperature_ev=1.0,
    delay=default_secondary_emission_delay_config(),
)

function secondary_emission_delay_config(see)
    defaults = default_secondary_emission_delay_config()
    hasproperty(see, :delay) || return defaults
    delay = merge(defaults, see.delay)
    surface = hasproperty(see.delay, :continuous_surface) ?
        merge(defaults.continuous_surface, see.delay.continuous_surface) :
        defaults.continuous_surface
    merge(delay, (continuous_surface=surface,))
end

secondary_emission_delay_state(delay) = delay.model === :continuous_surface ?
    ContinuousSurfaceEmissionState() : DelayedSecondaryEmissionState()

uses_expected_impact_flux(delay, state) = delay.enabled &&
    delay.model === :continuous_surface &&
    state isa ContinuousSurfaceEmissionState &&
    delay.continuous_surface.source_model === :expected_impact_flux

"Energy-dependent escaping SEE yield before self-consistent sheath recollection."
function secondary_electron_yield(energy_ev, see)
    energy_ev <= see.threshold_energy_ev && return 0.0
    span = see.peak_energy_ev - see.threshold_energy_ev
    x = (energy_ev - see.threshold_energy_ev) / span
    # Sternglass-type one-peak curve: delta=delta_max at x=1, with a
    # linear low-energy rise and exponential high-energy falloff.
    raw = see.peak_yield * x * exp(1 - x)
    return clamp(see.yield_scale * raw, 0.0, see.maximum_yield)
end

sample_secondary_energy_ev(rng, see) =
    -see.emitted_temperature_ev * log(max(rand(rng), eps(Float64)))

function emit_secondary_macroparticles!(electrons, debris, debris_cfg, see,
        number, face, center, radius, domain_length, rng)
    launch_offset = max(1.0e-9, 1.0e-6*radius)
    position = clamp(center + face*(radius+launch_offset), 0.0, domain_length)
    for _ in 1:number
        emitted_energy_ev = sample_secondary_energy_ev(rng, see)
        speed = sqrt(2*emitted_energy_ev/see.reference_temperature_ev)
        push!(electrons.x, position)
        push!(electrons.v, face*speed)
    end
    debris.total_charge_normalized += number * electrons.weight
    debris.total_charge_e +=
        number * debris_cfg.physical_particles_per_macroparticle
    return number
end

function deposit_delayed_excitation!(state::DelayedSecondaryEmissionState,
        amount, face, capacity_macroparticles)
    amount <= 0 && return 0.0
    trapped = face < 0 ? state.trapped_left_macroparticles :
                         state.trapped_right_macroparticles
    accepted = amount * clamp(1 - trapped/capacity_macroparticles, 0.0, 1.0)
    if face < 0
        state.trapped_left_macroparticles += accepted
    else
        state.trapped_right_macroparticles += accepted
    end
    state.deposited_excitation_macroparticles += accepted
    state.nonemissive_loss_macroparticles += amount-accepted
    return accepted
end

function deposit_delayed_excitation!(state::DelayedSecondaryEmissionState,
        amount, face, delay, debris_cfg, energy_ev)
    capacity_macroparticles = delay.trap_capacity_electrons /
                              debris_cfg.physical_particles_per_macroparticle
    deposit_delayed_excitation!(state, amount, face, capacity_macroparticles)
end

"Deposit a continuous, optionally energy-weighted surface excitation."
function deposit_delayed_excitation!(state::ContinuousSurfaceEmissionState,
        yield_macroparticles, face, delay, debris_cfg, energy_ev)
    yield_macroparticles <= 0 && return 0.0
    surface = delay.continuous_surface
    energy_factor = surface.energy_weighted_deposition ?
        energy_ev/surface.reference_impact_energy_ev : 1.0
    proposed = yield_macroparticles *
        debris_cfg.physical_particles_per_macroparticle * energy_factor
    excitation = face < 0 ? state.excitation_left_electrons :
                            state.excitation_right_electrons
    accepted = proposed * clamp(1-excitation/delay.trap_capacity_electrons,
                                0.0, 1.0)
    if face < 0
        state.excitation_left_electrons += accepted
    else
        state.excitation_right_electrons += accepted
    end
    state.deposited_excitation_electrons += accepted
    state.deposited_impact_energy_ev += energy_ev
    state.nonemissive_loss_electrons += proposed-accepted
    accepted
end

"Advance one face reservoir through cubic quenching and delayed release."
function advance_surface_reservoir(trapped, remainder, delay, dt,
                                   capacity_macroparticles)
    occupancy = clamp(trapped/capacity_macroparticles, 0.0, 1.0)
    # Exact update for dN/dt=-k*N^3/Nmax^2. This is a bounded physical
    # quenching channel and supplies the nonlinear saturation independently
    # of the unchanged reduced SL/vdPM surrogate.
    after_nonlinear = delay.nonlinear_saturation_rate == 0 ? trapped :
        trapped/sqrt(1 + 2delay.nonlinear_saturation_rate*dt*occupancy^2)
    nonlinear_loss = trapped-after_nonlinear
    release_fraction = -expm1(-dt/delay.relaxation_time_omega_pe)
    released = after_nonlinear*release_fraction
    remaining = after_nonlinear-released
    expected_emission = delay.release_efficiency*released + remainder
    emitted = floor(Int, expected_emission)
    new_remainder = expected_emission-emitted
    nonemissive_loss = nonlinear_loss + (1-delay.release_efficiency)*released
    return remaining, new_remainder, released, nonemissive_loss, emitted
end

"Release trap-mediated delayed secondaries after advancing the surface state."
function emit_delayed_secondary_electrons!(electrons, debris, debris_cfg, see,
        state::DelayedSecondaryEmissionState, dt, center, radius,
        domain_length, rng)
    delay = secondary_emission_delay_config(see)
    delay.enabled || return 0
    capacity_macroparticles = delay.trap_capacity_electrons /
                              debris_cfg.physical_particles_per_macroparticle
    total_emitted = 0
    total_released = total_loss = 0.0
    for face in (-1.0, 1.0)
        trapped = face < 0 ? state.trapped_left_macroparticles :
                             state.trapped_right_macroparticles
        remainder = face < 0 ? state.emission_remainder_left :
                               state.emission_remainder_right
        remaining, new_remainder, released, loss, emitted =
            advance_surface_reservoir(trapped, remainder, delay, dt,
                                      capacity_macroparticles)
        if face < 0
            state.trapped_left_macroparticles = remaining
            state.emission_remainder_left = new_remainder
        else
            state.trapped_right_macroparticles = remaining
            state.emission_remainder_right = new_remainder
        end
        emit_secondary_macroparticles!(electrons, debris, debris_cfg, see,
            emitted, face, center, radius, domain_length, rng)
        total_emitted += emitted
        total_released += released
        total_loss += loss
    end
    state.released_excitation_macroparticles += total_released
    state.nonemissive_loss_macroparticles += total_loss
    state.emitted_secondary_macroparticles += total_emitted
    return total_emitted
end


"Advance and emit from the continuous physical-electron surface state."
function emit_delayed_secondary_electrons!(electrons, debris, debris_cfg, see,
        state::ContinuousSurfaceEmissionState, dt, center, radius,
        domain_length, rng)
    delay = secondary_emission_delay_config(see)
    delay.enabled || return 0
    delay.model === :continuous_surface ||
        error("Continuous surface state requires delay.model=:continuous_surface")
    weight = debris_cfg.physical_particles_per_macroparticle
    total_emitted = 0
    total_released = total_loss = 0.0
    for face in (-1.0, 1.0)
        excitation = face < 0 ? state.excitation_left_electrons :
                                state.excitation_right_electrons
        remainder = face < 0 ? state.emission_remainder_left_electrons :
                               state.emission_remainder_right_electrons
        occupancy = clamp(excitation/delay.trap_capacity_electrons, 0.0, 1.0)
        after_nonlinear = delay.nonlinear_saturation_rate == 0 ? excitation :
            excitation/sqrt(1+2delay.nonlinear_saturation_rate*dt*occupancy^2)
        nonlinear_loss = excitation-after_nonlinear
        release_fraction = -expm1(-dt/delay.relaxation_time_omega_pe)
        released = after_nonlinear*release_fraction
        remaining = after_nonlinear-released
        expected_electrons = delay.release_efficiency*released + remainder
        emitted = floor(Int, expected_electrons/weight)
        new_remainder = expected_electrons-emitted*weight
        if face < 0
            state.excitation_left_electrons = remaining
            state.emission_remainder_left_electrons = new_remainder
        else
            state.excitation_right_electrons = remaining
            state.emission_remainder_right_electrons = new_remainder
        end
        emit_secondary_macroparticles!(electrons, debris, debris_cfg, see,
            emitted, face, center, radius, domain_length, rng)
        total_emitted += emitted
        total_released += released
        total_loss += nonlinear_loss+(1-delay.release_efficiency)*released
    end
    state.released_excitation_electrons += total_released
    state.nonemissive_loss_electrons += total_loss
    state.emitted_secondary_macroparticles += total_emitted
    total_emitted
end

"""Collect primary electrons and emit charge-conserving secondaries.

Capture decisions are completed before any emission random numbers are drawn.
This preserves the established MCC collection statistics. Emitted particles
are advanced by the ordinary PIC loop, so sheath escape or recollection is
resolved self-consistently rather than folded into an analytic escape factor.
"""
function collect_electrons_with_secondary_emission!(electrons, old_x, debris,
        debris_cfg, see, state::SecondaryEmissionState,
        delayed_state, center, radius,
        domain_length, capture_rng, prompt_emission_rng,
        polarization_state=nothing, polarization_cfg=nothing,
        see_probe_multiplier=1.0,
        charge_collection_cfg=default_charge_dependent_collection_config(),
        charge_collection_state=ChargeDependentCollectionState(),
        charge_collection_sample=ChargeDependentCollectionSample(0.0,0.0,0.0))
    see.enabled || error("SEE collector called while secondary emission is disabled")
    keep = trues(length(electrons.x))
    impacts = Tuple{Float64,Float64}[] # (incident energy eV, impacted face sign)
    left_edge, right_edge = center-radius, center+radius
    model = debris_cfg.electron_capture
    delay = secondary_emission_delay_config(see)
    delayed_fraction = delay.enabled ? delay.delayed_fraction : 0.0
    expected_surface_source = uses_expected_impact_flux(delay, delayed_state)
    for p in eachindex(electrons.x)
        crossed = max(old_x[p], electrons.x[p]) >= left_edge &&
                  min(old_x[p], electrons.x[p]) <= right_edge
        crossed || continue
        kinetic_energy = 0.5 * electrons.mass * electrons.v[p]^2
        face = old_x[p] <= center ? -1.0 : 1.0
        probability = charge_dependent_capture_probability(
            capture_probability(kinetic_energy, model),kinetic_energy,
            electrons.charge,face,charge_collection_sample,
            charge_collection_cfg,charge_collection_state)
        energy_ev = kinetic_energy * see.reference_temperature_ev
        polarization_multiplier = isnothing(polarization_state) ||
            isnothing(polarization_cfg) ? 1.0 :
            polarization_see_yield_multiplier(polarization_state,
                                              polarization_cfg,face)
        yield_multiplier=polarization_multiplier*see_probe_multiplier
        if !isnothing(polarization_state) && !isnothing(polarization_cfg)
            deposit_expected_polarization!(polarization_state,
                polarization_cfg, electrons.charge, probability, face,
                debris_cfg.physical_particles_per_macroparticle)
        end
        if expected_surface_source
            expected_delayed_yield = probability * delayed_fraction *
                clamp(secondary_electron_yield(energy_ev, see) *
                      yield_multiplier,0.0,see.maximum_yield)
            deposit_delayed_excitation!(delayed_state,
                expected_delayed_yield, face, delay, debris_cfg, energy_ev)
        end
        if rand(capture_rng) < probability
            keep[p] = false
            push!(impacts, (energy_ev, face))
        end
    end
    electrons.x = electrons.x[keep]
    electrons.v = electrons.v[keep]

    emitted = 0
    for (energy_ev, face) in impacts
        state.captured_primary_macroparticles += 1
        state.incident_energy_ev += energy_ev
        polarization_multiplier = isnothing(polarization_state) ||
            isnothing(polarization_cfg) ? 1.0 :
            polarization_see_yield_multiplier(polarization_state,
                                              polarization_cfg,face)
        yield_multiplier=polarization_multiplier*see_probe_multiplier
        total_yield = clamp(secondary_electron_yield(energy_ev, see)*
                            yield_multiplier,0.0,see.maximum_yield)
        prompt_yield = (1-delayed_fraction)*total_yield
        delayed_yield = delayed_fraction*total_yield
        state.macroparticle_remainder += prompt_yield
        delay.enabled && !expected_surface_source &&
            deposit_delayed_excitation!(delayed_state,
            delayed_yield, face, delay, debris_cfg, energy_ev)
        number = floor(Int, state.macroparticle_remainder)
        state.macroparticle_remainder -= number
        emit_secondary_macroparticles!(electrons, debris, debris_cfg, see,
            number, face, center, radius, domain_length, prompt_emission_rng)
        emitted += number
    end

    state.emitted_secondary_macroparticles += emitted
    collected = length(impacts)
    return collected*electrons.charge*electrons.weight, collected, emitted
end

# Backward-compatible entry point for focused callers and legacy tests. Passing
# one RNG retains the original draw order exactly; the production PIC loop uses
# independent capture and prompt-emission streams.
function collect_electrons_with_secondary_emission!(electrons, old_x, debris,
        debris_cfg, see, state::SecondaryEmissionState,
        delayed_state, center, radius,
        domain_length, rng)
    collect_electrons_with_secondary_emission!(electrons, old_x, debris,
        debris_cfg, see, state, delayed_state, center, radius, domain_length,
        rng, rng)
end

function append_secondary_emission_delay_history!(path, step, time,
        state::ContinuousSurfaceEmissionState;
        prompt_emitted=0, delayed_emitted=0, new_file=false)
    mkpath(dirname(path))
    open(path, new_file ? "w" : "a") do io
        new_file && println(io,
            "step\ttime_omega_pe\texcitation_left_electrons\t" *
            "excitation_right_electrons\texcitation_total_electrons\t" *
            "prompt_emitted_this_step\tdelayed_emitted_this_step\t" *
            "cumulative_deposited_excitation_electrons\t" *
            "cumulative_deposited_impact_energy_ev\t" *
            "cumulative_released_excitation_electrons\t" *
            "cumulative_nonemissive_loss_electrons\t" *
            "cumulative_delayed_emitted_macroparticles")
        total = state.excitation_left_electrons+
                state.excitation_right_electrons
        println(io, join((step,time,state.excitation_left_electrons,
            state.excitation_right_electrons,total,prompt_emitted,
            delayed_emitted,state.deposited_excitation_electrons,
            state.deposited_impact_energy_ev,
            state.released_excitation_electrons,
            state.nonemissive_loss_electrons,
            state.emitted_secondary_macroparticles),'\t'))
    end
end

delayed_surface_summary(state::DelayedSecondaryEmissionState) = (
    excitation=state.trapped_left_macroparticles+
               state.trapped_right_macroparticles,
    unit="macroparticles",
    emitted=state.emitted_secondary_macroparticles,
    loss=state.nonemissive_loss_macroparticles,
)

delayed_surface_summary(state::ContinuousSurfaceEmissionState) = (
    excitation=state.excitation_left_electrons+
               state.excitation_right_electrons,
    unit="electrons",
    emitted=state.emitted_secondary_macroparticles,
    loss=state.nonemissive_loss_electrons,
)

secondary_emission_delay_history_path(history_file) =
    joinpath(dirname(history_file), "secondary_emission_delay.tsv")

function append_secondary_emission_delay_history!(path, step, time, state;
        prompt_emitted=0, delayed_emitted=0, new_file=false)
    mkpath(dirname(path))
    open(path, new_file ? "w" : "a") do io
        new_file && println(io, "step\ttime_omega_pe\ttrapped_left_macroparticles\t" *
            "trapped_right_macroparticles\ttrapped_total_macroparticles\t" *
            "prompt_emitted_this_step\tdelayed_emitted_this_step\t" *
            "cumulative_deposited_excitation\tcumulative_released_excitation\t" *
            "cumulative_nonemissive_loss\tcumulative_delayed_emitted")
        total = state.trapped_left_macroparticles +
                state.trapped_right_macroparticles
        println(io, join((step,time,state.trapped_left_macroparticles,
            state.trapped_right_macroparticles,total,prompt_emitted,
            delayed_emitted,state.deposited_excitation_macroparticles,
            state.released_excitation_macroparticles,
            state.nonemissive_loss_macroparticles,
            state.emitted_secondary_macroparticles),'\t'))
    end
end
