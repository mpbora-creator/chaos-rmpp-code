# =============================================================================
# OPEN-RESERVOIR PLASMA WITH A COLLECTING DEBRIS
# =============================================================================
# This module-sized source file is intentionally isolated from the original
# periodic solver. It follows the supplied reservoir-PIC paper: particles exit
# without wrapping; fresh inflow is statistically independent of outflow;
# Poisson uses finite-domain Dirichlet boundaries; no mean charge is removed.

mutable struct OpenGrid1D
    nx::Int
    length::Float64
    dx::Float64
    x::Vector{Float64}
    ne::Vector{Float64}
    ni::Vector{Float64}
    debris_rho::Vector{Float64}
    rho::Vector{Float64}
    phi::Vector{Float64}
    electric_field::Vector{Float64}
end

function OpenGrid1D(nx, length; periodic=false)
    nx >= 3 || error("Open-grid nx must be at least 3")
    dx = periodic ? length/nx : length/(nx-1)
    coordinates = periodic ? collect(0:dx:length-dx) : collect(range(0.0,length;length=nx))
    z() = zeros(Float64, nx)
    return OpenGrid1D(nx, length, dx, coordinates,
                      z(), z(), z(), z(), z(), z())
end

mutable struct ReservoirAccumulators
    electron_left::Float64
    electron_right::Float64
    ion_left::Float64
    ion_right::Float64
end

ReservoirAccumulators() = ReservoirAccumulators(0.0, 0.0, 0.0, 0.0)

mutable struct DebrisMarkers
    positions::Vector{Float64}       # Fixed forever: these are static markers.
    charge_fractions::Vector{Float64} # Sum to one; define the spatial profile.
    total_charge_normalized::Float64
    total_charge_e::Float64
    initial_charge_normalized::Float64
end

mutable struct OpenDiagnostics
    time::Vector{Float64}
    mode_time::Vector{Float64}
    electron_mode::Vector{Float64}
    ion_mode::Vector{Float64}
    field_mode::Vector{Float64}
    electron_kinetic::Vector{Float64}
    ion_kinetic::Vector{Float64}
    field_energy::Vector{Float64}
    debris_charge_e::Vector{Float64}
    charge_error::Vector{Float64}
end

OpenDiagnostics() = OpenDiagnostics(Float64[], Float64[], Float64[], Float64[],
                                    Float64[], Float64[], Float64[], Float64[],
                                    Float64[], Float64[])

# Accurate elementary approximation to the standard-normal cumulative function.
function normal_cdf_open(x)
    ax = abs(x)
    t = 1 / (1 + 0.2316419 * ax)
    density = exp(-0.5 * ax^2) / sqrt(2pi)
    tail = density * t * (0.319381530 + t * (-0.356563782 +
           t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    positive_value = 1 - tail
    return x >= 0 ? positive_value : 1 - positive_value
end

"Incident number flux from a drifting 1D Maxwellian reservoir."
function reservoir_flux(density, drift, thermal_speed, side::Symbol)
    if thermal_speed == 0
        return side === :left ? density * max(drift, 0) : density * max(-drift, 0)
    end
    a = drift / thermal_speed
    gaussian = exp(-0.5 * a^2) / sqrt(2pi)
    # Avoid subtracting two nearly equal tiny numbers in the counter-streaming
    # tail.  For |a|>5, phi(a)-|a|Q(|a|) is evaluated with its Mills expansion.
    against_drift = (side === :right && a > 5) || (side === :left && a < -5)
    if against_drift
        z = abs(a)
        tail_flux = thermal_speed * gaussian *
                    max(0.0, z^-2 - 3z^-4 + 15z^-6 - 105z^-8)
        return density * tail_flux
    end
    if side === :left
        return max(0.0, density * (thermal_speed * gaussian + drift * normal_cdf_open(a)))
    else
        return max(0.0, density * (thermal_speed * gaussian - drift * normal_cdf_open(-a)))
    end
end

"Draw from |v| times a half drifting Maxwellian by bounded rejection sampling."
function sample_flux_velocity(rng, drift, thermal_speed, side::Symbol)
    if thermal_speed == 0
        return side === :left ? max(drift, 0.0) : min(drift, 0.0)
    end
    incoming_sign = side === :left ? 1.0 : -1.0

    # When the reservoir drift points away from this boundary, drawing an
    # ordinary drifting Maxwellian and waiting for a velocity in the opposite
    # tail is catastrophically slow (for example, right-injected ions when the
    # plasma flows to the right).  A Rayleigh proposal already contains the
    # required flux factor s=|v|.  The remaining counter-drift exponential is
    # then a well-behaved acceptance probability.
    if incoming_sign * drift < 0
        counter_drift = abs(drift)
        while true
            speed = thermal_speed * sqrt(-2log(max(rand(rng), eps(Float64))))
            rand(rng) < exp(-counter_drift * speed / thermal_speed^2) &&
                return incoming_sign * speed
        end
    end

    speed_cap = abs(drift) + 8thermal_speed
    while true
        velocity = drift + thermal_speed * randn(rng)
        incoming = side === :left ? velocity > 0 : velocity < 0
        incoming || continue
        abs(velocity) <= speed_cap || continue
        rand(rng) < abs(velocity) / speed_cap && return velocity
    end
end

function open_deposit_density!(density, species::Species, grid::OpenGrid1D)
    fill!(density, 0.0)
    for xp in species.x
        0 <= xp <= grid.length || continue
        coordinate = xp / grid.dx
        left = min(floor(Int, coordinate) + 1, grid.nx)
        if left == grid.nx
            density[end] += species.weight / grid.dx
        else
            fraction_right = coordinate - floor(coordinate)
            density[left] += species.weight * (1 - fraction_right) / grid.dx
            density[left + 1] += species.weight * fraction_right / grid.dx
        end
    end
    # The nonperiodic mesh contains the physical endpoints x=0 and x=L.
    # Their nodal control volumes have width dx/2, whereas interior nodes have
    # width dx. Dividing endpoint charge by a full dx creates an artificial
    # half-density dip in an otherwise uniform reservoir plasma. Correcting
    # by the half-cell volume also makes the trapezoidal density integral equal
    # the total particle weight exactly.
    density[1] *= 2
    density[end] *= 2
    return nothing
end

function open_interpolate(field, particle_x, grid::OpenGrid1D)
    values = similar(particle_x)
    @pic_threads for p in eachindex(particle_x)
        coordinate = clamp(particle_x[p], 0, grid.length) / grid.dx
        left = min(floor(Int, coordinate) + 1, grid.nx)
        if left == grid.nx
            values[p] = field[end]
        else
            fraction_right = coordinate - floor(coordinate)
            values[p] = field[left] * (1 - fraction_right) +
                        field[left + 1] * fraction_right
        end
    end
    return values
end

"Deposit net debris charge plus an optional zero-net dielectric dipole."
function deposit_debris_markers!(grid, debris::DebrisMarkers,
        polarization_state=nothing, polarization_cfg=nothing,
        normalized_charge_per_e=0.0,external_charge_normalized=0.0)
    fill!(grid.debris_rho, 0.0)
    isempty(debris.positions) && return nothing
    for (xp, fraction) in zip(debris.positions, debris.charge_fractions)
        marker_charge = (debris.total_charge_normalized+
                         external_charge_normalized)*fraction
        coordinate = xp / grid.dx
        left = min(floor(Int, coordinate) + 1, grid.nx - 1)
        fraction_right = coordinate - floor(coordinate)
        grid.debris_rho[left] += marker_charge * (1 - fraction_right) / grid.dx
        grid.debris_rho[left + 1] += marker_charge * fraction_right / grid.dx
    end
    if !isnothing(polarization_state) && !isnothing(polarization_cfg) &&
       polarization_cfg.enabled && polarization_state.polarization_charge_e != 0
        center = 0.5(first(debris.positions)+last(debris.positions))
        left_total = sum(f for (x,f) in zip(debris.positions,
                                            debris.charge_fractions) if x < center)
        right_total = sum(f for (x,f) in zip(debris.positions,
                                             debris.charge_fractions) if x >= center)
        left_total > 0 && right_total > 0 ||
            error("Dielectric polarization requires debris markers on both faces")
        polarization_normalized =
            polarization_state.polarization_charge_e*normalized_charge_per_e
        for (xp, fraction) in zip(debris.positions, debris.charge_fractions)
            marker_charge = xp < center ?
                -polarization_normalized*fraction/left_total :
                 polarization_normalized*fraction/right_total
            coordinate = xp/grid.dx
            left = min(floor(Int, coordinate)+1, grid.nx-1)
            fraction_right = coordinate-floor(coordinate)
            grid.debris_rho[left] += marker_charge*(1-fraction_right)/grid.dx
            grid.debris_rho[left+1] += marker_charge*fraction_right/grid.dx
        end
    end
    return nothing
end

"Solve -phi''=rho with phi(0),phi(L) prescribed. No mean is subtracted."
function solve_open_poisson!(grid, phi_left, phi_right)
    n = grid.nx - 2
    lower = fill(-1.0, n - 1)
    diagonal = fill(2.0, n)
    upper = fill(-1.0, n - 1)
    rhs = grid.dx^2 .* copy(grid.rho[2:end-1])
    rhs[1] += phi_left
    rhs[end] += phi_right

    # Thomas tridiagonal algorithm.
    for j in 2:n
        factor = lower[j - 1] / diagonal[j - 1]
        diagonal[j] -= factor * upper[j - 1]
        rhs[j] -= factor * rhs[j - 1]
    end
    interior = similar(rhs)
    interior[end] = rhs[end] / diagonal[end]
    for j in n-1:-1:1
        interior[j] = (rhs[j] - upper[j] * interior[j + 1]) / diagonal[j]
    end
    grid.phi[1], grid.phi[end] = phi_left, phi_right
    grid.phi[2:end-1] .= interior
    grid.electric_field[1] = -(grid.phi[2] - grid.phi[1]) / grid.dx
    grid.electric_field[end] = -(grid.phi[end] - grid.phi[end-1]) / grid.dx
    @. grid.electric_field[2:end-1] =
        -(grid.phi[3:end] - grid.phi[1:end-2]) / (2grid.dx)
    return nothing
end

function update_open_fields!(grid, electrons, ions, debris, cfg,
        polarization_state=nothing, polarization_cfg=nothing,
        normalized_charge_per_e=0.0,external_charge_normalized=0.0)
    open_deposit_density!(grid.ne, electrons, grid)
    open_deposit_density!(grid.ni, ions, grid)
    deposit_debris_markers!(grid, debris, polarization_state, polarization_cfg,
                            normalized_charge_per_e,external_charge_normalized)
    @. grid.rho = grid.ni - grid.ne + grid.debris_rho
    solve_open_poisson!(grid, cfg.phi_left, cfg.phi_right)
    return nothing
end

"Periodic CIC deposit used by the unified CPU engine."
function periodic_deposit_density!(density, species, grid::OpenGrid1D)
    fill!(density, 0.0)
    for xp in species.x
        coordinate = mod(xp, grid.length) / grid.length * grid.nx
        left = floor(Int, coordinate) + 1
        fraction_right = coordinate - floor(coordinate)
        right = mod1(left + 1, grid.nx)
        density[left] += species.weight * (1-fraction_right) * grid.nx/grid.length
        density[right] += species.weight * fraction_right * grid.nx/grid.length
    end
    return nothing
end

function periodic_interpolate(field, particle_x, grid::OpenGrid1D)
    values = similar(particle_x)
    @pic_threads for p in eachindex(particle_x)
        coordinate = mod(particle_x[p], grid.length) / grid.length * grid.nx
        left = floor(Int, coordinate) + 1
        fraction_right = coordinate - floor(coordinate)
        right = mod1(left + 1, grid.nx)
        values[p] = field[left]*(1-fraction_right) + field[right]*fraction_right
    end
    return values
end

"Periodic FFT Poisson solve on the unified engine grid."
function solve_periodic_poisson!(grid::OpenGrid1D)
    rho_hat = fft(grid.rho)
    modes = iseven(grid.nx) ? vcat(0:grid.nx÷2, (-grid.nx÷2+1):-1) :
            vcat(0:(grid.nx-1)÷2, (-(grid.nx-1)÷2):-1)
    k = (2pi/grid.length) .* modes
    phi_hat = zeros(ComplexF64, grid.nx)
    phi_hat[2:end] .= rho_hat[2:end] ./ k[2:end].^2
    electric_hat = -im .* k .* phi_hat
    grid.phi .= real(ifft(phi_hat))
    grid.electric_field .= real(ifft(electric_hat))
    return nothing
end

function update_unified_fields!(grid, electrons, ions, debris, cfg,
        polarization_state=nothing, polarization_cfg=nothing,
        normalized_charge_per_e=0.0,external_charge_normalized=0.0)
    if cfg.boundary_condition === :periodic
        periodic_deposit_density!(grid.ne, electrons, grid)
        periodic_deposit_density!(grid.ni, ions, grid)
        deposit_debris_markers!(grid, debris, polarization_state,
                                polarization_cfg, normalized_charge_per_e,
                                external_charge_normalized)
        @. grid.rho = grid.ni - grid.ne + grid.debris_rho
        grid.rho .-= sum(grid.rho)/grid.nx # required by periodic Poisson
        solve_periodic_poisson!(grid)
    else
        update_open_fields!(grid, electrons, ions, debris, cfg,
                            polarization_state, polarization_cfg,
                            normalized_charge_per_e,
                            external_charge_normalized)
    end
    return nothing
end

unified_interpolate(field, x, grid, cfg) =
    cfg.boundary_condition === :periodic ? periodic_interpolate(field, x, grid) :
                                           open_interpolate(field, x, grid)

function initial_open_species(name, charge, mass, density, temperature, drift,
                              number, grid, rng)
    x = ((collect(0:number-1) .+ rand(rng, number)) ./ number) .* grid.length
    v = drift .+ sqrt(temperature / mass) .* randn(rng, number)
    return Species(name, charge, mass, density, density * grid.length / number, x, v)
end

"Apply independent sinusoidal displacements that seed electron and ion branches."
function apply_open_perturbation!(electrons, ions, grid, cfg)
    cfg.perturbation.enabled || return nothing
    pe = cfg.perturbation
    ke = 2pi * pe.electron_mode / grid.length
    ki = 2pi * pe.ion_mode / grid.length
    @. electrons.x += pe.electron_amplitude * sin(ke * electrons.x)
    @. ions.x += pe.ion_amplitude * sin(ki * ions.x)
    if cfg.boundary_condition === :periodic
        @. electrons.x = mod(electrons.x, grid.length)
        @. ions.x = mod(ions.x, grid.length)
    else
        # Open runs cannot wrap a displaced particle to the opposite reservoir.
        @. electrons.x = clamp(electrons.x, 0.0, grid.length)
        @. ions.x = clamp(ions.x, 0.0, grid.length)
    end
    return nothing
end

"Absorb particles beyond the finite box and return their signed normalized charge."
function absorb_outer_boundaries!(species, grid)
    keep = @. (species.x >= 0) & (species.x <= grid.length)
    removed = length(keep) - count(keep)
    escaped_charge = removed * species.charge * species.weight
    species.x = species.x[keep]
    species.v = species.v[keep]
    return escaped_charge
end

function capture_probability(kinetic_energy, model)
    kinetic_energy <= model.threshold_energy && return 0.0
    excess = kinetic_energy - model.threshold_energy
    return model.maximum_probability * (1 - exp(-excess / model.energy_scale))
end

"Remove trajectories intersecting the debris according to the MCC energy law."
function collect_at_debris!(species, old_x, debris, debris_cfg, center, radius,
        rng, polarization_state=nothing, polarization_cfg=nothing,
        charge_collection_cfg=default_charge_dependent_collection_config(),
        charge_collection_state=ChargeDependentCollectionState(),
        charge_collection_sample=ChargeDependentCollectionSample(0.0,0.0,0.0))
    model = species.charge < 0 ? debris_cfg.electron_capture : debris_cfg.ion_capture
    keep = trues(length(species.x))
    collected = 0
    left_edge, right_edge = center - radius, center + radius
    for p in eachindex(species.x)
        crossed = max(old_x[p], species.x[p]) >= left_edge &&
                  min(old_x[p], species.x[p]) <= right_edge
        crossed || continue
        kinetic_energy = 0.5 * species.mass * species.v[p]^2
        face = old_x[p] <= center ? -1.0 : 1.0
        probability = charge_dependent_capture_probability(
            capture_probability(kinetic_energy, model),kinetic_energy,
            species.charge,face,charge_collection_sample,
            charge_collection_cfg,charge_collection_state)
        if !isnothing(polarization_state) && !isnothing(polarization_cfg)
            deposit_expected_polarization!(polarization_state,
                polarization_cfg, species.charge, probability, face,
                debris_cfg.physical_particles_per_macroparticle)
        end
        if rand(rng) < probability
            keep[p] = false
            collected += 1
        end
    end
    species.x = species.x[keep]
    species.v = species.v[keep]
    return collected * species.charge * species.weight, collected
end

function inject_from_reservoir!(species, side, cfg, dt, accumulator_value, rng)
    thermal_speed = sqrt(species.name === :electrons ?
                         cfg.electron_temperature : cfg.ion_temperature / cfg.mass_ratio)
    sound_speed = 1 / sqrt(cfg.mass_ratio)
    drift = cfg.flow_sound_speed * sound_speed
    flux = reservoir_flux(species.density, drift, thermal_speed, side)
    accumulator_value += flux * dt / species.weight
    number = floor(Int, accumulator_value)
    number = max(number, 0)
    accumulator_value -= number
    for _ in 1:number
        velocity = sample_flux_velocity(rng, drift, thermal_speed, side)
        position = side === :left ? rand(rng) * velocity * dt :
                   cfg.domain_length - rand(rng) * abs(velocity) * dt
        push!(species.x, clamp(position, 0, cfg.domain_length))
        push!(species.v, velocity)
    end
    return accumulator_value, number * species.charge * species.weight
end

function inject_all_reservoirs!(electrons, ions, cfg, accum, dt, rngs)
    injected_charge = 0.0
    accum.electron_left, q = inject_from_reservoir!(electrons, :left, cfg, dt,
        accum.electron_left, rngs.reservoir_electron_left)
    injected_charge += q
    accum.ion_left, q = inject_from_reservoir!(ions, :left, cfg, dt,
        accum.ion_left, rngs.reservoir_ion_left)
    injected_charge += q
    if cfg.two_sided_thermal_injection
        accum.electron_right, q = inject_from_reservoir!(electrons, :right, cfg,
            dt, accum.electron_right, rngs.reservoir_electron_right)
        injected_charge += q
        accum.ion_right, q = inject_from_reservoir!(ions, :right, cfg, dt,
            accum.ion_right, rngs.reservoir_ion_right)
        injected_charge += q
    end
    return injected_charge
end

# Backward-compatible helper: callers passing one RNG retain the historical
# interleaved draw order. Production runs pass isolated process streams.
function inject_all_reservoirs!(electrons, ions, cfg, accum, dt,
                                rng::AbstractRNG)
    shared = (reservoir_electron_left=rng, reservoir_electron_right=rng,
              reservoir_ion_left=rng, reservoir_ion_right=rng)
    inject_all_reservoirs!(electrons, ions, cfg, accum, dt, shared)
end

"Stochastically relax buffer velocities toward the undisturbed reservoir."
function apply_reservoir_sponge!(species, cfg, dt, rng)
    cfg.sponge.enabled || return nothing
    buffer = cfg.sponge.buffer_length
    buffer > 0 || return nothing
    thermal_speed = sqrt(species.name === :electrons ?
                         cfg.electron_temperature : cfg.ion_temperature / cfg.mass_ratio)
    drift = cfg.flow_sound_speed / sqrt(cfg.mass_ratio)
    for p in eachindex(species.x)
        edge_depth = species.x[p] < buffer ? buffer - species.x[p] :
                     species.x[p] > cfg.domain_length-buffer ?
                     species.x[p] - (cfg.domain_length-buffer) : 0.0
        edge_depth <= 0 && continue
        alpha = cfg.sponge.alpha_max * (edge_depth / buffer)^cfg.sponge.exponent
        if rand(rng) < 1 - exp(-alpha * dt)
            species.v[p] = drift + thermal_speed * randn(rng)
        end
    end
    return nothing
end

function open_kinetic_energy(species)
    return 0.5 * species.weight * species.mass * sum(abs2, species.v)
end

"Project central density/field structure onto the two seeded spatial modes."
function record_open_modes!(state, time, grid, cfg)
    cfg.frequency_diagnostics.enabled || return nothing
    fraction = cfg.frequency_diagnostics.central_fraction
    half_width = 0.5 * fraction * grid.length
    center = 0.5 * grid.length
    selected = findall(x -> abs(x-center) <= half_width, grid.x)
    pe = cfg.perturbation
    electron_basis = cos.(2pi * pe.electron_mode .* grid.x[selected] ./ grid.length)
    ion_basis = cos.(2pi * pe.ion_mode .* grid.x[selected] ./ grid.length)
    field_basis = sin.(2pi * pe.electron_mode .* grid.x[selected] ./ grid.length)
    push!(state.mode_time, time)
    push!(state.electron_mode,
          2sum((grid.ne[selected] .- cfg.ambient_density) .* electron_basis) /
          length(selected))
    push!(state.ion_mode,
          2sum((grid.ni[selected] .- cfg.ambient_density) .* ion_basis) /
          length(selected))
    push!(state.field_mode, 2sum(grid.electric_field[selected] .* field_basis) /
          length(selected))
    return nothing
end

function open_power_spectrum(times, signal, time_window)
    selected = findall(t -> time_window[1] <= t <= time_window[2], times)
    length(selected) >= 8 || return Float64[], Float64[]
    t = times[selected]
    values = signal[selected]
    # Open reservoirs can introduce a slow secular drift in particle number.
    # Remove its best-fit straight line before the FFT so it cannot masquerade
    # as an ion-frequency spectral peak.
    tmean = sum(t) / length(t)
    vmean = sum(values) / length(values)
    denominator = sum((t .- tmean).^2)
    slope = denominator > 0 ?
            sum((t .- tmean) .* (values .- vmean)) / denominator : 0.0
    values = values .- (vmean .+ slope .* (t .- tmean))
    window = 0.5 .- 0.5 .* cos.(2pi .* (0:length(values)-1) ./ (length(values)-1))
    transformed = fft(values .* window)
    last_positive = fld(length(values), 2) + 1
    omega = (2pi / ((t[2]-t[1]) * length(values))) .* collect(0:last_positive-1)
    power = abs2.(transformed[1:last_positive])
    maximum(power) > 0 && (power ./= maximum(power))
    return omega, power
end

function spectral_peak(omega, power, frequency_range)
    selected = findall(w -> frequency_range[1] <= w <= frequency_range[2], omega)
    isempty(selected) && return NaN
    return omega[selected[argmax(power[selected])]]
end

function late_debris_statistics(state, cfg)
    isempty(state.time) && return NaN, NaN
    start_time = (1-cfg.equilibrium_tail_fraction) * state.time[end]
    selected = findall(t -> t >= start_time, state.time)
    values = state.debris_charge_e[selected]
    mean_charge = sum(values) / length(values)
    sigma = sqrt(sum((values .- mean_charge).^2) / length(values))
    return mean_charge, sigma
end

function record_open_diagnostics!(state, time, grid, electrons, ions, debris,
                                  initial_total, cumulative_injected,
                                  cumulative_escaped, held_collection,
                                  cumulative_external_charge=0.0)
    push!(state.time, time)
    push!(state.electron_kinetic, open_kinetic_energy(electrons))
    push!(state.ion_kinetic, open_kinetic_energy(ions))
    push!(state.field_energy, 0.5 * grid.dx * sum(abs2, grid.electric_field))
    push!(state.debris_charge_e, debris.total_charge_e)
    actual = electrons.charge * electrons.weight * length(electrons.x) +
             ions.charge * ions.weight * length(ions.x) +
             debris.total_charge_normalized
    expected = initial_total + cumulative_injected - cumulative_escaped -
               held_collection + cumulative_external_charge
    push!(state.charge_error, actual - expected)
    return nothing
end

function plot_open_diagnostics(step, grid, electrons, ions, debris, state, cfg)
    panels = Plots.Plot[]
    if cfg.diagnostics.fields
        panel = plot(grid.x, grid.electric_field; label="E", xlabel="x/lambda_De",
                     ylabel="field", title="Fields")
        plot!(panel, grid.x, grid.phi; label="phi")
        cfg.debris.enabled && vline!(panel, [sum(debris.positions)/length(debris.positions)];
                                     label="debris", linestyle=:dash)
        push!(panels, panel)
    end
    if cfg.diagnostics.densities
        panel = plot(grid.x, grid.ne; label="n_e", xlabel="x/lambda_De",
                     ylabel="density", title="Plasma density")
        plot!(panel, grid.x, grid.ni; label="n_i")
        push!(panels, panel)
    end
    if cfg.diagnostics.phase_space
        stride, species = cfg.diagnostics.phase_space_stride, cfg.diagnostics.phase_space_species
        panel = plot(; xlabel="x/lambda_De", ylabel="v/v_te0", title="Phase space")
        species in (:electrons, :both) && scatter!(panel, electrons.x[1:stride:end],
            electrons.v[1:stride:end]; label="electrons", markersize=1,
            markerstrokewidth=0, alpha=0.5)
        species in (:ions, :both) && scatter!(panel, ions.x[1:stride:end],
            ions.v[1:stride:end]; label="ions", markersize=1,
            markerstrokewidth=0, alpha=0.5)
        push!(panels, panel)
    end
    if cfg.diagnostics.energies
        total = state.electron_kinetic .+ state.ion_kinetic .+ state.field_energy
        panel = plot(state.time, state.electron_kinetic; label="electron kinetic",
                     xlabel="t omega_pe", ylabel="energy", title="Energies")
        plot!(panel, state.time, state.ion_kinetic; label="ion kinetic")
        plot!(panel, state.time, state.field_energy; label="field")
        plot!(panel, state.time, total; label="total in box", linewidth=2)
        push!(panels, panel)
    end
    if cfg.diagnostics.debris_charge || cfg.diagnostics.charge_accounting
        panel = plot(; xlabel="t omega_pe", title="Charge diagnostics", legend=false)
        if cfg.diagnostics.debris_charge
            plot!(panel, state.time, state.debris_charge_e; label="Q_d/e",
                  ylabel="elementary charges")
            if length(state.time) >= 2
                mean_charge, sigma_charge = late_debris_statistics(state, cfg)
                hline!(panel, [mean_charge]; linestyle=:dash,
                       label="late mean=$(round(mean_charge,digits=2))e, sigma=$(round(sigma_charge,digits=2))e")
            end
        end
        cfg.diagnostics.charge_accounting && plot!(panel, state.time, state.charge_error;
            label="charge-accounting error", color=:red, linestyle=:dot)
        push!(panels, panel)
    end

    figure = isempty(panels) ? nothing : plot(panels...; layout=(length(panels),1),
             size=(950,300length(panels)), plot_title="$(cfg.boundary_condition) PIC, step $step")
    !isnothing(figure) && cfg.diagnostics.live_plots && display(figure)
    if !isnothing(figure) && cfg.diagnostics.save_plots
        mkpath(cfg.diagnostics.output_directory)
        savefig(figure, joinpath(cfg.diagnostics.output_directory,
                                "pic_$(lpad(step,6,'0')).png"))
    end
    return plot_open_frequency(step, state, cfg)
end

function plot_open_frequency(step, state, cfg)
    fd = cfg.frequency_diagnostics
    fd.enabled && cfg.perturbation.enabled || return nothing
    length(state.mode_time) >= 8 || return nothing

    omega_e, power_e = open_power_spectrum(state.mode_time, state.electron_mode,
                                            fd.electron_time_window)
    omega_i, power_i = open_power_spectrum(state.mode_time, state.ion_mode,
                                            fd.ion_time_window)
    if isempty(omega_e) || isempty(omega_i)
        return nothing
    end
    coherent_e = maximum(abs, state.electron_mode) >= fd.coherent_mode_threshold
    coherent_i = maximum(abs, state.ion_mode) >= fd.coherent_mode_threshold
    coherent_e || fill!(power_e, 0.0)
    coherent_i || fill!(power_i, 0.0)

    electron_peak = coherent_e ? spectral_peak(omega_e, power_e,
                                                fd.electron_frequency_range) : NaN
    ion_peak = coherent_i ? spectral_peak(omega_i, power_i,
                                          fd.ion_frequency_range) : NaN
    omega_pi = 1 / sqrt(cfg.mass_ratio)
    ki = 2pi * cfg.perturbation.ion_mode / cfg.domain_length
    warm_ion_prediction = sqrt((ki^2/(1+ki^2) +
                                3(cfg.ion_temperature/cfg.electron_temperature)*ki^2) /
                               cfg.mass_ratio)

    trace = plot(state.mode_time, state.electron_mode; label="electron density mode",
                 xlabel="t omega_pe", ylabel="mode amplitude",
                 title="Central mode traces")
    plot!(trace, state.mode_time, state.ion_mode; label="ion density mode")

    erange = (omega_e .>= fd.electron_frequency_range[1]) .&
             (omega_e .<= fd.electron_frequency_range[2])
    epanel = plot(omega_e[erange], power_e[erange]; label="electron-mode power",
                  xlabel="omega/omega_pe", ylabel="normalized power",
                  title="Electron branch peak=$(round(electron_peak,digits=4))")
    vline!(epanel, [1.0]; label="omega_pe=1", linestyle=:dash)

    irange = (omega_i .>= fd.ion_frequency_range[1]) .&
             (omega_i .<= fd.ion_frequency_range[2])
    ipanel = plot(omega_i[irange], power_i[irange]; label="ion-mode power",
                  xlabel="omega/omega_pe", ylabel="normalized power",
                  title="Ion branch peak=$(round(ion_peak,digits=5))")
    vline!(ipanel, [omega_pi]; label="omega_pi=$(round(omega_pi,digits=5))",
           linestyle=:dash)
    vline!(ipanel, [warm_ion_prediction];
           label="finite-k warm-ion=$(round(warm_ion_prediction,digits=5))",
           linestyle=:dot)

    figure = plot(trace, epanel, ipanel; layout=(3,1), size=(900,900),
                  plot_title="Open-reservoir frequency diagnostics, step $step")
    cfg.diagnostics.live_plots && display(figure)
    if cfg.diagnostics.save_plots
        mkpath(cfg.diagnostics.output_directory)
        savefig(figure, joinpath(cfg.diagnostics.output_directory,
                                "frequency_$(lpad(step,6,'0')).png"))
    end
    return (; electron_peak, ion_peak, omega_pi, warm_ion_prediction)
end

function instantaneous_charge_error(electrons, ions, debris, initial_total,
                                    cumulative_injected, cumulative_escaped,
                                    held_collection,cumulative_external_charge=0.0)
    actual = electrons.charge*electrons.weight*length(electrons.x) +
             ions.charge*ions.weight*length(ions.x) + debris.total_charge_normalized
    expected = initial_total + cumulative_injected - cumulative_escaped -
               held_collection + cumulative_external_charge
    return actual - expected
end

function append_time_history!(path, step, time, grid, electrons, ions, debris,
                              charge_error, cfg, sindy_options; new_file=false, collected_e=0,
                              collected_i=0, cumulative_injected=0.0,
                              cumulative_escaped=0.0)
    mkpath(dirname(path))
    open(path, new_file ? "w" : "a") do io
        new_file && println(io, "step\ttime_omega_pe\tNe\tNi\tcollected_e_this_step\tcollected_i_this_step\telectron_kinetic\tion_kinetic\tfield_energy\ttotal_energy\tQd_over_e\tcharge_error\tcumulative_injected_charge\tcumulative_escaped_charge\tmax_abs_E\tmean_ne\tmean_ni\tne_at_debris\tni_at_debris\tphi_at_debris\tE_at_debris\telectron_mean_v_local\tion_mean_v_local\telectron_temperature_local\tion_temperature_local\telectron_flux_local\tion_flux_local\tdne_dx_at_debris\tdni_dx_at_debris\tdphi_dx_at_debris\tdE_dx_at_debris\td2ne_dx2_at_debris\td2ni_dx2_at_debris\td2phi_dx2_at_debris\td2E_dx2_at_debris\tpotential_well_depth\tion_trapping_speed\tion_trapped_fraction\tion_trapped_density\tion_counterstreaming_index\telectron_inward_flux\tion_inward_flux\telectron_expected_capture_flux\tion_expected_capture_flux")
        ke = open_kinetic_energy(electrons)
        ki = open_kinetic_energy(ions)
        fe = 0.5*grid.dx*sum(abs2, grid.electric_field)
        center = cfg.debris.center_fraction*cfg.domain_length
        probe = sindy_options.local_probe_half_width
        em = local_species_moments(electrons,center,probe)
        im = local_species_moments(ions,center,probe)
        local_ne = scalar_at_position(grid.ne,center,grid,cfg)
        local_ni = scalar_at_position(grid.ni,center,grid,cfg)
        local_phi = scalar_at_position(grid.phi,center,grid,cfg)
        local_electric = scalar_at_position(grid.electric_field,center,grid,cfg)
        passes = sindy_options.spatial_derivative_smoothing_passes
        ne_derivatives = spatial_derivatives_at_position(grid.ne,center,grid,cfg,passes)
        ni_derivatives = spatial_derivatives_at_position(grid.ni,center,grid,cfg,passes)
        phi_derivatives = spatial_derivatives_at_position(grid.phi,center,grid,cfg,passes)
        electric_derivatives = spatial_derivatives_at_position(grid.electric_field,center,grid,cfg,passes)
        instability = local_instability_observables(electrons,ions,grid,cfg,center,
            sindy_options.instability_probe_half_width)
        println(io, join((step, time, length(electrons.x), length(ions.x),
                          collected_e, collected_i, ke, ki, fe, ke+ki+fe,
                          debris.total_charge_e, charge_error, cumulative_injected,
                          cumulative_escaped,maximum(abs,grid.electric_field),sum(grid.ne)/grid.nx,
                          sum(grid.ni)/grid.nx,local_ne,local_ni,local_phi,local_electric,
                          em.mean_velocity,im.mean_velocity,em.temperature,im.temperature,
                          em.density*em.mean_velocity,im.density*im.mean_velocity,
                          ne_derivatives.first,ni_derivatives.first,phi_derivatives.first,
                          electric_derivatives.first,ne_derivatives.second,ni_derivatives.second,
                          phi_derivatives.second,electric_derivatives.second,
                          instability.potential_well_depth,instability.ion_trapping_speed,
                          instability.ion_trapped_fraction,instability.ion_trapped_density,
                          instability.ion_counterstreaming_index,
                          instability.electron_inward_flux,instability.ion_inward_flux,
                          instability.electron_expected_capture_flux,
                          instability.ion_expected_capture_flux), '\t'))
    end
end

function save_checkpoint_atomic!(path, payload)
    mkpath(dirname(path))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        serialize(io, payload)
    end
    mv(temporary, path; force=true)
    return nothing
end

"Deterministic, checkpointable random streams isolated by physical process."
function simulation_rng_streams(seed)
    make_stream(id) = MersenneTwister(Int(seed) + 104_729*id)
    return (
        # Preserve the historical seed-to-initial-particle mapping exactly.
        initialization=MersenneTwister(Int(seed)),
        electron_capture=make_stream(2),
        ion_capture=make_stream(3),
        photoemission=make_stream(4),
        see_prompt=make_stream(5),
        see_delayed=make_stream(6),
        reservoir_electron_left=make_stream(7),
        reservoir_electron_right=make_stream(8),
        reservoir_ion_left=make_stream(9),
        reservoir_ion_right=make_stream(10),
        sponge_electron=make_stream(11),
        sponge_ion=make_stream(12),
        density_regulator_electron=make_stream(13),
        density_regulator_ion=make_stream(14),
    )
end

"Map a legacy single-RNG checkpoint without changing its continuation draws."
function legacy_shared_rng_streams(rng)
    return (
        initialization=rng,
        electron_capture=rng,
        ion_capture=rng,
        photoemission=rng,
        see_prompt=rng,
        see_delayed=rng,
        reservoir_electron_left=rng,
        reservoir_electron_right=rng,
        reservoir_ion_left=rng,
        reservoir_ion_right=rng,
        sponge_electron=rng,
        sponge_ion=rng,
        density_regulator_electron=rng,
        density_regulator_ion=rng,
    )
end

function run_boundary_plasma(parameters, cfg; restart=nothing, additional_steps=cfg.nsteps)
    boundary = make_boundary(cfg.boundary_condition) # Explicit separate BC object.
    photo = hasproperty(cfg, :photoemission) ? cfg.photoemission :
            default_photoemission_config()
    see = hasproperty(cfg, :secondary_electron_emission) ?
          cfg.secondary_electron_emission : default_secondary_emission_config()
    see_delay = secondary_emission_delay_config(see)
    polarization_cfg = dielectric_polarization_config(
        hasproperty(cfg, :dielectric_surface_polarization) ?
        cfg.dielectric_surface_polarization : nothing)
    see_probe_cfg=see_response_probe_config(
        hasproperty(cfg,:see_response_probe) ? cfg.see_response_probe : nothing)
    photo_bridge_cfg=photoemission_sl_vdpm_config(
        hasproperty(cfg,:photoemission_sl_vdpm_feedback) ?
        cfg.photoemission_sl_vdpm_feedback : nothing)
    periodic_photo_cfg=periodic_photoemission_driver_config(
        hasproperty(cfg,:periodic_photoemission_driver) ?
        cfg.periodic_photoemission_driver : nothing)
    photo_charge_coupling_cfg=photoemission_charge_coupling_config(
        hasproperty(cfg,:photoemission_charge_coupling) ?
        cfg.photoemission_charge_coupling : nothing)
    surface_photo_feedback_cfg=surface_state_photoemission_feedback_config(
        hasproperty(cfg,:surface_state_photoemission_feedback) ?
        cfg.surface_state_photoemission_feedback : nothing)
    photoelectron_sheath_cfg=photoelectron_sheath_feedback_config(
        hasproperty(cfg,:photoelectron_sheath_feedback) ?
        cfg.photoelectron_sheath_feedback : nothing)
    charge_collection_cfg=charge_dependent_collection_config(
        hasproperty(cfg,:charge_dependent_collection) ?
        cfg.charge_dependent_collection : nothing)
    external_charge_cfg=external_debris_charge_driver_config(
        hasproperty(cfg,:external_debris_charge_driver) ?
        cfg.external_debris_charge_driver : nothing)
    density_regulator_cfg=reservoir_density_regulator_config(
        hasproperty(cfg,:reservoir_density_regulator) ?
        cfg.reservoir_density_regulator : nothing)
    multimode_diagnostics_cfg=multimode_iicsi_diagnostics_config(
        hasproperty(cfg,:multimode_iicsi_diagnostics) ?
        cfg.multimode_iicsi_diagnostics : nothing)
    cfg.dt > 0 || error("time.dt must be positive")
    cfg.flow_sound_speed >= 0 || error("Use a nonnegative left-to-right flow_sound_speed")
    !cfg.debris.enabled || cfg.debris.size_debye > 0 ||
        error("Debris size_debye must be positive")
    !cfg.debris.enabled || cfg.debris.marker_count > 0 ||
        error("Debris marker_count must be positive")
    !cfg.debris.enabled || cfg.debris.charge_distribution === :gaussian ||
        error("Currently supported debris charge_distribution: :gaussian")
    !cfg.debris.enabled || cfg.debris.gaussian_sigma_fraction > 0 ||
        error("debris.gaussian_sigma_fraction must be positive")
    cfg.boundary_condition !== :open_reservoir || 2cfg.sponge.buffer_length < cfg.domain_length ||
        error("Two sponge buffers must fit inside the open domain")
    !cfg.debris.enabled || cfg.debris.collection_model === :mcc_energy ||
        error("Currently supported debris collection_model: :mcc_energy")
    if cfg.perturbation.enabled
        0 < cfg.perturbation.electron_mode < cfg.nx÷2 ||
            error("electron perturbation mode must lie below the grid Nyquist mode")
        0 < cfg.perturbation.ion_mode < cfg.nx÷2 ||
            error("ion perturbation mode must lie below the grid Nyquist mode")
    end
    parameters.compute_backend === :metal &&
        error("The unified debris/open engine currently requires compute_backend=:cpu")

    current_sindy_options=hasproperty(parameters,:sindy_export) ?
        parameters.sindy_export : default_sindy_export_options()
    saved_sindy_options=!isnothing(restart) && hasproperty(restart,:sindy_export) ?
        restart.sindy_export : current_sindy_options
    # Particle state, fields, random streams and cfg always come from the checkpoint.
    # Diagnostic export may be explicitly overridden because it does not alter
    # the resumed plasma dynamics. This permits recording a new sustained-state
    # E(x,t) interval without repeating the expensive transient.
    override_export=!isnothing(restart) &&
        hasproperty(current_sindy_options,:override_on_restart) &&
        current_sindy_options.override_on_restart
    sindy_options=override_export ? current_sindy_options : saved_sindy_options
    override_export && println("Restart diagnostic override enabled: using current sindy_export settings; saved plasma physics is unchanged.")
    sindy_options.local_probe_half_width > 0 ||
        error("sindy_export.local_probe_half_width must be positive")
    sindy_options.instability_probe_half_width > cfg.debris.size_debye ||
        error("sindy_export.instability_probe_half_width must exceed debris.size_debye")
    sindy_options.vlasov.every > 0 || error("sindy_export.vlasov.every must be positive")

    center = cfg.debris.center_fraction * cfg.domain_length
    radius = cfg.debris.size_debye
    if isnothing(restart)
        # Preserve the exact historical stochastic trajectory whenever the
        # delayed-SEE plug-in is inactive. Isolation is enabled only for the
        # delayed feedback whose event count would otherwise shift all draws.
        rng_streams = (see.enabled && see_delay.enabled) || photo_bridge_cfg.enabled ||
            density_regulator_cfg.enabled ?
            simulation_rng_streams(parameters.random_seed) :
            legacy_shared_rng_streams(MersenneTwister(parameters.random_seed))
        initialization_rng = rng_streams.initialization
        grid = OpenGrid1D(cfg.nx, cfg.domain_length;
                          periodic=cfg.boundary_condition === :periodic)
        sound_speed = 1 / sqrt(cfg.mass_ratio)
        drift = cfg.flow_sound_speed * sound_speed
        electrons = initial_open_species(:electrons, -1.0, 1.0, cfg.density_e,
            cfg.electron_temperature, drift, cfg.n_particles_e, grid,
            initialization_rng)
        ions = initial_open_species(:ions, 1.0, cfg.mass_ratio, cfg.density_i,
            cfg.ion_temperature, drift, cfg.n_particles_i, grid,
            initialization_rng)
        apply_open_perturbation!(electrons, ions, grid, cfg)

        if cfg.debris.enabled
            0 < center-radius && center+radius < cfg.domain_length ||
                error("Debris must lie strictly inside the domain")
            marker_positions = collect(range(center-radius, center+radius;
                                             length=cfg.debris.marker_count))
            sigma = cfg.debris.gaussian_sigma_fraction * radius
            marker_weights = @. exp(-0.5 * ((marker_positions-center)/sigma)^2)
            marker_weights ./= sum(marker_weights)
        else
            marker_positions = Float64[]
            marker_weights = Float64[]
        end
        charge_per_e = cfg.debris.enabled ?
                       electrons.weight/cfg.debris.physical_particles_per_macroparticle : 1.0
        initial_charge_e = cfg.debris.enabled ? cfg.debris.initial_charge_e : 0.0
        initial_q = initial_charge_e*charge_per_e
        debris = DebrisMarkers(marker_positions, marker_weights, initial_q,
                               initial_charge_e, initial_q)
        photoemission_state = PhotoemissionState()
        photoemission_sl_vdpm_state = PhotoemissionSLvdPMState(photo_bridge_cfg)
        photoemission_iicsi_filter_state = PhotoemissionIICSIFilterState()
        photoemission_charge_coupling_state=PhotoemissionChargeCouplingState()
        surface_state_photoemission_feedback_state=
            SurfaceStatePhotoemissionFeedbackState()
        photoelectron_sheath_feedback_state=PhotoelectronSheathFeedbackState()
        charge_dependent_collection_state=ChargeDependentCollectionState()
        secondary_emission_state = SecondaryEmissionState()
        delayed_secondary_emission_state =
            secondary_emission_delay_state(see_delay)
        dielectric_polarization_state = DielectricPolarizationState()
        external_debris_current_state=ExternalDebrisCurrentState()
        density_regulator_state=ReservoirDensityRegulatorState()
        accumulators = ReservoirAccumulators()
        diagnostics = OpenDiagnostics()
        update_unified_fields!(grid, electrons, ions, debris, cfg,
            dielectric_polarization_state, polarization_cfg, charge_per_e)
        record_open_modes!(diagnostics, 0.0, grid, cfg)
        initial_total = electrons.charge*electrons.weight*length(electrons.x) +
                        ions.charge*ions.weight*length(ions.x) + initial_q
        cumulative_injected = cumulative_escaped = held_collection = 0.0
        cumulative_external_charge=0.0
        completed_step = 0
        record_open_diagnostics!(diagnostics, 0.0, grid, electrons, ions, debris,
                                 initial_total, 0.0, 0.0, 0.0)
        frequency_summary = plot_open_diagnostics(0, grid, electrons, ions, debris,
                                                   diagnostics, cfg)
        if parameters.checkpoint.save_time_history
            append_time_history!(parameters.checkpoint.history_file, 0, 0.0, grid,
                electrons, ions, debris, 0.0,cfg,sindy_options;new_file=true)
        end
        if multimode_diagnostics_cfg.enabled
            append_multimode_iicsi_history!(
                multimode_iicsi_history_path(parameters.checkpoint.history_file),
                0,0.0,grid,multimode_diagnostics_cfg;new_file=true)
        end
        if see.enabled && see_delay.enabled
            append_secondary_emission_delay_history!(
                secondary_emission_delay_history_path(
                    parameters.checkpoint.history_file), 0, 0.0,
                delayed_secondary_emission_state; new_file=true)
        end
        if polarization_cfg.enabled
            append_dielectric_polarization_history!(
                dielectric_polarization_history_path(
                    parameters.checkpoint.history_file), 0, 0.0,
                dielectric_polarization_state, polarization_cfg; new_file=true)
        end
        if see_probe_cfg.enabled
            append_see_response_probe_history!(
                see_response_probe_history_path(parameters.checkpoint.history_file),
                0,0.0,see_response_probe_sample(see_probe_cfg,0.0);new_file=true)
        end
        if photo_bridge_cfg.enabled
            append_photoemission_sl_vdpm_history!(
                photoemission_sl_vdpm_history_path(
                    parameters.checkpoint.history_file), 0, 0.0,
                photoemission_sl_vdpm_state,photoemission_iicsi_filter_state;
                new_file=true)
        end
        if periodic_photo_cfg.enabled
            append_periodic_photoemission_history!(
                periodic_photoemission_history_path(
                    parameters.checkpoint.history_file),0,0.0,
                periodic_photoemission_sample(periodic_photo_cfg,0.0);
                new_file=true)
        end
        if photo_charge_coupling_cfg.enabled
            append_photoemission_charge_coupling_history!(
                photoemission_charge_coupling_history_path(
                    parameters.checkpoint.history_file),0,0.0,
                photoemission_charge_coupling_state;new_file=true)
        end
        if surface_photo_feedback_cfg.enabled
            append_surface_state_photoemission_feedback_history!(
                surface_state_photoemission_feedback_history_path(
                    parameters.checkpoint.history_file),0,0.0,
                surface_state_photoemission_feedback_state;new_file=true)
        end
        if photoelectron_sheath_cfg.enabled
            append_photoelectron_sheath_feedback_history!(
                photoelectron_sheath_feedback_history_path(
                    parameters.checkpoint.history_file),0,0.0,
                photoelectron_sheath_feedback_state;new_file=true)
        end
        if charge_collection_cfg.enabled
            initial_collection_sample=charge_dependent_collection_sample(
                grid,cfg,charge_collection_cfg)
            append_charge_dependent_collection_history!(
                charge_dependent_collection_history_path(
                    parameters.checkpoint.history_file),0,0.0,
                charge_dependent_collection_state,initial_collection_sample;
                new_file=true)
        end
        if external_charge_cfg.enabled
            sample=external_debris_charge_sample(external_charge_cfg,0.0,
                                                  charge_per_e)
            phi=scalar_at_position(grid.phi,center,grid,cfg)
            append_external_debris_charge_history!(
                external_debris_charge_history_path(
                    parameters.checkpoint.history_file),0,0.0,sample,
                debris.total_charge_e,phi;new_file=true)
        end
        if density_regulator_cfg.enabled
            append_reservoir_density_regulator_history!(
                reservoir_density_regulator_history_path(
                    parameters.checkpoint.history_file),0,0.0,electrons,ions,cfg,
                density_regulator_state,0,0.0;new_file=true)
        end
        if sindy_options.enabled && sindy_options.vlasov.enabled
            write_vlasov_snapshot!(sindy_options.vlasov.file,0,0.0,grid,electrons,
                                   ions,cfg,sindy_options.vlasov;new_file=true)
        end
    else
        rng_streams = hasproperty(restart, :rng_streams) ? restart.rng_streams :
            legacy_shared_rng_streams(restart.rng)
        grid = restart.grid
        electrons, ions, debris = restart.electrons, restart.ions, restart.debris
        photoemission_state = hasproperty(restart, :photoemission_state) ?
            restart.photoemission_state : PhotoemissionState()
        photoemission_sl_vdpm_state =
            hasproperty(restart, :photoemission_sl_vdpm_state) ?
            restart.photoemission_sl_vdpm_state :
            PhotoemissionSLvdPMState(photo_bridge_cfg)
        photoemission_iicsi_filter_state =
            hasproperty(restart,:photoemission_iicsi_filter_state) ?
            restart.photoemission_iicsi_filter_state :
            PhotoemissionIICSIFilterState()
        photoemission_charge_coupling_state=
            hasproperty(restart,:photoemission_charge_coupling_state) ?
            restart.photoemission_charge_coupling_state :
            PhotoemissionChargeCouplingState()
        surface_state_photoemission_feedback_state=
            hasproperty(restart,:surface_state_photoemission_feedback_state) ?
            restart.surface_state_photoemission_feedback_state :
            SurfaceStatePhotoemissionFeedbackState()
        photoelectron_sheath_feedback_state=
            hasproperty(restart,:photoelectron_sheath_feedback_state) ?
            restart.photoelectron_sheath_feedback_state :
            PhotoelectronSheathFeedbackState()
        charge_dependent_collection_state=
            hasproperty(restart,:charge_dependent_collection_state) ?
            restart.charge_dependent_collection_state :
            ChargeDependentCollectionState()
        secondary_emission_state = hasproperty(restart, :secondary_emission_state) ?
            restart.secondary_emission_state : SecondaryEmissionState()
        delayed_secondary_emission_state =
            hasproperty(restart, :delayed_secondary_emission_state) ?
            restart.delayed_secondary_emission_state :
            secondary_emission_delay_state(see_delay)
        dielectric_polarization_state =
            hasproperty(restart, :dielectric_polarization_state) ?
            restart.dielectric_polarization_state : DielectricPolarizationState()
        external_debris_current_state=
            hasproperty(restart,:external_debris_current_state) ?
            restart.external_debris_current_state : ExternalDebrisCurrentState()
        density_regulator_state=hasproperty(restart,:density_regulator_state) ?
            restart.density_regulator_state : ReservoirDensityRegulatorState()
        accumulators, diagnostics = restart.accumulators, restart.diagnostics
        initial_total = restart.initial_total
        cumulative_injected = restart.cumulative_injected
        cumulative_escaped = restart.cumulative_escaped
        held_collection = restart.held_collection
        cumulative_external_charge=
            hasproperty(restart,:cumulative_external_charge) ?
            restart.cumulative_external_charge : 0.0
        charge_per_e = restart.charge_per_e
        completed_step = restart.completed_step
        frequency_summary = nothing

        # Optional twin-trajectory diagnostic.  This acts only on a resumed
        # checkpoint at the explicitly requested step; particles, RNG streams,
        # accumulated surface states and driver phase remain bitwise identical.
        restart_charge_kick = hasproperty(parameters,
            :restart_debris_charge_perturbation) ?
            parameters.restart_debris_charge_perturbation :
            (enabled=false, delta_charge_e=0.0,
             apply_at_completed_step=completed_step)
        if restart_charge_kick.enabled
            restart_charge_kick.apply_at_completed_step == completed_step ||
                error("restart_debris_charge_perturbation was requested at step " *
                      "$(restart_charge_kick.apply_at_completed_step), but the " *
                      "loaded checkpoint is at step $completed_step")
            delta_charge_e = Float64(restart_charge_kick.delta_charge_e)
            isfinite(delta_charge_e) ||
                error("restart debris-charge perturbation must be finite")
            delta_charge_normalized = delta_charge_e * charge_per_e
            debris.total_charge_e += delta_charge_e
            debris.total_charge_normalized += delta_charge_normalized
            # Treat the deliberately imposed charge as part of the restarted
            # reference total so the ordinary conservation diagnostic remains
            # meaningful rather than reporting the diagnostic kick as an error.
            initial_total += delta_charge_normalized
            restart_time = completed_step * cfg.dt
            restart_external_sample = external_debris_charge_sample(
                external_charge_cfg, restart_time, charge_per_e)
            update_unified_fields!(grid, electrons, ions, debris, cfg,
                dielectric_polarization_state, polarization_cfg, charge_per_e,
                restart_external_sample.charge_normalized)
            println("Applied restart-only debris charge perturbation " *
                    "delta_Q/e=$delta_charge_e at step $completed_step, " *
                    "t=$restart_time")
        end
    end

    # When field export is first enabled on a continuation, record the exact
    # checkpoint state as the first frame before advancing it by one step.
    if !isnothing(restart) && sindy_options.enabled &&
       sindy_options.vlasov.enabled && !isfile(sindy_options.vlasov.file)
        write_vlasov_snapshot!(sindy_options.vlasov.file,completed_step,
            completed_step*cfg.dt,grid,electrons,ions,cfg,sindy_options.vlasov;
            new_file=true)
        println("Started sustained-state field snapshots at t=$(completed_step*cfg.dt): $(sindy_options.vlasov.file)")
    end
    delay_history_file = secondary_emission_delay_history_path(
        parameters.checkpoint.history_file)
    polarization_history_file = dielectric_polarization_history_path(
        parameters.checkpoint.history_file)
    see_probe_history_file=see_response_probe_history_path(
        parameters.checkpoint.history_file)
    photo_bridge_history_file=photoemission_sl_vdpm_history_path(
        parameters.checkpoint.history_file)
    periodic_photo_history_file=periodic_photoemission_history_path(
        parameters.checkpoint.history_file)
    photo_charge_coupling_history_file=
        photoemission_charge_coupling_history_path(
            parameters.checkpoint.history_file)
    surface_photo_feedback_history_file=
        surface_state_photoemission_feedback_history_path(
            parameters.checkpoint.history_file)
    photoelectron_sheath_history_file=
        photoelectron_sheath_feedback_history_path(
            parameters.checkpoint.history_file)
    charge_collection_history_file=charge_dependent_collection_history_path(
        parameters.checkpoint.history_file)
    external_charge_history_file=external_debris_charge_history_path(
        parameters.checkpoint.history_file)
    density_regulator_history_file=reservoir_density_regulator_history_path(
        parameters.checkpoint.history_file)
    multimode_diagnostics_history_file=multimode_iicsi_history_path(
        parameters.checkpoint.history_file)
    if !isnothing(restart) && see.enabled && see_delay.enabled &&
       !isfile(delay_history_file)
        append_secondary_emission_delay_history!(delay_history_file,
            completed_step, completed_step*cfg.dt,
            delayed_secondary_emission_state; new_file=true)
    end
    if !isnothing(restart) && polarization_cfg.enabled &&
       !isfile(polarization_history_file)
        append_dielectric_polarization_history!(polarization_history_file,
            completed_step, completed_step*cfg.dt,
            dielectric_polarization_state, polarization_cfg; new_file=true)
    end
    if !isnothing(restart) && see_probe_cfg.enabled &&
       !isfile(see_probe_history_file)
        time=completed_step*cfg.dt
        append_see_response_probe_history!(see_probe_history_file,
            completed_step,time,see_response_probe_sample(see_probe_cfg,time);
            new_file=true)
    end
    if !isnothing(restart) && photo_bridge_cfg.enabled &&
       !isfile(photo_bridge_history_file)
        append_photoemission_sl_vdpm_history!(photo_bridge_history_file,
            completed_step, completed_step*cfg.dt,
            photoemission_sl_vdpm_state,photoemission_iicsi_filter_state;
            new_file=true)
    end
    if !isnothing(restart) && periodic_photo_cfg.enabled &&
       !isfile(periodic_photo_history_file)
        time=completed_step*cfg.dt
        append_periodic_photoemission_history!(periodic_photo_history_file,
            completed_step,time,
            periodic_photoemission_sample(periodic_photo_cfg,time);
            new_file=true)
    end
    if !isnothing(restart) && external_charge_cfg.enabled &&
       !isfile(external_charge_history_file)
        time=completed_step*cfg.dt
        sample=external_debris_charge_sample(external_charge_cfg,time,
                                              charge_per_e)
        phi=scalar_at_position(grid.phi,center,grid,cfg)
        append_external_debris_charge_history!(external_charge_history_file,
            completed_step,time,sample,debris.total_charge_e,phi;new_file=true)
    end
    if !isnothing(restart) && surface_photo_feedback_cfg.enabled &&
       !isfile(surface_photo_feedback_history_file)
        append_surface_state_photoemission_feedback_history!(
            surface_photo_feedback_history_file,completed_step,
            completed_step*cfg.dt,surface_state_photoemission_feedback_state;
            new_file=true)
    end
    if !isnothing(restart) && photoelectron_sheath_cfg.enabled &&
       !isfile(photoelectron_sheath_history_file)
        append_photoelectron_sheath_feedback_history!(
            photoelectron_sheath_history_file,completed_step,
            completed_step*cfg.dt,photoelectron_sheath_feedback_state;
            new_file=true)
    end
    if !isnothing(restart) && charge_collection_cfg.enabled &&
       !isfile(charge_collection_history_file)
        sample=charge_dependent_collection_sample(grid,cfg,charge_collection_cfg)
        append_charge_dependent_collection_history!(
            charge_collection_history_file,completed_step,
            completed_step*cfg.dt,charge_dependent_collection_state,sample;
            new_file=true)
    end
    if !isnothing(restart) && multimode_diagnostics_cfg.enabled &&
       !isfile(multimode_diagnostics_history_file)
        append_multimode_iicsi_history!(multimode_diagnostics_history_file,
            completed_step,completed_step*cfg.dt,grid,multimode_diagnostics_cfg;
            new_file=true)
    end

    println("Starting $(cfg.boundary_condition) plasma: flow=$(cfg.flow_sound_speed) cs, " *
            "debris enabled=$(cfg.debris.enabled), Qd=$(debris.total_charge_e)e, " *
            "dynamic=$(cfg.debris.dynamic_charge), " *
            "photoemission=$(photo.enabled), " *
            "secondary_electron_emission=$(see.enabled), " *
            "delayed_secondary_emission=$(see.enabled && see_delay.enabled), " *
            "dielectric_polarization=$(polarization_cfg.enabled), " *
            "surface_state_photo_feedback=$(surface_photo_feedback_cfg.enabled), " *
            "photoelectron_sheath_feedback=$(photoelectron_sheath_cfg.enabled), " *
            "charge_dependent_collection=$(charge_collection_cfg.enabled), " *
            "multimode_iicsi_diagnostics=$(multimode_diagnostics_cfg.enabled), " *
            "SEE_response_probe=$(see_probe_cfg.enabled)")
    photo_bridge_cfg.enabled && println(
        "Photoemission SL-vdPM feedback enabled: measurement=$(photo_bridge_cfg.measurement), " *
        "tau rate=$(photo_bridge_cfg.tau_rate_per_omega_pe), " *
        "mean multiplier=$(photo_bridge_cfg.mean_rate_multiplier)")
    periodic_photo_cfg.enabled && println(
        "Periodic photoemission driver enabled: frequency=$(periodic_photo_cfg.frequency_omega_pe), " *
        "depth=$(periodic_photo_cfg.modulation_depth), " *
        "ramp=$(periodic_photo_cfg.ramp_start_time)+$(periodic_photo_cfg.ramp_duration)")
    photo_charge_coupling_cfg.enabled && println(
        "Photoemission-charge coupling enabled: frequency=$(photo_charge_coupling_cfg.frequency_omega_pe), " *
        "surface tau=$(photo_charge_coupling_cfg.surface_response_time_omega_pe), " *
        "phase lag=$(photo_charge_coupling_cfg.phase_lag)")
    if external_charge_cfg.enabled
        strength=external_charge_cfg.model===:sinusoidal_current ?
            "current amplitude=$(external_charge_cfg.current_amplitude_e_per_omega_pe)e/omega_pe" :
            "charge amplitude=$(external_charge_cfg.amplitude_e)e"
        println("External debris charge driver enabled: model=$(external_charge_cfg.model), " *
            "$strength, frequency=$(external_charge_cfg.frequency_omega_pe), " *
            "ramp=$(external_charge_cfg.ramp_start_time)+$(external_charge_cfg.ramp_duration)")
    end
    final_step = completed_step + additional_steps
    for step in completed_step+1:final_step
        begin_dielectric_polarization_step!(dielectric_polarization_state)
        time=step*cfg.dt
        see_probe_sample=see_response_probe_sample(see_probe_cfg,time)
        photo_rate_multiplier=advance_photoemission_sl_vdpm!(
            photoemission_sl_vdpm_state,photoemission_iicsi_filter_state,
            grid, cfg, photo_bridge_cfg, cfg.dt)
        periodic_photo_sample=periodic_photoemission_sample(periodic_photo_cfg,time)
        photo_charge_multiplier=advance_photoemission_charge_coupling!(
            photoemission_charge_coupling_state,grid,cfg,
            photo_charge_coupling_cfg,time,cfg.dt)
        total_photo_rate_multiplier=
            photo_rate_multiplier*periodic_photo_sample.multiplier*
            photo_charge_multiplier
        photoelectron_sheath_multiplier=
            advance_photoelectron_sheath_feedback!(
                photoelectron_sheath_feedback_state,grid,cfg,
                photoelectron_sheath_cfg,time,cfg.dt)
        total_photo_rate_multiplier*=photoelectron_sheath_multiplier
        old_e, old_i = copy(electrons.x), copy(ions.x)
        begin_charge_dependent_collection_step!(charge_dependent_collection_state)
        charge_collection_sample=charge_dependent_collection_sample(
            grid,cfg,charge_collection_cfg)
        ee = unified_interpolate(grid.electric_field, electrons.x, grid, cfg)
        ei = unified_interpolate(grid.electric_field, ions.x, grid, cfg)
        @. electrons.v += 0.5 * cfg.dt * electrons.charge / electrons.mass * ee
        @. ions.v += 0.5 * cfg.dt * ions.charge / ions.mass * ei
        @. electrons.x += cfg.dt * electrons.v
        @. ions.x += cfg.dt * ions.v

        # Measure the expected energy of the same debris-crossing electrons
        # seen by the MCC collector below. The resulting surface state changes
        # only the existing physical photoelectron rate.
        surface_photo_multiplier=
            advance_surface_state_photoemission_feedback!(
                surface_state_photoemission_feedback_state,electrons,old_e,
                debris,grid,cfg,surface_photo_feedback_cfg,cfg.dt)
        total_photo_rate_multiplier*=surface_photo_multiplier

        prompt_emitted = delayed_emitted = 0
        if cfg.debris.enabled
            if see.enabled
                qcoll_e, necoll, prompt_emitted =
                    collect_electrons_with_secondary_emission!(
                    electrons, old_e, debris, cfg.debris, see,
                    secondary_emission_state, delayed_secondary_emission_state,
                    center, radius,
                    cfg.domain_length, rng_streams.electron_capture,
                    rng_streams.see_prompt, dielectric_polarization_state,
                    polarization_cfg,see_probe_sample.multiplier,
                    charge_collection_cfg,charge_dependent_collection_state,
                    charge_collection_sample)
                delayed_emitted = emit_delayed_secondary_electrons!(electrons,
                    debris, cfg.debris, see, delayed_secondary_emission_state,
                    cfg.dt, center, radius, cfg.domain_length,
                    rng_streams.see_delayed)
            else
                qcoll_e, necoll = collect_at_debris!(electrons, old_e, debris,
                    cfg.debris, center, radius, rng_streams.electron_capture,
                    dielectric_polarization_state, polarization_cfg,
                    charge_collection_cfg,charge_dependent_collection_state,
                    charge_collection_sample)
            end
            qcoll_i, nicoll = collect_at_debris!(ions, old_i, debris,
                cfg.debris, center, radius, rng_streams.ion_capture,
                dielectric_polarization_state, polarization_cfg,
                charge_collection_cfg,charge_dependent_collection_state,
                charge_collection_sample)
        else
            qcoll_e, necoll = 0.0, 0
            qcoll_i, nicoll = 0.0, 0
        end
        collected_charge = qcoll_e + qcoll_i
        # Convert the actually transferred normalized charge back to physical
        # elementary charges. This remains correct if species counts differ.
        collected_charge_e = collected_charge / charge_per_e
        if cfg.debris.dynamic_charge
            debris.total_charge_normalized += collected_charge
            debris.total_charge_e += collected_charge_e
        else
            held_collection += collected_charge
            debris.total_charge_normalized = debris.initial_charge_normalized
        end
        external_charge_sample,external_delta_charge=
            advance_external_debris_current!(debris,
                external_debris_current_state,external_charge_cfg,time,cfg.dt,
                charge_per_e)
        cumulative_external_charge += external_delta_charge

        # Optional internal source. It runs before boundary absorption and
        # reservoir replenishment so emitted electrons obey the same open BC,
        # sponge, field solve and equilibrium-flow evolution as ambient ones.
        emitted_photoelectrons=emit_debris_photoelectrons!(electrons, debris,
            cfg.debris, photo, photoemission_state, cfg.dt, center,
            radius, cfg.domain_length, rng_streams.photoemission;
            rate_multiplier=total_photo_rate_multiplier)
        register_photoelectron_emission!(photoelectron_sheath_feedback_state,
            emitted_photoelectrons,cfg.debris,photoelectron_sheath_cfg)

        advance_dielectric_polarization!(dielectric_polarization_state,
                                         polarization_cfg, cfg.dt)

        if cfg.boundary_condition === :periodic
            apply_boundary!(electrons.x, grid.length, boundary)
            apply_boundary!(ions.x, grid.length, boundary)
        else
            cumulative_escaped += absorb_outer_boundaries!(electrons, grid)
            cumulative_escaped += absorb_outer_boundaries!(ions, grid)
            cumulative_injected += inject_all_reservoirs!(electrons, ions, cfg,
                accumulators, cfg.dt, rng_streams)
            regulator_charge,regulator_pairs,regulator_deficit =
                apply_reservoir_density_regulator!(electrons,ions,cfg,
                    density_regulator_cfg,density_regulator_state,cfg.dt,
                    rng_streams.density_regulator_electron,
                    rng_streams.density_regulator_ion)
            cumulative_injected += regulator_charge
            apply_reservoir_sponge!(electrons, cfg, cfg.dt,
                                    rng_streams.sponge_electron)
            apply_reservoir_sponge!(ions, cfg, cfg.dt,
                                    rng_streams.sponge_ion)
        end
        update_unified_fields!(grid, electrons, ions, debris, cfg,
            dielectric_polarization_state,polarization_cfg,charge_per_e,
            external_charge_sample.charge_normalized)
        record_open_modes!(diagnostics, step*cfg.dt, grid, cfg)
        ee = unified_interpolate(grid.electric_field, electrons.x, grid, cfg)
        ei = unified_interpolate(grid.electric_field, ions.x, grid, cfg)
        @. electrons.v += 0.5 * cfg.dt * electrons.charge / electrons.mass * ee
        @. ions.v += 0.5 * cfg.dt * ions.charge / ions.mass * ei

        charge_error = instantaneous_charge_error(electrons, ions, debris,
            initial_total, cumulative_injected, cumulative_escaped,
            held_collection,cumulative_external_charge)
        if parameters.checkpoint.save_time_history
            append_time_history!(parameters.checkpoint.history_file, step,
                step*cfg.dt,grid,electrons,ions,debris,charge_error,cfg,sindy_options;
                collected_e=necoll, collected_i=nicoll,
                cumulative_injected=cumulative_injected,
                cumulative_escaped=cumulative_escaped)
        end
        if see.enabled && see_delay.enabled
            append_secondary_emission_delay_history!(delay_history_file,
                step, step*cfg.dt, delayed_secondary_emission_state;
                prompt_emitted=prompt_emitted,
                delayed_emitted=delayed_emitted)
        end
        if polarization_cfg.enabled
            append_dielectric_polarization_history!(polarization_history_file,
                step, step*cfg.dt, dielectric_polarization_state,
                polarization_cfg)
        end
        if see_probe_cfg.enabled
            append_see_response_probe_history!(see_probe_history_file,
                step,time,see_probe_sample)
        end
        if photo_bridge_cfg.enabled
            append_photoemission_sl_vdpm_history!(photo_bridge_history_file,
                step,time,photoemission_sl_vdpm_state,
                photoemission_iicsi_filter_state)
        end
        if periodic_photo_cfg.enabled
            append_periodic_photoemission_history!(periodic_photo_history_file,
                step,time,periodic_photo_sample)
        end
        if photo_charge_coupling_cfg.enabled
            append_photoemission_charge_coupling_history!(
                photo_charge_coupling_history_file,step,time,
                photoemission_charge_coupling_state)
        end
        if surface_photo_feedback_cfg.enabled
            append_surface_state_photoemission_feedback_history!(
                surface_photo_feedback_history_file,step,time,
                surface_state_photoemission_feedback_state)
        end
        if photoelectron_sheath_cfg.enabled
            append_photoelectron_sheath_feedback_history!(
                photoelectron_sheath_history_file,step,time,
                photoelectron_sheath_feedback_state)
        end
        if charge_collection_cfg.enabled
            append_charge_dependent_collection_history!(
                charge_collection_history_file,step,time,
                charge_dependent_collection_state,charge_collection_sample)
        end
        if external_charge_cfg.enabled
            phi=scalar_at_position(grid.phi,center,grid,cfg)
            append_external_debris_charge_history!(external_charge_history_file,
                step,time,external_charge_sample,debris.total_charge_e,phi)
        end
        if density_regulator_cfg.enabled
            append_reservoir_density_regulator_history!(
                density_regulator_history_file,step,step*cfg.dt,electrons,ions,cfg,
                density_regulator_state,regulator_pairs,regulator_deficit)
        end
        if multimode_diagnostics_cfg.enabled
            append_multimode_iicsi_history!(multimode_diagnostics_history_file,
                step,time,grid,multimode_diagnostics_cfg)
        end
        if sindy_options.enabled && sindy_options.vlasov.enabled &&
           step % sindy_options.vlasov.every == 0
            write_vlasov_snapshot!(sindy_options.vlasov.file,step,step*cfg.dt,
                grid,electrons,ions,cfg,sindy_options.vlasov)
        end

        if step % cfg.diagnostics.every == 0 || step == final_step
            record_open_diagnostics!(diagnostics, step*cfg.dt, grid, electrons,
                ions, debris, initial_total, cumulative_injected,
                cumulative_escaped, held_collection,cumulative_external_charge)
            frequency_summary = plot_open_diagnostics(step, grid, electrons,
                ions, debris, diagnostics, cfg)
            println("  step $step: Ne=$(length(electrons.x)), " *
                    "Ni=$(length(ions.x)), Qd/e=$(debris.total_charge_e)")
        end

        # A full checkpoint is deliberately written only at the end of the
        # requested run segment. Per-step checkpointing would impose enormous
        # I/O overhead and is not needed for ordinary continuation runs.
        if parameters.checkpoint.enabled && step == final_step
            payload = (version=1, saved_at=Dates.now(), cfg=cfg,
                completed_step=step, completed_time=step*cfg.dt,
                rng=rng_streams.initialization, rng_streams=rng_streams,
                grid=grid, electrons=electrons, ions=ions, debris=debris,
                photoemission_state=photoemission_state,
                photoemission_sl_vdpm_state=photoemission_sl_vdpm_state,
                photoemission_iicsi_filter_state=photoemission_iicsi_filter_state,
                photoemission_charge_coupling_state=
                    photoemission_charge_coupling_state,
                surface_state_photoemission_feedback_state=
                    surface_state_photoemission_feedback_state,
                photoelectron_sheath_feedback_state=
                    photoelectron_sheath_feedback_state,
                charge_dependent_collection_state=
                    charge_dependent_collection_state,
                secondary_emission_state=secondary_emission_state,
                delayed_secondary_emission_state=
                    delayed_secondary_emission_state,
                dielectric_polarization_state=dielectric_polarization_state,
                external_debris_current_state=external_debris_current_state,
                density_regulator_state=density_regulator_state,
                accumulators=accumulators, diagnostics=diagnostics,
                initial_total=initial_total,
                cumulative_injected=cumulative_injected,
                cumulative_escaped=cumulative_escaped,
                held_collection=held_collection,charge_per_e=charge_per_e,
                cumulative_external_charge=cumulative_external_charge,
                sindy_export=sindy_options)
            save_checkpoint_atomic!(parameters.checkpoint.file, payload)
        end
    end
    mean_charge, sigma_charge = late_debris_statistics(diagnostics, cfg)
    if cfg.debris.enabled
        println("Late-time debris statistics: mean=$(mean_charge)e, sigma=$(sigma_charge)e")
    end
    if see.enabled
        captured = secondary_emission_state.captured_primary_macroparticles
        mean_impact = captured == 0 ? 0.0 :
            secondary_emission_state.incident_energy_ev/captured
        println("SEE statistics: captured primary macroparticles=$captured, " *
                "emitted secondary macroparticles=" *
                "$(secondary_emission_state.emitted_secondary_macroparticles), " *
                "mean incident energy=$(mean_impact) eV")
        if see_delay.enabled
            summary = delayed_surface_summary(delayed_secondary_emission_state)
            println("Delayed SEE statistics ($(see_delay.model)): surface excitation=" *
                "$(summary.excitation) $(summary.unit), emitted secondary " *
                "macroparticles=$(summary.emitted), nonemissive loss=" *
                "$(summary.loss) $(summary.unit)")
        end
    end
    if polarization_cfg.enabled
        println("Dielectric polarization statistics: P=" *
            "$(dielectric_polarization_state.polarization_charge_e)e, " *
            "capacity=$(polarization_cfg.polarization_capacity_electrons)e")
    end
    if surface_photo_feedback_cfg.enabled
        state=surface_state_photoemission_feedback_state
        println("Surface-state photoemission statistics: occupation=" *
            "$(state.occupation), multiplier=$(state.multiplier), " *
            "last expected impact energy=$(state.expected_impact_energy_ev) eV")
    end
    if photoelectron_sheath_cfg.enabled
        state=photoelectron_sheath_feedback_state
        println("Photoelectron-sheath statistics: occupation=$(state.occupation), " *
            "escape gate=$(state.escape_gate), space-charge factor=" *
            "$(state.space_charge_factor), multiplier=$(state.multiplier)")
    end
    if charge_collection_cfg.enabled
        state=charge_dependent_collection_state
        println("Charge-dependent collection statistics: electron candidates=" *
            "$(state.cumulative_electron_candidates), ion candidates=" *
            "$(state.cumulative_ion_candidates)")
    end
    !isnothing(frequency_summary) && println("Measured branches: electron=" *
        "$(frequency_summary.electron_peak), ion=$(frequency_summary.ion_peak); " *
        "references omega_pe=1, omega_pi=$(frequency_summary.omega_pi)")
    return (; grid, electrons, ions, debris, rng_streams,
            photoemission_state=photoemission_state,
            photoemission_sl_vdpm_state=photoemission_sl_vdpm_state,
            photoemission_iicsi_filter_state=photoemission_iicsi_filter_state,
            surface_state_photoemission_feedback_state=
                surface_state_photoemission_feedback_state,
            photoelectron_sheath_feedback_state=
                photoelectron_sheath_feedback_state,
            charge_dependent_collection_state=charge_dependent_collection_state,
            secondary_emission_state=secondary_emission_state,
            delayed_secondary_emission_state=
                delayed_secondary_emission_state,
            dielectric_polarization_state=dielectric_polarization_state,
            diagnostic_state=diagnostics)
end
