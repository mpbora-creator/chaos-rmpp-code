"Create uniformly positioned particles with a one-dimensional Maxwellian velocity."
function initialize_species(name::Symbol, charge, mass, density, temperature,
                            drift, number, grid::Grid1D, rng;
                            perturbation_amplitude=0.0, perturbation_mode=1)
    number > 0 || error("Each species needs at least one simulation particle")
    temperature >= 0 || error("Temperature cannot be negative")

    # Stratified positions reduce random density noise: one particle is placed
    # randomly inside each of N equal bins.
    x = ((collect(0:number-1) .+ rand(rng, number)) ./ number) .* grid.length
    k = 2pi * perturbation_mode / grid.length
    @. x += perturbation_amplitude * sin(k * x)
    @. x = mod(x, grid.length)

    # In normalized units a Maxwellian has standard deviation sqrt(T/m).
    thermal_speed = sqrt(temperature / mass)
    v = drift .+ thermal_speed .* randn(rng, number)
    weight = density * grid.length / number
    return Species(name, charge, mass, density, weight, x, v)
end

"Move one species through one full time step using kick-drift-kick."
function push_species!(species::Species, grid::Grid1D, dt,
                       boundary::AbstractBoundaryCondition)
    field_at_particles = interpolate_field(grid.electric_field, species.x, grid)
    acceleration = species.charge / species.mass
    @pic_threads for p in eachindex(species.x)
        species.v[p] += 0.5 * dt * acceleration * field_at_particles[p]
        species.x[p] += dt * species.v[p]
    end
    apply_boundary!(species.x, grid.length, boundary)
    return nothing
end

"Apply the second half of the velocity kick after the new field is known."
function finish_velocity_kick!(species::Species, grid::Grid1D, dt)
    field_at_particles = interpolate_field(grid.electric_field, species.x, grid)
    acceleration = species.charge / species.mass
    @pic_threads for p in eachindex(species.v)
        species.v[p] += 0.5 * dt * acceleration * field_at_particles[p]
    end
    return nothing
end
