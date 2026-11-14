# Apple-GPU implementation. Particle arrays remain on the Metal device between
# steps. Only the small nx-element density and electric-field grids cross the
# unified-memory boundary because the existing FFT solve is tiny and CPU based.

mutable struct MetalWorkspace
    density::Metal.MtlVector{Float32,Metal.PrivateStorage}
    electric_field::Metal.MtlVector{Float32,Metal.PrivateStorage}
end

function prepare_metal_workspace(grid::Grid1D)
    return MetalWorkspace(Metal.zeros(Float32, grid.nx),
                          Metal.zeros(Float32, grid.nx))
end

function to_metal(species::Species)
    return Species(species.name, species.charge, species.mass, species.density,
                   species.weight, Metal.mtl(Float32.(species.x)),
                   Metal.mtl(Float32.(species.v)))
end

function to_cpu(species::Species)
    return Species(species.name, species.charge, species.mass, species.density,
                   species.weight, Float64.(Array(species.x)),
                   Float64.(Array(species.v)))
end

function metal_deposit_kernel!(density, particle_x, weight_over_dx, dx, nx)
    p = Metal.thread_position_in_grid().x
    if p <= length(particle_x)
        coordinate = particle_x[p] / dx
        left = floor(Int32, coordinate) + Int32(1)
        fraction_right = coordinate - floor(coordinate)
        right = left == nx ? Int32(1) : left + Int32(1)
        @atomic density[left] += weight_over_dx * (Float32(1) - fraction_right)
        @atomic density[right] += weight_over_dx * fraction_right
    end
    return
end

function metal_push_kernel!(x, v, field, half_dt_q_over_m, dt, dx, length_x, nx)
    p = Metal.thread_position_in_grid().x
    if p <= length(x)
        coordinate = x[p] / dx
        left = floor(Int32, coordinate) + Int32(1)
        fraction_right = coordinate - floor(coordinate)
        right = left == nx ? Int32(1) : left + Int32(1)
        field_p = field[left] * (Float32(1) - fraction_right) +
                  field[right] * fraction_right
        new_v = v[p] + half_dt_q_over_m * field_p
        new_x = x[p] + dt * new_v
        new_x -= floor(new_x / length_x) * length_x # Periodic boundary.
        v[p] = new_v
        x[p] = new_x
    end
    return
end

function metal_kick_kernel!(x, v, field, half_dt_q_over_m, dx, nx)
    p = Metal.thread_position_in_grid().x
    if p <= length(x)
        coordinate = x[p] / dx
        left = floor(Int32, coordinate) + Int32(1)
        fraction_right = coordinate - floor(coordinate)
        right = left == nx ? Int32(1) : left + Int32(1)
        field_p = field[left] * (Float32(1) - fraction_right) +
                  field[right] * fraction_right
        v[p] += half_dt_q_over_m * field_p
    end
    return
end

function launch_geometry(number)
    threads = 256
    return threads, cld(number, threads)
end

function metal_deposit_density!(host_density, species, grid, workspace)
    fill!(workspace.density, 0f0)
    threads, groups = launch_geometry(length(species.x))
    Metal.@metal threads=threads groups=groups metal_deposit_kernel!(
        workspace.density, species.x, Float32(species.weight / grid.dx),
        Float32(grid.dx), Int32(grid.nx))
    copyto!(host_density, Array(workspace.density))
    return nothing
end

function update_charge_density_metal!(grid, electrons, ions, workspace;
                                      immobile_ion_background=false)
    metal_deposit_density!(grid.ne, electrons, grid, workspace)
    if immobile_ion_background
        fill!(grid.ni, ions.density)
    else
        metal_deposit_density!(grid.ni, ions, grid, workspace)
    end
    @. grid.rho = electrons.charge * grid.ne + ions.charge * grid.ni
    grid.rho .-= sum(grid.rho) / grid.nx
    return nothing
end

function copy_field_to_metal!(workspace, grid)
    copyto!(workspace.electric_field, Float32.(grid.electric_field))
    return nothing
end

function push_species_metal!(species, grid, dt, workspace)
    threads, groups = launch_geometry(length(species.x))
    Metal.@metal threads=threads groups=groups metal_push_kernel!(
        species.x, species.v, workspace.electric_field,
        Float32(0.5 * dt * species.charge / species.mass), Float32(dt),
        Float32(grid.dx), Float32(grid.length), Int32(grid.nx))
    return nothing
end

function finish_velocity_kick_metal!(species, grid, dt, workspace)
    threads, groups = launch_geometry(length(species.x))
    Metal.@metal threads=threads groups=groups metal_kick_kernel!(
        species.x, species.v, workspace.electric_field,
        Float32(0.5 * dt * species.charge / species.mass), Float32(grid.dx),
        Int32(grid.nx))
    return nothing
end
