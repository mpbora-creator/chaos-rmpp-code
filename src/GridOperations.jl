"Deposit particle number density on the grid with cloud-in-cell (CIC) weights."
function deposit_density!(density_grid, species::Species, grid::Grid1D)
    fill!(density_grid, 0.0)
    # Each CPU thread writes to its own grid column, avoiding races. The small
    # columns are summed only after all particles have been deposited.
    # maxthreadid() also covers Julia's interactive thread pool, whose numeric
    # ID may be one larger than nthreads() on recent Julia versions.
    thread_grids = zeros(Float64, grid.nx, Threads.maxthreadid())
    @pic_threads for p in eachindex(species.x)
        xp = species.x[p]
        coordinate = xp / grid.dx
        left = floor(Int, coordinate) + 1
        fraction_right = coordinate - floor(coordinate)
        right = mod1(left + 1, grid.nx) # Wrap the final cell periodically.
        thread = Threads.threadid()
        thread_grids[left, thread] += species.weight * (1 - fraction_right) / grid.dx
        thread_grids[right, thread] += species.weight * fraction_right / grid.dx
    end
    @pic_threads for cell in eachindex(density_grid)
        total = 0.0
        @inbounds for thread in axes(thread_grids, 2)
            total += thread_grids[cell, thread]
        end
        density_grid[cell] = total
    end
    return nothing
end

"Interpolate a grid field to particle positions using the same CIC shape."
function interpolate_field(field, particle_x, grid::Grid1D)
    result = similar(particle_x)
    @pic_threads for p in eachindex(particle_x)
        coordinate = particle_x[p] / grid.dx
        left = floor(Int, coordinate) + 1
        fraction_right = coordinate - floor(coordinate)
        right = mod1(left + 1, grid.nx)
        result[p] = field[left] * (1 - fraction_right) +
                    field[right] * fraction_right
    end
    return result
end

function update_charge_density!(grid::Grid1D, electrons::Species, ions::Species;
                                immobile_ion_background=false)
    deposit_density!(grid.ne, electrons, grid)
    if immobile_ion_background
        # A background is represented exactly on the grid: it does not move and
        # has no particle noise. This is the usual elementary two-stream model.
        fill!(grid.ni, ions.density)
    else
        deposit_density!(grid.ni, ions, grid)
    end
    @. grid.rho = electrons.charge * grid.ne + ions.charge * grid.ni
    # A periodic Poisson equation requires zero net charge. Removing the tiny
    # mean (normally sampling/roundoff error) is equivalent to setting k=0 to 0.
    grid.rho .-= sum(grid.rho) / grid.nx
    return nothing
end
