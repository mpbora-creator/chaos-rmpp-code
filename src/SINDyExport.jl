# Data export for the separate integral-SINDy project. The PIC remains only a
# data producer; equation discovery lives in its own sibling directory.

default_sindy_export_options() = (enabled=false,override_on_restart=false,
    local_probe_half_width=1.0,
    instability_probe_half_width=4.0,
    spatial_derivative_smoothing_passes=2,
    vlasov=(enabled=false,every=1,file="diagnostics/vlasov_snapshots.jls",
            nx_bins=64,nv_bins=64,electron_velocity_range=(-5.0,5.0),
            ion_velocity_range=(-0.1,0.1)))

function local_species_moments(species, center, half_width)
    selected = findall(x -> abs(x-center) <= half_width, species.x)
    isempty(selected) && return (density=0.0, mean_velocity=0.0, temperature=0.0)
    velocities = species.v[selected]
    mean_velocity = sum(velocities)/length(velocities)
    variance = sum((velocities .- mean_velocity).^2)/length(velocities)
    density = length(selected)*species.weight/(2half_width)
    return (; density, mean_velocity, temperature=species.mass*variance)
end

"""Compact observables of the debris-pinned ion phase-space structure.

The diagnostic region is `center ± half_width`. `ion_trapped_fraction` uses
the single-particle energy relative to the mean potential at the two edges of
that region. `ion_counterstreaming_index` is zero for a one-way population and
one for equal positive/negative ion flux. The inward and expected-capture
fluxes use |v|, not the signed bulk flux n<u>, so particles approaching the
debris from both sides contribute.
"""
function local_instability_observables(electrons,ions,grid,cfg,center,half_width)
    half_width>0 || error("instability_probe_half_width must be positive")
    left=max(0.0,center-half_width); right=min(grid.length,center+half_width)
    width=right-left
    width>0 || error("IICSI diagnostic region has zero width")
    phi_reference=0.5*(scalar_at_position(grid.phi,left,grid,cfg)+
                       scalar_at_position(grid.phi,right,grid,cfg))
    phi_center=scalar_at_position(grid.phi,center,grid,cfg)
    potential_well_depth=max(phi_reference-phi_center,0.0)

    function incident_fluxes(species,model)
        selected=findall(x->left<=x<=right,species.x)
        isempty(selected) && return (inward=0.0,expected=0.0,
                                     positive=0.0,negative=0.0)
        inward=expected=positive=negative=0.0
        for p in selected
            speed=abs(species.v[p])
            kinetic=0.5*species.mass*species.v[p]^2
            inward += speed
            expected += speed*capture_probability(kinetic,model)
            positive += max(species.v[p],0.0)
            negative += max(-species.v[p],0.0)
        end
        factor=species.weight/width
        return (inward=inward*factor,expected=expected*factor,
                positive=positive*factor,negative=negative*factor)
    end

    electron_flux=incident_fluxes(electrons,cfg.debris.electron_capture)
    ion_flux=incident_fluxes(ions,cfg.debris.ion_capture)
    flux_sum=ion_flux.positive+ion_flux.negative
    counterstreaming=flux_sum>0 ?
        2min(ion_flux.positive,ion_flux.negative)/flux_sum : 0.0

    selected_ions=findall(x->left<=x<=right,ions.x)
    trapped=0
    if !isempty(selected_ions) && potential_well_depth>0
        local_phi=unified_interpolate(grid.phi,ions.x[selected_ions],grid,cfg)
        for (j,p) in enumerate(selected_ions)
            relative_energy=0.5*ions.mass*ions.v[p]^2+
                            ions.charge*(local_phi[j]-phi_reference)
            relative_energy<0 && (trapped+=1)
        end
    end
    trapped_fraction=isempty(selected_ions) ? 0.0 : trapped/length(selected_ions)
    trapped_density=trapped*ions.weight/width
    trapping_speed=sqrt(2potential_well_depth/ions.mass)
    return (;potential_well_depth,ion_trapping_speed=trapping_speed,
        ion_trapped_fraction=trapped_fraction,ion_trapped_density=trapped_density,
        ion_counterstreaming_index=counterstreaming,
        electron_inward_flux=electron_flux.inward,ion_inward_flux=ion_flux.inward,
        electron_expected_capture_flux=electron_flux.expected,
        ion_expected_capture_flux=ion_flux.expected)
end

scalar_at_position(values, position, grid, cfg) =
    unified_interpolate(values, [position], grid, cfg)[1]

"""Return first and second spatial derivatives on the simulation grid."""
function grid_spatial_derivatives(values, grid, cfg)
    number = length(values)
    number >= 4 || error("At least four grid points are required for SINDy spatial derivatives")
    first_derivative = similar(values, Float64)
    second_derivative = similar(values, Float64)
    dx = grid.dx
    if cfg.boundary_condition === :periodic
        for j in eachindex(values)
            left, right = mod1(j-1,number), mod1(j+1,number)
            first_derivative[j] = (values[right]-values[left])/(2dx)
            second_derivative[j] = (values[right]-2values[j]+values[left])/dx^2
        end
    else
        @views first_derivative[2:end-1] .= (values[3:end].-values[1:end-2])./(2dx)
        @views second_derivative[2:end-1] .= (values[3:end].-2values[2:end-1].+values[1:end-2])./dx^2
        first_derivative[1] = (-3values[1]+4values[2]-values[3])/(2dx)
        first_derivative[end] = (3values[end]-4values[end-1]+values[end-2])/(2dx)
        second_derivative[1] = (2values[1]-5values[2]+4values[3]-values[4])/dx^2
        second_derivative[end] = (2values[end]-5values[end-1]+4values[end-2]-values[end-3])/dx^2
    end
    return first_derivative, second_derivative
end

function smooth_grid_signal(values, passes, cfg)
    passes >= 0 || error("spatial_derivative_smoothing_passes cannot be negative")
    result = Float64.(values)
    for _ in 1:passes
        old = copy(result)
        if cfg.boundary_condition === :periodic
            for j in eachindex(old)
                result[j] = (old[mod1(j-1,length(old))]+2old[j]+old[mod1(j+1,length(old))])/4
            end
        else
            @views result[2:end-1] .= (old[1:end-2].+2old[2:end-1].+old[3:end])./4
            result[1], result[end] = old[1], old[end]
        end
    end
    result
end

function spatial_derivatives_at_position(values, position, grid, cfg, smoothing_passes)
    smoothed = smooth_grid_signal(values,smoothing_passes,cfg)
    first_derivative, second_derivative = grid_spatial_derivatives(smoothed,grid,cfg)
    return (first=scalar_at_position(first_derivative,position,grid,cfg),
            second=scalar_at_position(second_derivative,position,grid,cfg))
end

"Histogram particles into f(x,v), with integral f dx dv equal to particle number."
function phase_space_histogram(species, grid, nx_bins, nv_bins, velocity_range)
    vmin, vmax = velocity_range
    vmin < vmax || error("SINDy velocity_range must be increasing")
    dx, dv = grid.length/nx_bins, (vmax-vmin)/nv_bins
    distribution = zeros(Float64, nx_bins, nv_bins)
    for p in eachindex(species.x)
        # Positions are already wrapped for PBC. For an open grid, x=L belongs
        # to the final bin rather than being incorrectly mapped back to x=0.
        position = clamp(species.x[p],0.0,prevfloat(grid.length))
        xbin = clamp(floor(Int,position/dx)+1,1,nx_bins)
        velocity = species.v[p]
        vmin <= velocity < vmax || continue
        vbin = clamp(floor(Int,(velocity-vmin)/dv)+1,1,nv_bins)
        distribution[xbin,vbin] += species.weight/(dx*dv)
    end
    x = (collect(1:nx_bins) .- 0.5).*dx
    v = vmin .+ (collect(1:nv_bins) .- 0.5).*dv
    return x, v, distribution
end

function write_vlasov_snapshot!(path, step, time, grid, electrons, ions, cfg,
                                 options; new_file=false)
    x, electron_v, fe = phase_space_histogram(electrons,grid,options.nx_bins,
        options.nv_bins,options.electron_velocity_range)
    _, ion_v, fi = phase_space_histogram(ions,grid,options.nx_bins,
        options.nv_bins,options.ion_velocity_range)
    frame = (kind=:vlasov_frame,version=1,step=step,time=time,
        boundary_condition=cfg.boundary_condition,x=x,electron_v=electron_v,
        ion_v=ion_v,electron_f=fe,ion_f=fi,
        electric_field=unified_interpolate(grid.electric_field,x,grid,cfg),
        potential=unified_interpolate(grid.phi,x,grid,cfg),
        electron_charge=-1.0,electron_mass=1.0,
        ion_charge=1.0,ion_mass=cfg.mass_ratio)
    mkpath(dirname(path))
    open(path,new_file ? "w" : "a") do io
        serialize(io,frame)
    end
    return nothing
end
