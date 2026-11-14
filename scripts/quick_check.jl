#!/usr/bin/env julia

# A tiny smoke test. It uses fewer particles and steps than the normal example.
project_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(project_root, "config", "user_parameters.jl"))
hasproperty(USER_PARAMETERS,:runtime_identifier) ||
    error("runtime_identifier is required in config/user_parameters.jl")
include(joinpath(project_root, "src", "BareBonesPIC.jl"))

quick_diagnostics = merge(USER_PARAMETERS.diagnostics,
    (every=2, live_plots=false, save_plots=false, delete_previous_plots=false,
     debris_charge=false,
     power_spectrum=false,
     output_directory=joinpath(project_root,"fig","quick_check")))
quick_parameters = merge(USER_PARAMETERS,
    (runtime_identifier="quick_check",
     control_run=merge(USER_PARAMETERS.control_run,(enabled=false,)),
     grid=merge(USER_PARAMETERS.grid,
                (nx=32, domain_length=32.0, scale_with_domain=false)),
     time=(dt=0.05, nsteps=4),
     plasma=merge(USER_PARAMETERS.plasma,(n_particles_e=500,n_particles_i=500)),
     boundary=merge(USER_PARAMETERS.boundary,(kind=:periodic,)),
     debris=merge(USER_PARAMETERS.debris,(enabled=false,)),
     photoemission=merge(USER_PARAMETERS.photoemission,(enabled=false,)),
     photoemission_charge_coupling=merge(
         USER_PARAMETERS.photoemission_charge_coupling,(enabled=false,)),
     surface_state_photoemission_feedback=merge(
         USER_PARAMETERS.surface_state_photoemission_feedback,(enabled=false,)),
     photoelectron_sheath_feedback=merge(
         USER_PARAMETERS.photoelectron_sheath_feedback,(enabled=false,)),
     secondary_electron_emission=merge(
         USER_PARAMETERS.secondary_electron_emission,(enabled=false,)),
     dielectric_surface_polarization=merge(
         USER_PARAMETERS.dielectric_surface_polarization,(enabled=false,)),
     reservoir_density_regulator=merge(
         USER_PARAMETERS.reservoir_density_regulator,(enabled=false,)),
     external_debris_charge_driver=merge(
         USER_PARAMETERS.external_debris_charge_driver,(enabled=false,)),
     checkpoint=merge(USER_PARAMETERS.checkpoint,
                      (enabled=false,resume_if_present=false,save_time_history=false)),
     sindy_export=merge(USER_PARAMETERS.sindy_export,(enabled=false,)),
     compute_backend=:cpu, diagnostics=quick_diagnostics))

result = BareBonesPIC.run_simulation(quick_parameters)
@assert all(isfinite, result.grid.electric_field)
@assert abs(sum(result.grid.rho)) < 1e-10

# A fresh run removes its configured histories and PNGs. Turning restart on
# bypasses this cleanup so that a continuation can append existing data.
mktempdir() do directory
    history=joinpath(directory,"history.tsv")
    snapshots=joinpath(directory,"snapshots.jls")
    plot_directory=joinpath(directory,"fig")
    mkpath(plot_directory)
    write(history,"old history")
    write(snapshots,"old snapshots")
    write(joinpath(plot_directory,"old_name.png"),"old png")
    cleanup_parameters=merge(quick_parameters,
        (time=(dt=0.05,nsteps=0),
         diagnostics=merge(quick_parameters.diagnostics,
             (save_plots=true,delete_previous_plots=true,
              output_directory=plot_directory)),
         checkpoint=merge(quick_parameters.checkpoint,
             (save_time_history=true,history_file=history,
              auto_separate_fresh_runs=false,
              delete_previous_history=true)),
         sindy_export=merge(USER_PARAMETERS.sindy_export,
             (enabled=true,vlasov=merge(USER_PARAMETERS.sindy_export.vlasov,
                 (enabled=true,file=snapshots,nx_bins=16,nv_bins=16))))))
    BareBonesPIC.run_simulation(cleanup_parameters)
    @assert !isfile(joinpath(plot_directory,"old_name.png"))
    @assert isfile(history) && !occursin("old history",read(history,String))
    @assert isfile(snapshots) && filesize(snapshots)>length("old snapshots")
end

# Reusing a runtime identifier for a fresh run must replace its existing
# products even when the scientific configuration has changed.
mktempdir() do directory
    history=joinpath(directory,"history.tsv")
    snapshots=joinpath(directory,"snapshots.jls")
    checkpoint=joinpath(directory,"checkpoint.jls")
    figures=joinpath(directory,"fig")
    base_parameters=merge(quick_parameters,
        (time=(dt=0.05,nsteps=1),
         diagnostics=merge(quick_parameters.diagnostics,
             (output_directory=figures,save_plots=false)),
         checkpoint=merge(USER_PARAMETERS.checkpoint,
             (enabled=true,file=checkpoint,resume_if_present=false,
              parameter_record_file=joinpath(directory,"restart.txt"),
              history_file=history,save_time_history=true,
              auto_separate_fresh_runs=true)),
         sindy_export=merge(USER_PARAMETERS.sindy_export,
             (enabled=true,vlasov=merge(USER_PARAMETERS.sindy_export.vlasov,
                 (enabled=true,every=1,file=snapshots,nx_bins=16,nv_bins=16))))))
    BareBonesPIC.run_simulation(base_parameters)
    changed_parameters=merge(base_parameters,
        (debris=merge(base_parameters.debris,(enabled=true,)),))
    BareBonesPIC.run_simulation(changed_parameters)
    @assert isfile(history)&&isfile(snapshots)&&isfile(checkpoint)
    routed_histories=filter(name->startswith(name,"history_")&&endswith(name,".tsv"),
                            readdir(directory))
    routed_snapshots=filter(name->startswith(name,"snapshots_")&&endswith(name,".jls"),
                            readdir(directory))
    routed_checkpoints=filter(name->startswith(name,"checkpoint_")&&endswith(name,".jls"),
                              readdir(directory))
    @assert isempty(routed_histories)
    @assert isempty(routed_snapshots)
    @assert isempty(routed_checkpoints)
end

# Verify the discovery-quality history contains the compact IICSI observables.
mktempdir() do directory
    history=joinpath(directory,"iicsi_history.tsv")
    feature_parameters=merge(quick_parameters,
        (time=(dt=0.05,nsteps=2),
         debris=merge(USER_PARAMETERS.debris,(enabled=true,)),
         checkpoint=merge(quick_parameters.checkpoint,
             (save_time_history=true,history_file=history)),
         sindy_export=merge(USER_PARAMETERS.sindy_export,
             (enabled=true,vlasov=merge(USER_PARAMETERS.sindy_export.vlasov,
                                        (enabled=false,))))))
    feature_result=BareBonesPIC.run_simulation(feature_parameters)
    header=split(readline(history),'\t')
    required=["potential_well_depth","ion_trapped_fraction",
              "ion_counterstreaming_index","electron_inward_flux",
              "ion_inward_flux","electron_expected_capture_flux",
              "ion_expected_capture_flux"]
    @assert all(name->name in header,required)
    values=split(readlines(history)[end],'\t')
    @assert length(values)==length(header)
    @assert all(isfinite,parse.(Float64,values))
end

# A diagnostic-only export override must preserve an exact restart while
# allowing a later sustained-state E(x,t) recording interval.
mktempdir() do directory
    checkpoint=joinpath(directory,"restart.jls")
    history=joinpath(directory,"history.tsv")
    snapshots=joinpath(directory,"sustained_fields.jls")
    first_export=merge(USER_PARAMETERS.sindy_export,
        (enabled=true,override_on_restart=false,
         vlasov=merge(USER_PARAMETERS.sindy_export.vlasov,
                      (enabled=false,file=snapshots,nx_bins=16,nv_bins=16))))
    first_parameters=merge(quick_parameters,
        (time=(dt=0.05,nsteps=2),
         debris=merge(USER_PARAMETERS.debris,(enabled=true,)),
         checkpoint=merge(USER_PARAMETERS.checkpoint,
             (enabled=true,file=checkpoint,resume_if_present=false,
              parameter_record_file=joinpath(directory,"restart.txt"),
              history_file=history,save_time_history=true)),
         sindy_export=first_export))
    first_result=BareBonesPIC.run_simulation(first_parameters)
    @assert isfile(checkpoint)
    @assert !isfile(snapshots)

    continued_export=merge(first_export,
        (override_on_restart=true,
         vlasov=merge(first_export.vlasov,(enabled=true,every=1))))
    continued_parameters=merge(first_parameters,
        (time=(dt=0.05,nsteps=2),
         checkpoint=merge(first_parameters.checkpoint,(resume_if_present=true,)),
         sindy_export=continued_export))
    continued_result=BareBonesPIC.run_simulation(continued_parameters)
    @assert isfile(snapshots)
    history_lines=readlines(history)
    @assert split(history_lines[end],'\t')[1]=="4"
    @assert all(isfinite,continued_result.grid.electric_field)
    @assert all(isfinite,first_result.grid.electric_field)
end
println("Quick check passed.")
