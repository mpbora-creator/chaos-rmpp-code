#!/usr/bin/env julia

# Small regression for the plug-in open-reservoir/debris mode. It deliberately
# disables figures and checks finite fields plus global charge accounting.
project_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(project_root, "config", "user_parameters.jl"))
hasproperty(USER_PARAMETERS,:runtime_identifier) ||
    error("runtime_identifier is required in config/user_parameters.jl")
include(joinpath(project_root, "src", "BareBonesPIC.jl"))

diagnostics = merge(USER_PARAMETERS.diagnostics,
    (every=5, live_plots=false, save_plots=false, delete_previous_plots=false,
     power_spectrum=false,
     output_directory=joinpath(project_root,"fig","quick_check_open")))
boundary = merge(USER_PARAMETERS.boundary,
    (kind=:open_reservoir,
     sponge=merge(USER_PARAMETERS.boundary.sponge,(buffer_length=4.0,))))
parameters = merge(USER_PARAMETERS,
    (runtime_identifier="quick_check_open",
     control_run=merge(USER_PARAMETERS.control_run,(enabled=false,)),
     grid=merge(USER_PARAMETERS.grid,
                (nx=64,domain_length=32.0,scale_with_domain=false)),
     time=(dt=0.05,nsteps=10),
     plasma=merge(USER_PARAMETERS.plasma,(n_particles_e=1000,n_particles_i=1000)),
     boundary=boundary, flow=(enabled=false,sound_speed=0.0),
     perturbation=merge(USER_PARAMETERS.perturbation,(enabled=false,)),
     checkpoint=merge(USER_PARAMETERS.checkpoint,
                      (enabled=false,resume_if_present=false,save_time_history=false)),
     sindy_export=merge(USER_PARAMETERS.sindy_export,(enabled=false,)),
     diagnostics=diagnostics, compute_backend=:cpu))

result = BareBonesPIC.run_simulation(parameters)
@assert all(isfinite, result.grid.phi)
@assert all(isfinite, result.grid.electric_field)
@assert maximum(abs, result.diagnostic_state.charge_error) < 1e-10
println("Open-reservoir quick check passed.")
