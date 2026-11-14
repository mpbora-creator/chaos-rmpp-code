#!/usr/bin/env julia

# This is the only file a beginner needs to run. It loads the editable input
# parameters, loads the simulation module, and starts the time loop.
project_root = normpath(joinpath(@__DIR__, ".."))
cd(project_root)

function usage()
    println("""
Usage:
  julia --project=. scripts/run_pic.jl
  julia --project=. scripts/run_pic.jl --use-saved-parameters RUNTIME_IDENTIFIER
  julia --project=. scripts/run_pic.jl --parameters-file FILE

The second form loads the exact archived configuration from:
  diagnostics/<RUNTIME_IDENTIFIER>/user_parameters.jl

The third form loads an explicitly named parameter file.
""")
end

if ARGS==["--help"]
    usage()
    exit(0)
end
length(ARGS) in (0,2)||begin
    usage()
    error("Expected no arguments, --use-saved-parameters RUNTIME_IDENTIFIER, "*
          "or --parameters-file FILE")
end

saved_identifier=nothing
explicit_parameters_file=nothing
if length(ARGS)==2
    if ARGS[1]=="--use-saved-parameters"
        saved_identifier=ARGS[2]
        !isempty(strip(saved_identifier))&&
            basename(saved_identifier)==saved_identifier&&
            saved_identifier∉(".","..")||
            error("RUNTIME_IDENTIFIER must be a non-empty directory name, not a path")
    elseif ARGS[1]=="--parameters-file"
        isempty(strip(ARGS[2]))&&error("FILE must be a non-empty path")
        explicit_parameters_file=abspath(ARGS[2])
    else
        error("Unknown option: $(ARGS[1])")
    end
end

user_parameters_source=!isnothing(explicit_parameters_file) ?
    explicit_parameters_file :
    (isnothing(saved_identifier) ?
        joinpath(project_root,"config","user_parameters.jl") :
        joinpath(project_root,"diagnostics",saved_identifier,
                 "user_parameters.jl"))
isfile(user_parameters_source)||
    error("User-parameters file not found: $user_parameters_source")
user_parameters_source_text=read(user_parameters_source,String)
include(user_parameters_source)
hasproperty(USER_PARAMETERS, :runtime_identifier) ||
    error("runtime_identifier is required in $user_parameters_source")
if !isnothing(saved_identifier)
    USER_PARAMETERS.runtime_identifier==saved_identifier||
        error("Saved file defines runtime_identifier="*
              "'$(USER_PARAMETERS.runtime_identifier)', expected '$saved_identifier'")
    println("Using archived user parameters: $user_parameters_source")
elseif !isnothing(explicit_parameters_file)
    println("Using explicit user-parameters file: $user_parameters_source")
end
RUN_PARAMETERS=merge(USER_PARAMETERS,(
    user_parameters_source_text=user_parameters_source_text,
    user_parameters_source_path=user_parameters_source,
))

# =============================================================================
# PLUGGABLE DEBRIS PHOTOEMISSION
# =============================================================================
# The physics implementation lives in src/DebrisPhotoemission.jl so the
# launcher and the established PIC paths stay independent. This subsection
# only reports the user-selected plug-in state; `enabled=false` is a strict
# no-op, including for archived parameter files that predate photoemission.
if hasproperty(USER_PARAMETERS, :photoemission) &&
   USER_PARAMETERS.photoemission.enabled
    photo = USER_PARAMETERS.photoemission
    println("Debris photoemission enabled: photon flux=" *
            "$(photo.photon_flux_m2_s) m^-2 s^-1, yield=$(photo.quantum_yield), " *
            "spectrum=$(photo.cold_temperature_ev)/$(photo.hot_temperature_ev) eV")
end
if hasproperty(USER_PARAMETERS, :photoemission_sl_vdpm_feedback) &&
   USER_PARAMETERS.photoemission_sl_vdpm_feedback.enabled
    bridge=USER_PARAMETERS.photoemission_sl_vdpm_feedback
    measurement=hasproperty(bridge,:measurement) ? bridge.measurement : :ion_density_mode
    measurement_detail=measurement===:local_electric_field_bandpass ?
        "band=$(bridge.band_frequency)" : "mode=$(bridge.mode_number)"
    println("Photoemission SL-vdPM feedback enabled: measurement=" *
            "$(measurement), $(measurement_detail), " *
            "tau rate=$(bridge.tau_rate_per_omega_pe), " *
            "mean multiplier=$(bridge.mean_rate_multiplier), " *
            "modulation depth=$(bridge.modulation_depth)")
end
if hasproperty(USER_PARAMETERS,:periodic_photoemission_driver) &&
   USER_PARAMETERS.periodic_photoemission_driver.enabled
    driver=USER_PARAMETERS.periodic_photoemission_driver
    println("Periodic photoemission driver enabled: waveform=$(driver.waveform), " *
            "depth=$(driver.modulation_depth), frequency=$(driver.frequency_omega_pe), " *
            "ramp=$(driver.ramp_start_time)+$(driver.ramp_duration)")
end
if hasproperty(USER_PARAMETERS,:external_debris_charge_driver) &&
   USER_PARAMETERS.external_debris_charge_driver.enabled
    driver=USER_PARAMETERS.external_debris_charge_driver
    println("External debris charge driver enabled: model=$(driver.model), " *
            "amplitude=$(driver.amplitude_e)e, frequency=$(driver.frequency_omega_pe), " *
            "ramp=$(driver.ramp_start_time)+$(driver.ramp_duration)")
end
if hasproperty(USER_PARAMETERS,:surface_state_photoemission_feedback) &&
   USER_PARAMETERS.surface_state_photoemission_feedback.enabled
    surface=USER_PARAMETERS.surface_state_photoemission_feedback
    println("Impact-driven surface-state photoemission feedback enabled: " *
        "tau=$(surface.relaxation_time_omega_pe)/omega_pe, " *
        "capacity=$(surface.capacity_electrons)e, " *
        "deposition=$(surface.deposition_efficiency), " *
        "photo gain=$(surface.photoemission_gain), " *
        "cubic quenching=$(surface.nonlinear_quenching_rate)")
end
if hasproperty(USER_PARAMETERS,:photoelectron_sheath_feedback) &&
   USER_PARAMETERS.photoelectron_sheath_feedback.enabled
    sheath=USER_PARAMETERS.photoelectron_sheath_feedback
    println("Photoelectron-sheath escape feedback enabled: " *
        "activation=$(sheath.activation_time_omega_pe)/omega_pe, " *
        "residence=$(sheath.residence_time_omega_pe)/omega_pe, " *
        "return=$(sheath.return_time_omega_pe)/omega_pe, " *
        "capacity=$(sheath.capacity_electrons)e, trapped fraction=" *
        "$(sheath.trapped_fraction), space-charge strength=" *
        "$(sheath.space_charge_strength)")
end
if hasproperty(USER_PARAMETERS,:charge_dependent_collection) &&
   USER_PARAMETERS.charge_dependent_collection.enabled
    collection=USER_PARAMETERS.charge_dependent_collection
    println("Charge-dependent collection enabled: model=$(collection.model), " *
        "electron coupling=$(collection.electron_coupling_strength), " *
        "ion coupling=$(collection.ion_coupling_strength), multiplier bounds=" *
        "$(collection.minimum_cross_section_multiplier).." *
        "$(collection.maximum_cross_section_multiplier)")
end
if hasproperty(USER_PARAMETERS,:multimode_iicsi_diagnostics) &&
   USER_PARAMETERS.multimode_iicsi_diagnostics.enabled
    modes=USER_PARAMETERS.multimode_iicsi_diagnostics
    println("Read-only multimode IICSI diagnostics enabled: modes=$(modes.modes), " *
        "window=$(modes.window_fraction), every=$(modes.every) steps")
end

# Independent electron-impact SEE plug-in. Archived configurations without
# this block use the disabled default supplied by the simulation module.
if hasproperty(USER_PARAMETERS, :secondary_electron_emission) &&
   USER_PARAMETERS.secondary_electron_emission.enabled
    see = USER_PARAMETERS.secondary_electron_emission
    println("Debris secondary electron emission enabled: model=$(see.model), " *
            "yield scale=$(see.yield_scale), peak yield=$(see.peak_yield) " *
            "at $(see.peak_energy_ev) eV")
    # The simulation module is loaded below, after the optional threaded
    # relaunch. Read the user-facing delay block directly at this stage.
    delay = hasproperty(see, :delay) ? see.delay : (enabled=false,)
    if delay.enabled
        println("  delayed SEE enabled: fraction=$(delay.delayed_fraction), " *
                "relaxation=$(delay.relaxation_time_omega_pe)/omega_pe, " *
                "trap capacity=$(delay.trap_capacity_electrons) electrons, " *
                "release efficiency=$(delay.release_efficiency), " *
                "nonlinear saturation rate=$(delay.nonlinear_saturation_rate)")
    end
end
if hasproperty(USER_PARAMETERS, :dielectric_surface_polarization) &&
   USER_PARAMETERS.dielectric_surface_polarization.enabled
    pol=USER_PARAMETERS.dielectric_surface_polarization
    println("Dielectric surface polarization enabled: source=$(pol.source_model), " *
            "fraction=$(pol.deposition_fraction), tau=$(pol.relaxation_time_omega_pe), " *
            "capacity=$(pol.polarization_capacity_electrons)e, " *
            "cubic saturation=$(pol.nonlinear_saturation_rate), " *
            "SEE modulation=$(pol.modulate_see_yield)")
end
if hasproperty(USER_PARAMETERS,:see_response_probe) &&
   USER_PARAMETERS.see_response_probe.enabled
    probe=USER_PARAMETERS.see_response_probe
    println("SEE response probe enabled: $(probe.waveform), " *
            "amplitude=$(probe.amplitude), frequency=$(probe.start_frequency)" *
            "..$(probe.end_frequency), time=$(probe.start_time)..$(probe.end_time)")
end

# Julia fixes its thread count at process startup. On the configured MacBook,
# transparently relaunch this script once with the requested number of threads.
# Other machines proceed in the original serial process.
thread_policy=USER_PARAMETERS.machine_threading

# A relaunched child must not survive if its one-purpose launcher is closed or
# interrupted by the terminal/application. Check the parent relationship in a
# lightweight asynchronous watchdog while the simulation runs.
if get(ENV,"PIC_PARALLEL_RELAUNCH","0")=="1" &&
   haskey(ENV,"PIC_LAUNCHER_PID")
    expected_parent=parse(Int,ENV["PIC_LAUNCHER_PID"])
    @async begin
        while true
            sleep(0.5)
            if ccall(:getppid,Cint,())!=expected_parent
                println("PIC launcher exited; stopping the threaded child process.")
                exit(130)
            end
        end
    end
end

if thread_policy.automatic && thread_policy.enable_on_matching_machine &&
   get(ENV,"PIC_PARALLEL_RELAUNCH","0")!="1"
    models=unique(info.model for info in Sys.cpu_info())
    model=isempty(models) ? "unknown" : join(models,", ")
    matches=Sys.isapple() && occursin("apple",lowercase(Sys.MACHINE)) &&
        gethostname()==thread_policy.required_hostname &&
        model==thread_policy.required_cpu_model &&
        length(Sys.cpu_info())==thread_policy.required_hardware_cores
    if matches && Threads.nthreads()!=thread_policy.threads
        println("Detected configured MacBook Pro; relaunching PIC with " *
                "$(thread_policy.threads) Julia threads.")
        command=`$(Base.julia_cmd()) --threads=$(thread_policy.threads) --project=$(project_root) $(abspath(@__FILE__)) $ARGS`
        child_command=setenv(command,
            "PIC_PARALLEL_RELAUNCH"=>"1",
            "PIC_LAUNCHER_PID"=>string(getpid()))
        # An asynchronously launched process otherwise receives disconnected
        # standard streams in some terminals. Explicitly attach the threaded
        # child to this terminal so progress and errors remain visible.
        child=run(pipeline(child_command,
                           stdin=stdin,stdout=stdout,stderr=stderr);wait=false)
        atexit(() -> process_running(child) && kill(child,Base.SIGTERM))
        try
            wait(child)
        catch error
            process_running(child) && kill(child,Base.SIGTERM)
            rethrow(error)
        end
        success(child) || error("Threaded PIC child exited with code $(child.exitcode)")
        exit()
    end
end

include(joinpath(project_root, "src", "BareBonesPIC.jl"))

BareBonesPIC.run_simulation(RUN_PARAMETERS)

# Optional, read-only post-run projection onto Eqs. (15)-(17) of the
# Stuart-Landau charging model.  This does not participate in the PIC time
# loop and is an exact no-op when sl_charging_reduction.enabled=false or absent.
if hasproperty(RUN_PARAMETERS,:sl_charging_reduction) &&
   RUN_PARAMETERS.sl_charging_reduction.enabled
    include(joinpath(project_root,"src","SLChargingReduction.jl"))
    SLChargingReduction.fit_sl_charging_reduction(RUN_PARAMETERS)
end
