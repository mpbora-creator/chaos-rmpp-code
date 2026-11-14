"Return effective grid/particle resolution, optionally scaled with box length."
function domain_scaled_resolution(p)
    if !p.grid.scale_with_domain
        return (; nx=p.grid.nx, n_particles_e=p.plasma.n_particles_e,
                n_particles_i=p.plasma.n_particles_i, scale_factor=1.0)
    end
    p.grid.reference_domain_length > 0 ||
        error("grid.reference_domain_length must be positive")
    p.grid.reference_nx >= 3 || error("grid.reference_nx must be at least 3")
    p.plasma.reference_n_particles_e > 0 ||
        error("plasma.reference_n_particles_e must be positive")
    p.plasma.reference_n_particles_i > 0 ||
        error("plasma.reference_n_particles_i must be positive")
    factor = p.grid.domain_length / p.grid.reference_domain_length
    reference_dx = p.boundary.kind === :periodic ?
                   p.grid.reference_domain_length/p.grid.reference_nx :
                   p.grid.reference_domain_length/(p.grid.reference_nx-1)
    scaled_nx = p.boundary.kind === :periodic ?
                round(Int, p.grid.domain_length/reference_dx) :
                round(Int, p.grid.domain_length/reference_dx) + 1
    return (; nx=max(3, scaled_nx),
            n_particles_e=max(1, round(Int, p.plasma.reference_n_particles_e*factor)),
            n_particles_i=max(1, round(Int, p.plasma.reference_n_particles_i*factor)),
            scale_factor=factor)
end

function validate_runtime_identifier(p)
    hasproperty(p,:runtime_identifier) ||
        error("runtime_identifier is required in config/user_parameters.jl")
    identifier=p.runtime_identifier
    identifier isa AbstractString&&!isempty(strip(identifier)) ||
        error("runtime_identifier must be a non-empty string")
    basename(identifier)==identifier&&identifier∉(".","..") ||
        error("runtime_identifier must be a directory name, not a path")
    return nothing
end

"Check the independent user choices before allocating a large simulation."
function validate_user_parameters(p, resolution)
    validate_runtime_identifier(p)
    is_control=endswith(p.runtime_identifier,"_control")
    if hasproperty(p,:control_run)
        p.control_run.end_time>=0||
            error("control_run.end_time must be nonnegative")
        is_control==p.control_run.enabled||
            error("A runtime ending in _control must have control_run.enabled=true")
    end
    is_control&&p.debris.enabled&&
        error("A runtime ending in _control must have debris.enabled=false")
    photo = hasproperty(p, :photoemission) ? p.photoemission :
            default_photoemission_config()
    see = hasproperty(p, :secondary_electron_emission) ?
          p.secondary_electron_emission : default_secondary_emission_config()
    see_delay = secondary_emission_delay_config(see)
    polarization = dielectric_polarization_config(
        hasproperty(p, :dielectric_surface_polarization) ?
        p.dielectric_surface_polarization : nothing)
    see_probe=see_response_probe_config(hasproperty(p,:see_response_probe) ?
        p.see_response_probe : nothing)
    photo_bridge=photoemission_sl_vdpm_config(
        hasproperty(p,:photoemission_sl_vdpm_feedback) ?
        p.photoemission_sl_vdpm_feedback : nothing)
    periodic_photo=periodic_photoemission_driver_config(
        hasproperty(p,:periodic_photoemission_driver) ?
        p.periodic_photoemission_driver : nothing)
    photo_charge_coupling=photoemission_charge_coupling_config(
        hasproperty(p,:photoemission_charge_coupling) ?
        p.photoemission_charge_coupling : nothing)
    surface_photo_feedback=surface_state_photoemission_feedback_config(
        hasproperty(p,:surface_state_photoemission_feedback) ?
        p.surface_state_photoemission_feedback : nothing)
    photoelectron_sheath=photoelectron_sheath_feedback_config(
        hasproperty(p,:photoelectron_sheath_feedback) ?
        p.photoelectron_sheath_feedback : nothing)
    charge_collection=charge_dependent_collection_config(
        hasproperty(p,:charge_dependent_collection) ?
        p.charge_dependent_collection : nothing)
    external_charge=external_debris_charge_driver_config(
        hasproperty(p,:external_debris_charge_driver) ?
        p.external_debris_charge_driver : nothing)
    density_regulator=reservoir_density_regulator_config(
        hasproperty(p,:reservoir_density_regulator) ?
        p.reservoir_density_regulator : nothing)
    multimode_diagnostics=multimode_iicsi_diagnostics_config(
        hasproperty(p,:multimode_iicsi_diagnostics) ?
        p.multimode_iicsi_diagnostics : nothing)
    photo.enabled && !p.debris.enabled &&
        error("photoemission.enabled=true requires debris.enabled=true")
    if multimode_diagnostics.enabled
        !isempty(multimode_diagnostics.modes) ||
            error("multimode_iicsi_diagnostics.modes must not be empty")
        all(m->m isa Integer && 0<m<resolution.nx÷2,
            multimode_diagnostics.modes) ||
            error("multimode IICSI modes must lie between 1 and spatial Nyquist")
        left,right=multimode_diagnostics.window_fraction
        0<=left<right<=1 ||
            error("multimode diagnostic window must lie in [0,1]")
        multimode_diagnostics.every>0 ||
            error("multimode_iicsi_diagnostics.every must be positive")
    end
    photo.enabled && !p.debris.dynamic_charge &&
        error("photoemission.enabled=true requires debris.dynamic_charge=true")
    photo.enabled && p.debris.physical_particles_per_macroparticle <= 0 &&
        error("photoemission requires positive debris.physical_particles_per_macroparticle")
    photo.photon_flux_m2_s >= 0 || error("photoemission.photon_flux_m2_s must be nonnegative")
    0 <= photo.quantum_yield <= 1 || error("photoemission.quantum_yield must lie in [0,1]")
    0 <= photo.illumination_fraction <= 1 || error("photoemission.illumination_fraction must lie in [0,1]")
    photo.photon_energy_ev > 0 || error("photoemission.photon_energy_ev must be positive")
    photo.reference_density_m3 > 0 || error("photoemission.reference_density_m3 must be positive")
    photo.reference_temperature_ev > 0 || error("photoemission.reference_temperature_ev must be positive")
    photo.cold_temperature_ev > 0 && photo.hot_temperature_ev > 0 ||
        error("photoelectron temperatures must be positive")
    0 <= photo.hot_fraction <= 1 || error("photoemission.hot_fraction must lie in [0,1]")
    photo.sunward_side in (:left,:right,:both) ||
        error("photoemission.sunward_side must be :left, :right, or :both")
    (isnothing(photo.effective_emitting_area_m2) || photo.effective_emitting_area_m2 > 0) ||
        error("photoemission.effective_emitting_area_m2 must be nothing or positive")
    photo_bridge.enabled && !photo.enabled &&
        error("photoemission_sl_vdpm_feedback.enabled=true requires photoemission.enabled=true")
    periodic_photo.enabled && !photo.enabled &&
        error("periodic_photoemission_driver.enabled=true requires photoemission.enabled=true")
    photo_charge_coupling.enabled && !photo.enabled &&
        error("photoemission_charge_coupling.enabled=true requires photoemission.enabled=true")
    surface_photo_feedback.enabled && !photo.enabled &&
        error("surface_state_photoemission_feedback.enabled=true requires photoemission.enabled=true")
    surface_photo_feedback.enabled && !see.enabled &&
        error("surface-state photoemission feedback requires secondary_electron_emission.enabled=true")
    photoelectron_sheath.enabled && !photo.enabled &&
        error("photoelectron_sheath_feedback.enabled=true requires photoemission.enabled=true")
    if photoelectron_sheath.enabled
        photoelectron_sheath.model===:space_charge_return ||
            error("photoelectron sheath model must be :space_charge_return")
        photoelectron_sheath.activation_time_omega_pe>=0 ||
            error("photoelectron sheath activation time must be nonnegative")
        photoelectron_sheath.residence_time_omega_pe>0 ||
            error("photoelectron sheath residence time must be positive")
        photoelectron_sheath.return_time_omega_pe>0 ||
            error("photoelectron sheath return time must be positive")
        photoelectron_sheath.capacity_electrons>0 ||
            error("photoelectron sheath capacity must be positive")
        0<=photoelectron_sheath.trapped_fraction<=1 ||
            error("photoelectron sheath trapped_fraction must lie in [0,1]")
        photoelectron_sheath.potential_baseline_time_omega_pe>0 ||
            error("photoelectron sheath potential baseline time must be positive")
        photoelectron_sheath.potential_scale_normalized>0 ||
            error("photoelectron sheath potential scale must be positive")
        photoelectron_sheath.potential_return_coupling>=0 ||
            error("photoelectron sheath potential coupling must be nonnegative")
        photoelectron_sheath.space_charge_strength>=0 ||
            error("photoelectron sheath space-charge strength must be nonnegative")
        0<photoelectron_sheath.minimum_escape_multiplier<=1 ||
            error("photoelectron sheath minimum multiplier must lie in (0,1]")
    end
    charge_collection.enabled && !p.debris.enabled &&
        error("charge_dependent_collection.enabled=true requires debris.enabled=true")
    if charge_collection.enabled
        charge_collection.model===:flow_corrected_oml ||
            error("charge-dependent collection model must be :flow_corrected_oml")
        charge_collection.sampling_offset_debye>=0 ||
            error("collection sampling offset must be nonnegative")
        charge_collection.barrier_smoothing_normalized>0 ||
            error("collection barrier smoothing must be positive")
        charge_collection.electron_coupling_strength>=0 &&
            charge_collection.ion_coupling_strength>=0 ||
            error("collection coupling strengths must be nonnegative")
        0<charge_collection.minimum_cross_section_multiplier<=1 ||
            error("minimum collection multiplier must lie in (0,1]")
        charge_collection.maximum_cross_section_multiplier>=1 ||
            error("maximum collection multiplier must be at least one")
        charge_collection.minimum_cross_section_multiplier<=
            charge_collection.maximum_cross_section_multiplier ||
            error("invalid collection multiplier bounds")
    end
    if surface_photo_feedback.enabled
        surface_photo_feedback.source_model===:expected_incident_energy ||
            error("surface-state source_model must be :expected_incident_energy")
        surface_photo_feedback.relaxation_time_omega_pe>0 ||
            error("surface-state relaxation time must be positive")
        surface_photo_feedback.capacity_electrons>0 ||
            error("surface-state capacity must be positive")
        0<=surface_photo_feedback.deposition_efficiency<=1 ||
            error("surface-state deposition_efficiency must lie in [0,1]")
        surface_photo_feedback.reference_impact_energy_ev>0 ||
            error("surface-state reference impact energy must be positive")
        surface_photo_feedback.nonlinear_quenching_rate>=0 ||
            error("surface-state nonlinear quenching must be nonnegative")
        surface_photo_feedback.photoemission_gain>=0 ||
            error("surface-state photoemission gain must be nonnegative")
        surface_photo_feedback.potential_baseline_time_omega_pe>0 ||
            error("surface-state potential baseline time must be positive")
        surface_photo_feedback.potential_scale_normalized>0 ||
            error("surface-state potential scale must be positive")
        surface_photo_feedback.potential_escape_coupling>=0 ||
            error("surface-state potential escape coupling must be nonnegative")
        0<=surface_photo_feedback.minimum_rate_multiplier<=
            surface_photo_feedback.maximum_rate_multiplier ||
            error("invalid surface-state photo-rate multiplier bounds")
    end
    if photo_charge_coupling.enabled
        photo_charge_coupling.frequency_omega_pe>0 ||
            error("photoemission charge-coupling frequency must be positive")
        0<=photo_charge_coupling.modulation_depth<1 ||
            error("photoemission charge-coupling modulation_depth must lie in [0,1)")
        photo_charge_coupling.surface_response_time_omega_pe>0 ||
            error("photoemission surface response time must be positive")
        photo_charge_coupling.nonlinear_saturation>=0 ||
            error("photoemission nonlinear saturation must be nonnegative")
        photo_charge_coupling.potential_baseline_time_omega_pe>0 ||
            error("photoemission potential baseline time must be positive")
        photo_charge_coupling.potential_scale_normalized>0 ||
            error("photoemission potential scale must be positive")
        photo_charge_coupling.potential_coupling_strength>=0 ||
            error("photoemission potential coupling must be nonnegative")
        0<=photo_charge_coupling.minimum_rate_multiplier<=
            photo_charge_coupling.maximum_rate_multiplier ||
            error("invalid photoemission charge-coupling multiplier bounds")
    end
    if periodic_photo.enabled
        periodic_photo.waveform===:sinusoidal ||
            error("periodic photoemission waveform must be :sinusoidal")
        0<=periodic_photo.modulation_depth<1 ||
            error("periodic photoemission modulation_depth must lie in [0,1)")
        periodic_photo.frequency_omega_pe>0 ||
            error("periodic photoemission frequency_omega_pe must be positive")
        periodic_photo.ramp_start_time>=0 ||
            error("periodic photoemission ramp_start_time must be nonnegative")
        periodic_photo.ramp_duration>=0 ||
            error("periodic photoemission ramp_duration must be nonnegative")
    end
    external_charge.enabled && !p.debris.enabled &&
        error("external_debris_charge_driver.enabled=true requires debris.enabled=true")
    if external_charge.enabled
        external_charge.model in (:prescribed_charge_offset,:sinusoidal_current) ||
            error("external debris charge model must be :prescribed_charge_offset or :sinusoidal_current")
        external_charge.model!==:prescribed_charge_offset ||
            external_charge.amplitude_e>0 ||
            error("external debris charge amplitude_e must be positive")
        external_charge.model!==:sinusoidal_current ||
            external_charge.current_amplitude_e_per_omega_pe>0 ||
            error("external current amplitude must be positive")
        external_charge.model!==:sinusoidal_current || p.debris.dynamic_charge ||
            error("sinusoidal external current requires debris.dynamic_charge=true")
        external_charge.model!==:sinusoidal_current ||
            external_charge.maximum_absolute_driven_charge_e>0 ||
            error("maximum_absolute_driven_charge_e must be positive")
        external_charge.frequency_omega_pe>0 ||
            error("external debris charge frequency_omega_pe must be positive")
        external_charge.ramp_start_time>=0 ||
            error("external debris charge ramp_start_time must be nonnegative")
        external_charge.ramp_duration>=0 ||
            error("external debris charge ramp_duration must be nonnegative")
        p.boundary.kind===:open_reservoir ||
            error("external debris charge driver requires open-reservoir Poisson boundaries")
    end
    if density_regulator.enabled
        p.boundary.kind===:open_reservoir ||
            error("reservoir density regulator requires boundary.kind=:open_reservoir")
        density_regulator.target_density_e>0 &&
            density_regulator.target_density_i>0 ||
            error("reservoir density-regulator target densities must be positive")
        density_regulator.relaxation_time_omega_pe>0 ||
            error("reservoir density-regulator relaxation time must be positive")
        0<density_regulator.maximum_fractional_correction_per_step<=1 ||
            error("density-regulator maximum fractional correction must lie in (0,1]")
        density_regulator.injection_side in (:upstream,:left,:right) ||
            error("density-regulator injection_side must be :upstream, :left, or :right")
        density_regulator.inject_charge_neutral_pairs===true ||
            error("the regulator currently requires inject_charge_neutral_pairs=true")
        isapprox(p.plasma.density_e,p.plasma.density_i;rtol=0,atol=0) ||
            error("charge-neutral pair regulation requires equal reference electron and ion densities")
    end
    photo_bridge.enabled && !p.debris.enabled &&
        error("photoemission SL-vdPM feedback requires debris.enabled=true")
    if photo_bridge.enabled
        photo_bridge.measurement in (:local_electric_field_bandpass,:ion_density_mode) ||
            error("unsupported photoemission feedback measurement")
        if photo_bridge.measurement === :ion_density_mode
            0 < photo_bridge.mode_number < resolution.nx÷2 ||
                error("photoemission feedback mode_number must lie below grid Nyquist")
        end
        0 <= photo_bridge.window_fraction[1] < photo_bridge.window_fraction[2] <= 1 ||
            error("photoemission feedback window_fraction must lie inside [0,1]")
        photo_bridge.electric_field_reference > 0 ||
            error("photoemission feedback electric_field_reference must be positive")
        0 < photo_bridge.band_frequency[1] < photo_bridge.band_frequency[2] ||
            error("photoemission feedback band_frequency must be positive and increasing")
        photo_bridge.baseline_relaxation_time_omega_pe > 0 ||
            error("photoemission feedback baseline relaxation must be positive")
        photo_bridge.mode_reference_amplitude > 0 ||
            error("photoemission feedback mode_reference_amplitude must be positive")
        photo_bridge.quadrature_clip > 0 ||
            error("photoemission feedback quadrature_clip must be positive")
        photo_bridge.tau_rate_per_omega_pe > 0 ||
            error("photoemission feedback tau_rate_per_omega_pe must be positive")
        photo_bridge.output_scale > 0 ||
            error("photoemission feedback output_scale must be positive")
        0 <= photo_bridge.modulation_depth <= 1 ||
            error("photoemission feedback modulation_depth must lie in [0,1]")
        0 <= photo_bridge.minimum_rate_multiplier <= photo_bridge.maximum_rate_multiplier ||
            error("invalid photoemission feedback multiplier bounds")
        photo_bridge.mean_rate_multiplier >= 0 ||
            error("photoemission feedback mean_rate_multiplier must be nonnegative")
    end
    see.enabled && !p.debris.enabled &&
        error("secondary_electron_emission.enabled=true requires debris.enabled=true")
    see.enabled && !p.debris.dynamic_charge &&
        error("secondary_electron_emission.enabled=true requires debris.dynamic_charge=true")
    see.model === :electron_impact ||
        error("secondary_electron_emission.model must be :electron_impact")
    see.yield_curve === :sternglass ||
        error("secondary_electron_emission.yield_curve must be :sternglass")
    see.yield_scale >= 0 || error("secondary_electron_emission.yield_scale must be nonnegative")
    see.peak_yield >= 0 || error("secondary_electron_emission.peak_yield must be nonnegative")
    see.threshold_energy_ev >= 0 || error("SEE threshold_energy_ev must be nonnegative")
    see.peak_energy_ev > see.threshold_energy_ev ||
        error("SEE peak_energy_ev must exceed threshold_energy_ev")
    see.maximum_yield >= 0 || error("SEE maximum_yield must be nonnegative")
    see.emitted_temperature_ev > 0 || error("SEE emitted_temperature_ev must be positive")
    see.reference_temperature_ev > 0 || error("SEE reference_temperature_ev must be positive")
    see_delay.model in (:discrete_trap, :continuous_surface) ||
        error("SEE delay.model must be :discrete_trap or :continuous_surface")
    0 <= see_delay.delayed_fraction <= 1 ||
        error("SEE delay.delayed_fraction must lie in [0,1]")
    see_delay.relaxation_time_omega_pe > 0 ||
        error("SEE delay.relaxation_time_omega_pe must be positive")
    see_delay.trap_capacity_electrons > 0 ||
        error("SEE delay.trap_capacity_electrons must be positive")
    0 <= see_delay.release_efficiency <= 1 ||
        error("SEE delay.release_efficiency must lie in [0,1]")
    see_delay.nonlinear_saturation_rate >= 0 ||
        error("SEE delay.nonlinear_saturation_rate must be nonnegative")
    surface = see_delay.continuous_surface
    surface.source_model in (:captured_events, :expected_impact_flux) ||
        error("SEE continuous_surface.source_model must be :captured_events or :expected_impact_flux")
    surface.reference_impact_energy_ev > 0 ||
        error("SEE continuous_surface.reference_impact_energy_ev must be positive")
    surface.energy_weighted_deposition isa Bool ||
        error("SEE continuous_surface.energy_weighted_deposition must be Bool")
    surface.conservative_emission_remainder === true ||
        error("SEE continuous surface requires conservative_emission_remainder=true")
    see.enabled && p.debris.physical_particles_per_macroparticle <= 0 &&
        error("secondary electron emission requires positive physical_particles_per_macroparticle")
    polarization.enabled && !p.debris.enabled &&
        error("dielectric_surface_polarization.enabled=true requires debris.enabled=true")
    polarization.source_model === :expected_impact_asymmetry ||
        error("dielectric polarization source_model must be :expected_impact_asymmetry")
    polarization.relaxation_time_omega_pe > 0 ||
        error("dielectric polarization relaxation_time_omega_pe must be positive")
    polarization.polarization_capacity_electrons > 0 ||
        error("dielectric polarization capacity must be positive")
    0 <= polarization.deposition_fraction <= 1 ||
        error("dielectric polarization deposition_fraction must lie in [0,1]")
    polarization.nonlinear_saturation_rate >= 0 ||
        error("dielectric polarization nonlinear_saturation_rate must be nonnegative")
    polarization.modulate_see_yield isa Bool ||
        error("dielectric polarization modulate_see_yield must be Bool")
    polarization.see_yield_coupling >= 0 ||
        error("dielectric polarization see_yield_coupling must be nonnegative")
    polarization.minimum_yield_multiplier >= 0 ||
        error("dielectric polarization minimum_yield_multiplier must be nonnegative")
    polarization.maximum_yield_multiplier >=
        polarization.minimum_yield_multiplier ||
        error("maximum_yield_multiplier must not be below minimum_yield_multiplier")
    polarization.modulate_see_yield && !see.enabled &&
        error("modulate_see_yield=true requires secondary electron emission")
    see_probe.enabled && !see.enabled &&
        error("see_response_probe.enabled=true requires secondary electron emission")
    see_probe.waveform === :log_chirp ||
        error("see_response_probe.waveform must be :log_chirp")
    0 <= see_probe.amplitude < 1 ||
        error("see_response_probe.amplitude must lie in [0,1)")
    see_probe.start_frequency > 0 ||
        error("see_response_probe.start_frequency must be positive")
    see_probe.end_frequency > see_probe.start_frequency ||
        error("see_response_probe.end_frequency must exceed start_frequency")
    see_probe.start_time >= 0 ||
        error("see_response_probe.start_time must be nonnegative")
    see_probe.end_time > see_probe.start_time ||
        error("see_response_probe.end_time must exceed start_time")
    p.boundary.kind in (:periodic, :open_reservoir) ||
        error("boundary.kind must be :periodic or :open_reservoir")
    p.compute_backend in (:cpu, :metal) || error("compute_backend must be :cpu or :metal")
    resolution.nx >= 3 || error("effective grid nx must be at least 3")
    p.grid.domain_length > 0 || error("grid.domain_length must be positive")
    p.time.dt > 0 || error("time.dt must be positive")
    p.time.nsteps >= 0 || error("time.nsteps cannot be negative")
    p.diagnostics.every > 0 || error("diagnostics.every must be positive")
    p.diagnostics.phase_space_stride > 0 || error("phase_space_stride must be positive")
    p.diagnostics.phase_space_species in (:electrons, :ions, :both) ||
        error("phase_space_species must be :electrons, :ions, or :both")
    p.flow.sound_speed >= 0 || error("flow.sound_speed must be nonnegative")
    0 < p.debris.center_fraction < 1 ||
        error("debris.center_fraction must lie strictly between 0 and 1")
    p.plasma.density_e > 0 && p.plasma.density_i > 0 ||
        error("Both plasma densities must be positive")
    if p.perturbation.enabled
        0 < p.perturbation.electron_mode < resolution.nx÷2 ||
            error("electron perturbation mode must be below the grid Nyquist mode")
        0 < p.perturbation.ion_mode < resolution.nx÷2 ||
            error("ion perturbation mode must be below the grid Nyquist mode")
    end
    p.diagnostics.power_spectrum && !p.perturbation.enabled &&
        error("power_spectrum=true requires perturbation.enabled=true; no artificial peak is generated for an unperturbed run")
    # With debris disabled this diagnostic is a useful control: it records the
    # identically zero auxiliary charge rather than making a matched no-debris
    # run change an otherwise identical diagnostic configuration.
    p.compute_backend === :metal &&
        error("The unified boundary/flow/debris engine currently uses threaded CPU arrays. Set compute_backend=:cpu.")
    return nothing
end

"Number of steps needed to reach an absolute control end time."
function requested_additional_steps(p,saved=nothing)
    is_control=endswith(p.runtime_identifier,"_control")&&
               hasproperty(p,:control_run)
    is_control||return p.time.nsteps
    target=ceil(Int,p.control_run.end_time/p.time.dt)
    completed=isnothing(saved) ? 0 : saved.completed_step
    return max(0,target-completed)
end

"Force every `_control` invocation to be a clean run from t=0."
function enforce_fresh_control_run(p)
    endswith(p.runtime_identifier,"_control")||return p
    hasproperty(p,:control_run)||
        error("A _control runtime requires the control_run configuration block")
    println("Control runtime detected: ignoring checkpoint.resume_if_present, "*
            "disabling debris, and starting at t=0.")
    return merge(p,(
        debris=merge(p.debris,(enabled=false,)),
        photoemission=merge(hasproperty(p,:photoemission) ? p.photoemission :
                            default_photoemission_config(),(enabled=false,)),
        secondary_electron_emission=merge(
            hasproperty(p,:secondary_electron_emission) ? p.secondary_electron_emission :
            default_secondary_emission_config(),(enabled=false,)),
        dielectric_surface_polarization=merge(
            hasproperty(p,:dielectric_surface_polarization) ?
            p.dielectric_surface_polarization :
            default_dielectric_polarization_config(),(enabled=false,)),
        see_response_probe=merge(hasproperty(p,:see_response_probe) ?
            p.see_response_probe : default_see_response_probe_config(),
            (enabled=false,)),
        photoemission_sl_vdpm_feedback=merge(
            hasproperty(p,:photoemission_sl_vdpm_feedback) ?
            p.photoemission_sl_vdpm_feedback :
            default_photoemission_sl_vdpm_config(), (enabled=false,)),
        periodic_photoemission_driver=merge(
            hasproperty(p,:periodic_photoemission_driver) ?
            p.periodic_photoemission_driver :
            default_periodic_photoemission_driver_config(),(enabled=false,)),
        photoemission_charge_coupling=merge(
            hasproperty(p,:photoemission_charge_coupling) ?
            p.photoemission_charge_coupling :
            default_photoemission_charge_coupling_config(),(enabled=false,)),
        surface_state_photoemission_feedback=merge(
            hasproperty(p,:surface_state_photoemission_feedback) ?
            p.surface_state_photoemission_feedback :
            default_surface_state_photoemission_feedback_config(),
            (enabled=false,)),
        photoelectron_sheath_feedback=merge(
            hasproperty(p,:photoelectron_sheath_feedback) ?
            p.photoelectron_sheath_feedback :
            default_photoelectron_sheath_feedback_config(),(enabled=false,)),
        charge_dependent_collection=merge(
            hasproperty(p,:charge_dependent_collection) ?
            p.charge_dependent_collection :
            default_charge_dependent_collection_config(),(enabled=false,)),
        external_debris_charge_driver=merge(
            hasproperty(p,:external_debris_charge_driver) ?
            p.external_debris_charge_driver :
            default_external_debris_charge_driver_config(),(enabled=false,)),
        reservoir_density_regulator=merge(
            hasproperty(p,:reservoir_density_regulator) ?
            p.reservoir_density_regulator :
            default_reservoir_density_regulator_config(),(enabled=false,)),
        multimode_iicsi_diagnostics=merge(
            hasproperty(p,:multimode_iicsi_diagnostics) ?
            p.multimode_iicsi_diagnostics :
            default_multimode_iicsi_diagnostics_config(),(enabled=false,)),
        checkpoint=merge(p.checkpoint,(
            resume_if_present=false,
            auto_separate_fresh_runs=false,
            delete_previous_history=true,
        )),
    ))
end

"Preserve the exact editable input file used to start a fresh runtime."
function preserve_user_parameters!(p)
    hasproperty(p.checkpoint,:user_parameters_copy)||return nothing
    destination=abspath(p.checkpoint.user_parameters_copy)
    mkpath(dirname(destination))
    if hasproperty(p,:user_parameters_source_text)
        open(destination,"w") do io
            write(io,p.user_parameters_source_text)
        end
    else
        source=normpath(joinpath(@__DIR__,"..","config","user_parameters.jl"))
        cp(source,destination;force=true)
    end
    println("Saved user parameters: $destination")
    return nothing
end

function empty_directory!(directory,label)
    isdir(directory)||return 0
    removed=0
    for name in readdir(directory)
        rm(joinpath(directory,name);recursive=true)
        removed+=1
    end
    println("Fresh run requested: removed $removed existing $label item(s) from $directory/")
    return removed
end

"Clear one runtime identifier before a fresh run, independent of config identity."
function clean_existing_runtime_outputs!(p)
    cp=p.checkpoint
    cp.resume_if_present&&return 0
    validate_runtime_identifier(p)
    identifier=p.runtime_identifier
    export_options=hasproperty(p,:sindy_export) ?
        p.sindy_export : default_sindy_export_options()
    removed=0

    # The identifier owns the entire diagnostics/<runtime_identifier>
    # directory, including checkpoint, histories, SINDy data, and old analysis
    # products. Any additional custom active paths are cleaned file-by-file.
    active_paths=String[]
    if cp.enabled
        push!(active_paths,cp.file,cp.parameter_record_file)
    end
    if cp.save_time_history
        push!(active_paths,cp.history_file,run_identity_record_path(cp.history_file))
    end
    export_options.enabled&&export_options.vlasov.enabled&&
        push!(active_paths,export_options.vlasov.file)

    expected_diagnostic=joinpath("diagnostics",identifier)
    removed+=empty_directory!(expected_diagnostic,"diagnostic")
    for path in unique(active_paths)
        # The standard directory has already been cleared wholesale. Custom
        # active output paths are removed individually.
        abspath(dirname(path))==abspath(expected_diagnostic)&&continue
        if isfile(path)
            rm(path)
            removed+=1
            println("Fresh run requested: removed existing file $path")
        end
    end

    # Figures belong to the configured run output directory. Remove every
    # file/subdirectory, not only PNGs or names known to the current config.
    removed+=empty_directory!(p.diagnostics.output_directory,"figure")
    return removed
end

"Translate the clean public configuration into the reservoir/debris engine inputs."
function unified_engine_config(p, resolution)
    flow_speed = p.flow.enabled ? p.flow.sound_speed : 0.0
    frequency = merge(p.diagnostics.frequency, (enabled=p.diagnostics.power_spectrum,))
    diagnostics = (
        every=p.diagnostics.every,
        fields=p.diagnostics.fields,
        densities=p.diagnostics.densities,
        phase_space=p.diagnostics.phase_space,
        phase_space_species=p.diagnostics.phase_space_species,
        phase_space_stride=p.diagnostics.phase_space_stride,
        energies=p.diagnostics.energies,
        debris_charge=p.diagnostics.debris_charge,
        charge_accounting=p.diagnostics.charge_accounting,
        live_plots=p.diagnostics.live_plots,
        save_plots=p.diagnostics.save_plots,
        output_directory=p.diagnostics.output_directory,
    )
    return (
        boundary_condition=p.boundary.kind,
        nx=resolution.nx, domain_length=p.grid.domain_length,
        dt=p.time.dt, nsteps=p.time.nsteps,
        n_particles_e=resolution.n_particles_e,
        n_particles_i=resolution.n_particles_i,
        density_e=p.plasma.density_e, density_i=p.plasma.density_i,
        ambient_density=p.plasma.density_e,
        electron_temperature=p.plasma.temperature_e,
        ion_temperature=p.plasma.temperature_i,
        mass_ratio=p.plasma.mass_ratio,
        flow_sound_speed=flow_speed,
        two_sided_thermal_injection=p.boundary.two_sided_thermal_injection,
        phi_left=p.boundary.phi_left, phi_right=p.boundary.phi_right,
        sponge=p.boundary.sponge,
        perturbation=p.perturbation,
        frequency_diagnostics=frequency,
        equilibrium_tail_fraction=p.diagnostics.equilibrium_tail_fraction,
        debris=p.debris,
        photoemission=hasproperty(p,:photoemission) ? p.photoemission :
                      default_photoemission_config(),
        secondary_electron_emission=hasproperty(p,:secondary_electron_emission) ?
                      p.secondary_electron_emission : default_secondary_emission_config(),
        dielectric_surface_polarization=
            hasproperty(p,:dielectric_surface_polarization) ?
            p.dielectric_surface_polarization :
            default_dielectric_polarization_config(),
        see_response_probe=hasproperty(p,:see_response_probe) ?
            p.see_response_probe : default_see_response_probe_config(),
        photoemission_sl_vdpm_feedback=
            hasproperty(p,:photoemission_sl_vdpm_feedback) ?
            p.photoemission_sl_vdpm_feedback :
            default_photoemission_sl_vdpm_config(),
        periodic_photoemission_driver=
            hasproperty(p,:periodic_photoemission_driver) ?
            p.periodic_photoemission_driver :
            default_periodic_photoemission_driver_config(),
        photoemission_charge_coupling=
            hasproperty(p,:photoemission_charge_coupling) ?
            p.photoemission_charge_coupling :
            default_photoemission_charge_coupling_config(),
        surface_state_photoemission_feedback=
            hasproperty(p,:surface_state_photoemission_feedback) ?
            p.surface_state_photoemission_feedback :
            default_surface_state_photoemission_feedback_config(),
        photoelectron_sheath_feedback=
            hasproperty(p,:photoelectron_sheath_feedback) ?
            p.photoelectron_sheath_feedback :
            default_photoelectron_sheath_feedback_config(),
        charge_dependent_collection=
            hasproperty(p,:charge_dependent_collection) ?
            p.charge_dependent_collection :
            default_charge_dependent_collection_config(),
        external_debris_charge_driver=
            hasproperty(p,:external_debris_charge_driver) ?
            p.external_debris_charge_driver :
            default_external_debris_charge_driver_config(),
        reservoir_density_regulator=
            hasproperty(p,:reservoir_density_regulator) ?
            p.reservoir_density_regulator :
            default_reservoir_density_regulator_config(),
        multimode_iicsi_diagnostics=
            hasproperty(p,:multimode_iicsi_diagnostics) ?
            p.multimode_iicsi_diagnostics :
            default_multimode_iicsi_diagnostics_config(),
        diagnostics=diagnostics,
    )
end

"Configuration identity used to keep physically different fresh runs separate."
function fresh_run_identity(p,cfg)
    # File destinations, live display choices, and machine threading do not
    # alter the simulated physics. Everything physical, numerical, diagnostic,
    # and export-related remains in the identity.
    diagnostics=merge(cfg.diagnostics,
        (live_plots=false,save_plots=false,output_directory=""))
    clean_cfg=merge(cfg,(diagnostics=diagnostics,))
    # SEE is a pluggable debris-current model and must not rename an existing
    # runtime selected by runtime_identifier. Its complete configuration is
    # still archived in user_parameters.jl and in the checkpoint cfg.
    if hasproperty(clean_cfg,:secondary_electron_emission)
        clean_cfg=Base.structdiff(clean_cfg,
            (secondary_electron_emission=nothing,))
    end
    if hasproperty(clean_cfg,:dielectric_surface_polarization)
        clean_cfg=Base.structdiff(clean_cfg,
            (dielectric_surface_polarization=nothing,))
    end
    if hasproperty(clean_cfg,:see_response_probe)
        clean_cfg=Base.structdiff(clean_cfg,(see_response_probe=nothing,))
    end
    if hasproperty(clean_cfg,:photoemission_sl_vdpm_feedback)
        clean_cfg=Base.structdiff(clean_cfg,
            (photoemission_sl_vdpm_feedback=nothing,))
    end
    if hasproperty(clean_cfg,:periodic_photoemission_driver)
        clean_cfg=Base.structdiff(clean_cfg,
            (periodic_photoemission_driver=nothing,))
    end
    if hasproperty(clean_cfg,:photoemission_charge_coupling)
        clean_cfg=Base.structdiff(clean_cfg,
            (photoemission_charge_coupling=nothing,))
    end
    if hasproperty(clean_cfg,:surface_state_photoemission_feedback)
        clean_cfg=Base.structdiff(clean_cfg,
            (surface_state_photoemission_feedback=nothing,))
    end
    if hasproperty(clean_cfg,:photoelectron_sheath_feedback)
        clean_cfg=Base.structdiff(clean_cfg,
            (photoelectron_sheath_feedback=nothing,))
    end
    if hasproperty(clean_cfg,:charge_dependent_collection)
        clean_cfg=Base.structdiff(clean_cfg,
            (charge_dependent_collection=nothing,))
    end
    if hasproperty(clean_cfg,:external_debris_charge_driver)
        clean_cfg=Base.structdiff(clean_cfg,
            (external_debris_charge_driver=nothing,))
    end
    if hasproperty(clean_cfg,:reservoir_density_regulator)
        clean_cfg=Base.structdiff(clean_cfg,
            (reservoir_density_regulator=nothing,))
    end
    if hasproperty(clean_cfg,:multimode_iicsi_diagnostics)
        clean_cfg=Base.structdiff(clean_cfg,
            (multimode_iicsi_diagnostics=nothing,))
    end
    export_options=hasproperty(p,:sindy_export) ?
        p.sindy_export : default_sindy_export_options()
    clean_vlasov=merge(export_options.vlasov,(file="",))
    clean_export=merge(export_options,
        (override_on_restart=false,vlasov=clean_vlasov,))
    return (
        random_seed=p.random_seed,
        compute_backend=p.compute_backend,
        grid=p.grid,
        time=p.time,
        plasma=p.plasma,
        boundary=p.boundary,
        flow=p.flow,
        perturbation=p.perturbation,
        debris=p.debris,
        photoemission=hasproperty(p,:photoemission) ? p.photoemission :
                      default_photoemission_config(),
        cfg=clean_cfg,
        sindy_export=clean_export,
        save_time_history=p.checkpoint.save_time_history,
    )
end

"Small deterministic FNV-1a fingerprint; unlike Base.hash it is stable across sessions."
function stable_run_fingerprint(identity)
    hash=UInt64(0xcbf29ce484222325)
    for byte in codeunits(repr(identity))
        hash=(hash⊻UInt64(byte))*UInt64(0x100000001b3)
    end
    string(hash;base=16,pad=16)
end

function suffixed_file(path,suffix)
    stem,extension=splitext(path)
    stem*"_"*suffix*extension
end

run_identity_record_path(history_file)=splitext(history_file)[1]*"_run_identity.txt"

function recorded_fingerprint(path)
    isfile(path)||return nothing
    open(path,"r") do io
        eof(io)&&return nothing
        fields=split(readline(io),'\t';limit=2)
        length(fields)==2&&fields[1]=="fingerprint" ? fields[2] : nothing
    end
end

"Compare a legacy checkpoint when no run-identity record existed yet."
function legacy_checkpoint_matches(checkpoint_file,cfg)
    isfile(checkpoint_file)||return false
    try
        saved=open(deserialize,checkpoint_file)
        hasproperty(saved,:cfg)||return false
        sanitize(value)=merge(value,(diagnostics=merge(value.diagnostics,
            (live_plots=false,save_plots=false,output_directory="")),))
        return repr(sanitize(saved.cfg))==repr(sanitize(cfg))
    catch
        # An unreadable/foreign checkpoint must never authorize overwriting
        # existing histories.
        return false
    end
end

function write_fresh_run_record(path,fingerprint,tag,identity,p)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"fingerprint\t",fingerprint)
        println(io,"tag\t",tag)
        println(io,"created\t",Dates.now())
        println(io,"\nScientific/numerical run identity:")
        show(io,MIME("text/plain"),identity)
        println(io,"\n\nEffective output paths:")
        println(io,"history\t",p.checkpoint.history_file)
        println(io,"checkpoint\t",p.checkpoint.file)
        println(io,"restart_record\t",p.checkpoint.parameter_record_file)
        if hasproperty(p,:sindy_export)&&p.sindy_export.enabled&&
           p.sindy_export.vlasov.enabled
            println(io,"vlasov\t",p.sindy_export.vlasov.file)
        end
        println(io,"figures\t",p.diagnostics.output_directory)
    end
end

"Route a changed fresh configuration away from an existing run's products."
function route_fresh_run_outputs(p,cfg,resolution)
    cp=p.checkpoint
    automatic=!hasproperty(cp,:auto_separate_fresh_runs)||
              cp.auto_separate_fresh_runs
    identity=fresh_run_identity(p,cfg)
    fingerprint=stable_run_fingerprint(identity)
    base_record=run_identity_record_path(cp.history_file)
    export_options=hasproperty(p,:sindy_export) ?
        p.sindy_export : default_sindy_export_options()
    possible_paths=String[cp.history_file,cp.file]
    export_options.enabled&&export_options.vlasov.enabled&&
        push!(possible_paths,export_options.vlasov.file)
    collision=any(isfile,possible_paths)||
        (isdir(p.diagnostics.output_directory)&&
         any(name->endswith(lowercase(name),".png"),
             readdir(p.diagnostics.output_directory)))

    if !automatic||!collision
        write_fresh_run_record(base_record,fingerprint,"base",identity,p)
        return p,cfg
    end

    previous=recorded_fingerprint(base_record)
    base_physics_files=any(isfile,possible_paths)
    # A legacy figure directory by itself is not evidence of a conflicting
    # physics run. The configured delete_previous_plots policy owns it.
    same_configuration=previous==fingerprint||
        (isnothing(previous)&&legacy_checkpoint_matches(cp.file,cfg))||
        (isnothing(previous)&&!base_physics_files)
    if same_configuration
        write_fresh_run_record(base_record,fingerprint,"base",identity,p)
        return p,cfg
    end

    debris_label=cfg.debris.enabled ? "debris" : "no_debris"
    flow_label=replace(replace(string(cfg.flow_sound_speed),"."=>"p"),"-"=>"m")
    tag="$(cfg.boundary_condition)_$(debris_label)_flow$(flow_label)_$(fingerprint[1:8])"
    routed_checkpoint=merge(cp,(
        file=suffixed_file(cp.file,tag),
        parameter_record_file=suffixed_file(cp.parameter_record_file,tag),
        history_file=suffixed_file(cp.history_file,tag),
    ))
    routed_diagnostics=merge(p.diagnostics,(
        output_directory=p.diagnostics.output_directory*"_"*tag,))
    routed_export=export_options
    if export_options.enabled&&export_options.vlasov.enabled
        routed_export=merge(export_options,(
            vlasov=merge(export_options.vlasov,(
                file=suffixed_file(export_options.vlasov.file,tag),)),))
    end
    routed=merge(p,(checkpoint=routed_checkpoint,
                    diagnostics=routed_diagnostics,
                    sindy_export=routed_export))
    routed_cfg=unified_engine_config(routed,resolution)
    routed_record=run_identity_record_path(routed_checkpoint.history_file)
    write_fresh_run_record(routed_record,fingerprint,tag,identity,routed)
    println("Fresh-run configuration differs from the existing dataset.")
    println("Preserving all previous histories, checkpoint, and figures.")
    println("New run tag: $tag")
    println("  history:   $(routed_checkpoint.history_file)")
    println("  snapshots: $(routed_export.vlasov.file)")
    println("  checkpoint: $(routed_checkpoint.file)")
    println("  figures:   $(routed_diagnostics.output_directory)/")
    return routed,routed_cfg
end

"Delete PNG figures from this run's dedicated output directory."
function clean_previous_plots!(options)
    # A run that is not saving figures must never alter an output directory.
    (options.save_plots && options.delete_previous_plots) || return 0
    isdir(options.output_directory) || return 0
    removed = 0
    for name in readdir(options.output_directory)
        if endswith(lowercase(name), ".png")
            rm(joinpath(options.output_directory, name))
            removed += 1
        end
    end
    println("Removed $removed previous diagnostic plot(s) from $(options.output_directory)/")
    return removed
end

"Remove histories belonging to a fresh run, but never during checkpoint restart."
function clean_previous_histories!(p, export_options)
    cp = p.checkpoint
    delete_history = !hasproperty(cp, :delete_previous_history) ||
                     cp.delete_previous_history
    delete_history || return 0

    paths = String[]
    cp.save_time_history && push!(paths, cp.history_file)
    if export_options.enabled && export_options.vlasov.enabled
        push!(paths, export_options.vlasov.file)
    end

    removed = 0
    for path in unique(paths)
        if isfile(path)
            rm(path)
            println("Removed previous history file $path")
            removed += 1
        end
    end
    return removed
end

function load_checkpoint_if_requested(p)
    cp = p.checkpoint
    cp.resume_if_present && isfile(cp.file) || return nothing
    open(cp.file, "r") do io
        saved = deserialize(io)
        saved.version == 1 || error("Unsupported checkpoint version $(saved.version)")
        return saved
    end
end

function write_restart_parameter_record(saved, p)
    path = p.checkpoint.parameter_record_file
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "PIC checkpoint restart record")
        println(io, "Restarted at: ", Dates.now())
        println(io, "Checkpoint file: ", p.checkpoint.file)
        println(io, "Saved completed step: ", saved.completed_step)
        println(io, "Saved normalized time: ", saved.completed_time)
        println(io, "Additional steps requested now: ", p.time.nsteps)
        println(io, "\nSaved effective configuration (this controls the resumed physics):")
        show(io, MIME("text/plain"), saved.cfg)
        println(io, "\n\nCurrent user parameters (physical values ignored during restart):")
        show(io, MIME("text/plain"), p)
        println(io)
    end
    println("Wrote restart parameter record to $path")
end

function validate_sindy_restart_files(p)
    options = hasproperty(p,:sindy_export) ? p.sindy_export : default_sindy_export_options()
    if options.enabled && p.checkpoint.save_time_history && isfile(p.checkpoint.history_file)
        header = split(readline(p.checkpoint.history_file),'\t')
        required = ["ne_at_debris","dne_dx_at_debris","d2ne_dx2_at_debris",
                    "ion_trapped_fraction","ion_counterstreaming_index",
                    "electron_expected_capture_flux","ion_expected_capture_flux"]
        all(column -> column in header, required) || error("The existing checkpoint/history predates the current SINDy spatial-derivative export. " *
            "A fresh PIC run is required for the IICSI phase-space observables. " *
            "Set checkpoint.resume_if_present=false once to start a fresh discovery-quality dataset; " *
            "the old checkpoint and results are otherwise left untouched.")
    end
    return nothing
end

"Run a thermal electron-ion PIC simulation from independent user choices."
function run_simulation(parameters)
    parameters=enforce_fresh_control_run(parameters)
    # Older user-parameter and checkpoint-era files remain runnable even if
    # they predate the machine_threading input block.
    threading_policy = hasproperty(parameters, :machine_threading) ?
                       parameters.machine_threading : default_machine_threading_policy()
    configure_machine_threading!(threading_policy)
    resolution = domain_scaled_resolution(parameters)
    validate_user_parameters(parameters, resolution)
    !parameters.checkpoint.resume_if_present&&
        clean_existing_runtime_outputs!(parameters)
    saved = load_checkpoint_if_requested(parameters)
    if !isnothing(saved)
        validate_sindy_restart_files(parameters)
        write_restart_parameter_record(saved, parameters)
        additional_steps=requested_additional_steps(parameters,saved)
        println("Restarting exactly from checkpoint step $(saved.completed_step), " *
                "t=$(saved.completed_time); running $additional_steps additional steps.")
        return run_boundary_plasma(parameters, saved.cfg;
                                   restart=saved,
                                   additional_steps=additional_steps)
    end
    cfg = unified_engine_config(parameters, resolution)
    parameters,cfg=route_fresh_run_outputs(parameters,cfg,resolution)
    clean_previous_plots!(parameters.diagnostics)
    export_options = hasproperty(parameters,:sindy_export) ?
                     parameters.sindy_export : default_sindy_export_options()
    clean_previous_histories!(parameters, export_options)
    preserve_user_parameters!(parameters)
    if export_options.enabled && export_options.vlasov.enabled
        export_options.vlasov.nx_bins >= 3 || error("sindy_export.vlasov.nx_bins must be at least 3")
        export_options.vlasov.nv_bins >= 3 || error("sindy_export.vlasov.nv_bins must be at least 3")
    end
    if parameters.grid.scale_with_domain
        println("Domain scaling: L=$(parameters.grid.domain_length), nx=$(cfg.nx), " *
                "Ne=$(cfg.n_particles_e), Ni=$(cfg.n_particles_i), " *
                "dx=$(cfg.domain_length / (cfg.boundary_condition === :periodic ? cfg.nx : cfg.nx-1))")
    end
    return run_boundary_plasma(parameters, cfg;
                               additional_steps=requested_additional_steps(parameters))
end
