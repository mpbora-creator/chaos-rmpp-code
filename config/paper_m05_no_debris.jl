# Reconstructed on 2026-07-30 from the saved run identity, restart record,
# checkpoint, and history for diagnostics/debris_flow_3.0.
#
# The run used 50,000-step increments. Its history reached step 100,000
# (t*omega_pe=5000) through one checkpoint continuation.

runtime_identifier="equilibrium_flow_M05_no_debris"

control_run=(enabled=false,end_time=1000.0)
base_runtime_identifier=endswith(runtime_identifier,"_control") ?
    chop(runtime_identifier;tail=length("_control")) : runtime_identifier
is_control_run=control_run.enabled||endswith(runtime_identifier,"_control")
runtime_identifier=is_control_run ?
    base_runtime_identifier*"_control" : base_runtime_identifier
simulation_dt=0.05
normal_run_nsteps=20_000 # matched debris-free control: t_end*omega_pe = 1000
control_nsteps=ceil(Int,control_run.end_time/simulation_dt)
runtime_identifier isa AbstractString&&!isempty(strip(runtime_identifier))||
    error("runtime_identifier must be a non-empty string")
basename(runtime_identifier)==runtime_identifier&&
    runtime_identifier∉(".","..")||
    error("runtime_identifier must be a directory name, not a path")
diagnostic_directory=joinpath(@__DIR__,"..","diagnostics",runtime_identifier)

const USER_PARAMETERS=(
    runtime_identifier=runtime_identifier,
    base_runtime_identifier=base_runtime_identifier,
    control_runtime_identifier=base_runtime_identifier*"_control",
    control_run=(enabled=is_control_run,end_time=control_run.end_time),
    random_seed=42,
    compute_backend=:cpu,
    machine_threading=(
        automatic=false,
        enable_on_matching_machine=true,
        #required_hostname="Scarlett-Z.local",
        #required_cpu_model="Apple M5",
        required_hostname="",
        required_cpu_model="",
        #required_hardware_cores=10,
        required_hardware_cores=0,
        threads=10,
        manual_enabled=true,
    ),
    grid=(
        domain_length=350.0,
        nx=1024,
        scale_with_domain=true,
        reference_domain_length=136.0,
        reference_nx=1024,
    ),
    time=(
        dt=simulation_dt,
        nsteps=is_control_run ? control_nsteps : normal_run_nsteps,
    ),
    checkpoint=(
        enabled=true,
        file=joinpath(diagnostic_directory,"checkpoint.jls"),
        resume_if_present=false,
        auto_separate_fresh_runs=true,
        parameter_record_file=
            joinpath(diagnostic_directory,"restarted_state_parameters.txt"),
        user_parameters_copy=
            joinpath(diagnostic_directory,"user_parameters.jl"),
        history_file=joinpath(diagnostic_directory,"time_history.tsv"),
        save_time_history=true,
        delete_previous_history=true,
    ),
    # Disabled for this fresh forcing-parameter scan.  This switch is used
    # only for deliberately perturbed checkpoint twins.
    restart_debris_charge_perturbation=(
        enabled=false,
        delta_charge_e=0.0,
        apply_at_completed_step=20_000,
    ),
    sindy_export=(
        enabled=true,
        override_on_restart=false,
        local_probe_half_width=1.0,
        instability_probe_half_width=4.0,
        spatial_derivative_smoothing_passes=2,
        vlasov=(
            enabled=true, # Required for the coherent ion-density mode fit.
            every=25,
            file=joinpath(diagnostic_directory,"vlasov_snapshots.jls"),
            nx_bins=64,
            nv_bins=64,
            electron_velocity_range=(-5.0,5.0),
            ion_velocity_range=(-0.08,0.08),
        ),
    ),
    # Optional post-run reduction onto the autonomous theoretical model
    # xi'=y, y'=eps(1-xi^2)y-(1+H*z)xi+K*z,
    # z'=C*xi+D*y-Lambda*z.  This block reads the completed history only; it
    # never modifies the PIC evolution.  Set enabled=false for a strict no-op.
    sl_charging_reduction=(
        enabled=false,
        start_time_omega_pe=500.0,
        end_time_omega_pe=1800.0,
        observables=(:ion_density_mode,:ion_mean_v_local,:phi_at_debris,
                     :potential_well_depth,:ne_at_debris,:ni_at_debris),
        frequency_band=(0.015,0.05),
        carrier_halfwidth_fraction=0.35,
        charge_lowpass_frequency=0.10,
        # Fixed Q_* for comparisons between runs; do not renormalize it to
        # each trajectory's standard deviation.
        charge_scale_e=1.0e8,
        training_fraction=0.70,
        edge_fraction=0.05,
        # Integrate the equations over this physical-time window before
        # regression. This suppresses PIC derivative noise without changing
        # the theoretical SL->vdPM equation library.
        integral_window_omega_pe=2.0,
        ridge=1.0e-6,
    ),
    plasma=(
        n_particles_e=167_000,
        n_particles_i=167_000,
        reference_n_particles_e=167_000,
        reference_n_particles_i=167_000,
        density_e=1.0,
        density_i=1.0,
        temperature_e=1.0,
        temperature_i=0.01,
        mass_ratio=1836.0,
    ),
    boundary=(
        kind=:open_reservoir,
        phi_left=0.0,
        phi_right=0.0,
        two_sided_thermal_injection=true,
        sponge=(
            enabled=true,
            buffer_length=20.0,
            alpha_max=0.20,
            exponent=2.0,
        ),
    ),
    # Optional slow feedback that restores only the common electron/ion
    # deficit. It leaves the ordinary open reservoir and sponge unchanged and
    # is an exact no-op when enabled=false.
    reservoir_density_regulator=(
        enabled=true,
        target_density_e=1.0,
        target_density_i=1.0,
        relaxation_time_omega_pe=100.0,
        maximum_fractional_correction_per_step=1.0e-3,
        injection_side=:upstream,
        inject_charge_neutral_pairs=true,
    ),
    flow=(enabled=true,sound_speed=0.5),
    
    perturbation=(
        enabled=false,
        electron_amplitude=0.02,
        electron_mode=1,
        ion_amplitude=0.02,
        ion_mode=40,
    ),
    debris=(
        enabled=false,
        center_fraction=0.3,
        size_debye=1.0,
        marker_count=64,
        charge_distribution=:gaussian,
        gaussian_sigma_fraction=1/3,
        # Start near the late self-consistent negative floating charge measured
        # in the established M=1.25, size_debye=1 trajectory.  This suppresses
        # the large Qd=0 charging overshoot; the plasma/sheath still relaxes
        # self-consistently from the initialized particle distribution.
        initial_charge_e=-1.35e9,
        dynamic_charge=true,
        physical_particles_per_macroparticle=1.0e6,
        collection_model=:mcc_energy,
        electron_capture=(
            maximum_probability=0.020,
            threshold_energy=0.0,
            energy_scale=0.20,
        ),
        ion_capture=(
            maximum_probability=0.010,
            threshold_energy=0.0,
            energy_scale=0.05,
        ),
    ),
    # Existing physical photoelectron source. The reference temperature is the
    # same 100 eV dimensionalization used by SEE, so size_debye=1 maps to the
    # corresponding physical projected area. Disable this and the bridge below
    # to recover the original collecting-debris IICSI.
    photoemission=(
        # Base physical emitter. Set this false, together with the two switches
        # documented below, to recover the original collecting-debris IICSI.
        enabled=false,
        photon_flux_m2_s=2.5e15,
        photon_energy_ev=20.0,
        quantum_yield=0.10,
        illumination_fraction=0.50,
        sunward_side=:left,
        reference_density_m3=1.0e11,
        reference_temperature_ev=100.0,
        effective_emitting_area_m2=nothing,
        cold_temperature_ev=2.7,
        hot_temperature_ev=10.0,
        hot_fraction=0.05,
    ),
    # Closed-loop interface to the unchanged theoretical nonlinear charging
    # equation z'=Cx+Dy-Lambda*z+E*x*z-G*z^3. It measures causal quadratures of
    # the local debris electric field in the IICSI band and only modulates the
    # physical photoelectron rate.
    photoemission_sl_vdpm_feedback=(
        enabled=false,
        measurement=:local_electric_field_bandpass,
        electric_field_reference=0.01,
        band_frequency=(0.015,0.05),
        baseline_relaxation_time_omega_pe=12.0,
        mode_number=40,
        window_fraction=(0.10,0.90),
        mode_reference_amplitude=0.01,
        quadrature_clip=2.5,
        tau_rate_per_omega_pe=0.20,
        coefficients=(
            Lambda=0.19072799832134626,
            C=-2.376110373607035,
            D=-0.9267200175037225,
            E=-0.7870886774439758,
            G=1.183328923259103,
        ),
        initial_z=0.0,
        mean_rate_multiplier=1.0,
        modulation_depth=0.80,
        output_scale=1.0,
        minimum_rate_multiplier=0.20,
        maximum_rate_multiplier=1.80,
    ),
    # Independent external optical driver. This multiplies the existing
    # physical photoemission rate after a smooth turn-on; it is exactly unity
    # and consumes no random numbers when disabled.
    periodic_photoemission_driver=(
        enabled=false,
        waveform=:sinusoidal,
        modulation_depth=0.60,
        frequency_omega_pe=0.019,
        phase=0.0,
        ramp_start_time=100.0,
        ramp_duration=100.0,
    ),
    # NEW PLUG-IN: delayed, saturating optical response coupled through the
    # actual debris potential. It only multiplies the physical emitter above.
    # Set enabled=false to remove this coupling without changing base emission.
    photoemission_charge_coupling=(
        enabled=false,
        frequency_omega_pe=0.022,
        # Restore the best earlier near-quadrature photoemission coupling.
        phase_lag=pi/2,
        modulation_depth=0.80,
        ramp_start_time=50.0,
        ramp_duration=100.0,
        surface_response_time_omega_pe=20.0,
        nonlinear_saturation=0.3,
        potential_baseline_time_omega_pe=100.0,
        potential_scale_normalized=0.50,
        potential_coupling_strength=0.35,
        minimum_rate_multiplier=0.20,
        maximum_rate_multiplier=1.80,
    ),
    # PLUG-IN: impact-driven finite surface-state reservoir. Incident electron
    # energy populates the state; it relaxes and saturates, then modulates only
    # the existing physical photoemission rate through the resolved debris
    # potential. Set enabled=false for an exact no-op and the preceding model.
    surface_state_photoemission_feedback=(
        enabled=false,
        source_model=:expected_incident_energy,
        # The measured branch has omega0 ~= 0.203, so tau_s=5 gives
        # omega0*tau_s ~= 1: appreciable causal phase lag without a long run.
        relaxation_time_omega_pe=5.0,
        # About twice the equilibrium debris charge. With the measured
        # captured-electron rate this keeps the optical correction moderate
        # while retaining the finite surface-state dynamics.
        capacity_electrons=2.0e9,
        deposition_efficiency=0.50,
        reference_impact_energy_ev=100.0,
        nonlinear_quenching_rate=5.0,
        # Preserve the gain-8 feedback branch. The restored 2.5e15 photon
        # flux and moderate quenching avoid the overdriven, strongly clipped
        # behavior of the stationary_sat20 trial.
        photoemission_gain=8.0,
        potential_baseline_time_omega_pe=50.0,
        potential_scale_normalized=0.50,
        potential_escape_coupling=0.30,
        minimum_rate_multiplier=0.50,
        maximum_rate_multiplier=1.80,
    ),
    # PLUG-IN: unresolved sub-cell photoelectron sheath and return-current
    # response. The ordinary photoemitter still launches real PIC electrons;
    # those particles retain their self-consistent fields, boundary escape and
    # debris recollection. This state only supplies finite surface memory and
    # smooth space-charge-limited escape. It has no hard upper multiplier.
    #
    # Set enabled=false to recover the preceding physical photoemission model
    # exactly. Set photoemission.enabled=false as well (and disable SEE) to
    # recover the original collecting-debris IICSI.
    photoelectron_sheath_feedback=(
        enabled=false,
        model=:space_charge_return,
        # Let the imposed initial charge and open reservoir establish their
        # ordinary floating-potential baseline before closing this feedback.
        activation_time_omega_pe=300.0,
        # Comparable to the observed IICSI phase time, but short enough for a
        # 50,000-step screening run to contain many response times.
        residence_time_omega_pe=8.0,
        return_time_omega_pe=4.0,
        # The configured physical photon flux supplies about 1.22e6 e per
        # omega_pe^-1. With the residence/return times below the expected
        # unresolved cloud is about 1.6e6 e, i.e. occupation ~0.16: strong
        # enough to measure but still perturbative relative to |Qd|~1.35e9 e.
        capacity_electrons=1.0e7,
        trapped_fraction=0.50,
        potential_baseline_time_omega_pe=100.0,
        potential_scale_normalized=0.50,
        potential_return_coupling=0.75,
        # Smooth (1+s*N/Ncap)^(-3/2) space-charge suppression.
        space_charge_strength=1.0,
        minimum_escape_multiplier=0.10,
    ),
    # Pluggable external RF current. Unlike :prescribed_charge_offset, this
    # enters the dynamic debris charge and therefore feeds back on collection.
    # Set enabled=false to recover the original self-charging IICSI exactly.
    external_debris_charge_driver=(
        enabled=false,
        model=:sinusoidal_current,
        amplitude_e=0.0, # Used only by model=:prescribed_charge_offset.
        current_amplitude_e_per_omega_pe=6.9e6,
        maximum_absolute_driven_charge_e=2.5e8,
        frequency_omega_pe=0.022,
        phase=0.0,
        ramp_start_time=50.0,
        ramp_duration=100.0,
    ),
    # Optional electron-impact secondary electron emission. This is independent
    # of photoemission and is an exact no-op when enabled=false.
    secondary_electron_emission=(
        enabled=false,
        model=:electron_impact,
        yield_curve=:sternglass,
        # Keep SEE perturbative so electron/ion collection still sets a
        # negative debris potential and sustains the debris-driven IICSI.
        # The 1.0 trial produced delta_eff ~= 0.45 and charge reversal.
        # Weak prompt-SEE screening value.  This is half of the earlier
        # perturbative 0.20 case and far below the charge-reversing 1.0 case.
        yield_scale=0.10,
        # Example graphite-like curve; replace with measured material values.
        peak_yield=1.0,
        peak_energy_ev=300.0,
        threshold_energy_ev=10.0,
        maximum_yield=1.5,
        emitted_temperature_ev=2.0,
        # Physical eV represented by one normalized electron-temperature unit.
        reference_temperature_ev=100.0,
        # Optional trap-mediated delayed channel. The Sternglass yield above
        # remains the total available yield: delayed_fraction reallocates part
        # of it from prompt emission into a finite surface reservoir.
        delay=(
            # Prompt SEE only: no trap population, delay, or extra memory
            # variable is active in this reduction experiment.
            enabled=false,
            # Pluggable choices: :discrete_trap preserves the earlier delayed
            # SEE implementation; :continuous_surface uses a physical-electron
            # surface excitation with conservative sub-macroparticle release.
            model=:continuous_surface,
            delayed_fraction=0.90,
            relaxation_time_omega_pe=20.0,
            trap_capacity_electrons=1.0e7,
            release_efficiency=1.0,
            # Cubic occupancy-dependent non-emissive quenching rate.
            nonlinear_saturation_rate=0.30,
            continuous_surface=(
                # :captured_events restores the original stochastic-impact
                # source; :expected_impact_flux uses the mean captured fraction
                # of every crossing macroparticle to drive surface excitation.
                source_model=:expected_impact_flux,
                energy_weighted_deposition=true,
                reference_impact_energy_ev=300.0,
                conservative_emission_remainder=true,
            ),
        ),
    ),
    # Optional dielectric surface-memory feedback. P is deposited as -P on
    # the left face and +P on the right face, so it changes the local sheath
    # while adding exactly zero net charge. Set enabled=false for the previous
    # SEE+IICSI model; disabling SEE as well recovers pure collecting IICSI.
    dielectric_surface_polarization=(
        enabled=false,
        source_model=:expected_impact_asymmetry,
        relaxation_time_omega_pe=20.0,
        polarization_capacity_electrons=2.0e7,
        deposition_fraction=0.05,
        nonlinear_saturation_rate=0.10,
        # Face-local dielectric trapped charge modifies the existing
        # Sternglass yield but cannot exceed these conservative bounds.
        modulate_see_yield=false,
        see_yield_coupling=0.15,
        minimum_yield_multiplier=0.80,
        maximum_yield_multiplier=1.20,
    ),
    # Small deterministic open-loop identification input. This multiplies the
    # existing SEE yield, consumes no random numbers, and is exactly unity when
    # disabled. It is a diagnostic probe, not a proposed plasma process.
    see_response_probe=(
        enabled=false,
        waveform=:log_chirp,
        amplitude=0.02,
        start_frequency=0.002,
        end_frequency=0.20,
        start_time=100.0,
        end_time=1000.0,
    ),
    diagnostics=(
        every=5000,
        fields=true,
        densities=true,
        phase_space=true,
        phase_space_species=:ions,
        phase_space_stride=20,
        energies=true,
        debris_charge=true,
        charge_accounting=true,
        power_spectrum=false,
        live_plots=false,
        save_plots=true,
        output_directory=joinpath(@__DIR__,"..","fig",runtime_identifier),
        delete_previous_plots=true,
        equilibrium_tail_fraction=0.50,
        frequency=(
            central_fraction=0.50,
            electron_time_window=(0.0,100.0),
            ion_time_window=(0.0,2500.0),
            electron_frequency_range=(0.5,1.5),
            ion_frequency_range=(0.015,0.05),
            coherent_mode_threshold=1.0e-5,
        ),
    ),
)
