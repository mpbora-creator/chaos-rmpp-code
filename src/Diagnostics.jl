"Values retained in memory to draw the energy history without CSV files."
mutable struct DiagnosticsState
    time::Vector{Float64}
    energy_time::Vector{Float64}
    electron_kinetic::Vector{Float64}
    ion_kinetic::Vector{Float64}
    field_energy::Vector{Float64}
    mode_signal::Vector{Float64}
    mode_amplitude::Vector{Float64}
end

"Prepare the figure directory and an initially empty diagnostic history."
function prepare_diagnostics(options)
    options.save_plots && mkpath(options.output_directory)
    return DiagnosticsState(Float64[], Float64[], Float64[], Float64[],
                            Float64[], Float64[], Float64[])
end

function kinetic_energy(species::Species)
    return 0.5 * species.weight * species.mass * sum(abs2, species.v)
end

"Record the selected electric-field mode every step for spectral diagnostics."
function record_mode!(state, time, grid, options)
    push!(state.time, time)
    electric_modes = fft(grid.electric_field) ./ grid.nx
    selected_mode = electric_modes[options.spectrum_mode + 1]
    push!(state.mode_signal, imag(selected_mode))
    push!(state.mode_amplitude, abs(selected_mode))
    return nothing
end

"Record energies only at the plotting cadence to avoid needless synchronization."
function record_energies!(state, time, grid, electrons, ions)
    push!(state.energy_time, time)
    push!(state.electron_kinetic, kinetic_energy(electrons))
    push!(state.ion_kinetic, kinetic_energy(ions))
    push!(state.field_energy, 0.5 * grid.dx * sum(abs2, grid.electric_field))

    # For the sinusoidal position perturbation used here, the imaginary part of
    # E_k is the signed oscillating amplitude. Recording it at every simulation
    # step gives an evenly sampled signal suitable for a temporal FFT.
    return nothing
end

"Fit log(mode amplitude) = constant + gamma*time inside a chosen time window."
function fitted_growth_rate(state, window)
    indices = findall(t -> window[1] <= t <= window[2], state.time)
    length(indices) >= 2 || return NaN
    x = state.time[indices]
    y = log.(max.(state.mode_amplitude[indices], eps(Float64)))
    xmean, ymean = sum(x) / length(x), sum(y) / length(y)
    denominator = sum((x .- xmean) .^ 2)
    denominator > 0 || return NaN
    return sum((x .- xmean) .* (y .- ymean)) / denominator
end

"Plot mode growth and straight-line fits for the selected linear/nonlinear windows."
function growth_rate_plot(state, options)
    log_amplitude = log.(max.(state.mode_amplitude, eps(Float64)))
    panel = plot(state.time, log_amplitude; label="log |E_k|",
        xlabel="t omega_pe", ylabel="log mode amplitude",
        title="Two-stream mode growth")
    for (name, window, style) in
        (("linear", options.linear_growth_window, :dash),
         ("nonlinear", options.nonlinear_growth_window, :dot))
        gamma = fitted_growth_rate(state, window)
        if isfinite(gamma)
            indices = findall(t -> window[1] <= t <= window[2], state.time)
            x = state.time[indices]
            y = log_amplitude[indices]
            intercept = sum(y .- gamma .* x) / length(x)
            plot!(panel, x, intercept .+ gamma .* x; linestyle=style, linewidth=2,
                  label="$name gamma=$(round(gamma, digits=3))")
        end
    end
    return panel
end

"Return angular frequency and normalized power for the recorded mode signal."
function power_spectrum(state)
    number = length(state.mode_signal)
    number >= 3 || return Float64[], Float64[]
    dt = state.time[2] - state.time[1]
    signal = state.mode_signal .- sum(state.mode_signal) / number
    # A Hann window reduces spectral leakage caused by the finite run length.
    window = 0.5 .- 0.5 .* cos.(2pi .* (0:number-1) ./ (number - 1))
    transformed = fft(signal .* window)
    last_positive = fld(number, 2) + 1
    omega = (2pi / (number * dt)) .* collect(0:last_positive-1)
    power = abs2.(transformed[1:last_positive])
    maximum(power) > 0 && (power ./= maximum(power))
    return omega, power
end

"Build one diagnostic dashboard from whichever plots the user enabled."
function make_diagnostic_plot(step, time, grid, electrons, ions, options, state)
    panels = Plots.Plot[]
    if options.fields
        push!(panels, plot(grid.x, grid.electric_field; label="E",
             xlabel="x / lambda_De", ylabel="Field", title="Fields: t=$(round(time, digits=3))"))
        plot!(panels[end], grid.x, grid.phi; label="phi")
    end
    if options.densities
        push!(panels, plot(grid.x, grid.ne; label="n_e", xlabel="x / lambda_De",
             ylabel="Density", title="Particle densities"))
        plot!(panels[end], grid.x, grid.ni; label="n_i")
    end
    if options.phase_space
        stride = options.phase_space_stride
        electron_sample = 1:stride:length(electrons.x)
        phase = scatter(electrons.x[electron_sample], electrons.v[electron_sample];
            markersize=1.5, markerstrokewidth=0, alpha=0.55, label="electrons",
            xlabel="x / lambda_De", ylabel="v / v_te0", title="Phase space")
        if options.plot_ion_phase_space
            ion_sample = 1:stride:length(ions.x)
            scatter!(phase, ions.x[ion_sample], ions.v[ion_sample];
                markersize=1.5, markerstrokewidth=0, alpha=0.55, label="ions")
        end
        push!(panels, phase)
    end
    if options.energies
        total = state.electron_kinetic .+ state.ion_kinetic .+ state.field_energy
        energy = plot(state.energy_time, state.electron_kinetic; label="electron kinetic",
            xlabel="t omega_pe", ylabel="Energy", title="Energy history")
        plot!(energy, state.energy_time, state.ion_kinetic; label="ion kinetic")
        plot!(energy, state.energy_time, state.field_energy; label="field")
        plot!(energy, state.energy_time, total; label="total", linewidth=2)
        push!(panels, energy)
    end
    if options.power_spectrum
        omega, power = power_spectrum(state)
        if !isempty(omega)
            keep = omega .<= 2.5 # Focus on the electron-frequency range.
            coherent = maximum(state.mode_amplitude) >=
                       options.coherent_mode_threshold
            displayed_power = coherent ? power[keep] : zeros(sum(keep))
            title = coherent ? "Power spectrum" :
                    "No coherent oscillation (particle noise only)"
            label = coherent ? "mode $(options.spectrum_mode)" : "no resolved signal"
            spectrum = plot(omega[keep], displayed_power; label=label,
                xlabel="angular frequency / omega_pe", ylabel="Normalized coherent power",
                title=title, ylims=(0, 1.08))
            vline!(spectrum, [1.0]; label="omega_pe = 1", linestyle=:dash,
                   linewidth=2)
            push!(panels, spectrum)
        end
    end
    if options.growth_rates
        push!(panels, growth_rate_plot(state, options))
    end
    isempty(panels) && return nothing
    return plot(panels...; layout=(length(panels), 1),
                size=(900, 260 * length(panels)), plot_title="PIC step $step")
end

"Update the live display and/or save the same dashboard as a PNG image."
function update_diagnostics!(step, time, grid, electrons, ions, options, state)
    figure = make_diagnostic_plot(step, time, grid, electrons, ions, options, state)
    isnothing(figure) && return nothing
    # On a headless computer, set live_plots=false and retain save_plots=true.
    options.live_plots && display(figure)
    if options.save_plots
        tag = lpad(step, 6, '0')
        savefig(figure, joinpath(options.output_directory, "pic_$(tag).png"))
    end
    return nothing
end
