"Solve periodic -phi'' = rho spectrally, then calculate E = -phi'."
function solve_poisson_fft!(grid::Grid1D)
    rho_hat = fft(grid.rho)

    # FFT ordering is 0,1,...,N/2,-N/2+1,...,-1 for an even grid.
    integer_modes = if iseven(grid.nx)
        vcat(0:(grid.nx ÷ 2), (-grid.nx ÷ 2 + 1):-1)
    else
        vcat(0:((grid.nx - 1) ÷ 2), (-(grid.nx - 1) ÷ 2):-1)
    end
    wave_numbers = (2pi / grid.length) .* integer_modes

    phi_hat = zeros(ComplexF64, grid.nx)
    phi_hat[2:end] .= rho_hat[2:end] ./ (wave_numbers[2:end] .^ 2)
    phi_hat[1] = 0.0 # Potential is defined only up to an arbitrary constant.

    electric_hat = -im .* wave_numbers .* phi_hat
    grid.phi .= real(ifft(phi_hat))
    grid.electric_field .= real(ifft(electric_hat))
    return nothing
end

