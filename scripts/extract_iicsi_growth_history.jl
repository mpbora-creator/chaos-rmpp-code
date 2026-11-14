using Serialization
using LinearAlgebra
using Printf
using Statistics
using FFTW

const root = normpath(joinpath(@__DIR__, ".."))
const source = get(ENV, "CHAOS_VLASOV_FILE",
    joinpath(root, "data", "raw", "vlasov_snapshots.jls"))
const output = get(ENV, "CHAOS_GROWTH_OUTPUT",
    joinpath(root, "data", "derived", "iicsi_growth_history_M05.tsv"))
isfile(source) || error("Vlasov archive not found: $source. Set CHAOS_VLASOV_FILE to the raw .jls archive.")

times = Float64[]
profiles = Vector{Float64}[]
open(source, "r") do io
    while !eof(io)
        frame = deserialize(io)
        push!(times, Float64(frame.time))
        push!(profiles, vec(sum(frame.ion_f; dims=2)))
    end
end

order = sortperm(times)
times = times[order]
profiles = profiles[order]
# A restart records its initial state at the join time; retain one copy.
keep = vcat(true, diff(times) .> 0)
times = times[keep]
profiles = profiles[keep]
density = reduce(hcat, profiles)'
density .-= mean(density; dims=1)
factor = svd(density; full=false)
score = factor.U[:, 1] .* factor.S[1]
sum(score .* collect(1.0:length(score))) < 0 && (score .*= -1)

dt = median(diff(times))
# Remove the slow secular component, mirror-pad to suppress FFT edge wrap,
# and retain only the documented IICSI temporal band.
time_centered = times .- mean(times)
score .-= mean(score) .+ time_centered .* (
    dot(time_centered, score .- mean(score)) / dot(time_centered, time_centered))
n = length(score)
padded = vcat(reverse(score), score, reverse(score))
np = length(padded)
spectrum = fft(padded)
analytic_spectrum = zeros(ComplexF64, np)
for k in 0:np-1
    frequency = k <= np÷2 ? k/(np*dt) : (k-np)/(np*dt)
    if 0.015 <= frequency <= 0.05
        analytic_spectrum[k+1] = 2spectrum[k+1]
    end
end
analytic = ifft(analytic_spectrum)
amplitude = abs.(analytic[n+1:2n])

# Smooth only the carrier-scale ripple in the analytic envelope.
smooth_halfwidth = max(2, round(Int, 12.5 / dt))
smoothed = similar(amplitude)
for i in eachindex(amplitude)
    lo, hi = max(1, i-smooth_halfwidth), min(length(amplitude), i+smooth_halfwidth)
    smoothed[i] = mean(@view amplitude[lo:hi])
end
amplitude = smoothed ./ maximum(smoothed)

# Local least-squares slope of log(A) over a 100/omega_pe window.
slope_halfwidth = max(3, round(Int, 50.0 / dt))
growth = similar(amplitude)
for i in eachindex(amplitude)
    lo, hi = max(1, i-slope_halfwidth), min(length(amplitude), i+slope_halfwidth)
    tw = @view times[lo:hi]
    yw = log.(max.(@view(amplitude[lo:hi]), 1.0e-12))
    local_time = tw .- mean(tw)
    growth[i] = dot(local_time, yw .- mean(yw)) / max(dot(local_time, local_time), eps())
end

mkpath(dirname(output))
open(output, "w") do io
    println(io, "# leading ion-density mode from full Vlasov archive")
    println(io, "# envelope: analytic amplitude in 0.015--0.05 omega_pe band, smoothed over 25/omega_pe")
    println(io, "# growth: 100/omega_pe local log-slope of the band-limited envelope")
    println(io, "# time_omega_pe\tnormalized_envelope\teffective_growth_rate_over_omega_pe")
    for i in eachindex(times)
        @printf(io, "%.12g\t%.12g\t%.12g\n", times[i], amplitude[i], growth[i])
    end
end

@printf("wrote %d samples, dt=%.6g, leading-mode variance fraction=%.6g\n",
    length(times), dt, factor.S[1]^2 / sum(abs2, factor.S))
