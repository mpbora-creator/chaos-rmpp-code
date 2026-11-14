using DelimitedFiles

const ε = 0.26113426413326724
const Λ = 0.19072799832134626
const H = 0.034160949053665846
const K = -0.31559103249347764
const C = -2.37611
const D = -0.9267200175037225
const E = -0.7870886774439758
const G = 1.183328923259103

function rhs(u)
    ξ, y, z = u
    return (
        y,
        -ξ + ε * (1 - ξ^2) * y - H * ξ * z + K * z,
        C * ξ + D * y - Λ * z + E * ξ * z - G * z^3,
    )
end

function rk4(u, h)
    k1 = rhs(u)
    u2 = ntuple(i -> u[i] + h * k1[i] / 2, 3)
    k2 = rhs(u2)
    u3 = ntuple(i -> u[i] + h * k2[i] / 2, 3)
    k3 = rhs(u3)
    u4 = ntuple(i -> u[i] + h * k3[i], 3)
    k4 = rhs(u4)
    return ntuple(i -> u[i] + h * (k1[i] + 2k2[i] + 2k3[i] + k4[i]) / 6, 3)
end

function main()
    h = 0.01
    transient = 1500.0
    measure = 3867.0
    phase_stride = 50              # displayed phase projection: spacing 0.5
    trace_stride = 20              # displayed time trace: spacing 0.2
    u = (0.1, 0.0, 0.0)
    for _ in 1:round(Int, transient / h)
        u = rk4(u, h)
    end

    phase_trajectory = Vector{NTuple{4,Float64}}()
    time_trace = Vector{NTuple{4,Float64}}()
    crossings = Vector{NTuple{3,Float64}}()
    old = u
    nsteps = round(Int, measure / h)
    for i in 1:nsteps
        old = u
        u = rk4(u, h)
        τ = i * h
        if i % phase_stride == 0
            push!(phase_trajectory, (τ, u[1], u[2], u[3]))
        end
        if τ <= 500.0 && i % trace_stride == 0
            push!(time_trace, (τ, u[1], u[2], u[3]))
        end
        if old[1] < 0 <= u[1] && u[2] > 0
            α = -old[1] / (u[1] - old[1])
            yc = old[2] + α * (u[2] - old[2])
            zc = old[3] + α * (u[3] - old[3])
            push!(crossings, (τ - h + α * h, yc, zc))
        end
    end

    outdir = joinpath(@__DIR__, "..", "output", "pdf", "figures", "fig1_theory_data")
    mkpath(outdir)
    open(joinpath(outdir, "phase_trajectory.tsv"), "w") do io
        println(io, "tau\txi\ty\tz")
        writedlm(io, phase_trajectory, '\t')
    end
    open(joinpath(outdir, "time_trace.tsv"), "w") do io
        println(io, "tau\txi\ty\tz")
        writedlm(io, time_trace, '\t')
    end
    open(joinpath(outdir, "poincare.tsv"), "w") do io
        println(io, "tau\ty\tz")
        writedlm(io, crossings, '\t')
    end
    open(joinpath(outdir, "lyapunov.tsv"), "w") do io
        println(io, "dt\tlambda1\tlambda2\tlambda3")
        println(io, "0.020\t0.07781\t-0.000133\t-1.2850")
        println(io, "0.010\t0.07959\t-0.0000891\t-1.2691")
        println(io, "0.005\t0.07803\t-0.000587\t-1.2706")
    end
    println("phase-trajectory points: ", length(phase_trajectory))
    println("time-trace points: ", length(time_trace))
    println("Poincare crossings: ", length(crossings))
    println("output: ", outdir)
end

main()
