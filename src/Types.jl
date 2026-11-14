"A collection of numerical particles belonging to one physical species."
mutable struct Species{TX<:AbstractVector,TV<:AbstractVector}
    name::Symbol
    charge::Float64
    mass::Float64
    density::Float64
    weight::Float64
    x::TX
    v::TV
end

"Uniform one-dimensional mesh. Arrays contain values at x_j = (j-1) dx."
mutable struct Grid1D
    nx::Int
    length::Float64
    dx::Float64
    x::Vector{Float64}
    ne::Vector{Float64}
    ni::Vector{Float64}
    rho::Vector{Float64}
    phi::Vector{Float64}
    electric_field::Vector{Float64}
end

function Grid1D(nx::Integer, length::Real)
    nx > 1 || error("nx must be greater than one")
    length > 0 || error("domain_length must be positive")
    dx = Float64(length) / nx
    x = collect(0:dx:(Float64(length) - dx))
    zeros_on_grid = () -> zeros(Float64, nx)
    return Grid1D(nx, Float64(length), dx, x, zeros_on_grid(),
                  zeros_on_grid(), zeros_on_grid(), zeros_on_grid(),
                  zeros_on_grid())
end
