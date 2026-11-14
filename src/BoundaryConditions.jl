# Boundary behavior is kept behind this small interface. Particle-pushing code
# therefore does not need to know which boundary condition is selected.
abstract type AbstractBoundaryCondition end

"Particles leaving one side immediately re-enter through the opposite side."
struct PeriodicBC <: AbstractBoundaryCondition end

"Particles may leave and are replenished by independent external reservoirs."
struct OpenReservoirBC <: AbstractBoundaryCondition end

function make_boundary(name::Symbol)
    name === :periodic && return PeriodicBC()
    name === :open_reservoir && return OpenReservoirBC()
    error("Boundary $(name) is not implemented. Add a subtype of " *
          "AbstractBoundaryCondition and its apply_boundary! method.")
end

function apply_boundary!(x::AbstractVector, domain_length::Real, ::PeriodicBC)
    # mod maps every coordinate, including a negative one, into [0, L).
    @. x = mod(x, domain_length)
    return nothing
end
