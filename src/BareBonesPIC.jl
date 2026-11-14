module BareBonesPIC

using FFTW
using Metal
using Atomix: @atomic
using Plots
using Random
using Serialization
using Dates

include("HardwarePolicy.jl")
include("Types.jl")
include("BoundaryConditions.jl")
include("Particles.jl")
include("GridOperations.jl")
include("PoissonFFT.jl")
include("MetalBackend.jl")
include("Diagnostics.jl")
include("SINDyExport.jl")
include("DebrisPhotoemission.jl")
include("PhotoemissionSLvdPMFeedback.jl")
include("PeriodicPhotoemissionDriver.jl")
include("PhotoemissionChargeCoupling.jl")
include("SurfaceStatePhotoemissionFeedback.jl")
include("PhotoelectronSheathFeedback.jl")
include("ChargeDependentCollection.jl")
include("ExternalDebrisChargeDriver.jl")
include("DebrisSEEResponseProbe.jl")
include("DebrisSecondaryEmission.jl")
include("DielectricSurfacePolarization.jl")
include("ReservoirDensityRegulator.jl")
include("MultimodeIICSIDiagnostics.jl")
include("OpenReservoirDebris.jl")
include("Simulation.jl")

export run_simulation

end # module
