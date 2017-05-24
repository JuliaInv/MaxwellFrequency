using Base.Test

println(" test module MaxwellFrequency")
include("testMaxwellFwd.jl")
include("testMaxwellSolvers.jl")
println(" MaxwellFrequency: All tests passed!")