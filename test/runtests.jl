if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
println(" test module MaxwellFrequency")
include("testMaxwellFwd.jl")
include("testMaxwellSolvers.jl")
println(" MaxwellFrequency: All tests passed!")