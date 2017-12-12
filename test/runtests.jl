using Base.Test

println(" test module MaxwellFrequency")
@testset "Maxwell Frequency" begin
include("testMaxwellFwd.jl")
include("testMaxwellSolvers.jl")
include("testStorage.jl")
include("testMaxwellJOcTree.jl")
include("testMultiParameterGlobalToLocal.jl")
println(" MaxwellFrequency: All tests passed!")
end
