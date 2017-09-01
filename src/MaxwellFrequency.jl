module MaxwellFrequency

using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers 
using KrylovMethods

import jInv.ForwardShare.ForwardProbType


const mu0 = 4*pi*1e-7
#export use_iw
#const use_iw = true  # false # = true for i*omega, = false for -i*omega

include("MaxwellFreqParam.jl")
include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("solveMaxFreq.jl")
include("getMaxwellFreqMatrix.jl")

hasJOcTree = false
try
    using JOcTree
    hasJOcTree = true
catch
end

include("Utils/getMaxwellFreqParamOT.jl")
include("Utils/getOTMeshFromTxRx.jl")

end
