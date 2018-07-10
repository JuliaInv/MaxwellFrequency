module MaxwellFrequency

using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
using KrylovMethods
using JOcTree

import jInv.ForwardShare.ForwardProbType

const mu0 = 4*pi*1e-7
export mu0

include("MaxwellFreqParam.jl")
include("MaxwellFreqModel.jl")
include("getData.jl")
include("getSensMat.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("solveMaxFreq.jl")
include("getMaxwellFreqMatrix.jl")

end
