module MaxwellFrequency

using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers 
using KrylovMethods
hasJOcTree = false
try
  using JOcTree
  hasJOcTree = true
catch
end

import jInv.ForwardShare.getData
import jInv.ForwardShare.getSensTMatVec
import jInv.ForwardShare.getSensMatVec

# ========== Setting up the parameter structure for the forward problem ========
export MaxwellFreqParam
import jInv.ForwardShare.ForwardProbType
import jInv.ForwardShare.prepareMesh2Mesh
export prepareMesh2Mesh
export ForwardProbType
type MaxwellFreqParam{T} <: ForwardProbType
	Mesh::T     # mesh                       
	Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}  # used for calculating rhs   
	Obs::Union{Array{Complex128},SparseMatrixCSC}  # transpose interpolation matrix from fields to receivers
	Fields::Array{Complex128}    # solution to the fwd problem
	freq::Float64  # includes 2*pi term
	Ainv::AbstractSolver   # stores MUMPSfactorization
	fname::AbstractString  # filename where pfor is stored
end  # type MaxwellParam


Base.copy(P::MaxwellFreqParam) = MaxwellFreqParam(P.M, P.Sources, P.Obs, P.Fields, P.freq, P.Ainv, P.fname)

#Dummy edge constraints function for non-Octree meshes
getEdgeConstraints(M::AbstractMesh) = 1.0

export MaxwellFreqParamSE
type MaxwellFreqParamSE{T} <: ForwardProbType
	Mesh::T     # mesh                       
	Sources::Union{Array{Complex128},SparseMatrixCSC}  # used for calculating rhs   
	Obs::Union{Array{Complex128},SparseMatrixCSC}  # transpose interpolation matrix from fields to receivers
	freq::Float64  # includes 2*pi term
	Ainv::AbstractSolver   # stores MUMPSfactorization
	Sens::Array{Complex128} # sensitivity matrix
	fname::AbstractString  # filename where pfor is stored
end  # type MaxwellParam



include("getMaxwellFreqParam.jl")
include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("solveMaxFreq.jl")
if hasJOcTree
  include("getOTMeshFromTxRx.jl")
end
end  # module MaxwellFrequency
