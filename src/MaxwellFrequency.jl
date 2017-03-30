module MaxwellFrequency

using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers 
using KrylovMethods

import jInv.ForwardShare.ForwardProbType
export MaxwellFreqParam, getMaxwellFreqParam

"""
type MaxwellFrequency.MaxwellFreqParam <: ForwardProbType

defines one MaxwellFrequency problem

Fields:

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    Fields::Array{Complex128}
        - solution to the fwd problem
    freq:Float64
        - angular frequency (includes 2*pi term)
    Ainv::AbstractSolver
    fname::AbstractString     -> Not used??
"""
type MaxwellFreqParam{T} <: ForwardProbType
    Mesh::T
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    Fields::Array{Complex128}
    freq::Float64
    Ainv::AbstractSolver
    fname::AbstractString
end

Base.copy(P::MaxwellFreqParam) = MaxwellFreqParam(P.M, P.Sources, P.Obs, P.Fields, P.freq, P.Ainv, P.fname)

"""
function getMaxwellFreqParam
    
constructor for MaxwellFreqParam

Required Inputs

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    Fields::Array{Complex128}
        - solution to the fwd problem
    freq:Float64
        - angular frequency (includes 2*pi term)
    linSolParam::AbstractSolver
    fname::AbstractString
        - filename where pfor is stored. Not used??
"""
function getMaxwellFreqParam(Mesh::AbstractMesh, 
                             Sources, 
                             Obs, 
                             fields, 
                             freq, 
                             linSolParam::AbstractSolver; 
                             fname="")
    
    if isempty(fields)
        fields = Array(Complex128,0,0)
    end

    return MaxwellFreqParam(Mesh, Sources, Obs, fields, freq, linSolParam, fname)
end

export MaxwellFreqParamSE, getMaxwellFreqParamSE

"""
type MaxwellFrequency.MaxwellFreqParam <: ForwardProbType

defines one MaxwellFrequency problem

Fields:

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
        - used for calculating rhs   
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    Fields::Array{Complex128}
        - solution to the fwd problem
    freq:Float64
        - angular frequency (includes 2*pi term)
    Ainv::AbstractSolver
    Sens::Array{Complex128} 
        - sensitivity matrix
    fname::AbstractString
        - filename where pfor is stored. Not used??
"""
type MaxwellFreqParamSE{T} <: ForwardProbType
    Mesh::T
    Sources::Union{Array{Complex128},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    freq::Float64
    Ainv::AbstractSolver
    Sens::Array{Complex128}
    fname::AbstractString
end

"""
function getMaxwellFreqParamSE
    
constructor for MaxwellFreqParamSE

Required Inputs

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    Fields::Array{Complex128}
        - solution to the fwd problem
    freq:Float64
        - angular frequency (includes 2*pi term)
    linSolParam::AbstractSolver
    fname::AbstractString
        - filename where pfor is stored. Not used??
"""
function getMaxwellFreqParamSE(M::AbstractMesh, 
                               Sources,
                               Obs,
                               freq,
                               linSolParam::AbstractSolver;
                               fname="")
    
    Sens = Array(Complex128,0,0)
    
    return MaxwellFreqParamSE(M, Sources, Obs, freq, linSolParam, Sens, fname)
end

include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("solveMaxFreq.jl")

hasJOcTree = false
try
    using JOcTree
    hasJOcTree = true
catch
end

include("Utils/getMaxwellFreqParamOT.jl")
include("Utils/getOTMeshFromTxRx.jl")

end
