
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
"""
type MaxwellFreqParam{T} <: ForwardProbType
    Mesh::T
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    Fields::Array{Complex128}
    freq::Float64
    Ainv::AbstractSolver
end

Base.copy(P::MaxwellFreqParam) = MaxwellFreqParam(P.M, P.Sources, P.Obs, P.Fields, P.freq, P.Ainv)

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
"""
function getMaxwellFreqParam(Mesh::AbstractMesh, 
                             Sources, 
                             Obs, 
                             fields, 
                             freq, 
                             linSolParam::AbstractSolver)
    
    if isempty(fields)
        fields = Array(Complex128,0,0)
    end

    return MaxwellFreqParam(Mesh, Sources, Obs, fields, freq, linSolParam)
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
"""
type MaxwellFreqParamSE{T} <: ForwardProbType
    Mesh::T
    Sources::Union{Array{Complex128},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    freq::Float64
    Ainv::AbstractSolver
    Sens::Array{Complex128}
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
"""
function getMaxwellFreqParamSE(M::AbstractMesh, 
                               Sources,
                               Obs,
                               freq,
                               linSolParam::AbstractSolver)
    
    Sens = Array(Complex128,0,0)
    
    return MaxwellFreqParamSE(M, Sources, Obs, freq, linSolParam, Sens)
end
