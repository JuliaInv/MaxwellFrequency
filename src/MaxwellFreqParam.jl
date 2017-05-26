
export MaxwellFreqParam, getMaxwellFreqParam

# Supported options for various settings
supportedSensitivityMethods = [:Implicit; :Explicit]

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
    Sens::Array{Complex128}
        - field for storage of explicit sensitivities (if required)
    freq:Float64
        - angular frequency (includes 2*pi term)
    Ainv::AbstractSolver
"""
type MaxwellFreqParam{T} <: ForwardProbType
    Mesh::T
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    Fields::Array{Complex128}
    Sens::Array{Complex128}
    freq::Float64
    Ainv::AbstractSolver
    sensitivityMethod::Symbol
end

Base.copy(P::MaxwellFreqParam) = MaxwellFreqParam(P.M, P.Sources, P.Obs, P.Fields, P.freq, P.Ainv, P.sensitivityMethod)

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
                             linSolParam::AbstractSolver;
                             sensitivityMethod::Symbol=:Implicit)
    
    # Check that user has chosen valid settings for categorical options
    in(sensitivityMethod,supportedSensitivityMethods) || error("Invalid sensitivity method")

    if isempty(fields)
        fields = Array(Complex128,0,0)
    end

    return MaxwellFreqParam(Mesh, Sources, Obs, fields, Array(Complex128,0,0), freq, linSolParam, sensitivityMethod)
end
