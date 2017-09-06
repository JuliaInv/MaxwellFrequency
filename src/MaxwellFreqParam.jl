
export MaxwellFreqParam, getMaxwellFreqParam

# Supported options for various settings
supportedSensitivityMethods = [:Implicit; :Explicit]
supportedStorageLevels      = [:Factors; :Matrices; :None]

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
    Matrices::Array{SparseMatrixCSC}
        - field for storing system matrix
    freq:Float64
        - angular frequency (includes 2*pi term)
    Ainv::AbstractSolver
    sensitivityMethod::Symbol
        - Options are :Implicit, :Explicit. Default value is :Implicit.
    storageLevel::Symbol
        - Options are :Factors, :Matrices, :None.
"""
type MaxwellFreqParam <: ForwardProbType
    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    Fields::Array{Complex128}
    Sens::Array{Complex128}
    Matrices::Array{SparseMatrixCSC}
    freq::Float64
    Ainv::AbstractSolver
    sensitivityMethod::Symbol
    storageLevel::Symbol
end

Base.copy(P::MaxwellFreqParam) = MaxwellFreqParam(P.M, P.Sources, P.Obs, P.Fields, P.Sens, P.Matrices, P.freq, P.Ainv,
                                                  P.sensitivityMethod, P.storageLevel)

"""
function getMaxwellFreqParam

constructor for MaxwellFreqParam

Required Inputs

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    freq:Float64
        - angular frequency (includes 2*pi term)
    linSolParam::AbstractSolver

Optional keyword arguments:

    Fields::Array{Complex128}
        - solution to the fwd problem
    sensitivityMethod::Symbol
        - Options are :Implicit, :Explicit. Default value is :Implicit.
    storageLevel::Symbol
        - Options are :Factors, :Matrices, :None.

"""
function getMaxwellFreqParam(Mesh::AbstractMesh,
                             Sources,
                             Obs,
                             freq,
                             linSolParam::AbstractSolver;
                             fields::Array{Complex128}=Array{Complex128}(0,0),
                             sensitivityMethod::Symbol=:Implicit,
                             storageLevel::Symbol=:Factors)

    # Check that user has chosen valid settings for categorical options
    in(sensitivityMethod,supportedSensitivityMethods) || error("Invalid sensitivity method")
    in(storageLevel,supportedStorageLevels) || error("Unknown storageLevel selection")

    return MaxwellFreqParam(Mesh, Sources, Obs, fields,
                            Array{Complex128}(0, 0), Array{SparseMatrixCSC}(0, 0),
                            freq, linSolParam, sensitivityMethod, storageLevel)
end
