
export MaxwellFreqParam, getMaxwellFreqParam

# Supported options for various settings
supportedSensitivityMethods = [:Implicit; :Explicit]
supportedStorageLevels      = [:Factors; :Matrices; :None]
supportedTimeConventions    = [:ExpMinusImOmegaT, :ExpPlusImOmegaT]

"""
type MaxwellFrequency.MaxwellFreqParam <: ForwardProbType

defines one MaxwellFrequency problem

Fields:

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    Fields::Array{Complex128}
        - solution to the forward problem
    Sens::Array{Complex128}
        - field for storage of explicit sensitivities (if required)
    Matrices::Array{SparseMatrixCSC}
        - field for storing system matrix
    frequency:Float64
        - frequency
    Ainv::AbstractSolver
    sensitivityMethod::Symbol
        - Options are :Implicit, :Explicit. Default value is :Implicit
    storageLevel::Symbol
        - Options are :Factors, :Matrices, :None
    timeConvention::Symbol
        - Options are :ExpMinusImOmegaT, :ExpPlusImOmegaT
    metaData::Dict{Any,Any}
        - additional optional information
"""
type MaxwellFreqParam <: ForwardProbType
    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
    Fields::Array{Complex128}
    Sens::Array{Complex128}
    Matrices::Array{SparseMatrixCSC}
    frequency::Float64
    Ainv::AbstractSolver
    sensitivityMethod::Symbol
    storageLevel::Symbol
    timeConvention::Symbol
    metaData::Dict
end

Base.copy(P::MaxwellFreqParam) = MaxwellFreqParam(P.M, P.Sources, P.Obs, P.Fields, P.Sens, P.Matrices, P.frequency, P.Ainv,
                                                  P.sensitivityMethod, P.storageLevel, P.timeConvention, P.metaData)

"""
function getMaxwellFreqParam

constructor for MaxwellFreqParam

Required Inputs

    Mesh::AbstractMesh
    Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC}
    Obs::Union{Array{Complex128},SparseMatrixCSC}
        - transpose interpolation matrix from fields to receivers
    frequency:Float64
        - frequency
    linSolParam::AbstractSolver

Optional keyword arguments:

    Fields::Array{Complex128}
        - solution to the forward problem
    sensitivityMethod::Symbol
        - Options are :Implicit, :Explicit. Default value is :Implicit
    storageLevel::Symbol
        - Options are :Factors, :Matrices, :None
    timeConvention::Symbol
        - Options are :ExpMinusImOmegaT, :ExpPlusImOmegaT

"""
function getMaxwellFreqParam(Mesh::AbstractMesh,
                             Sources::Union{Array{Complex128},Array{Float64},SparseMatrixCSC},
                             Obs::Union{Array{Complex128},SparseMatrixCSC},
                             frequency::Float64,
                             linSolParam::AbstractSolver;
                             Fields::Array{Complex128}=Array{Complex128}(0,0),
                             sensitivityMethod::Symbol=:Implicit,
                             storageLevel::Symbol=:Factors,
                             timeConvention::Symbol=:ExpMinusImOmegaT,
                             metaData::Dict=Dict())

    # Check that user has chosen valid settings for categorical options
    in(sensitivityMethod,supportedSensitivityMethods) || error("Invalid sensitivity method")
    in(storageLevel,supportedStorageLevels) || error("Unknown storageLevel selection")
    in(timeConvention,supportedTimeConventions) || error("Unknown timeConvention selection")

    Sens = Array{Complex128}(0,0)
    Matrices = Array{SparseMatrixCSC}(0)

    return MaxwellFreqParam(Mesh, Sources, Obs, Fields, Sens, Matrices, frequency, linSolParam,
                            sensitivityMethod, storageLevel, timeConvention, metaData)
end

import Base.clear!
function clear!(param::MaxwellFreqParam)
    param.Fields = Array{Complex128}(0,0)
    param.Sens = Array{Complex128}(0,0)
    param.Matrices = Array{SparseMatrixCSC}(0)
    clear!(param.Ainv)
    return param
end
