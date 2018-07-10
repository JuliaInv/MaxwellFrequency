export MaxwellFreqModel, MaxwellFreqModelDerivative, MaxwellFreqModelTransform

import jInv.InverseSolve: AbstractModel,AbstractModelDerivative,
                          AbstractModelTransform

# struct MaxwellFreqModel{T<:Real} <: AbstractModel
#   sigmaCell::Vector{T}
#   sigmaEdge::Vector{T}
#   muCell::Vector{T}
#   muFace::Vector{T}
#   activeInversionProperties::Vector{String}
# end
validProperties = ["sigmaCell","muCell"] # Add later:"sigmaEdge","muFace"]
struct MaxwellFreqModel{T<:Real} <: AbstractModel
    values::Dict{String,Vector{T}}
    activeInversionProperties::Vector{String}
end


# struct MaxwellFreqModelDerivative{A1<:AbstractArray{T,N}   where {T,N},
#                                   A2<:AbstractArray{T2,N2} where {T2,N2},
#                                   A3<:AbstractArray{T3,N3} where {T3,N3},
#                                   A4<:AbstractArray{T4,N4} where {T4,N4}} <:
#                                   AbstractModelDerivative
#     DSigmaCell::A1
#     DSigmaEdge::A2
#     DMuCell::A3
#     DMuFace::A4
# end
struct MaxwellFreqModelDerivative <: AbstractModelDerivative
    mats::Dict{String,AbstractArray}
    activeInversionProperties::Vector{String}
end

struct MaxwellFreqModelTransform{T<:Real,N<:Integer} <: AbstractModelTransform
    PCell::SparseMatrixCSC{T,N}
    PEdge::SparseMatrixCSC{T,N}
end

function MaxwellFreqModel()
    return MaxwellFreqModel(Dict{String,Vector{Float64}}(),Vector{String}())
end

function MaxwellFreqModel{T<:Real}(sigma::Vector{T})
    return MaxwellFreqModel(Dict("sigmaCell"=>sigma),["sigmaCell"])
end

import Base.isempty
isempty(m::MaxwellFreqModel) = isempty(m.values)
# import Base.zeros
# zeros(sigma::MaxwellFreqModel) = MaxwellFreqModel(zeros(sigma.sigmaCell), zeros(sigma.sigmaEdge))

# Import stuff from base in order to define multiplication and addition
# for MaxwellFreqModels.
import Base.+
import Base.-
import Base.*
import Base.Ac_mul_B # A' * B
# export A_mul_B # not defined in Base

# Define addition and subtraction of two models
function +(m1::MaxwellFreqModel, m2::MaxwellFreqModel)
    values = merge(+,m1.values,m2.values)
    active = unique([m1.activeInversionProperties;m2.activeInversionProperties])
    return MaxwellFreqModel(values,active)
end

function -(m1::MaxwellFreqModel, m2::MaxwellFreqModel)
    values = merge(-,m1.values,m2.values)
    active = unique([m1.activeInversionProperties;m2.activeInversionProperties])
    return MaxwellFreqModel(values,active)
end

# Define multiplication of model by scalar
# Should this be in place?
function *(a::Union{Real,UniformScaling}, m::MaxwellFreqModel)
    for key in keys(m.values)
        m.values[key] = a*m.values[key]
    end
  return m
end

# Define multiplication by model derivative
# Vector{T} → MaxwellFreqModelDerivative
function *{T<:Real}(D::MaxwellFreqModelDerivative, x::Vector{T})
    values = Dict{String,Vector{T}}()
    for key in D.activeInversionProperties
        values[key] = D.mats[key]*x
    end
    MaxwellFreqModel(values,D.activeInversionProperties)
end

function Ac_mul_B(D::MaxwellFreqModelDerivative, y::MaxwellFreqModel)
    key1 = first(keys(D.mats))
    T    = eltype(D.mats[key1])
    xOut = zeros(T,size(D.mats[key1],2))
    for key in D.activeInversionProperties
        xOut += D.mats[key]'*y.values[key]
    end
    return xOut
end

# Define multiplication by interpolation matrix
function *{T}(P::MaxwellFreqModelTransform, x::MaxwellFreqModel{T})
    values = Dict{String,Vector{T}}()
    for key in x.activeInversionProperties
        if in(key,["sigmaCell","muCell"])
            values[key] = P.PCell*x.values[key]
        elseif in(key,["sigmaEdge","muEdge"])
            values[key] = P.PEdge*x.values[key]
        end
    end
    return MaxwellFreqModel(values, x.activeInversionProperties)
end

function Ac_mul_B{T,N}(P::MaxwellFreqModelTransform{T,N}, y::MaxwellFreqModel{T})
    values = Dict{String,Vector{T}}()
    for key in y.activeInversionProperties
        if in(key,["sigmaCell","muCell"])
            values[key] = P.PCell'*y.values[key]
        elseif in(key,["sigmaEdge","muEdge"])
            values[key] = P.PEdge'*y.values[key]
        end
    end
    return MaxwellFreqModel(values, y.activeInversionProperties)
end
