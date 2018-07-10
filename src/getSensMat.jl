
import jInv.ForwardShare.getSensMat
export getSensMat

"""
function data, param = getSensMat(m, param)

Computes the data observed by survey defined in param for MaxwellFreqModelmodel
m. Used by getData for when sensitivityMethod = :Explicit.

inputs:

        m::MaxwellFreqModel     - A conductivity and susceptibility model
        param::MaxwellFreqParam - MaxwellFreqParam

output:
        sensMat:Array{Complex{Float64}} - Sensitivity Matrix

"""
function getSensMat(m::MaxwellFreqModel, param::MaxwellFreqParam)
    @assert m.activeInversionProperties == ["sigmaCell"] "Explicit sensitivities only supported for cell conductivity inversions"
    sigma = m.values["sigmaCell"]

    if param.sensitivityMethod == :Implicit
        x = eye(size(param.Sources,2)*size(param.Obs,2))
        if isempty(param.Fields)
            d, param = getData(sigma, param)
        end
        U = param.Fields
        P = param.Obs

        Ne, = getEdgeConstraints(param.Mesh)

        iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * param.frequency

        X = reshape(complex(x),size(P,2),size(U,2),size(x,2))
        sensMat = zeros(Complex128,length(sigma),size(x,2))

        for i=1:size(U,2)
            Z = -Ne'*(conj(P)*X[:,i,:])
            Z, = solveMaxFreq(Z,m,param,1)
            u = U[:,i]
            dAdm = getdEdgeMassMatrix(param.Mesh,sigma,u)
            dAdm = iw*Ne'*dAdm
            sensMat += (dAdm'*Z)
        end
    elseif param.sensitivityMethod == :Explicit
        d, param = getData(sigma, param)
        sensMat = param.Sens'
    end

    return sensMat
end

function getSensMat(sigma::Vector{Float64}, param::MaxwellFreqParam)
    values = Dict{String,Vector{Float64}}("sigmaCell"=>sigma,"muCell"=>fill(mu0,param.Mesh.nc))
    activeInversionProperties = ["sigmaCell"]
    m = MaxwellFreqModel(values,activeInversionProperties)
    getSensMat(m,param)
end
