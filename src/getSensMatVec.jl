import jInv.ForwardShare.getSensMatVec
export getSensMatVec

"""
function matv = getSensMatVec(x, sigma, param)

Multiplies the sensitvity matrix by a vector

inputs:
        x::Vector               - A vector
        sigma::Vector{Float64}  - A conductivity model
        param::MaxwellFreqParam - MaxwellFreqParam

output:
        matv:Array{Complex{Float64}}  - Solution (J*x)

"""
function getSensMatVec(x::Vector, sigma::Vector{Float64}, param::MaxwellFreqParam)

    if param.sensitivityMethod == :Implicit

        U = param.Fields
        w = param.freq
        P = param.Obs
        Ne, = getEdgeConstraints(param.Mesh)

        if param.useIw
           iw =  complex(0., w)
        else
           iw = -complex(0., w)
        end

        matv = zeros(Complex128, size(P, 2), size(U, 2))

        for i=1:size(U, 2)
           # u = U[:, i]
           # dAdm = getdEdgeMassMatrix(param.Mesh, u)
            dAdmx  = DerivativeTimesVector(param.Mesh, sigma, U[:,i], x)
           # z = -im*w*Ne'*dAdm*x
            z = iw * (Ne' * dAdmx)
            z, = solveMaxFreq(z, sigma, param, 0)
            z = vec(Ne*z)
            matv[:, i] = -P'*z
        end

    elseif param.sensitivityMethod == :Explicit

        if isempty(param.Sens)
            warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
            Dc,param = getData(sigma,param)
        end
        matv = param.Sens*x

    else
        error("getSensTMatVec: Invalid sensitivity method")
    end
    return vec(matv)
end
