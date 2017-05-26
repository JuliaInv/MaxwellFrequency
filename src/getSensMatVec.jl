import jInv.ForwardShare.getSensMatVec
export getSensMatVec

function getSensMatVec(x::Vector, sigma::Vector, param::MaxwellFreqParam)
    # Sens Mat Vec for FV disctretization
    # matv = J*x

    if param.sensitivityMethod == :Implicit

        U = param.Fields
        w = param.freq
        P = param.Obs
        Ne, = getEdgeConstraints(param.Mesh)
        A = getMaxwellFreqMatrix(sigma, param)
        
        matv = zeros(Complex128, size(P, 2), size(U, 2))

        for i=1:size(U, 2)
            u = U[:, i] 
            dAdm = getdEdgeMassMatrix(param.Mesh, u)
            z = -im*w*Ne'*dAdm*x
            z, = solveMaxFreq(A, z, sigma, param.Mesh, w, param.Ainv, 0)
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
