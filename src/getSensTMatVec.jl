import jInv.ForwardShare.getSensTMatVec
export getSensTMatVec

function getSensTMatVec(x::Vector, sigma::Vector, param::MaxwellFreqParam)
    # SensT Mat Vec for FV disctretization

    if param.sensitivityMethod == :Implicit
        U = param.Fields
        w = param.freq
        P = param.Obs
        Ne, = getEdgeConstraints(param.Mesh)

        A = getMaxwellFreqMatrix(sigma, param)

        X    = reshape(complex(x), size(P,2), size(U,2))
        matv = zeros(size(sigma))
      
        for i=1:size(U, 2)
            u = U[:, i] 
            dAdm = getdEdgeMassMatrix(param.Mesh, u)
            dAdm = -im*w*Ne'*dAdm
            z = -Ne'*(P*X[:, i])
            z, = solveMaxFreq(A, z, sigma, param.Mesh, w, param.Ainv, 1)
            z = vec(z)
            matv += real(dAdm'*z)
        end

    elseif param.sensitivityMethod == :Explicit
        if isempty(param.Sens)
            warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
            Dc, param = getData(sigma, param)
        end
        matv = real(param.Sens'*x)
    else
        error("getSensMatVec: Invalid sensitivity method")
    end

    return matv
end
