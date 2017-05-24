import jInv.ForwardShare.getSensMatVec
export getSensMatVec

# = ========= The forward sensitivity ===========================

function getSensMatVec(x::Vector,
				   	     sigma::Vector,  # conductivity on fwd mesh
					        param::MaxwellFreqParam)
    # Sens Mat Vec for FV disctretization on OcTreeMesh
    # matv = J*x
	U    = param.Fields
	w    = param.freq
	P    = param.Obs
	
    Ne,   = getEdgeConstraints(param.Mesh)
    Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
    Msig = Ne' * Msig * Ne
    
    A = getMaxwellFreqMatrix(sigma, param)
	
	matv = zeros(Complex128,size(P,2),size(U,2))

	for i=1:size(U,2)
		u = U[:,i] 
		dAdm    = getdEdgeMassMatrix(param.Mesh,u)
		z       = -im*w*Ne'*dAdm*x
		z,      = solveMaxFreq(A,z,Msig,param.Mesh,w,param.Ainv,0)
		z       = vec(Ne*z)
		matv[:,i] = -P'*z	
	end
	return vec(matv)
end # function getSensMatVec for FV OcTreeMesh

function getSensMatVec(x::Vector,sigma::Vector,param::MaxwellFreqParamSE)
	if isempty(param.Sens)
		warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
		Dc,param = getData(sigma,param)
	end
	return param.Sens*x
end
