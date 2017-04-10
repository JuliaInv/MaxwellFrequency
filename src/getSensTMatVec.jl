import jInv.ForwardShare.getSensTMatVec
export getSensTMatVec

# = ========= The Transpose sensitivity ===========================

function getSensTMatVec(x::Vector,sigma::Vector,param::MaxwellFreqParam)
# SensT Mat Vec for FV disctretization on OcTree mesh
	U    = param.Fields
	w    = param.freq
	P    = param.Obs
  
    Ne, = getEdgeConstraints(param.Mesh)
    Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
    Msig = Ne' * Msig * Ne
    A = getMaxwellFreqMatrix(sigma, param)

	X    = reshape(complex(x),size(P,2),size(U,2))
	matv = zeros(size(sigma))
  
	for i=1:size(U,2)
		u     = U[:,i] 
		dAdm  = getdEdgeMassMatrix(param.Mesh,u)
		dAdm  = -im*w*Ne'*dAdm
		z     = -Ne'*(P*X[:,i])
		z,    = solveMaxFreq(A,z,Msig,param.Mesh,w,param.Ainv,1)
		z     = vec(z)
		matv += real(dAdm'*z)
	end
	
	return matv
end

function getSensTMatVec(x::SparseMatrixCSC,sigma::Vector,param::MaxwellFreqParam)
	# SensT Mat Vec for FV disctretization on OcTree mesh
	x = full(x)
	U    = param.Fields
	w    = param.freq
	P    = param.Obs
	
    Ne, = getEdgeConstraints(param.Mesh)
    Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
    Msig = Ne' * Msig * Ne
    A = getMaxwellFreqMatrix(sigma, param)

	X    = reshape(complex(x),size(P,2),size(U,2),size(x,2))
	matv = zeros(Complex128,length(sigma),size(x,2))	

	for i=1:size(U,2)
		if !all(X.==0)
			
			Z     = -Ne'*(P*X[:,i,:])
			Z,    = solveMaxFreq(A,Z,Msig,param.Mesh,w,param.Ainv,1)
			u     = U[:,i] 
			dAdm  = getdEdgeMassMatrix(param.Mesh,u)
			dAdm  = -im*w*Ne'*dAdm
			matv +=  (dAdm'*Z)
		end
	end
	
	return matv
end

function getSensTMatVec(x::Vector,sigma::Vector,param::MaxwellFreqParamSE)
	if isempty(param.Sens)
		warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
		Dc,param = getData(sigma,param)
	end
	return real(param.Sens'*x)
end