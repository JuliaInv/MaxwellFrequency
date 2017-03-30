import jInv.ForwardShare.getSensTMatVec
export getSensTMatVec

# = ========= The Transpose sensitivity ===========================

function getSensTMatVec(x::Vector,sigma::Vector,param::MaxwellFreqParam)
# SensT Mat Vec for FV disctretization on OcTree mesh
	mu   = 4*pi*1e-7
	U    = param.Fields
	w    = param.freq
	P    = param.Obs
  
	Curl = getCurlMatrix(param.Mesh)
	Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
	Mmu  = getFaceMassMatrix(param.Mesh,vec(zeros(size(sigma)).+1/mu))
  
  # eliminate hanging edges and faces
	Ne,   = getEdgeConstraints(param.Mesh)
  Nf,Qf = getFaceConstraints(param.Mesh)
  
  Curl = Qf  * Curl * Ne
  Msig = Ne' * Msig * Ne
  Mmu  = Nf' * Mmu  * Nf

	A   = Curl' * Mmu * Curl - (im * w) * Msig

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
	mu   = 4*pi*1e-7
	U    = param.Fields
	w    = param.freq
	P    = param.Obs
	
	Curl = getCurlMatrix(param.Mesh)
	Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
	Mmu  = getFaceMassMatrix(param.Mesh,fill(1/mu,length(sigma)))
  
  # eliminate hanging edges and faces
	Ne,   = getEdgeConstraints(param.Mesh)
  Nf,Qf = getFaceConstraints(param.Mesh)
  
  Curl = Qf  * Curl * Ne
  Msig = Ne' * Msig * Ne
  Mmu  = Nf' * Mmu  * Nf

	A   = Curl' * Mmu * Curl - (im * w) * Msig
	
	X    = reshape(complex(x),size(P,2),size(U,2),size(x,2))
	#matv = zeros(Complex128,length(sigma),size(P,2),size(U,2))
	matv = zeros(Complex128,length(sigma),size(x,2))	

	for i=1:size(U,2)
		if !all(X.==0)
			
			Z     = -Ne'*(P*X[:,i,:])
			Z,    = solveMaxFreq(A,Z,Msig,param.Mesh,w,param.Ainv,1)
			u     = U[:,i] 
			dAdm  = getdEdgeMassMatrix(param.Mesh,u)
			dAdm  = -im*w*Ne'*dAdm
			#matv[:,:,i] = (dAdm'*Z)
			matv +=  (dAdm'*Z)
		end
	end
	#matv = reshape(matv,length(sigma),size(P,2)*size(U,2))
	
	return matv
end

function getSensTMatVec(x::Vector,sigma::Vector,param::MaxwellFreqParamSE)
	if isempty(param.Sens)
		warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
		Dc,param = getData(sigma,param)
	end
	return real(param.Sens'*x)
end