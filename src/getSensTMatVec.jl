export getSensTMatVec

# = ========= The Transpose sensitivity ===========================

function getSensTMatVec(x::Vector,
                        sigma::Vector,
                        param::MaxwellFreqParam)
# SensT Mat Vec for FV disctretization on OcTree mesh
  # mu   = 4*pi*1e-7
   U    = param.Fields
   w    = param.freq
   P    = param.Obs
  
   # eliminate hanging edges and faces
   Ne,   = getEdgeConstraints(param.Mesh)

   A = spzeros(1,1)    # not needed
   Msig = spzeros(1,1) # not needed

   X    = reshape(complex(x),size(P,2),size(U,2))
   matv = zeros(size(sigma))
  
   for i=1:size(U,2)
     # u     = U[:,i] 
     # dAdm  = getdEdgeMassMatrix(param.Mesh,u)
     # dAdm  = -im*w*Ne'*dAdm
      z     = -Ne'*(P*X[:,i])
      z,    = solveMaxFreq(A,z,Msig,param.Mesh,w,param.Ainv,1)  # Solve A'*z = z

     # matv += real(  -(im*w)*dAdm' * (Ne*z) )
      z     = Ne * vec(z)
      dAdmTz  = DerivativeTrTimesVector(param.Mesh, U[:,i], z)
      matv += w*imag(dAdmTz)  # real(  -(im*w)* dAdmTz )
   end  # i

   return matv
end  # function getSensTMatVec


function getSensTMatVec(x::SparseMatrixCSC,sigma::Vector,param::MaxwellFreqParam)
   # SensT Mat Vec for FV disctretization on OcTree mesh
   x = full(x)
 #  mu   = 4*pi*1e-7
   U    = param.Fields
   w    = param.freq
   P    = param.Obs
   
  
  # eliminate hanging edges and faces
   Ne,   = getEdgeConstraints(param.Mesh)
  
   A = spzeros(1,1)    # not needed
   Msig = spzeros(1,1) # not needed
   
   X    = reshape(complex(x),size(P,2),size(U,2),size(x,2))
   #matv = zeros(Complex128,length(sigma),size(P,2),size(U,2))
   matv = zeros(Complex128,length(sigma),size(x,2))   

   for i=1:size(U,2)
      if !all(X.==0)
         
         Z     = -Ne'*(P*squeeze(X[:,i,:],2))
         Z,    = solveMaxFreq(A,Z,Msig,param.Mesh,w,param.Ainv,1)  # Solve A'*Z = Z
         u     = U[:,i] 
         dAdm  = getdEdgeMassMatrix(param.Mesh,u)
         dAdm  = -(im*w)*(Ne'*dAdm)
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