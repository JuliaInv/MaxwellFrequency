export getSensMatVec

# = ========= The forward sensitivity ===========================

function getSensMatVec(x::Vector,
				   	     sigma::Vector,  # conductivity on fwd mesh
					        param::MaxwellFreqParam)
   # Sens Mat Vec for FV disctretization on OcTreeMesh
   # matv = J*x
 #  mu   = 4*pi*1e-7
   U    = param.Fields
   w    = param.freq
   P    = param.Obs
   
   # eliminate hanging edges and faces
   Ne,   = getEdgeConstraints(param.Mesh)
   iw = complex(0., w)

   A = spzeros(1,1)    # not needed
   Msig = spzeros(1,1) # not needed
   
   matv = zeros(Complex128,size(P,2),size(U,2))

   for i=1:size(U,2)
     # u = U[:,i] 
     # dAdm    = getdEdgeMassMatrix(param.Mesh,u)
      dAdmx    = DerivativeTimesVector(param.Mesh, U[:,i], x)
     # z       = -(im*w) * Ne' * (dAdm*x)
      z       = -iw * (Ne' * dAdmx)
      z,      = solveMaxFreq(A,z,Msig,param.Mesh,w,param.Ainv,0)  # Solve A*z = z
      z       = vec(Ne*z)
      matv[:,i] = -P'*z 
   end  # i
   
   return vec(matv)
end # function getSensMatVec for FV OcTreeMesh

function getSensMatVec(x::Vector,sigma::Vector,param::MaxwellFreqParamSE)
	if isempty(param.Sens)
		warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
		Dc,param = getData(sigma,param)
	end
	return param.Sens*x
end
