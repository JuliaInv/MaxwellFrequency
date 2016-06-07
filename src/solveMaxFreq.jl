export solveMaxFreq

function solveMaxFreq(A, rhs, Msig,
                      M::AbstractMesh,
                      w::Real,
                      linSolParam::AbstractSolver,
                      flag=0)
	if linSolParam.doClear == 1
		iterTol = 1e-10
		mu = 4*pi*1e-7
		# setup preconditioner using Aphi system

		Grad = getNodalGradientMatrix(M)
		Ace  = getEdgeAverageMatrix(M)
		Ane  = abs(Grad)
		V    = getVolume(M)
		v    = diag(V)

		#println(size(Ane'),"   ",size(Ace'),"    ",size(v))
		muInvCells = Ane'*(Ace'*v)
		Mmuinvn    = sdiag(muInvCells)

		STBa = Grad*Mmuinvn*Grad'
		Aap = [A + STBa          -(im*w)*Msig*Grad; 
				 -(im*w)*Grad'*Msig  -(im*w)*Grad'*Msig*Grad];
 
		Map(x) = ssor(Aap,x,out=-1)[1]

		P  = [speye(size(A,2)); Grad']
		MM(x) = P'*Map(P*x)
		linSolParam.Ainv = MM
	end
	en, = solveLinearSystem(A,rhs,linSolParam,flag)

	return en, linSolParam
end  # function solveMaxFreq



function solveMaxFreq(A, rhs, Msig,
                      M::AbstractMesh,
                      w::Real,
                      linSolParam::MUMPSsolver,
                      flag=0)
# 
# 	Solve the maxwell system using MUMPS
#
# 	en, = solveMaxFreq(A,rhs,Msig,M::AbstractMesh,w::Real,linSolParam::MUMPSsolver,flag=0)

#

sym = 2 # structurally symmetric
en, = solveLinearSystem(A,rhs,linSolParam,flag)
#en, = solveLinearSystem(A,rhs,linSolParam,sym,flag)

return en, linSolParam
end # function solveMaxFreq
