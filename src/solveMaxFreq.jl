export solveMaxFreq

function solveMaxFreq(A, rhs, Msig,
                      M::AbstractMesh,
                      w::Real,  # omega (=0 when called from solveMTsystem)
                      linSolParam::AbstractSolver,
                      flag=0)    # flag=0 for solving Ax=rhs, =1 for A'x=rhs
# A, Msig, M, and w are only needed when initializing (doClear).

   if linSolParam.doClear == 1
      if w != 0.
         println("getting preconditioner in solveMaxFreq")
         iterTol = 1e-10
         mu = 4*pi*1e-7
         # setup preconditioner using Aphi system

         Grad = getNodalGradientMatrix(M)
        # Ace  = getEdgeAverageMatrix(M)
        # Ane  = abs(Grad)
         V    = getVolume(M)
         v    = diag(V)

         An2c = getNodalAverageMatrix(M)
         Nn,Qn = getNodalConstraints(M)
         Ne,Qe = getEdgeConstraints(M)
         Grad = Qe * Grad * Nn
         #println(size(Ane'),"   ",size(Ace'),"    ",size(v))
         #muInvCells = Ane'*(Ace'*v)
         muInvCells = An2c' * (v/mu)
         Mmuinvn    = sdiag(muInvCells)
         iw = complex(0., w)

         #STBa = Grad*Mmuinvn*Grad'
         STBa = Grad * Nn' * Mmuinvn * Nn * Grad'
         Nn=[] ; Qn=[] ; Ne=[] ; Qe=[] ; An2c=[]

         Aap = [ A + STBa         -iw*Msig*Grad ; 
                -iw*Grad'*Msig    -iw*Grad'*Msig*Grad]

         Map(x) = ssor(Aap,x,out=-1)[1]

         P  = [speye(size(A,2)); Grad']
         MM(x) = P'*Map(P*x)
         linSolParam.Ainv = MM
         
         linSolParam.doClear = false  # make sure the preconditioner isn't cleared in iterativeWrapper
         linSolParam.AA = A  # save for later

#     else
#        Map(x) = ssor(A,x,out=-1)[1]
#        MM(x) = Map(x)
#        linSolParam.Ainv = MM
      end
   end  # if linSolParam.doClear == 1
   
   en, = solveLinearSystem(A,rhs,linSolParam,flag)

   return en, linSolParam
end  # function solveMaxFreq



function solveMaxFreq(A, rhs, Msig, 
                      M::AbstractMesh,
                      w::Real,
                      linSolParam::MUMPSsolver,
                      flag=0)  # flag=0 for solving Ax=rhs, =1 for A'x=rhs

# A is only needed when calling factorMUMPS
# Msig, M, and w are never needed.

# 
# 	Solve the maxwell system using MUMPS
#
# 	en, = solveMaxFreq(A,rhs,Msig,M::AbstractMesh,w::Real,linSolParam::MUMPSsolver,flag=0)

#

linSolParam.sym = 2 # structurally symmetric
en, = solveLinearSystem(A,rhs,linSolParam,flag)

return en, linSolParam
end # function solveMaxFreq
