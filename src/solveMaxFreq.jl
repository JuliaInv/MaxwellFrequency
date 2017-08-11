export solveMaxFreq

"""
function en, linSolParam = solveMaxFreq(A, rhs, sigma, M, w, linSolParam, flag)

Solves the Maxwell frequency domain forward problem.

inputs:
      A::SparseMatrixCSC          - Stiffness matrix
      rhs::Array                  - Right hand side
      sigma::Vector{Float64}      - Conductivity model
      mesh::AbstractMesh          - The mesh
      w::Real                     - The angular frequency being solved
      linSolParam::AbstractSolver - Linear solver parameters
      doTranspose::Int=0          - Solve the transposed system

output:
      en:Array{Complex{Float64}}  - Solution (electric fields on cell edges)
      linSolParam::AbstractSolver - Linear solver parameters

"""
function solveMaxFreq(A::SparseMatrixCSC,
                      rhs::Union{Array,SparseMatrixCSC},
                      sigma::Vector{Float64},
                  #    mesh::AbstractMesh,
                  #    w::Real,
                      linSolParam::AbstractSolver,
                      doTranspose::Int=0)

    linSolParam.sym = 2 # structurally symmetric
    en, linSolParam = solveLinearSystem(A,rhs,linSolParam,doTranspose)

    return en, linSolParam
end

function solveMaxFreq(A::SparseMatrixCSC,
                      rhs::Union{Array,SparseMatrixCSC},
                      Msig,
                      mesh::AbstractMesh,
                      w::Real,  # omega (=0 when called from solveMTsystem)
                      linSolParam::IterativeSolver,
                      doTranspose::Int=0)    # doTranspose=0 for solving Ax=rhs, =1 for A'x=rhs
# A, Msig, M, and w are only needed when initializing (doClear).

    # setup preconditioner using Aphi system
   if linSolParam.doClear == 1
      if w != 0.
         println("getting preconditioner in solveMaxFreq")
         iterTol = 1e-10
         mu = 4*pi*1e-7
         # setup preconditioner using Aphi system

         Grad = getNodalGradientMatrix(mesh)
        # Ace  = getEdgeAverageMatrix(M)
        # Ane  = abs(Grad)
         V    = getVolume(mesh)
         v    = diag(V)

         An2c = getNodalAverageMatrix(mesh)
         Nn,Qn = getNodalConstraints(mesh)
         Ne,Qe = getEdgeConstraints(mesh)
         Grad = Qe * Grad * Nn
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

      end
   end  # if linSolParam.doClear == 1
   
   en, = solveLinearSystem(A,rhs,linSolParam,doTranspose)

   return en, linSolParam
end  # function solveMaxFreq


function solveMaxFreq(rhs::Union{Array,SparseMatrixCSC},
                      sigma::Vector{Float64},
                      param::MaxwellFreqParam,
                      doTranspose::Int=0)

    if param.storageLevel == :Matrices
        if isempty(param.Matrices)
            A = getMaxwellFreqMatrix(sigma, param)
            push!(param.Matrices, A)
        else
            A = param.Matrices[1]
        end
    elseif param.storageLevel == :Factors
        if param.Ainv.Ainv == []
            A = getMaxwellFreqMatrix(sigma, param)
        else
            A = speye(Complex{Float64}, 0)
        end
    else
        A = getMaxwellFreqMatrix(sigma, param)
    end

  #  en, param.Ainv = solveMaxFreq(A,rhs,sigma,param.Mesh,param.freq,param.Ainv,doTranspose)
    en, param.Ainv = solveMaxFreq(A,rhs,sigma, param.Ainv,doTranspose)
    if (param.storageLevel != :Factors) & ~isa(param.Ainv, IterativeSolver)
        clear!(param.Ainv)
    end

    return en, param.Ainv
end
