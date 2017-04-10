
export solveMaxFreq

"""
function en, linSolParam = solveMaxFreq(A, rhs, Msig, M, w,linSolParam, flag)

Solves the Maxwell frequency domain forward problem.

inputs: 
        A::SparseMatrixCSC          - Stiffness matrix
        rhs::Array                  - Right hand side
        MsigE::SparseMatrixCSC      - The edge conductivity mass matrix
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
                      MsigE::SparseMatrixCSC,
                      mesh::AbstractMesh,
                      w::Real,
                      linSolParam::AbstractSolver,
                      doTranspose::Int=0)

    linSolParam.sym = 2 # structurally symmetric
    en, linSolParam = solveLinearSystem(A,rhs,linSolParam,doTranspose)

    return en, linSolParam
end

function solveMaxFreq(A::SparseMatrixCSC, 
                      rhs::Union{Array,SparseMatrixCSC}, 
                      MsigE::SparseMatrixCSC,
                      mesh::AbstractMesh,
                      w::Real,
                      linSolParam::IterativeSolver,
                      doTranspose::Int=0)

    # setup preconditioner using Aphi system
    if linSolParam.doClear == 1
        MmuN = getNodalMassMatrix(mesh,vec(zeros(mesh.nc).+1/mu0))
        Grad = getNodalGradientMatrix(mesh)

        STBa = Grad*MmuN*Grad'
        Aap = [A + STBa          -(im*w)*MsigE*Grad; 
             -(im*w)*Grad'*MsigE  -(im*w)*Grad'*MsigE*Grad]

        Map(x) = ssor(Aap,x,out=-1)[1]
        P  = [speye(size(A,2)); Grad']
        MM(x) = P'*Map(P*x)

        linSolParam.Ainv = MM
        linSolParam.doClear = 0
    end

    linSolParam.sym = 2 # structurally symmetric
    en, linSolParam = solveLinearSystem(A,rhs,linSolParam,doTranspose)

    return en, linSolParam    

end
