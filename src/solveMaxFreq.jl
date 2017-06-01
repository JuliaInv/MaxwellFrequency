
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
                      sigma::Vector{Float64},
                      mesh::AbstractMesh,
                      w::Real,
                      linSolParam::IterativeSolver,
                      doTranspose::Int=0)

    # setup preconditioner using Aphi system
    if linSolParam.doClear == 1
        MmuN = getNodalMassMatrix(mesh,  fill(1./mu0, mesh.nc))
        Grad = getNodalGradientMatrix(mesh)

        Ne, = getEdgeConstraints(mesh)
        MsigE = getEdgeMassMatrix(mesh, sigma)
        MsigE = Ne' * MsigE * Ne

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

    en, param.Ainv = solveMaxFreq(A,rhs,sigma,param.Mesh,param.freq,param.Ainv,doTranspose)
    if (param.storageLevel != :Factors) & ~isa(param.Ainv, IterativeSolver)
        clear!(param.Ainv)
    end

    return en, param.Ainv
end
