
export solveMaxFreq

"""
function en, linSolParam = solveMaxFreq(A, rhs, m, M, w, linSolParam, flag)

Solves the Maxwell frequency domain forward problem.

inputs:
      A::SparseMatrixCSC          - Stiffness matrix
      rhs::Array                  - Right hand side
      m::MaxwellFreqModel         - Model
      mesh::AbstractMesh          - The mesh
      f::Float64                  - The frequency
      linSolParam::AbstractSolver - Linear solver parameters
      doTranspose::Int=0          - Solve the transposed system

output:
      en:Array{Complex{Float64}}  - Solution (electric fields on cell edges)
      linSolParam::AbstractSolver - Linear solver parameters

"""
function solveMaxFreq(A::SparseMatrixCSC,
                      rhs::Union{Array,SparseMatrixCSC},
                      m::MaxwellFreqModel,
                      mesh::AbstractMesh,
                      f::Float64,
                      linSolParam::AbstractSolver,
                      doTranspose::Int=0)

    en, linSolParam = solveLinearSystem(A,rhs,linSolParam,doTranspose)

    return en, linSolParam
end

function solveMaxFreq(A::SparseMatrixCSC,
                      rhs::Union{Array,SparseMatrixCSC},
                      m::MaxwellFreqModel,
                      mesh::AbstractMesh,
                      f::Float64,
                      linSolParam::IterativeSolver,
                      doTranspose::Int=0)

    sigma = m.values["sigmaCell"]
    muInv = haskey(m.values,"muCell") ? 1./m.values["muCell"] : fill(1./mu0, mesh.nc)

    # setup preconditioner using Aphi system
    if linSolParam.doClear == 1
        MmuN = getNodalMassMatrix(mesh, muInv)
        Grad = getNodalGradientMatrix(mesh)

        Ne, = getEdgeConstraints(mesh)
        MsigE = getEdgeMassMatrix(mesh, sigma)
        MsigE = Ne' * MsigE * Ne

        iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * f

        STBa = Grad*MmuN*Grad'
        Aap = [A + STBa         iw*MsigE*Grad
               iw*Grad'*MsigE   iw*Grad'*MsigE*Grad]

        Map(x) = ssor(Aap,x,out=-1)[1]
        P  = [speye(size(A,2)); Grad']
        MM(x) = P'*Map(P*x)

        linSolParam.Ainv = MM
        linSolParam.doClear = 0
    end

    en, linSolParam = solveLinearSystem(A,rhs,linSolParam,doTranspose)

    return en, linSolParam

end

function solveMaxFreq(rhs::Union{Array,SparseMatrixCSC},
                      m::MaxwellFreqModel,
                      param::MaxwellFreqParam,
                      doTranspose::Int=0)

    if param.storageLevel == :Matrices
        if isempty(param.Matrices)
            A = getMaxwellFreqMatrix(m, param)
            push!(param.Matrices, A)
        else
            A = param.Matrices[1]
        end
    elseif param.storageLevel == :Factors
        if param.Ainv.Ainv == []
            A = getMaxwellFreqMatrix(m, param)
        else
            A = speye(Complex{Float64}, 0)
        end
    else
        A = getMaxwellFreqMatrix(m, param)
    end

    en, param.Ainv = solveMaxFreq(A,rhs,m,param.Mesh,param.frequency,param.Ainv,doTranspose)
    if (param.storageLevel != :Factors) & ~isa(param.Ainv, IterativeSolver)
        clear!(param.Ainv)
    end

    return en, param.Ainv
end
