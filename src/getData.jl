import jInv.ForwardShare.getData
export getData

"""
function data, param = getData(m, param[, doClear=false])

Computes the data observed by survey defined in param for model m

inputs:

        m::Vector or MaxwellFreqModel - A conductivity model
        param::MaxwellFreqParam           - MaxwellFreqParam
        doClear::Bool                     - ???

output:
        data:Array{Complex{Float64}} - Computed data
        param::MaxwellFreqParam      - MaxwellFreqParam

"""
function getData(sigma::Vector{Float64}, param::MaxwellFreqParam, doClear::Bool=false)
    values = Dict{String,Vector{Float64}}("sigmaCell"=>sigma,"muCell"=>fill(mu0,param.Mesh.nc))
    activeInversionProperties = ["sigmaCell"]
    m = MaxwellFreqModel(values,activeInversionProperties)
    getData(m,param,doClear)
end

function getData(m::MaxwellFreqModel, param::MaxwellFreqParam, doClear::Bool=false)

    if param.sensitivityMethod == :Implicit

        S = param.Sources
        P = param.Obs
        Ne, = getEdgeConstraints(param.Mesh)

        iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * param.frequency

        if param.storageLevel == :Matrices
        # Initialize matrix storage if needed. Note that if a direct solver is
        # being used and factorizatons are being stored, then matrices don't
        # need to be stored
            param.Matrices = Vector{SparseMatrixCSC}()
        elseif param.storageLevel == :Factors
        # Clear the solver
            clear!(param.Ainv)
        end

        # A = getMaxwellFreqMatrix(sigma, param)
        if size(S,1) == size(Ne,1)
            rhs = -iw * (Ne' * S)
        elseif size(S,1) == size(Ne,2)
            rhs = -iw * S
        else
            error("Invalid source term")
        end

        param.Ainv.doClear = 1
        U, param.Ainv = solveMaxFreq(rhs, m, param, 0)
        param.Ainv.doClear = 0
        if size(P,1) == size(Ne,2) # if P = Ne'*P, then P' = P'*Ne
            D = P.' * U
            U = Ne * U
        elseif size(P,1) == size(Ne,1)
            U = Ne * U
            D = P.' * U
        else
            error("Invalid size of P")
        end

        param.Fields = U
        if doClear
            # clear fields and factorization
            clear!(param,clearAll=false)
        end

    elseif param.sensitivityMethod == :Explicit

        tempParam = getMaxwellFreqParam(param.Mesh, param.Sources, param.Obs, param.frequency, param.Ainv,
          sensitivityMethod = :Implicit, storageLevel = param.storageLevel, timeConvention = param.timeConvention)
        param.Sens=[]
        D,tempParam = getData(m,tempParam)
        J = getSensMat(m, tempParam)
        param.Sens = J'
        clear!(tempParam)

    else
        error("getData: Invalid sensitivity method")

    end

	return D, param
end

"""
function data, param = getSensMat(m, param)

Computes the data observed by survey defined in param for MaxwellFreqModelmodel
m. Used by getData for when sensitivityMethod = :Explicit.

inputs:

        m::MaxwellFreqModel     - A conductivity and susceptibility model
        param::MaxwellFreqParam - MaxwellFreqParam

output:
        sensMat:Array{Complex{Float64}} - Sensitivity Matrix

"""
function getSensMat(m::MaxwellFreqModel, param::MaxwellFreqParam)
    @assert m.activeInversionProperties == ["sigmaCell"] "Explicit sensitivities only supported for cell conductivity inversions"
    sigma = m.values["sigmaCell"]
    x = eye(size(param.Sources,2)*size(param.Obs,2))
    U = param.Fields
    P = param.Obs

    Ne, = getEdgeConstraints(param.Mesh)

    iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * param.frequency

    X = reshape(complex(x),size(P,2),size(U,2),size(x,2))
    sensMat = zeros(Complex128,length(sigma),size(x,2))

    for i=1:size(U,2)
        Z = -Ne'*(conj(P)*X[:,i,:])
        Z, = solveMaxFreq(Z,m,param,1)
        u = U[:,i]
        dAdm = getdEdgeMassMatrix(param.Mesh,sigma,u)
        dAdm = iw*Ne'*dAdm
        sensMat += (dAdm'*Z)
    end

    return sensMat
end
