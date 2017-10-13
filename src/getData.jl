import jInv.ForwardShare.getData
export getData

"""
function data, param = getData(sigma, param[, doClear=false])

Computes the data observed by survey defined in param for conductivity model sigma

inputs:

        sigma::Vector           - A conductivity model
        param::MaxwellFreqParam - MaxwellFreqParam
        doClear::Bool           - ???

output:
        data:Array{Complex{Float64}} - Computed data
        param::MaxwellFreqParam      - MaxwellFreqParam

"""
function getData(sigma::Vector{Float64}, param::MaxwellFreqParam, doClear::Bool=false)

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
        rhs = -iw * (Ne' * S)

        param.Ainv.doClear = 1
        U, param.Ainv = solveMaxFreq(rhs, sigma, param, 0)
        param.Ainv.doClear = 0

        U = Ne * U
        D = P' * U

        param.Fields = U
        if doClear
            # clear fields and factorization
            clear!(param,clearAll=false)
        end

    elseif param.sensitivityMethod == :Explicit

        tempParam = getMaxwellFreqParam(param.Mesh, param.Sources, param.Obs, param.frequency, param.Ainv,
          sensitivityMethod = :Implicit, storageLevel = param.storageLevel, timeConvention = param.timeConvention)
        param.Sens=[]
        D,tempParam = getData(sigma,tempParam)
        J = getSensMat(sigma, tempParam)
        param.Sens = J'
        clear!(tempParam)

    else
        error("getData: Invalid sensitivity method")

    end

	return D, param
end

"""
function data, param = getSensMat(sigma, param)

Computes the data observed by survey defined in param for conductivity model sigma.
Used by getData for when sensitivityMethod = :Explicit.

inputs:

        sigma::Vector           - A conductivity model
        param::MaxwellFreqParam - MaxwellFreqParam

output:
        sensMat:Array{Complex{Float64}} - Sensitivity Matrix

"""
function getSensMat(sigma::Vector{Float64}, param::MaxwellFreqParam)
    x = eye(size(param.Sources,2)*size(param.Obs,2))
    U = param.Fields
    P = param.Obs

    Ne, = getEdgeConstraints(param.Mesh)

    iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * param.frequency

    X = reshape(complex(x),size(P,2),size(U,2),size(x,2))
    sensMat = zeros(Complex128,length(sigma),size(x,2))

    for i=1:size(U,2)
        Z = -Ne'*(P*X[:,i,:])
        Z, = solveMaxFreq(Z,sigma,param,1)
        u = U[:,i]
        dAdm = getdEdgeMassMatrix(param.Mesh,sigma,u)
        dAdm = iw*Ne'*dAdm
        sensMat += (dAdm'*Z)
    end

    return sensMat
end
