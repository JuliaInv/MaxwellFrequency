import jInv.ForwardShare.getSensTMatVec
export getSensTMatVec

"""
function matv = getSensTMatVec(x, sigma, param)

Multiplies the transposed sensitvity matrix by a vector

inputs:
        x::Vector               - A vector
        sigma::T                - A model. T can be a vector or a MaxwellFreqModel
        param::MaxwellFreqParam - MaxwellFreqParam

output:
        MaxwellFreqModel containing (dData/dm_i).T*x for each physical property
        i active in the inversion.

"""
function getSensTMatVec(x::Vector{Complex{Float64}},sigma::Vector{Float64},
                                 param::MaxwellFreqParam)
    activeInversionProperties = ["sigmaCell"]
    mvalues = Dict{String,Vector{Float64}}("sigmaCell"=>sigma,"muCell"=>fill(mu0,param.Mesh.nc))
    m  = MaxwellFreqModel(mvalues,activeInversionProperties)
    JTx = getSensTMatVec(x,m,param)
    return JTx.values["sigmaCell"]
end
function getSensTMatVec(x::Vector, m::MaxwellFreqModel, param::MaxwellFreqParam)

    Mesh  = param.Mesh
    sigma = m.values["sigmaCell"]
    mu    = m.values["muCell"]
    JTx   = Dict{String,Vector{Float64}}()

    if param.sensitivityMethod == :Implicit
        U   = param.Fields
        P   = param.Obs
        Ne, = getEdgeConstraints(Mesh)
        invertSigma = false; invertMu = false
        if in("sigmaCell",m.activeInversionProperties)
            invertSigma      = true
            JTx["sigmaCell"] = zeros(size(sigma))
        end
        if in("muCell",m.activeInversionProperties)
            invertMu      = true
            DmuinvDmu     = spdiagm(-1./(mu.^2))
            Nf,Qf,        = getFaceConstraints(Mesh)
            Curl          = getCurlMatrix(Mesh)
            Curl          = Qf*Curl
            JTx["muCell"] = zeros(size(mu))
        end

        iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * param.frequency

        X    = reshape(complex(x), size(P,2), size(U,2))

        for i=1:size(U, 2)
            u  = U[:, i]
            if size(P,1) == size(Ne,2)
                z = -conj(P)*X[:,i]                
            elseif size(P,1) == size(Ne,1)
                z  = -Ne'*(conj(P)*X[:, i])
            end
            z, = solveMaxFreq(z, m, param, 1)
            z  = vec(z)
            if invertSigma
                dAdsigTz = dEdgeMassMatrixTrTimesVector(Mesh, sigma, u,
                                                           conj(iw)*Ne*z)
                JTx["sigmaCell"] += real(dAdsigTz)
            end
            if invertMu
                NfCurlu = Nf*(Curl*u)
                dAdmuTz = dFaceMassMatrixTrTimesVector(Mesh, mu, NfCurlu,
                                                          Nf*(Curl*(Ne*z)))
                JTx["muCell"] += real(DmuinvDmu*dAdmuTz)
            end
        end
    elseif param.sensitivityMethod == :Explicit
        if isempty(param.Sens)
            warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
            Dc, param = getData(m, param)
        end
        JTx = Dict{String,Vector{Float64}}("sigmaCell"=>real(param.Sens'*x))
    else
        error("getSensMatVec: Invalid sensitivity method")
    end

    return MaxwellFreqModel(JTx,m.activeInversionProperties)
end
