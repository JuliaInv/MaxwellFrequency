import jInv.ForwardShare.getSensMatVec
export getSensMatVec

"""
function matv = getSensMatVec(x, sigma, param)

Multiplies the sensitvity matrix by a vector

inputs:
        x::T                    - A vector or MaxwellFreqModel
        sigma::T                - A conductivity (and possibly
                                  susceptibility) model
        param::MaxwellFreqParam - MaxwellFreqParam
        where T can be a vector or MaxwellFreqModel

output:
        matv:Array{Complex{Float64}}  - Solution (J*x)

"""
function getSensMatVec(x::Vector, sigma::Vector{Float64}, param::MaxwellFreqParam)
    activeInversionProperties = ["sigmaCell"]
    mvalues = Dict{String,Vector{Float64}}("sigmaCell"=>sigma,"muCell"=>fill(mu0,param.Mesh.nc))
    xvalues = Dict{String,Vector{Float64}}("sigmaCell"=>x,"muCell"=>fill(mu0,param.Mesh.nc))
    m  = MaxwellFreqModel(mvalues,activeInversionProperties)
    xm = MaxwellFreqModel(xvalues,activeInversionProperties)
    getSensMatVec(xm,m,param)
end

function getSensMatVec(x::MaxwellFreqModel, m::MaxwellFreqModel, param::MaxwellFreqParam)

    # Recompute fields if necessary
    if isempty(param.Fields)
        dDummy,param = getData(m,param)
    end


    Mesh = param.Mesh
    sigma = m.values["sigmaCell"]
    mu    = m.values["muCell"]

    if param.sensitivityMethod == :Implicit

        U   = param.Fields
        P   = param.Obs
        Ne, = getEdgeConstraints(Mesh)

        invertSigma = false; invertMu = false
        if in("sigmaCell",x.activeInversionProperties)
            invertSigma = true
            DsigDmx     = x.values["sigmaCell"]
        end
        if in("muCell",x.activeInversionProperties)
            invertMu  = true
            DmuDmx    = x.values["muCell"]
            DmuinvDmu = spdiagm(-1./(mu.^2))
            Nf,Qf,    = getFaceConstraints(Mesh)
            Curl      = getCurlMatrix(Mesh)
            Curl      = Qf*Curl
        end

        iw = (param.timeConvention == :ExpMinusImOmegaT ? -im : im) * 2 * pi * param.frequency

        matv = zeros(Complex128, size(P, 2), size(U, 2))
        for i=1:size(U, 2)
            u = U[:, i]
            nActiveEdges = 0
            try
                nActiveEdges = size(Ne,2)
            catch # size(Ne,2) doesn't work if Ne is uniformScaling
                nActiveEdges = sum(Mesh.ne)
            end
            z = zeros(Complex128,nActiveEdges)
            if invertSigma
                dAdsigx = dEdgeMassMatrixTimesVector(Mesh, sigma, u, DsigDmx)
                z += iw*Ne'*dAdsigx
            end
            if invertMu
                NfCurlu     = Nf*(Curl*u)
                dMmuinvdmx  = dFaceMassMatrixTimesVector(Mesh, sigma, NfCurlu,
                                                         DmuinvDmu*DmuDmx)
                z          += Ne'*(Curl'*(Nf'*dMmuinvdmx))
            end
            z, = solveMaxFreq(z, m, param, 0)
            if size(P,1) == size(Ne,2)
                z = vec(z)
            elseif size(P,1) == size(Ne,1)
                z = vec(Ne*z)
            else
                error("Invalid size of P")
            end
            matv[:, i] = -P.'*z
        end

    elseif param.sensitivityMethod == :Explicit

        if isempty(param.Sens)
            warn("getSensTMatVec: Recomputing data to get sensitvity matrix. This should be avoided.")
            Dc,param = getData(m,param)
        end
        matv = param.Sens*x.values["sigmaCell"]

    else
        error("getSensTMatVec: Invalid sensitivity method")
    end
    return vec(matv)
end
