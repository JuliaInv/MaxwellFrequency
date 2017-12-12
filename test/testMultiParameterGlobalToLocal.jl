using MaxwellFrequency
using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
using JOcTree
using jInv.InverseSolve

include("Maxwell-derivative-test.jl")

# Create OcTree mesh
h  = [1/32,1/32,1/64]
n  = [128,64,32]
x0 = [-2.0,-1.0,0.0]
S  = createOcTreeFromBox(x0[1],x0[2],x0[3],n[1],n[2],n[3],h[1],h[2],h[3],0.0,0.0,0.0,0.0,0.0,1.0,2,2)
M  = getOcTreeMeshFV(S,h,x0=x0)

# Setup problem with random sources (hanging edge elimination is essential)
Curl    = getCurlMatrix(M)
Ne,     = getEdgeConstraints(M)
Nf,Qf   = getFaceConstraints(M)
Curl    = Qf * Curl * Ne
Sources = Ne * (Curl' * complex.(randn(size(Curl,1),2),randn(size(Curl,1),2)))
Obs     = Ne * (Curl' * complex.(sprandn(size(Curl,1),2,0.001),sprandn(size(Curl,1),2,0.001)))
freq    = 1e2
Ainv    = getMUMPSsolver([],1)
pFor    = getMaxwellFreqParam(M,Sources,Obs,freq,Ainv)

# Isotropic conductivity and susceptibility model
#sigma = exp.(randn(M.nc))
sigma = randn(M.nc)
pFor.sensitivityMethod = :Implicit
activeProps = ["sigmaCell","muCell"]
scaleFac = 100
n    = M.nc
chi  = scaleFac*rand(n)
mu   = mu0*(1+chi)
dsig = randn(n)
# dsig[sigma+dsig .< 1e-8] = 1e-8
dmu  = mu0*(scaleFac/2)*randn(n)
dmu[mu+dmu .< mu0] = 0.0

values = Dict{String,Vector{Float64}}("sigmaCell"=>exp.(sigma),"muCell"=>mu)
m0 = MaxwellFreqModel(values,activeProps)

# Call getData
D,pFor = getData(m0,pFor)

# Modfun to be passed to MisfitParam
function modFun(sig)
    n = M.nc
    values = Dict{String,Vector{Float64}}(
               "sigmaCell"=>exp.(sig[1:n]), "muCell"=>sig[n+1:end])
    activeProps = ["sigmaCell", "muCell"]
    m = MaxwellFreqModel(values,activeProps)
    dmVals = Dict{String,SparseMatrixCSC{Float64,Int64}}(
                  "sigmaCell"=>[spdiagm(exp.(sig[1:n])) spzeros(n,n)],
                  "muCell"=>[spzeros(n,n) speye(n)])
    dm = MaxwellFreqModelDerivative(dmVals,activeProps)
    return m,dm
end

#Global to local
P = MaxwellFreqModelTransform(speye(M.nc),spzeros(0,0))
mref = MaxwellFreqModel()
gloc = GlobalToLocal(P,mref)

# Construct MisfitParam
pMis = MisfitParam(pFor,ones(Complex{Float64},length(D)),D,SSDFun,modFun,gloc)

function f(mdum)
    d, = computeMisfit(mdum,pMis,false,true)
    return d
end

function df(x,sig)
    sigma,dsigma = pMis.modelfun(sig)

    sigmaloc = interpGlobalToLocal(sigma,pMis.gloc.PForInv,pMis.gloc.sigmaBackground)
    xloc     = interpGlobalToLocal(dsigma*x,pMis.gloc.PForInv)
    Jx       = vec(getSensMatVec(xloc,sigmaloc,pMis.pFor))
end


# Derivative test
checkDerivativeMaxwellOcTreeFV,error,order =
    checkDerivativeMax(f,df,[sigma;mu],v=[dsig;dmu],base=2.0)
@test checkDerivativeMaxwellOcTreeFV

# Adjoint test
x   = randn(2*M.nc)
u   = complex.(randn(size(D)),randn(size(D)))
v   = df(x,[sigma;mu])

m,dm = modFun([sigma;mu])
w    = dm'*interpLocalToGlobal(
       getSensTMatVec(vec(u),m,pMis.pFor),
       pMis.gloc.PForInv)
uv   = real(dot(vec(u),vec(v)))
wx   = dot(w,x)
tol  = 1e-13*max(abs(uv),abs(wx))
@test uv ≈ wx atol=tol
