
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.LinearSolvers
using jInv.Utils

using JOcTree

using MaxwellFrequency

# Create a Mesh
nx = 128
ny = 128
nz = 128
hx, hy, hz = 10., 10., 10.
x0 = -hx*nx/2
y0 = -hy*ny/2
z0 = -hz*nz/2
xa, xb = -125., 125.
ya, yb = -125., 125.
za, zb = -100., 50.
nf = 2
nc = 2

S = createOcTreeFromBox(x0, y0, z0, nx, ny, nz, hx, hy, hz, xa, xb, ya, yb, za, zb, nf, nc)
mesh = getOcTreeMeshFV(S,[hx,hy,hz]; x0=[x0,y0,z0])
gcc = getCellCenteredGrid(mesh)
exportUBCOcTreeMesh("mesh.msh", mesh)

# Create the true model
model = fill(1e-8, mesh.nc)
model[gcc[:,3].<0.] = 1e-3
blkInd = (abs.(gcc[:,1]) .< 25.) .& (abs.(gcc[:,2]) .< 25.) .& (gcc[:,3].<-30.) .& (gcc[:,3].>-70.)
model[blkInd] = 0.1
exportUBCOcTreeModel("trueModel.con", mesh, model)

# Define the survey
frq = logspace(2,4,4)
txSp = 100.
xc = [x for x in -100.:txSp:100., y in -100.:txSp:100., z in [10.]][:]
yc = [y for x in -100.:txSp:100., y in -100.:txSp:100., z in [10.]][:]
zc = [z for x in -100.:txSp:100., y in -100.:txSp:100., z in [10.]][:]
txCP = [xc yc zc]

txLoop = 20*[-1. -1. 0.; 
              1. -1. 0.;
              1.  1. 0.;
             -1.  1. 0.;
             -1. -1. 0.]

rxLoop = 2*[-1. -1. 0.; 
             1. -1. 0.;
             1.  1. 0.;
            -1.  1. 0.;
            -1. -1. 0.]
            
Sources = spzeros(Complex128, sum(mesh.ne), size(txCP,1))
for i = 1:size(txCP,1)
    Sources[:,i] = complex(getEdgeIntegralOfPolygonalChain(mesh,txCP[i,:]'.+txLoop))
end

Receivers = spzeros(Complex128, sum(mesh.ne), size(txCP,1))
for i = 1:size(txCP,1)
    Receivers[:,i] = complex(getEdgeIntegralOfPolygonalChain(mesh,txCP[i,:]'.+rxLoop))
end
Obs = typeof(Receivers)[copy(Receivers) for i = 1:length(frq)]
for i = 1:size(txCP,1)
    for j = 1:length(frq)
        Obs[j][:,i] *= complex(0.0, 1.0 / (mu0 * frq[j]))
    end
end

solver = getMUMPSsolver()

# Setup pFor
nFreqs = length(frq)
pFor   = Array{RemoteChannel}(nFreqs)
workerList = workers()
nw         = length(workerList)
for i = 1:nFreqs
    fields = Array{Complex128}(0, 0)
    pFor[i] = initRemoteChannel(getMaxwellFreqParam,workerList[i%nw+1],
                                mesh,Sources,Obs[i],frq[i],solver)
end

# Generate synthetic data
dobs, pFor = getData(model, pFor)

# Setup initial & reference model
isactive = gcc[:,3].<0.
sigmaBackground = zeros(mesh.nc)
sigmaBackground[.~isactive] = 1e-8

Iact = speye(Bool,mesh.nc)
Iact = Iact[:,find(isactive)]

m0 = fill(log(1e-3), size(Iact,2))
mref = fill(log(1e-3), size(Iact,2))

boundsLow = fill(log(1e-5),size(Iact,2)) 
boundsHigh = fill(log(1e0),size(Iact,2))

# Set up misfit function
Dobs = [fetch(dobs[i]) for i in 1:length(frq)]
Wd = [0.1*abs.(Dobs[i])+1e-6 for i in 1:length(frq)]
Wd = Wd + im*Wd
pMis = getMisfitParam(pFor, Wd, Dobs, SSDFun, Iact, sigmaBackground)

# Setup regularization
regfun = wdiffusionReg
regparams = [sqrt(1.0), sqrt(1.0), sqrt(1.0), 5e-7]
regfunw(m,mreff,Mm) = wdiffusionReg(m,mreff,Mm,Iact=Iact,C=regparams)

# Setup inversion
modfun = expMod
beta = 1e-32
cgit = 20 
maxit = 3

pInv = getInverseParam(mesh, modfun, regfunw, beta, mref,
                       boundsLow, boundsHigh,
                       pcgMaxIter=cgit, maxIter=maxit)

# Run it
mc,Dc,flag = projGNCG(m0,pInv,pMis)

exportUBCOcTreeModel("inv.con", mesh, Iact*modfun(mc)[1] + sigmaBackground)
