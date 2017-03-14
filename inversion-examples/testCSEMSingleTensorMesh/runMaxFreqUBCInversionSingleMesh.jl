#Load parameters
include("parametersForInversion.jl")

# ------- Read data ---------------------------------------

if length(datafile) != 4
	error("Expect data in new UBC data format")
end
only_loc  = false
datainput = readDataFile(datafile[1], only_loc)
src       = readTrxRcvFile(datafile[2])
rcv       = readTrxRcvFile(datafile[3])
frq       = readFrqFile(datafile[4])

# extract trx, rcv, frq that actually used and convert to simple arrays
i    = Int64[datainput[i].trx_idx for i = 1:length(datainput)]
j    = sort(unique(i))
src  = src[j]
itrx = indexin(i,j)

i    = Int64[datainput[i].rcv_idx for i = 1:length(datainput)]
j    = sort(unique(i))
rcv  = rcv[j]
ircv = indexin(i,j)

i    = Int64[datainput[i].frq_idx for i = 1:length(datainput)]
j    = sort(unique(i))
frq  = frq[j]
ifrq = indexin(i,j)

# copy data and weights to dense 2D arrays, one for each frequency
Dobs = Array{Complex128,2}[zeros(Complex128, length(rcv), length(src)) for i = 1:length(frq)]
Wd   = Array{Complex128,2}[zeros(Complex128, length(rcv), length(src)) for i = 1:length(frq)]
for i = 1:length(datainput)
	Dobs[ifrq[i]][ircv[i],itrx[i]] = complex(
		datainput[i].dobs[1],
		datainput[i].dobs[2])
	Wd[ifrq[i]][ircv[i],itrx[i]]   = complex(
		datainput[i].sd[1] == 0.0 ? 0.0 : 1.0 / datainput[i].sd[1],
		datainput[i].sd[2] == 0.0 ? 0.0 : 1.0 / datainput[i].sd[2])
end

src = Array{Float64,2}[src[i].trxpts' for i = 1:length(src)]
rcv = Array{Float64,2}[rcv[i].trxpts' for i = 1:length(rcv)]
frq = Float64[frq[i].omega for i = 1:length(frq)]

x1 =  Inf
x2 = -Inf
y1 =  Inf
y2 = -Inf
z1 =  Inf
z2 = -Inf
for i = 1:length(src)
	x1 = min(x1, minimum(src[i][:,1]))
	x2 = max(x2, maximum(src[i][:,1]))
	y1 = min(y1, minimum(src[i][:,2]))
	y2 = max(y2, maximum(src[i][:,2]))
	z1 = min(z1, minimum(src[i][:,3]))
	z2 = max(z2, maximum(src[i][:,3]))
end
for i = 1:length(rcv)
	x1 = min(x1, minimum(rcv[i][:,1]))
	x2 = max(x2, maximum(rcv[i][:,1]))
	y1 = min(y1, minimum(rcv[i][:,2]))
	y2 = max(y2, maximum(rcv[i][:,2]))
	z1 = min(z1, minimum(rcv[i][:,3]))
	z2 = max(z2, maximum(rcv[i][:,3]))
end

# ------- Generate mesh --------------------------------------------------

println("Generating mesh")

z1 = -1041.79 - 150
z2 = -1041.79 + 150

ncore = round(Int64,[x2-x1;y2-y1;z2-z1]./h0)

h1 = [h0[1]*(expFac.^collect(nPadxy:-1:1))
      h0[1]*ones(ncore[1])
      h0[1]*(expFac.^collect(1:nPadxy))]

h2 = [h0[2]*(expFac.^collect(nPadxy:-1:1))
      h0[2]*ones(ncore[2])
      h0[2]*(expFac.^collect(1:nPadxy))]

h3 = [h0[3]*(expFac.^collect(nPadz:-1:1))
      h0[3]*ones(ncore[3])
      h0[3]*(expFac.^collect(1:nPadz))]

x0    = zeros(3)
x0[1] = x1-sum(h1[1:nPadxy])
x0[2] = y1-sum(h2[1:nPadxy])
x0[3] = z1-sum(h3[1:nPadz])

M = getTensorMesh3D(h1,h2,h3,x0)
display(M)

meshL = [sum(h1);sum(h2);sum(h3)]
xn    = x0 + meshL  # opposite corner
if x1 < x0[1] || x2 > xn[1] ||
   y1 < x0[2] || y2 > xn[2] ||   
   z1 < x0[3] || z2 > xn[3]
   error("Data outside of mesh.")
end

# # ----- Generate initial model -------------------------------------------
# 
println("Generating initial model")

#Mean elevation of topo is -1041.79
zSurf   = -1041.79
Xc      = getCellCenteredGrid(M)
IBck    = find(Xc[:,3] .>  zSurf)
IGrnd   = find(Xc[:,3] .<= zSurf)
Iact    = speye(Bool,M.nc)
Iact    = Iact[:,IGrnd]
IactBck = speye(Bool,M.nc)
IactBck = IactBck[:,IBck]

# model a half space
sigma    = ones(length(IGrnd))*halfSpaceCond
sigmaBck = ones(length(IBck))*backCond

sigmaBackground = IactBck * sigmaBck
 
#------------ Set up forward problem -------------------------------------

println("Setting up forward problem")
tic()

# transmitters
Sources = zeros(Complex128, sum(M.ne), length(src))
for i = 1:length(src)
        ply          = [src[i][:,1:2] zSurf*ones(size(src[i],1))]
	Sources[:,i] = complex(getEdgeIntegralOfPolygonalChain(M,ply,normalize=false))
end

# receivers
Receivers = spzeros(Complex128, sum(M.ne), length(rcv))
for i = 1:length(rcv)
        ply            = [rcv[i][:,1:2] zSurf*ones(size(rcv[i],1))]
        if sum(ply[1,:]-ply[2,:]) == 0.0
          ply[1,3] = zSurf-10
          ply[2,3] = zSurf+10
        end
	Receivers[:,i] = complex(getEdgeIntegralOfPolygonalChain(M,ply,normalize=true))
end
# for closed loops, scale by i / (omega * mu0)
Obs = typeof(Receivers)[copy(Receivers) for i = 1:length(frq)]
for i = 1:length(rcv)
	if norm(rcv[i][1,:] - rcv[i][end,:]) < 1e-16
		for j = 1:length(frq)
			Obs[j][:,i] *= complex(0.0, 1.0 / (mu0 * frq[j]))
		end
	end
end

#linear solver
linSolParam = getMUMPSsolver([],1,0,2)

# create forward solver for each frequency
nFreqs = length(frq)
pFor   = Array(RemoteChannel,nFreqs)
workerList = workers()
nw         = length(workerList)
for i = 1:nFreqs
  if doSE
    pFor[i] = initRemoteChannel(getMaxwellFreqParamSE,workerList[i%nw+1],
                                    M,Sources,Obs[i],fields,frq[i],linSolParam)
  else
    fields = Array(Complex128, 0, 0)
    pFor[i] = initRemoteChannel(getMaxwellFreqParam,workerList[i%nw+1],
                                M,Sources,Obs[i],fields,frq[i],linSolParam)
  end
end
toc()

# setup inverse param
println("Setup Inverse Param")
mref = fill(log(halfSpaceCond), size(Iact,2))

boundsLow  = fill(log(BL),size(Iact,2)) 
boundsHigh = fill(log(BH),size(Iact,2))    

pMisRF = getMisfitParam(pFor, Wd, Dobs, misfun,Iact,sigmaBackground)

regfunw(m,mreff,Mm) = wdiffusionReg(m,mreff,Mm,Iact=Iact,C=regparams)

pInv = getInverseParam(M,modfun,regfunw,beta,mref,
                       boundsLow,boundsHigh,
                       pcgMaxIter=cgit,maxIter=maxit)

print("=======  Start inversion =========\n")
tic()
m0 = fill(log(halfSpaceCond), size(Iact,2))
mc,Dc,flag = projGNCG(m0,pInv,pMisRF)
toc()

