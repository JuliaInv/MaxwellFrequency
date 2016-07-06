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
xn = x0 + meshL  # opposite corner
if x1 < x0[1] || x2 > xn[1] ||
   y1 < x0[2] || y2 > xn[2] ||   
   z1 < x0[3] || z2 > xn[3]
   error("Data outside of mesh.")
end

# ------- Read topography -----------

# cell size
# h = meshL ./ n
# 
# topogrid = readTopo( topofile, n, x0, h )
# 
# t1 = minimum(topogrid)
# t2 = maximum(topogrid)
# 
# if t1 < x0[3] || t2 > xn[3]
#    error("Topography outside of mesh.")
# end   
# 
# # Figure out the number of surface cells for each point in
# # the x,y grid.
# itopo = getItopo(h,n,x0, topogrid)

# ------- Generate mesh --------------------------------------------------

println("Generating mesh")
tic()

z1 = min(t1, z1)
z2 = max(t2, z2)
# nf = 2
# nc = 2

# S = createOcTreeFromBox(
# 	x0[1], x0[2], x0[3],
# 	n[1], n[2], n[3],
# 	h[1], h[2], h[3],
# 	x1, x2, y1, y2, z1, z2,
# 	nf, nc)
# if doFV
# 	M = getOcTreeMeshFV(S, h; x0 = x0)
# else
# 	M = getOcTreeMeshFEM(S, h; x0 = x0)
# end
h1    = [x1-h0[1]*(expFac.^collect(nPadxy:-1:1))
         collect(x1:h0[1]:x2)
         x2+h0[1]*(expFac.^collect(1:nPadxy))

h2    = [y1-h0[2]*(expFac.^collect(nPadxy:-1:1))
         collect(y1:h0[2]:y2)
         y2+h0[2]*(expFac.^collect(1:nPadxy))

h1    = [z1-h0[3]*(expFac.^collect(nPadz:-1:1))
         collect(z1:h0[3]:z2)
         z2+h0[3]*(expFac.^collect(1:nPadz))

M = getTensorMesh3D(h1,h2,h3,x0)

toc()

display(M)

# exportOcTreeMeshRoman("mesh.txt",M)
# 
# # ----- Generate initial model -------------------------------------------
# 
# println("Generating initial model")
# tic()
# 
# sigma, sigmaBck, isactive = getInitialmodel(M, itopo, halfSpaceCond, backCond)
# 
# Iact    = speye(Bool,M.nc)
# Iact    = Iact[:,find(isactive)]
# IactBck = speye(Bool,M.nc)
# IactBck = IactBck[:,find(!isactive)]
# 

#Mean elevation of topo is -1041.79

# sigmaBackground = IactBck * sigmaBck
# 
# toc()
# 
#------------ Set up forward problem -------------------------------------

println("Setting up forward problem")
tic()

nEx,nEy,nEz = getEdgeNumbering(S)

# transmitters
Sources = spzeros(Complex128, sum(M.ne), length(src))
for i = 1:length(src)
	Sources[:,i] = complex(getEdgeIntegralOfPolygonalChain(M,src[i],nEx,nEy,nEz,normalize=false))
end

# receivers
Receivers = spzeros(Complex128, sum(M.ne), length(rcv))
for i = 1:length(rcv)
	Receivers[:,i] = complex(getEdgeIntegralOfPolygonalChain(M,rcv[i],nEx,nEy,nEz,normalize=true))
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

# linear solver
linSolParam = getMUMPSsolver([],1,0,2)

# create forward solver for each frequency
nFreqs = length(frq)
pFor   = Array(RemoteRef{Channel{Any}},nFreqs)
for i = 1:nFreqs
  if doSE
  	pFor[i] = @spawn getMaxwellFreqParamSE(M,Sources,Obs[i],fields,frq[i],linSolParam)
  else
  	fields = Array(Complex128, 0, 0)
  	pFor[i] = @spawn getMaxwellFreqParam(M,Sources,Obs[i],fields,frq[i],linSolParam)
  end
end
toc()

# setup inverse param
println("Setup Inverse Param")
mref = fill(log(halfSpaceCond), size(Iact,2))

boundsLow  = fill(log(BL),size(Iact,2)) 
boundsHigh = fill(log(BH),size(Iact,2))    

# pMod = prepareGlobalToLocal(speye(M.nc),Iact,sigmaBackground,"")
# pMod = GlobalToLocal[pMod for i = 1:length(pFor)]

pMisRF = getMisfitParam(pFor, Wd, Dobs, misfun,Iact,sigmaBackground)


if regfun == wdiffusionReg
#   surfweight =
#      try
#         surfweight
#      catch
#         [1.0]  # default value
#      end

   

   if !isdefined(:surfweight)
      surfweight = [1.0]
   end

   if length(surfweight) >= 1 && any( surfweight .!= 1 )
      Weights = getInterfaceWeights( M, itopo, surfweight, regparams[1:3] )
      regparams = vcat( regparams[4], Weights )
   else
      println("No interface weights.")
   end
end  # regfun == wdiffusionReg
regfunw(m,mreff,Mm) = wdiffusionReg(m,mreff,Mm,Iact=Iact,C=regparams)

pInv = getInverseParam(M,modfun,regfunw,beta,mref,
                       boundsLow,boundsHigh,
                       pcgMaxIter=cgit,maxIter=maxit)

## only fwd model
# print("=======  Start forward model =========\n")
# m = log(sigma)
# tic()
# Dc,F,dF,d2F,pFor = computeMisfit(m,pInv.modelfun,pInv.model,pInv.misfit,pFor,Dobs,Wd,false)
# toc()

print("=======  Start inversion =========\n")
tic()
m0 = fill(log(halfSpaceCond), size(Iact,2))
mc,Dc,flag = projGNCG(m0,pInv,pMisRF)
toc()

