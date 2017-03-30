# get inversion parameters
include("parametersForInversion.jl")

# ------- SETUP PARAMETERS FOR THE MODEL AND DATA

trx, h, itopo = setupMeshParam(datafile, topofile, n,x0,meshL )

#trx = trx[1:3]


# ------- Preparing inverse mesh -----------------------------------------

println("Prepare Inverse Mesh")
tic()

#Minv = setupBigOctreeMesh(trx,h,n,x0,itopo,depth_core_inv,mincellfactor);
Minv = setupBigOctreeMeshPolygon(trx,h,n,x0,itopo,depth_core_inv,mincellfactor);

toc()
println("Inverse mesh has ", Minv.nc, " cells")

exportUBCOcTreeMesh("meshInv.txt",Minv)

# ----- Prepare the model to be run -------------------------------------
# get conductivity model
println("Generating initial model")
tic()


sigma, sigmaBck, isactive = getInitialmodel(Minv,itopo, halfSpaceCond, backCond)

Iact    = speye(Bool,Minv.nc)
Iact    = Iact[:,find(isactive)]
IactBck = speye(Bool,Minv.nc)
IactBck = IactBck[:,find(!isactive)]
toc()

sigmaBackground = IactBck * sigmaBck




############################### Setup to remote forward and misfit param  ######################################

#------------ Preparing forward mesh -------------------------------------
println("Generating small forward modelling meshes in parallel")
tic()

#getMUMPSsolver(Ainv=[],doClear=1,ooc=0,sym=0)
#linSolParam = getMUMPSsolver([],1,0,2)
linSolParam = getIterativeSolver(KrylovMethods.bicgstb)
  
meshingParam    = Array{Any}(5)
meshingParam[1] = nsmallcells 
meshingParam[2] = mincellsize
meshingParam[3] = itopo
meshingParam[4] = depth_core
meshingParam[5] = mincellfactor
  
pFor = getMaxwellFreqParam(x0,n,h,trx,itopo,meshingParam,linSolParam,fname,doFV,doSE)


toc()
println("Prepare interpolation matrices")
tic()

if isa(ninterp, Integer)
	Mesh2Mesh  = prepareMesh2MeshOT(pFor,Minv,ninterp,compact)
	#Mesh2Mesh  = prepareMesh2Mesh(pFor,Minv,compact)
	
else
	Mesh2Mesh  = prepareMesh2Mesh(pFor,Minv,compact)	
end
toc()

# ----- Compute data!!! ------------------------------------------------
# println("Compute Data")
# tic()
# 	sg = Iact*sigma + sigmaBackground
# 	Dcomp,pFor = getData(sg,pFor,Mesh2Mesh,true)
# toc()
#
#
# # ----- Prepare for inversion ------------------------------------------------
# #
# #
# Dobs = Array(Any,length(pFor))
# for k=1:length(pFor)
# 	Dobs[k] = fetch(Dcomp[k])
# end
#
# Wd   = Array(Array{Complex128},length(pFor))
# for i=1:length(pFor)
# 	Wd[i]       = complex(ones(size(Dobs[i])))./(mean(abs(real(Dobs[i])))/100) +
# 			1im * complex(ones(size(Dobs[i])))./(mean(abs(imag(Dobs[i])))/100);
# end

Wd    = Array(Array{Complex128},length(pFor))
Dobs  = Array(Array{Complex128},length(pFor))
for i=1:length(trx)
   
   if typeof(trx[i]) == MaxwellUtils.Transmitter
      dobsi, Wdi = getDobsWdFromTrx(trx[i])
   elseif typeof(trx[i]) == MaxwellUtils.TransmitterOmega
      dobsi = trx[i].Dobs
      Wdi   = trx[i].Wd
   else
      error("bad trx type")
   end

   Wd[i] = Wdi; Dobs[i] = dobsi
end  # i



# Write misfit parameters to remote workers. The following function clears the contents of Mesh2Mesh and pFor!
pMisRF = getMisfitParam(pFor, Wd, Dobs, misfun,Iact,sigmaBackground,Mesh2Mesh)



############################### Setup to inverse param  ######################################
println("Setup Inverse Param")



mref = fill(log(halfSpaceCond), size(Iact,2))

boundsLow  = fill(log(BL),size(Iact,2)) 
boundsHigh = fill(log(BH),size(Iact,2))    


if regtype == wdiffusionReg
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
      Weights = getInterfaceWeights( Minv, itopo, surfweight, regparams[1:3] )
      regparams = vcat( regparams[4], Weights );
      regfun(m, mref, M ) = wdiffusionReg(m, mref, M , Iact = Iact,C = regparams)
# 	  regfun2(m, mref, M , Iact) = wdiffusionReg(m, mref, M , Iact,regparams)
# 	  regfun = regfun2;
   else
      println("No interface weights.")
   end
end  # regfun == wdiffusionReg

							  
pInv = getInverseParam(Minv,modfun,regfun,beta,mref,boundsLow,boundsHigh,
                       maxStep=1.0,pcgMaxIter=cgit, pcgTol =0.1 ,minUpdate=1e-4,maxIter=maxit,
					   HesPrec=getSSORRegularizationPreconditioner(1.0,1e-15,10))

## only fwd model
# print("=======  Start forward model =========\n")
# m = log(sigma)
# tic()
# Dc,F,dF,d2F,pFor = computeMisfit(m,pMisRF,false)
# toc()

print("=======  Start inversion =========\n")
tic()
m0 = fill(log(halfSpaceCond), size(Iact,2))

sig,dsig = pInv.modelfun(m0)
Dc,F,dF,d2F,pMisRF,tMis = computeMisfit(sig,pMisRF,true)
mc,Dc,flag = projGNCG(m0,pInv,pMisRF);
toc()
