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


