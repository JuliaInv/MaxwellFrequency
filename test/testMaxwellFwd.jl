using MaxwellFrequency
using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
using KrylovMethods
hasJOcTree = false
try
  using JOcTree
  hasJOcTree = true
catch
  hasJOcTree = false
end

omega = [0 12 0 4 0 3]
n     = [32;32;16]
Mr     = getRegularMesh(omega,n)

# setup random sources
Curl     = getCurlMatrix(Mr)
Sources  = map(Complex{Float64},Curl'*randn(sum(Mr.nf),5))
freq     = 1e2
Obs      = Curl'*sprandn(sum(Mr.nf),20,0.001)
fields   = []
Ainv     = getMUMPSsolver([],1)
#Ainv    = getIterativeSolver(100,1e-7,[],[],1)
paramr   = getMaxwellFreqParam(Mr,Sources,Obs,fields,freq,Ainv)
paramrSE = getMaxwellFreqParamSE(Mr,Sources,Obs,freq,Ainv)

# conductivitie
p = rand(Mr.nc)
m = p

# call getData
Dr, paramr = getData(m,paramr);

@testset "Mesh consistency tests" begin

# Test implicit sensitivities
# derivative test
println("Testing regular mesh implicit sensitivities")
function f(sigdum)
  d, = getData(sigdum,paramr)
  return d
end
  
# df(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramr)
# checkDerivativeMaxwellRegular,errorr,orderr = checkDerivative(f,df,m;nSuccess=4)
# @test checkDerivativeMaxwellRegular
# 
# # adjoint test
# x = randn(length(m))
# u = randn(100)
# v = getSensMatVec(x,m,paramr)
# w = getSensTMatVec(vec(u),m,paramr)
# @test_approx_eq_eps norm(vec(u)'*vec(v))  norm(w'*x) 1e-13
# 
# # #Test explicit sensitivities
# function testFunSE(sigma,param,paramSE,v=[])
# 	if isempty(paramSE.Sens)
# 		dc,paramSE = getData(sigma,paramSE)
# 	else
# 		dc,param = getData(sigma,param)
# 	end
# 	if !isempty(v)
# 		dDv = getSensMatVec(vec(v),sigma,paramSE)
# 		return vec(dc),vec(dDv)
# 	else
# 		return vec(dc)
# 	end
# end
# println("Testing regular mesh explicit sensitivities")
# paramrSE  = getMaxwellFreqParamSE(Mr,Sources,Obs,freq,Ainv)
# #Derivative test
# f(x,v=[]) = testFunSE(x,paramr,paramrSE,v)
# checkDerivativeMaxwellRegularSE,errorrSE,orderrSE = checkDerivative(f,m)
# @test checkDerivativeMaxwellRegularSE
# 
# #Adjoint test
# x = randn(length(m))
# u = randn(size(Dr))
# v = getSensMatVec(x,m,paramrSE)
# w = getSensTMatVec(vec(u),m,paramrSE)
# @test_approx_eq_eps norm(vec(u)'*vec(v))  norm(w'*x) 1e-13

# create corresponding tensor mesh
h1 = Mr.h[1]*ones(n[1])
h2 = Mr.h[2]*ones(n[2])
h3 = Mr.h[3]*ones(n[3])
Mt = getTensorMesh3D(h1,h2,h3)

# call getData
paramt     = getMaxwellFreqParam(Mt,Sources,Obs,fields,freq,Ainv)
Dt, paramt = getData(m,paramt);

#Test that regular mesh data agree with tensor mesh data
RE1 = norm(Dt-Dr)/norm(Dt)
@test RE1 < 1e-9

#Test implicit sensitivities
# println("Testing tensor mesh implicit sensitivities")
# #Derivative test
# function f(sigdum)
#   d, = getData(sigdum,paramt)
#   return d
# end
#   
# df(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramt)
# checkDerivativeMaxwellTensor,errort,ordert = checkDerivative(f,df,m)
# @test checkDerivativeMaxwellTensor
# 
# # adjoint test
# x = randn(length(m))
# u = randn(size(Dt))
# v = getSensMatVec(x,m,paramt)
# w = getSensTMatVec(vec(u),m,paramt)
# @test_approx_eq_eps norm(vec(u)'*vec(v)) norm(w'*x) 1e-13
# 
# #Test explicit sensitivities
# # #Test explicit sensitivities
# println("Testing tensor mesh explicit sensitivities")
# paramtSE  = getMaxwellFreqParamSE(Mt,Sources,Obs,freq,Ainv)
# #Derivative test
# f(x,v=[]) = testFunSE(x,paramt,paramtSE,v)
# checkDerivativeMaxwellTensorSE,errortSE,ordertSE = checkDerivative(f,m)
# @test checkDerivativeMaxwellTensorSE
# 
# #Adjoint test
# x = randn(length(m))
# u = randn(size(Dt))
# v = getSensMatVec(x,m,paramtSE)
# w = getSensTMatVec(vec(u),m,paramtSE)
# @test_approx_eq_eps norm(vec(u)'*vec(v)) norm(w'*x) 1e-13

if hasJOcTree
  #Create corresponding OcTree mesh
  i,j,k = ndgrid(1:n[1],1:n[2],1:n[3])
  bsz   = ones(Int,prod(n))
  S     = sparse3(vec(i),vec(j),vec(k),bsz,n)
  h     = Mr.h
  Mofv  = getOcTreeMeshFV(S,h)

  # call getData
  paramofv       = getMaxwellFreqParam(Mofv,Sources,Obs,fields,freq,Ainv)
  Dofv, paramofv = getData(m,paramofv)
  
  #Test that tensor data agree with OcTreeFV data
  RE2 = norm(Dt-Dofv)/norm(Dt)
  @test RE2 < 1e-2
  
#   #Test implicit sensitivities
#   println("Testing OcTree mesh FV implicit sensitivities")
#   #Derivative test
#   function f(sigdum)
#   d, = getData(sigdum,paramofv)
#   return d
# end
#   
# df(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramofv)
#   checkDerivativeMaxwellOcTreeFV,errorofv,orderofv = checkDerivative(f,df,m)
#   @test checkDerivativeMaxwellOcTreeFV
#   
#   # adjoint test
#   x = randn(length(m))
#   u = randn(size(Dofv))
#   v = getSensMatVec(x,m,paramofv)
#   w = getSensTMatVec(vec(u),m,paramofv)
#   @test_approx_eq_eps norm(vec(u)'*vec(v)) norm(w'*x) 1e-13
#   
#   #Test explicit sensitivities
#   # #Test explicit sensitivities
#   println("Testing OcTree mesh FV explicit sensitivities")
#   paramofvSE  = getMaxwellFreqParamSE(Mofv,Sources,Obs,freq,Ainv)
#   #Derivative test
#   f(x,v=[]) = testFunSE(x,paramofv,paramofvSE,v)
#   checkDerivativeMaxwellOcTreeFVSE,errorofvSE,orderofvSE = checkDerivative(f,m)
#   @test checkDerivativeMaxwellOcTreeFVSE
#   
#   # adjoint test
#   x = randn(length(m))
#   u = randn(size(Dofv))
#   v = getSensMatVec(x,m,paramofvSE)
#   w = getSensTMatVec(vec(u),m,paramofvSE)
#   @test_approx_eq_eps norm(vec(u)'*vec(v)) norm(w'*x) 1e-13
end
end #End of mesh consistency test set

warn("MaxwellFrequency: testMaxwellFwd: OcTree mesh FEM not currently tested")
# # create corresponding OcTree mesh FEM
# Mofem  = getOcTreeMeshFEM(S,h)
# 
# # call getData
# paramofem        = getMaxwellFreqParam(Mofem,Sources,Obs,fields,freq,Ainv)
# Dofem, paramofem = getData(m,paramofem)
# RE3 = norm(Dofv-Dofem)/norm(Dofv)
# @test RE3 < 1e-2
# 
# #Test implicit sensitivities
# println("Testing OcTree mesh FEM implicit sensitivities")
# #Derivative test
# f(x,v=[]) = testFun(x,paramofem,v)
# checkDerivativeMaxwellOcTreeFEM,errorofem,orderofem = checkDerivative(f,m)
# @test checkDerivativeMaxwellOcTreeFEM
# 
# # adjoint test
# x = randn(length(m))
# u = randn(size(Dofem))
# v = getSensMatVec(x,m,paramofem)
# w = getSensTMatVec(vec(u),m,paramofem)
# @test norm(vec(u)'*vec(v) - w'*x)<1e-4
# 
# #Test explicit sensitivities
# # #Test explicit sensitivities
# println("Testing OcTree mesh FEM explicit sensitivities")
# paramofemSE  = getMaxwellFreqParamSE(Mofem,Sources,Obs,freq,Ainv)
# #Derivative test
# f(x,v=[]) = testFunSE(x,paramofem,paramofemSE,v)
# checkDerivativeMaxwellOcTreeFEMSE,errorofemSE,orderofemSE = checkDerivative(f,m)
# @test checkDerivativeMaxwellOcTreeFEMSE
# 
# # adjoint test
# x = randn(length(m))
# u = randn(size(Dofem))
# v = getSensMatVec(x,m,paramofemSE)
# w = getSensTMatVec(vec(u),m,paramofemSE)
# @test norm(vec(u)'*vec(v) - w'*x)<1e-4
