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
include("Maxwell-derivative-test.jl")

omega = [0 12 0 4 0 3]
n     = [16;16;14]
Mr     = getRegularMesh(omega,n)

# setup random sources
Curl     = getCurlMatrix(Mr)
Sources  = map(Complex{Float64},Curl'*randn(sum(Mr.nf),5))
freq     = 1e2
Obs      = Curl'*sprandn(sum(Mr.nf),20,0.001)
Ainv     = getMUMPSsolver([],1)
paramr   = getMaxwellFreqParam(Mr,Sources,Obs,[],freq,Ainv)
paramrSE = getMaxwellFreqParam(Mr,Sources,Obs,[],freq,Ainv,sensitivityMethod=:Explicit)

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
df(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramr)
checkDerivativeMaxwellRegular,errorr,orderr = checkDerivativeMax(f,df,m;nSuccess=4)
@test checkDerivativeMaxwellRegular

# adjoint test
x = randn(length(m))
u = randn(100)
v = getSensMatVec(x,m,paramr)
w = getSensTMatVec(vec(u),m,paramr)
@test norm((vec(u))' * vec(v)) ≈ norm(w' * x) atol=1.0e-13

# Test explicit sensitivities
println("Testing regular mesh explicit sensitivities")
paramrSE  = getMaxwellFreqParam(Mr,Sources,Obs,[],freq,Ainv,sensitivityMethod=:Explicit)

# Derivative test
function f1(sigdum)
  d, = getData(sigdum,paramrSE)
  return d
end
df1(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramrSE)
checkDerivativeMaxwellRegularSE,errorrSE,orderrSE = checkDerivativeMax(f1,df1,m;nSuccess=4)
@test checkDerivativeMaxwellRegularSE

# Adjoint test
x = randn(length(m))
u = randn(size(Dr))
v = getSensMatVec(x,m,paramrSE)
w = getSensTMatVec(vec(u),m,paramrSE)
@test norm((vec(u))' * vec(v)) ≈ norm(w' * x) atol=1.0e-13

# create corresponding tensor mesh
h1 = Mr.h[1]*ones(n[1])
h2 = Mr.h[2]*ones(n[2])
h3 = Mr.h[3]*ones(n[3])
Mt = getTensorMesh3D(h1,h2,h3)

# call getData
paramt     = getMaxwellFreqParam(Mt,Sources,Obs,[],freq,Ainv)
Dt, paramt = getData(m,paramt);

#Test that regular mesh data agree with tensor mesh data
RE1 = norm(Dt-Dr)/norm(Dt)
@test RE1 < 1e-9

# Test implicit sensitivities
println("Testing tensor mesh implicit sensitivities")
#Derivative test
function f2(sigdum)
  d, = getData(sigdum,paramt)
  return d
end
  
df2(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramt)
checkDerivativeMaxwellTensor,errort,ordert = checkDerivativeMax(f2,df2,m)
@test checkDerivativeMaxwellTensor

# adjoint test
x = randn(length(m))
u = randn(size(Dt))
v = getSensMatVec(x,m,paramt)
w = getSensTMatVec(vec(u),m,paramt)
@test norm((vec(u))' * vec(v)) ≈ norm(w' * x) atol=1.0e-13

# Test explicit sensitivities
println("Testing tensor mesh explicit sensitivities")
paramtSE  = getMaxwellFreqParam(Mr,Sources,Obs,[],freq,Ainv,sensitivityMethod=:Explicit)
# Derivative test
function f3(sigdum)
  d, = getData(sigdum,paramtSE)
  return d
end
  
df3(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramtSE)
checkDerivativeMaxwellTensorSE,errortSE,ordertSE = checkDerivativeMax(f3,df3,m)
@test checkDerivativeMaxwellTensorSE

# Adjoint test
x = randn(length(m))
u = randn(size(Dt))
v = getSensMatVec(x,m,paramtSE)
w = getSensTMatVec(vec(u),m,paramtSE)
@test norm((vec(u))' * vec(v)) ≈ norm(w' * x) atol=1.0e-13

if hasJOcTree
  # Create corresponding OcTree mesh
  i,j,k = ndgrid(1:n[1],1:n[2],1:n[3])
  bsz   = ones(Int,prod(n))
  S     = sparse3(vec(i),vec(j),vec(k),bsz,n)
  h     = Mr.h
  Mofv  = getOcTreeMeshFV(S,h)

  # call getData
  paramofv       = getMaxwellFreqParam(Mofv,Sources,Obs,[],freq,Ainv)
  Dofv, paramofv = getData(m,paramofv)
  
  #Test that tensor data agree with OcTreeFV data
  RE2 = norm(Dt-Dofv)/norm(Dt)
  @test RE2 < 1e-2
  
  # Test implicit sensitivities
  println("Testing OcTree mesh FV implicit sensitivities")
  # Derivative test
  function f4(sigdum)
    d, = getData(sigdum,paramofv)
    return d
  end
  df4(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramofv)
  checkDerivativeMaxwellOcTreeFV,errorofv,orderofv = checkDerivativeMax(f4,df4,m)
  @test checkDerivativeMaxwellOcTreeFV
  
  # adjoint test
  x = randn(length(m))
  u = randn(size(Dofv))
  v = getSensMatVec(x,m,paramofv)
  w = getSensTMatVec(vec(u),m,paramofv)
  @test norm((vec(u))' * vec(v)) ≈ norm(w' * x) atol=1.0e-13
  
  #Test explicit sensitivities
  println("Testing OcTree mesh FV explicit sensitivities")
  paramofvSE  = getMaxwellFreqParam(Mofv,Sources,Obs,[],freq,Ainv,sensitivityMethod=:Explicit)
                                    
  #Derivative test
  function f(sigdum)
    d, = getData(sigdum,paramofvSE)
    return d
  end
  df(zdum,sigdum) = getSensMatVec(zdum,sigdum,paramofvSE)
  checkDerivativeMaxwellOcTreeFVSE,errorofvSE,orderofvSE = checkDerivativeMax(f,df,m)
  @test checkDerivativeMaxwellOcTreeFVSE
  
  # adjoint test
  x = randn(length(m))
  u = randn(size(Dofv))
  v = getSensMatVec(x,m,paramofvSE)
  w = getSensTMatVec(vec(u),m,paramofvSE)
  @test norm((vec(u))' * vec(v)) ≈ norm(w' * x) atol=1.0e-13
end
end #End of mesh consistency test set
