using MaxwellFrequency
using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
using MUMPS
using KrylovMethods

omega = [0 200 0 200 0 200]
n = [16; 20; 32]
mesh = getRegularMesh(omega,n)

sigma = rand(mesh.nc)*1e-1
mu0 = 4*pi*1e-7
w = 1.0*2*pi

# get Maxwell system
Curl = getCurlMatrix(mesh)
Msig = getEdgeMassMatrix(mesh,vec(sigma))
Mmu  = getFaceMassMatrix(mesh,vec(zeros(mesh.nc).+1/mu0))
A    = Curl'*Mmu*Curl + 1im*w*Msig

rhs = randn(size(A,1)) + 1im*randn(size(A,1))
rhs = Curl'*(Curl*rhs)

@testset "Solver tests" begin

bicgstbIM(A,b;M=[],tol=1e-1,maxIter=10,out=-1) = bicgstb(A,b,M1=M,tol=tol,maxIter=maxIter,out=out)

# get solver types
solvers  = [ getJuliaSolver(), 
             getMUMPSsolver(), 
             getIterativeSolver(bicgstbIM, maxIter=1000, tol=1e-7) ]

# solvers  = [ getIterativeSolver(bicgstbIM, maxIter=1000, tol=1e-7) ]

for s = solvers
	println("   testing Maxwell solver for ", typeof(s), "...")
	tic()
	xi, = solveMaxFreq(A,rhs,Msig,mesh,w,s)
	solveTime = toq()
	@test norm(A*xi-rhs)/norm(rhs) < 1e-6
	print("   passed in ",solveTime ," sec!\n")
end

end