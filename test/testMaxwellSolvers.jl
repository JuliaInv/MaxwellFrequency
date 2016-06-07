using MaxwellFrequency
using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
using MUMPS
using KrylovMethods

omega = [0 200 0 200 0 200]
n     = [16; 20; 32]
Mr     = getRegularMesh(omega,n)

sigma = rand(Mr.nc)*1e-1
mu    = 4*pi*1e-7
w     = 1.0*2*pi

# get Maxwell system
Curl = getCurlMatrix(Mr)
Msig = getEdgeMassMatrix(Mr,vec(sigma))
Mmu  = getFaceMassMatrix(Mr,vec(zeros(size(sigma)).+1/mu))
A    = Curl'*Mmu*Curl + 1im*w*Msig

rhs = randn(size(A,1)) + 1im*randn(size(A,1))
#rhs = Curl'*(Curl*rhs)

@testset "Solver tests" begin

# get solver types
solvers  = ( getMUMPSsolver(),getIterativeSolver(:bicgstb))
solvers[2].maxIter=1000
xi = []
for s = solvers
	print("   testing Maxwell solver for ", typeof(s), "...")
	tic()
	xi, = solveMaxFreq(A,rhs,Msig,Mr,w,s)
	solveTime = toq()
	@test norm(A*xi-rhs)/norm(rhs) < 1e-7
	print("passed in ",solveTime ," sec!\n")
end
end #End of solver test set