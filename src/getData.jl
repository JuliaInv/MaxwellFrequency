export getData

# = ========= The forward problem ===========================

function getData(sigma,   # conductivity
	              param::MaxwellFreqParam,
	              doClear::Bool=false)
	# Maxwell Forward problem
	
	doClearAll=false

	mu   = 4*pi*1e-7
	w    = param.freq
	S    = param.Sources
	P    = param.Obs
	
	Curl = getCurlMatrix(param.Mesh)
	
	Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
	Mmu  = getFaceMassMatrix(param.Mesh,fill(1/mu,length(sigma)))
  
  # eliminate hanging edges and faces
	Ne,   = getEdgeConstraints(param.Mesh)
  Nf,Qf = getFaceConstraints(param.Mesh)
  
  Curl = Qf  * Curl * Ne
  Msig = Ne' * Msig * Ne
  Mmu  = Nf' * Mmu  * Nf

	A   = Curl' * Mmu * Curl - (im * w) * Msig
	rhs = (im * w) * full(Ne' * S)
	
	param.Ainv.doClear = 1
	U, param.Ainv = solveMaxFreq(A, rhs, Msig,
	                             param.Mesh, w, param.Ainv,0)
	param.Ainv.doClear = 0
  
	U = Ne * U
	D = P' * U
	
	param.Fields = U
	if doClear
		# clear fields and factorization
		clear!(param,clearAll=doClearAll)
	end
	
	return D, param # data, MaxwellParam
end # function getData

if hasJOcTree
function getData(sigma,  # conductivity
	              param::MaxwellFreqParam{OcTreeMeshFEM},
	              doClear::Bool=false)
	# Maxwell Forward problem for FEM discretization on OcTreeMesh
	doClearAll=false

	mu   = 4*pi*1e-7
	w    = param.freq
	S    = param.Sources
	P    = param.Obs
	
	K,M,Msig = getMatricesFEM(param.Mesh,sigma)
	
	# get the null space matrix
	N, = getEdgeConstraints(param.Mesh)   
	A  = N'*(K/mu - im*w*Msig)*N
	
	param.Ainv.doClear = 1
	rhs = im*w*full(N'*S)
	U, param.Ainv = solveMaxFreq(A, rhs, Msig,
	                             param.Mesh, w, param.Ainv,0)

	param.Ainv.doClear = 0
	U = N*U
	D = P'*U	
	param.Fields =U
	if doClear
		# clear fields and factorization
		clear!(param,clearAll=doClearAll)
	end
	
	return D, param  # data, MaxwellParam
end # function getData
end


function getData(sigma::Array{Float64,1},param::MaxwellFreqParamSE,doClear::Bool=false)
	tempParam = getMaxwellFreqParam(param.Mesh,param.Sources,param.Obs,[],param.freq, param.Ainv)
	param.Sens=[]
	D,tempParam = getData(sigma,tempParam)
	ns = size(param.Sources,2)
	nr = size(param.Obs,2)	
	v      = zeros(ns*nr)
	J    = getSensTMatVec(speye(ns*nr),sigma,tempParam)
	param.Sens = J'
	clear!(tempParam)	
	return D, param
end