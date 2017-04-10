import jInv.ForwardShare.getData
export getData

# = ========= The forward problem ===========================

function getData(sigma,   # conductivity
	              param::MaxwellFreqParam,
	              doClear::Bool=false)
	# Maxwell Forward problem
	
	doClearAll=false

	w    = param.freq
	S    = param.Sources
	P    = param.Obs
	
	A = getMaxwellFreqMatrix(sigma, param)

	Ne,   = getEdgeConstraints(param.Mesh)
	Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
  	Msig = Ne' * Msig * Ne
  
  	rhs = (im * w) * (Ne' * S)
	
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