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
	J = getSensMat(sigma, tempParam)
	param.Sens = J'
	clear!(tempParam)	
	return D, param
end

function getSensMat(sigma, param)
    x = eye(size(param.Sources,2)*size(param.Obs,2))
    U    = param.Fields
    w    = param.freq
    P    = param.Obs
    
    Ne, = getEdgeConstraints(param.Mesh)
    Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
    Msig = Ne' * Msig * Ne
    A = getMaxwellFreqMatrix(sigma, param)

    X    = reshape(complex(x),size(P,2),size(U,2),size(x,2))
    matv = zeros(Complex128,length(sigma),size(x,2))    

    for i=1:size(U,2)
        Z     = -Ne'*(P*X[:,i,:])
        Z,    = solveMaxFreq(A,Z,Msig,param.Mesh,w,param.Ainv,1)
        u     = U[:,i]
        dAdm  = getdEdgeMassMatrix(param.Mesh,u)
        dAdm  = -im*w*Ne'*dAdm
        matv +=  (dAdm'*Z)
    end
    
    return matv
end
