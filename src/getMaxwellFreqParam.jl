export getMaxwellFreqParam, getMaxwellFreqParamSE

using MaxwellUtils

#  prepare the param where all the fields are known
function getMaxwellFreqParam(M::AbstractMesh, Sources, Obs, fields, freq, linSolParam::AbstractSolver; fname="")
    
    # Construct a new MaxwellParam object.
    if isempty(fields)
        fields = Array(Complex128,0,0)
    end

    return MaxwellFreqParam(M, Sources, Obs, fields, freq, linSolParam, fname)
end


function getMaxwellFreqParamSE(M::AbstractMesh, Sources, Obs, freq, linSolParam::AbstractSolver; fname="")
	
    # Construct a new MaxwellParam object.
	Sens = Array(Complex128,0,0)
	
    return MaxwellFreqParamSE(M, Sources, Obs, freq, linSolParam, Sens, fname)
end
