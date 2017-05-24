
if hasJOcTree
	
export getOTMeshFromTxRx

function getOTMeshFromTxRx(x0::Array{Float64,1},
						     		n::Array{Int64,1},
						     		h::Array{Float64,1},	
							 		Srcs::Array{Array{Float64},1},
							 		Recs::Array{Array{Float64},1},
							 		nf=2,
							 		nc=2,
									Box::Array{Float64,2}=[2  2; 2  2; 2  4]  # 3x2 of how many h's to pad the txrx 
									)
										 
	x1 = Inf;  x2 = -Inf
	y1 = Inf;  y2 = -Inf
	z1 = Inf;  z2 = -Inf
	
	for k = 1:length(Srcs) 
		Src=Srcs[k]
		x1 = min(minimum(Src[:,1]) - Box[1,1]*h[1],x1)
		x2 = max(maximum(Src[:,1]) + Box[1,2]*h[1],x2)
		y1 = min(minimum(Src[:,2]) - Box[2,1]*h[2],y1)
		y2 = max(maximum(Src[:,2]) + Box[2,2]*h[2],y2)
		z1 = min(minimum(Src[:,3]) - Box[3,1]*h[3],z1)
		z2 = max(maximum(Src[:,3]) + Box[3,2]*h[3],z2)
	end

	for k=1:length(Recs) 
		Rec=Recs[k]
		x1 = min(minimum(Rec[:,1]) - Box[1,1]*h[1],x1)
		x2 = max(maximum(Rec[:,1]) + Box[1,2]*h[1],x2)
		y1 = min(minimum(Rec[:,2]) - Box[2,1]*h[2],y1)
		y2 = max(maximum(Rec[:,2]) + Box[2,2]*h[2],y2)
		z1 = min(minimum(Rec[:,3]) - Box[3,1]*h[3],z1)
		z2 = max(maximum(Rec[:,3]) + Box[3,2]*h[3],z2)
	end
	L = n.*h
	if (x1<x0[1]) || (x2>x0[1]+L[1]) || (y1<x0[2]) || (y2>x0[2]+L[2]) || (z1<x0[3]) || (z2>x0[3]+L[3])
		error("sources/receivers are outside to domain") 
	end

	S  = createOcTreeFromBox(x0[1],x0[2],x0[3],n[1],n[2],n[3],h[1],h[2],h[3],x1,x2,y1,y2,z1,z2,nf,nc)

   return S
end	

end