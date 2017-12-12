
using MaxwellFrequency
using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
hasJOcTree = false
try
  using JOcTree
  hasJOcTree = true
catch
  hasJOcTree = false
end
include("Maxwell-derivative-test.jl")

if hasJOcTree
    @testset "MaxwellFrequency Octree tests" begin
        # Create OcTree mesh
        h  = [1/32,1/32,1/64]
        n  = [128,64,32]
        x0 = [-2.0,-1.0,0.0]
        S  = createOcTreeFromBox(x0[1],x0[2],x0[3],n[1],n[2],n[3],h[1],h[2],h[3],0.0,0.0,0.0,0.0,0.0,1.0,2,2)
        M  = getOcTreeMeshFV(S,h,x0=x0)

        # Setup problem with random sources (hanging edge elimination is essential)
        Curl    = getCurlMatrix(M)
        Ne,     = getEdgeConstraints(M)
        Nf,Qf   = getFaceConstraints(M)
        Curl    = Qf * Curl * Ne
        Sources = Ne * (Curl' * complex.(randn(size(Curl,1),2),randn(size(Curl,1),2)))
        Obs     = Ne * (Curl' * complex.(sprandn(size(Curl,1),2,0.001),sprandn(size(Curl,1),2,0.001)))
        freq    = 1e2
        Ainv    = getMUMPSsolver([],1)
        param   = getMaxwellFreqParam(M,Sources,Obs,freq,Ainv)

        # Isotropic conductivity
        m1 = exp.(randn(M.nc))
        # Diagonally anisotropic conductivity
        m3 = exp.(rand(3*M.nc))
        # Generally anisotropic conductivity
        m6 = zeros(6*M.nc)
        for i = 1:M.nc
            A   = randn(3,3)
            u,V = eig(A + A')
            B   = V' * diagm(exp.(u)) * V
            m6[       i] = B[1,1]
            m6[  M.nc+i] = B[2,2]
            m6[2*M.nc+i] = B[3,3]
            m6[3*M.nc+i] = B[1,2]
            m6[4*M.nc+i] = B[1,3]
            m6[5*M.nc+i] = B[2,3]
        end

        @testset "Conductivity only" begin
            for (sensMethod, sensLabel) in ((:Explicit, "explicit"), (:Implicit, "implicit"))
                param.sensitivityMethod = sensMethod

                # Test implicit sensitivities
                println("Testing OcTree mesh FV ", sensLabel, " sensitivities")

                for (m,mLabel) in ((m1,"isotropic"), (m3,"diagonally anisotropic"), (m6, "generally anisotropic"))
                    println("Testing ", mLabel, " conductivity")

                    # Call getData
                    D,param = getData(m,param)

                    # Derivative test
                    function f(sigdum)
                        d, = getData(sigdum,param)
                        return d
                    end
                    df(zdum,sigdum) = getSensMatVec(zdum,sigdum,param)
                    checkDerivativeMaxwellOcTreeFV,error,order = checkDerivativeMax(f,df,m)
                    @test checkDerivativeMaxwellOcTreeFV

                    # adjoint test
                    x   = randn(length(m))
                    u   = complex.(randn(size(D)),randn(size(D)))
                    v   = getSensMatVec(x,m,param)
                    w   = getSensTMatVec(vec(u),m,param)
                    uv  = real(dot(vec(u),vec(v)))
                    wx  = dot(w,x)
                    tol = 1e-13*max(abs(uv),abs(wx))
                    @test uv ≈ wx atol=tol
                end
            end
        end  # End conductivity only testset

        @testset "Susceptibility only" begin
            param.sensitivityMethod = :Implicit
            scaleFac = 100
            for (m,mLabel) in ((m1,"isotropic"), (m3,"diagonally anisotropic"))
                println("Testing ", mLabel, " Susceptibility")

                n   = mLabel == "isotropic" ? M.nc : 3*M.nc
                chi = scaleFac*rand(n)
                mu  = mu0*(1+chi)
                dmu = mu0*(scaleFac/2)*randn(n)
                dmu[mu+dmu .< mu0] = 0.0

                values = Dict{String,Vector{Float64}}("sigmaCell"=>m,"muCell"=>mu)
                m0 = MaxwellFreqModel(values,["muCell"])

                # Call getData
                D,param = getData(m0,param)

                # Derivative test
                function f(mudum)
                    valsIn = Dict{String,Vector{Float64}}("sigmaCell"=>m,"muCell"=>mudum)
                    mIn = MaxwellFreqModel(valsIn,["muCell"])
                    d, = getData(mIn,param)
                    return d
                end
                function df(zdum,mudum)
                    zvals  = Dict{String,Vector{Float64}}("muCell"=>zdum)
                    muvals = Dict{String,Vector{Float64}}("sigmaCell"=>m,"muCell"=>mudum)
                    zIn = MaxwellFreqModel(zvals,["muCell"])
                    muIn = MaxwellFreqModel(muvals,["muCell"])
                    return getSensMatVec(zIn,muIn,param)
                end
                checkDerivativeMaxwellOcTreeFV,error,order =
                    checkDerivativeMax(f,df,mu,v=dmu,base=2.0)
                @test checkDerivativeMaxwellOcTreeFV

                # adjoint test
                x   = randn(length(m))
                xm  = MaxwellFreqModel(Dict{String,Vector{Float64}}("muCell"=>x),["muCell"])
                u   = complex.(randn(size(D)),randn(size(D)))
                v   = getSensMatVec(xm,m0,param)
                w   = getSensTMatVec(vec(u),m0,param).values["muCell"]
                uv  = real(dot(vec(u),vec(v)))
                wx  = dot(w,x)
                tol = 1e-13*max(abs(uv),abs(wx))
                @test uv ≈ wx atol=tol
            end
        end  # End sus only testset

        @testset "Conductivity and susceptibility" begin
            param.sensitivityMethod = :Implicit
            activeProps = ["sigmaCell","muCell"]
            scaleFac = 100
            for (m,mLabel) in ((m1,"isotropic"), (m3,"diagonally anisotropic"))
                println("Testing ", mLabel, " Con and sus")

                n    = mLabel == "isotropic" ? M.nc : 3*M.nc
                chi  = scaleFac*rand(n)
                mu   = mu0*(1+chi)
                dsig = randn(n)
                dsig[m+dsig .< 1e-8] = 1e-8
                dmu  = mu0*(scaleFac/2)*randn(n)
                dmu[mu+dmu .< mu0] = 0.0

                values = Dict{String,Vector{Float64}}("sigmaCell"=>m,"muCell"=>mu)
                m0 = MaxwellFreqModel(values,activeProps)

                # Call getData
                D,param = getData(m0,param)

                # Derivative test
                function f(mdum)
                    valsIn = Dict{String,Vector{Float64}}(
                               "sigmaCell"=>mdum[1:n],"muCell"=>mdum[n+1:end])
                    mIn = MaxwellFreqModel(valsIn,activeProps)
                    d, = getData(mIn,param)
                    return d
                end
                function df(zdum,mdum)
                    zvals  = Dict{String,Vector{Float64}}(
                               "sigmaCell"=>zdum[1:n],"muCell"=>zdum[n+1:end])
                    muvals = Dict{String,Vector{Float64}}(
                               "sigmaCell"=>mdum[1:n],"muCell"=>mdum[n+1:end])
                    zIn = MaxwellFreqModel(zvals,activeProps)
                    muIn = MaxwellFreqModel(muvals,activeProps)
                    return getSensMatVec(zIn,muIn,param)
                end
                checkDerivativeMaxwellOcTreeFV,error,order =
                    checkDerivativeMax(f,df,[m;mu],v=[dsig;dmu],base=2.0)
                @test checkDerivativeMaxwellOcTreeFV

                # adjoint test
                x   = randn(2*n)
                xvals  = Dict{String,Vector{Float64}}(
                           "sigmaCell"=>x[1:n],"muCell"=>x[n+1:end])
                xm  = MaxwellFreqModel(xvals,activeProps)
                u   = complex.(randn(size(D)),randn(size(D)))
                v   = getSensMatVec(xm,m0,param)
                wm  = getSensTMatVec(vec(u),m0,param)
                w   = [wm.values["sigmaCell"]; wm.values["muCell"]]
                uv  = real(dot(vec(u),vec(v)))
                wx  = dot(w,x)
                tol = 1e-13*max(abs(uv),abs(wx))
                @test uv ≈ wx atol=tol
            end
        end  # End sus and con testset

        @testset "Conductivity and susceptibility with Ne'* P proj mtx" begin
            Sources = Curl' * complex.(randn(size(Curl,1),2),randn(size(Curl,1),2))
            Obs     = Curl' * complex.(sprandn(size(Curl,1),2,0.001),sprandn(size(Curl,1),2,0.001))
            freq    = 1e2
            Ainv    = getMUMPSsolver([],1)
            param   = getMaxwellFreqParam(M,Sources,Obs,freq,Ainv)
            param.sensitivityMethod = :Implicit
            activeProps = ["sigmaCell","muCell"]
            scaleFac = 100
            for (m,mLabel) in ((m1,"isotropic"), (m3,"diagonally anisotropic"))
                println("Testing ", mLabel, " Con and sus")

                n    = mLabel == "isotropic" ? M.nc : 3*M.nc
                chi  = scaleFac*rand(n)
                mu   = mu0*(1+chi)
                dsig = randn(n)
                dsig[m+dsig .< 1e-8] = 1e-8
                dmu  = mu0*(scaleFac/2)*randn(n)
                dmu[mu+dmu .< mu0] = 0.0

                values = Dict{String,Vector{Float64}}("sigmaCell"=>m,"muCell"=>mu)
                m0 = MaxwellFreqModel(values,activeProps)

                # Call getData
                D,param = getData(m0,param)

                # Derivative test
                function f(mdum)
                    valsIn = Dict{String,Vector{Float64}}(
                               "sigmaCell"=>mdum[1:n],"muCell"=>mdum[n+1:end])
                    mIn = MaxwellFreqModel(valsIn,activeProps)
                    d, = getData(mIn,param)
                    return d
                end
                function df(zdum,mdum)
                    zvals  = Dict{String,Vector{Float64}}(
                               "sigmaCell"=>zdum[1:n],"muCell"=>zdum[n+1:end])
                    muvals = Dict{String,Vector{Float64}}(
                               "sigmaCell"=>mdum[1:n],"muCell"=>mdum[n+1:end])
                    zIn = MaxwellFreqModel(zvals,activeProps)
                    muIn = MaxwellFreqModel(muvals,activeProps)
                    return getSensMatVec(zIn,muIn,param)
                end
                checkDerivativeMaxwellOcTreeFV,error,order =
                    checkDerivativeMax(f,df,[m;mu],v=[dsig;dmu],base=2.0)
                @test checkDerivativeMaxwellOcTreeFV

                # adjoint test
                x   = randn(2*n)
                xvals  = Dict{String,Vector{Float64}}(
                           "sigmaCell"=>x[1:n],"muCell"=>x[n+1:end])
                xm  = MaxwellFreqModel(xvals,activeProps)
                u   = complex.(randn(size(D)),randn(size(D)))
                v   = getSensMatVec(xm,m0,param)
                wm  = getSensTMatVec(vec(u),m0,param)
                w   = [wm.values["sigmaCell"]; wm.values["muCell"]]
                uv  = real(dot(vec(u),vec(v)))
                wx  = dot(w,x)
                tol = 1e-13*max(abs(uv),abs(wx))
                @test uv ≈ wx atol=tol
            end
        end  # End sus and con with Ne'*P proj matrix testset
    end
end
