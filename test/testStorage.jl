
using Base.Test

using MaxwellFrequency
using jInv.Mesh
using jInv.LinearSolvers

# Setup problem
mesh = getRegularMesh([0 10 0 10 0 10], [10; 10; 10])

# Tx & Rx
Curl = getCurlMatrix(mesh)
Sources = map(Complex{Float64}, Curl'*randn(sum(mesh.nf), 5))
Obs = Curl'*sprandn(sum(mesh.nf), 20, 0.001)

# Solver
solver = getMUMPSsolver([], 0)

# Build test vectors
sigma = rand(mesh.nc)
x = rand(mesh.nc)
y = rand(size(Obs, 2)*size(Sources, 2)) + im*rand(size(Obs, 2)*size(Sources, 2))

@testset "Implicit Sensitivity Storage Tests" begin
    storageLevel = :Implicit
    @testset "No Storage" begin

        param = getMaxwellFreqParam(mesh, Sources, Obs, 100., solver)
        param.sensitivityMethod = storageLevel
        param.storageLevel = :None

        D, param = getData(sigma, param)
     
        # No storage of A-matrix
        @test isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        Jx = getSensMatVec(x, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        Jty = getSensTMatVec(y, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        # Adjoint test
        @test real(y' * Jx) ≈ real(x' * Jty) atol=1.0e-14

    end

    @testset "Matrix Storage" begin

        param = getMaxwellFreqParam(mesh, Sources, Obs, 100., solver)
        param.sensitivityMethod = :Implicit
        param.storageLevel = :Matrices

        D, param = getData(sigma, param)

        # Storing of A-matrix
        @test ~isempty(param.Matrices)
        # No storage of factors    
        @test param.Ainv.Ainv == []

        Jx = getSensMatVec(x, sigma, param)

        # Storing of A-matrix
        @test ~isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        Jty = getSensTMatVec(y, sigma, param)

        # Storing of A-matrix
        @test ~isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        # Adjoint test
        @test real(y' * Jx) ≈ real(x' * Jty) atol=1.0e-14
    end

    @testset "Factor Storage" begin

        param = getMaxwellFreqParam(mesh, Sources, Obs, 100., solver)
        param.storageLevel=:Factors

        D, param = getData(sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # Storing factors
        @test param.Ainv.Ainv != []

        Jx = getSensMatVec(x, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # Storing factors
        @test param.Ainv.Ainv != []

        Jty = getSensTMatVec(y, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # Storing factors
        @test param.Ainv.Ainv != []

        # Adjoint test
        @test real(y' * Jx) ≈ real(x' * Jty) atol=1.0e-14
    end
end

@testset "Explicit Sensitivity Storage Tests" begin
    storageLevel = :Explicit
    @testset "No Storage" begin

        param = getMaxwellFreqParam(mesh, Sources, Obs, 100., solver)
        param.sensitivityMethod = storageLevel
        param.storageLevel = :None

        D, param = getData(sigma, param)
     
        # No storage of A-matrix
        @test isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        Jx = getSensMatVec(x, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        Jty = getSensTMatVec(y, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        # Adjoint test
        @test real(y' * Jx) ≈ real(x' * Jty) atol=1.0e-14

    end

    @testset "Matrix Storage" begin

        param = getMaxwellFreqParam(mesh, Sources, Obs, 100., solver)
        param.sensitivityMethod = :Implicit
        param.storageLevel = :Matrices

        D, param = getData(sigma, param)

        # Storing of A-matrix
        @test ~isempty(param.Matrices)
        # No storage of factors    
        @test param.Ainv.Ainv == []

        Jx = getSensMatVec(x, sigma, param)

        # Storing of A-matrix
        @test ~isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        Jty = getSensTMatVec(y, sigma, param)

        # Storing of A-matrix
        @test ~isempty(param.Matrices)
        # No storage of factors
        @test param.Ainv.Ainv == []

        # Adjoint test
        @test real(y' * Jx) ≈ real(x' * Jty) atol=1.0e-14
    end

    @testset "Factor Storage" begin

        param = getMaxwellFreqParam(mesh, Sources, Obs, 100., solver)
        param.storageLevel=:Factors

        D, param = getData(sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # Storing factors
        @test param.Ainv.Ainv != []

        Jx = getSensMatVec(x, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # Storing factors
        @test param.Ainv.Ainv != []

        Jty = getSensTMatVec(y, sigma, param)

        # No storage of A-matrix
        @test isempty(param.Matrices)
        # Storing factors
        @test param.Ainv.Ainv != []

        # Adjoint test
        @test real(y' * Jx) ≈ real(x' * Jty) atol=1.0e-14
    end
end
