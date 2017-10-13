export getMaxwellFreqMatrix

"""
function MaxwellFreq.getMaxwellFreqMatrix

builds finite volume discretization of the Maxweel Frequency domain operator

    A = Curl * mu0 * Curl ± im * 2 * pi * f * sigma(x)

Inputs:

    sigma::Vector{Float64} - conductivities (cell-centered)
    f::Float64             - frequency
    mesh::AbstractMesh     - mesh from jInv.Mesh

Output:

    A  - PDE operator (sparse matrix)

"""
function getMaxwellFreqMatrix(sigma::Vector{Float64}, f::Float64,  mesh::AbstractMesh, minusIOmega::Bool)

    Curl = getCurlMatrix(mesh)

    Msig = getEdgeMassMatrix(mesh, sigma)
    Mmu  = getFaceMassMatrix(mesh, fill(1/mu0, mesh.nc))

    # eliminate hanging edges and faces
    Ne, = getEdgeConstraints(mesh)
    Nf,Qf = getFaceConstraints(mesh)

    Curl = Qf  * Curl * Ne
    Msig = Ne' * Msig * Ne
    Mmu  = Nf' * Mmu  * Nf

    iw = (minusIOmega ? -im : im) * 2 * pi * f

    A   = Curl' * Mmu * Curl + iw * Msig

    return A
end

function getMaxwellFreqMatrix(sigma::Vector{Float64}, pFor::MaxwellFreqParam, doClear::Bool=true)
    A = getMaxwellFreqMatrix(sigma, pFor.frequency, pFor.Mesh, pFor.timeConvention == :ExpMinusImOmegaT)
    return A
end
