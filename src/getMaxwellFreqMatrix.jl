export getMaxwellFreqMatrix

"""
function MaxwellFreq.getMaxwellFreqMatrix

builds finite volume discretization of the Maxweel Frequency domain operator

    A = Curl * mu0 * Curl - im * w * sigma(x)

Inputs:

    sigma::Vector{Float64} - conductivities (cell-centered)
    w::Float64             - Angular Frequency
    mesh::AbstractMesh     - mesh from jInv.Mesh

Output:

    A  - PDE operator (sparse matrix)

"""
function getMaxwellFreqMatrix(sigma::Vector{Float64}, w::Float64,  mesh::AbstractMesh)

    Curl = getCurlMatrix(mesh)

    Msig = getEdgeMassMatrix(mesh, sigma)
    Mmu  = getFaceMassMatrix(mesh, fill(1/mu0, mesh.nc))

    # eliminate hanging edges and faces
    Ne, = getEdgeConstraints(mesh)
    Nf,Qf = getFaceConstraints(mesh)

    Curl = Qf  * Curl * Ne
    Msig = Ne' * Msig * Ne
    Mmu  = Nf' * Mmu  * Nf

    A   = Curl' * Mmu * Curl - (im * w) * Msig

    return A
end

function getMaxwellFreqMatrix(sigma::Vector{Float64},pFor::MaxwellFreqParam,doClear::Bool=true)
    A = getMaxwellFreqMatrix(sigma,pFor.freq,pFor.Mesh)
    return A
end
