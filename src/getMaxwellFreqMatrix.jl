export getMaxwellFreqMatrix

"""
function MaxwellFreq.getMaxwellFreqMatrix

builds finite volume discretization of the Maxwell Frequency domain operator

    A = Curl * mu * Curl ± im * 2 * pi * f * sigma(x)

Inputs:

    m::MaxwellFreqModel - conductivities (cell-centered)
    f::Float64             - frequency
    mesh::AbstractMesh     - mesh from jInv.Mesh

Output:

    A  - PDE operator (sparse matrix)

"""
function getMaxwellFreqMatrix(m::MaxwellFreqModel, f::Float64,  mesh::AbstractMesh, minusIOmega::Bool)

    sigma = m.values["sigmaCell"]
    muInv = haskey(m.values,"muCell") ? 1./m.values["muCell"] : fill(1./mu0, mesh.nc)

    Curl = getCurlMatrix(mesh)

    Msig = getEdgeMassMatrix(mesh, sigma)
    Mmu  = getFaceMassMatrix(mesh, muInv)

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

function getMaxwellFreqMatrix(m::MaxwellFreqModel, pFor::MaxwellFreqParam, doClear::Bool=true)
    A = getMaxwellFreqMatrix(m, pFor.frequency, pFor.Mesh, pFor.timeConvention == :ExpMinusImOmegaT)
    return A
end
