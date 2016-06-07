using JOcTree
#using Utils
using MaxwellFrequency
using MaxwellUtils
using jInv.LinearSolvers
using jInv.InverseSolve

# ------- SETUP PARAMETERS FOR THE MODEL AND DATA

# data and topo files
#datafile = ["data_inv.txt"]
datafile = [ "data_inv.txt_data",
             "data_inv.txt_trx",
             "data_inv.txt_rcv",
             "data_inv.txt_frq" ]

topofile = "topo.txt"

# of cells in base mesh
n     = vec([ 256  256 1024 ])
# corner of the mesh
x0    = vec([ 7.56912000E+05  4.15763700E+06  -7852.178 ])
# total mesh lengths  
meshL = vec([ 6400. 6400. 10240. ])


# parameters for meshing
nsmallcells    = vec([1 1 1])  #  # of small cells around each point.
mincellsize    = 1  #  minimum cell size in the data area
depth_core     = vec([50. 100. 200.])  # how far to go down in the core region for the fwd meshes
depth_core_inv = vec([200. 300. 500.])  # how far to go down in the core region for the inv meshes

mincellfactor = 1    # minimum cellsize below topo


# parameters for the forward problem
fname = ""    # leave empty for now
doFV  = true  # use finite volume (other option FEM for finite elements)
doSE  = false # SE = sensitivity explicit - store sensitivities and not factorizations
 
# reference conductivity 
halfSpaceCond = 0.98
backCond = 3.0


# lower bounds
BL = 1e-6
# Higher bounds
BH = 1e+4

# Regularization function 
regtype = wdiffusionReg
# parameters for the regularization function
regparams = [sqrt(1.0), sqrt(1.0), sqrt(1.0), 5e-7]  # alphax  alphay  alphaz  alphas

# For TVp regularization
#regfun = wPTV
#regparams = [1e-2, 1.0, 2.0, 3.0, 0.1] # epsilon, alphax  alphay  alphaz,  p


beta = 1e-32
# misfit function
misfun = SSDFun
#  inner CG iter
cgit = 20 
# maximum iter for the inversion
maxit = 5


# approximate mesh interpolation matrix (inv -> fwd) using [2^ninterp]^3 quadrature points
# (set ninterp = [] to use full interpolation matrix)
ninterp = 3

# store global mesh to local mesh interpolation matrix in compact form
compact = true

# model parameter is log conductivity
modfun = expMod

surfweight = vec([1000.  500.  100. 50. 25.])  # surface interface weights
#surfweight = []

return