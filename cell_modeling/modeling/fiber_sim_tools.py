from dolfin import *
#import mshr
import ufl
import meshio
import numpy as np
from fiber_tools import *

mesh_dir = "../mesh/0609G1/"

# likely want to allow inputs for this to allow user to easily input files.
def setup():
    '''
    Create simple cell geometry, define different material domains,
    initiate finite element function spaces for scalar and vector variables.
    '''
    mesh = Mesh()
    #with XDMFFile('e10_dense_pruned.xdmf') as infile:
    with XDMFFile(mesh_dir + "new100.xdmf") as infile:
        infile.read(mesh)
    domains = MeshFunction('size_t',mesh,mesh.topology().dim())
    with XDMFFile(mesh_dir + "new100.xdmf") as infile:
        infile.read(domains)
    '''
    100: gel
    200: cytoplasm
    300: nucleus
    '''
    # degree of the finite element space
    degree = 1
    V = VectorFunctionSpace(mesh, 'P', degree)
    V0 = FunctionSpace(mesh, 'P', degree)
    return mesh, domains, V0, V

