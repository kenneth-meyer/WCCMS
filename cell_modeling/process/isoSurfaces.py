from sim_tools import level_sets
import numpy as np
import meshio

# this is me using the functions John wrote to create isosurface plots for my ellipsoidal contraction data

path = "./"
out_path = "./ellipsoid_txt_data/"
filename = ""
mesh_path = "e10_dense_surface.xdmf"

#loading displacements from what I outputted earlier; going backwards/kinda unnecessary. but don't want to have to 
# convert vtu to numpy every single time.
u = np.loadtxt(out_path + "displacement" + ".txt", delimiter=" ")

# need to load the surface mesh; can only have the surface mesh
surf_mesh = meshio.read(mesh_path)
surf_vert = np.array(surf_mesh.points)
surf_conn = np.array(surf_mesh.cells[0].data)

## Isosurfaces
sets = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]
level_sets(sets, surf_vert, surf_conn, u, out_path)