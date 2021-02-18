#This code is from John's github, I edited it to my needs.
#utilizing this script to better visualize my results in paraview
import meshio
import numpy as np

#should create separate files for meshes, models, data, etc. will help me. should do this tomorrow.
path = "./"
out_path = "./circular_txt_data/"
filename = "e_dense_10"
physical_num = 1  # physical number marking surface

# Read mesh
msh = meshio.read(path + filename + ".msh")

#not really used in this file I don't think

# Get triangle and tet connectivity
for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data  # 2D connectivity

    elif  cell.type == "tetra":
        tetra_cells = cell.data     # 3D connectivity

# Get physical labels
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]
    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]

# Get surface cells
surf_cells = np.column_stack((triangle_cells, triangle_data)) 
surf_cells = surf_cells[surf_cells[:,-1] == physical_num]    # Extract cells on surface
surf_cells = surf_cells[:,0:-1]

# Get Nodes
nodes = msh.points # all nodes
nodes_new = []
surf_cells_new = surf_cells
vert_map = {}      # maps old nodes to new nodes
num_vert = 0

# Rewrite connectivity
for i, face in enumerate(surf_cells_new):
    for j, vert in enumerate(face):

        if vert in vert_map:
            surf_cells_new[i][j] = vert_map[vert]
        else:
            vert_map[vert] = num_vert          # add key to map
            surf_cells_new[i][j] = num_vert    
            nodes_new.append(nodes[int(vert)]) # add node
            num_vert += 1

nodes_new = np.array(nodes_new)

# below code saves the output in the desired format so my visualizations can be improved.

# All vertices
np.savetxt(out_path + "vertices.txt", nodes, delimiter=" ")

# Surface vertices and faces
np.savetxt(out_path + "surf_vertices.txt", nodes_new, delimiter=" ")
np.savetxt(out_path + "surf_faces.txt", surf_cells_new, delimiter=" ", fmt='%d')
