# This is a script that converts .msh files to .xdmf files

# author: Kenneth L Meyer
# date: 12/21/2020

# first draft of this script was simply taken from fenics discussion
# post expaining how to do this, slide 21 link in Xinzeng's ppt
# altered to make it easier to read

# basically the same as John's code; copied his to see if i had a typo in mine/something I was missing

import meshio
import os

#path = os.path.dirname(os.path.realpath(__file__))

msh = meshio.read("./e_dense_10.msh")

for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data

    elif  cell.type == "tetra":
        tetra_cells = cell.data

# Get physical labels
for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]
    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]
        tetra_data2 = msh.cell_data_dict["gmsh:geometrical"][key]


# Save separate, trying to add extra data missing in fenicsMesh.xdmf
tetra_mesh = meshio.Mesh(points=msh.points, 
                        cells=[("tetra", tetra_cells)],
                        cell_data={"gmsh:physical": [tetra_data],"gmsh:geometrical": [tetra_data2]})
                        #field_data={)

triangle_mesh =meshio.Mesh(points=msh.points,
                           cells=[("triangle", triangle_cells)],
                           cell_data={"triangle":[triangle_data]})

# write to xdmf
meshio.write("e10_dense_pruned.xdmf", tetra_mesh)
# writing the surface mesh; need surface mesh only to create isosurfaces
meshio.write("e10_dense_surface.xdmf", triangle_mesh)

#I need to also use the file john and jp have to output this as .txt data to help my paraview visualizations


#meshio.write("mesh.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]}))
#meshio.write("mf.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                        #cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
#meshio.write("cf.xdmf", meshio.Mesh(
    #points=msh.points, cells={"tetra": msh.cells["tetra"]},
    #cell_data={"tetra": {"name_to_read": msh.cell_data["tetra"]["gmsh:physical"]}}))



#old, deleted code:

#for cell in msh.cells: 
#   if cell.type == "triangle":
#       tri_cells = cell.data
#   elif cell.type == "tetra":
#       tetra_cells = cell.data

#don't fully understand the part where we "rewrite connectivity"

#for key in msh.cell_data_dict["gmsh:physical"].keys():
#   if key == "triangle":
#       tri_data = msh.cell_data_dict["gmsh:physical"][key]
#   elif key == "tetra":
#       tetra_data = msh.cell_data_dict["gmsh:physical"][key] #not sure why this is here

# generating meshes after the 2D and 3D meshes have been seperated
#tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})
#tri_mesh = meshio.Mesh(points=msh.points,
#                       cells=[("triangle", tri_cells)],
#                       cell_data={"triangle":[tri_data]})



# writing 3d and 2d meshes to be used in fenics and visualization, respectively
#meshio.write("fenicsMesh.xdmf", tetra_mesh)
#meshio.write("surfaceMesh.xdmf", tri_mesh)
