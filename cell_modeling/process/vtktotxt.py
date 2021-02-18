# This script will be used to output displacement data as a txt file, had to copy the data previously.

import meshio
import numpy as np

path = "./"
out_path = "./circular_txt_data/"
filename = ""
physical_num = 1  # physical number marking surface (not sure what this does)
'''
likely want to read all of the vtk files, but I think just using the final displacement is ok for now.
I think I can just view the sf orienation data through the .pvd files; can use pvd and csv for different things
within the same paraview state file
'''
# copying vtk filename and location, going to save it in txt data folder
solution_final_x = meshio.read("result_circular/solution_x000010.vtu")
solution_final_y = meshio.read("result_circular/solution_y000010.vtu")
solution_final_z = meshio.read("result_circular/solution_z000010.vtu")

'''
	This was taken from John's code, it is basically a way to convert the data after running the vtk sim.
	Should do this for the stress fibers and the displacements.

	Should also incorporate this into the code itself; running this as a script is somewhat unnecessary.

	Would be cool to maybe create a big dataset and combine the vtk files; would allow for a time-series
	csv file if possible. not entirely sure how to format this though, would need to go into paraview documentation
'''

#this will be scrappy, definitely want to be able to reformat
#also...just suck it up and get better at paraview. I just gotta do that at the end of the day. I live .txt better though.

#extract point data (the displacement data) for x,y,z coordiantes
ux = solution_final_x.point_data
uy = solution_final_y.point_data
uz = solution_final_z.point_data

#extract displacements from dict and organize into numpy array
x = ux["ux"]
y = uy["uy"]
z = uz["uz"]

u_arr = np.column_stack((x,y,z))

# Tabulate dof coordinates
#u_arr = u.compute_vertex_values()  # 1-d numpy array
#length = np.shape(u_arr)[0]
#u_arr = np.reshape(u_arr, (length//3, 3), order="F") # Fortran ordering (IDK WHAT THIS MEANS)

# Save txt solution
np.savetxt(out_path + "displacement" + ".txt", u_arr, delimiter=" ")
np.savetxt(out_path + "xdisplacement" + ".txt", x, delimiter=" ")
np.savetxt(out_path + "ydisplacement" + ".txt", y, delimiter=" ")
np.savetxt(out_path + "zdisplacement" + ".txt", z, delimiter=" ")

#also saving x,y,z displacement individually if I messed up. can fix later and use current matlab post_process script.
