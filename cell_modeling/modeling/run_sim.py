# this file will handle running simulations.

#Ideas:
'''
    1. Accept txtfile or json file input data to determine what sim we are going to run
    2. create functions to output data in a pretty way
        2a. use paraview functions to make super cool paraview visualizations
        2b. create an option to run viz off completely
    3. 

'''

from dolfin import *
#import mshr
import ufl
import meshio
import numpy as np
from fiber_tools import *
from fiber_sim_tools import *

# this data can be inputed through a txt file
mesh_dir = "../mesh/0609G1/"
intput = ""
output_dir = "../output_data/test2"
mesh_file = "new100.xdmf"

# don't touch this yet, idk what it means.
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 5
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
set_log_level(20)


def run(output_dir,ffc_options):
    '''
    Define the solution, and external fields.
    Run solver to compute the solution and save it periodically.
    '''

    mesh, domains, V0, V = setup(mesh_dir,mesh_file)
    u = Function(V)
    # time-dependent field of body force
    B = Expression(('0.','0.','t*0.'),t=0.,element=V.ufl_element())

    # spatially defined contractile strength function, likely want to define this within the run class
    #f = Expression(('t*5'),t=0.,element=V0.ufl_element())

    # time is not advancing for some reason; likely there's an error in my subclass
    f = ContractileStrength(t=0.,element=V0.ufl_element())

    # no traction boundary condition
    T = Constant((0.,0.,0.))

    # time increment
    step = 20
    dt = 1./step
    freq_checkout = 2 
    # file names to save the solutions
    vtkfile_x = File(output_dir + '/solution_x.pvd')
    vtkfile_y = File(output_dir + '/solution_y.pvd')
    vtkfile_z = File(output_dir + '/solution_z.pvd')
    vtkfile_material = File(output_dir + '/domains.pvd')
    vtkfile_material << domains    
    vtkfile_mx = File(output_dir + '/deformed_fiber_x.pvd')
    vtkfile_my = File(output_dir + '/deformed_fiber_y.pvd')
    vtkfile_mz = File(output_dir + '/deformed_fiber_z.pvd')

    # within this, we can add certain info/things
    for n in range(step+1):
        print('n = %d'%(n))
        t = n*dt
        u, B, m = solver(u,mesh,domains,V0,V,B,T,f,ffc_options)

        if n%freq_checkout is 0:
            # split the solution to get displacement in x and y directions
            ux = u.sub(0)
            uy = u.sub(1)
            uz = u.sub(2)
            print(norm(ux))
            print(norm(uy))
            print(norm(uz))
            ux.rename('ux','x disp')
            uy.rename('uy','y disp')
            uz.rename('uz','z disp')
            proj_m = project(m,V)
            mx = proj_m.sub(0)
            my = proj_m.sub(1)
            mz = proj_m.sub(2)
            mx.rename('mx','x fiber')
            my.rename('my','y fiber')
            mz.rename('mz','z fiber')
            # save the solution
            vtkfile_x << (ux,t)
            vtkfile_y << (uy,t)
            vtkfile_z << (uz,t)
            vtkfile_mx << (mx,t)
            vtkfile_my << (my,t)
            vtkfile_mz << (mz,t)

            # at each timestep, output the displacement data as a txt file; NOT WORKING AS OF 3/4/2021 12:19pm
            #vtktotxt_disp(vtkfile_x,vtkfile_y,vtkfile_z,output_dir)

        # advance the fields
        B.t = B.t+dt
        f.t = f.t+dt 
        print(f.t)

if __name__ == '__main__':
    run(output_dir,ffc_options)
