from dolfin import *
#import mshr
import ufl
import meshio
import numpy as np
from fiber_tools import *

#mesh_dir = "../mesh/0609G1/"
#mesh_file = "new100.xdmf"

# likely want to allow inputs for this to allow user to easily input files.
def setup(mesh_dir,mesh_file):
    '''
    Create simple cell geometry, define different material domains,
    initiate finite element function spaces for scalar and vector variables.
    '''
    mesh = Mesh()
    #with XDMFFile('e10_dense_pruned.xdmf') as infile:
    with XDMFFile(mesh_dir + mesh_file) as infile:
        infile.read(mesh)
    domains = MeshFunction('size_t',mesh,mesh.topology().dim())
    with XDMFFile(mesh_dir + mesh_file) as infile:
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

def solver(u, mesh, domains, V0, V, B, T, f,ffc_options):
    '''
    Solve the boundary value problem with body force B and boundary traction T
    and active fiber contraction f
    u is the solution to be computed
    '''
    # boundary of the full domain where Dirichlet condition is prescribed
    def boundary(x, on_boundary):
        return on_boundary
    u_D = Expression(('0.','0.','0.'), degree=2)
    bc = DirichletBC(V, u_D, boundary)

    # Define variational problem
    du = TrialFunction(V)
    v = TestFunction(V)

    d = u.geometric_dimension()
    I = Identity(d)
    F = I+grad(u)
    C = F.T*F
    Ic = tr(C)
    J = det(F)
    C_bar = C/J**(2./3)
    Ic_bar = tr(C_bar)

    # prescribe material properites for different regions
    E_0 = 1.
    E_c = 10.
    E_n = 100.
    nu_0 = 0.49
    nu_c = 0.2
    nu_n = 0.4999
    mu_0 = E_0/2/(1+nu_0) 
    mu_c = E_c/2/(1+nu_c) 
    mu_n = E_n/2/(1+nu_n) 
    K_0 = E_0/3/(1-2*nu_0)
    K_c = E_c/3/(1-2*nu_c)
    K_n = E_n/3/(1-2*nu_n)
    dx = Measure('dx',domain=mesh,subdomain_data=domains)
    psi_0 = mu_0/2*(Ic_bar-3) + K_0/2*(ln(J))**2
    psi_c = mu_c/2*(Ic_bar-3) + K_c/2*(ln(J))**2
    psi_n = mu_n/2*(Ic_bar-3) + K_n/2*(ln(J))**2
 
    # assemble the total potential energy
    Pi = psi_0*dx(100)+psi_c*dx(200)+psi_n*dx(300) - dot(B,u)*dx('everywhere') - dot(T,u)*ds

    #Scheme 1
    #ef = as_vector([0,0,1]) # fiber orientation (undeformed)
    #Scheme 2

    # could list the distributions as scheme numbers, could be an easy/quick way.
    
    # ef is the initial, undeformed fiber orientation...is this the desired orientaion???
    ef = FiberOrientation() # class that defines sf orientation foudn in fiber_tools.py.
    I4 = sqrt(dot(ef,C*ef))
    m = F*ef/I4 # fiber orientation (deformed)

    # calcualtes cauchy stress in the sf
    Tsf = f*I4/J*as_matrix([[m[0]*m[0], m[0]*m[1], m[0]*m[2]],
            [m[1]*m[0], m[1]*m[1], m[1]*m[2]],
            [m[2]*m[0], m[2]*m[1], m[2]*m[2]]
            ])
    # could try to split f into 2 components to put sf density into the model
    # for now, assume that everything is uniformly distributed
    # still assuming fiber only concentrates along mean direction.

    # pilo-kirchoff stress in the sf
    Psf = J*Tsf*inv(F.T)

    # take Gateaux derivative of Pi
    A = derivative(Pi, u, v) + inner(Psf,grad(v))*dx(200)
    # calculate Jacobian
    J = derivative(A, u, du)

    # Compute solution
    solve(A == 0, u, bc, J=J, form_compiler_parameters=ffc_options)
    return u, B, m

'''
def run(output_dir,ffc_options):
    
    #Define the solution, and external fields.
    #Run solver to compute the solution and save it periodically.
    

    mesh, domains, V0, V = setup()
    u = Function(V)
    # time-dependent field of body force
    B = Expression(('0.','0.','t*0.'),t=0.,element=V.ufl_element())

    # spatially defined contractile strength function, likely want to define this within the run class
    f = Expression(('t*5'),t=0.,element=V0.ufl_element())

    # time is not advancing for some reason; likely there's an error in my subclass
    #f = ContractileStrength(t=0.,element=V0.ufl_element())

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
        # advance the fields
        B.t = B.t+dt
        f.t = f.t+dt 
        print(f.t)
'''

# converts the vtk output data to txt data, which can be used to compute surface normals
# likely also want to translate the post_process.m file into a python function so I don't have to go back and forth
def vtktotxt_disp(x_file,y_file,z_file,out_path):
    # copying vtk filename and location, going to save it in txt data folder
    solution_final_x = meshio.read(x_file)
    solution_final_y = meshio.read(y_file)
    solution_final_z = meshio.read(z_file)

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

    # Save txt solution, which is the most important thing
    np.savetxt(out_path + "displacement" + ".txt", u_arr, delimiter=" ")
    #np.savetxt(out_path + "xdisplacement" + ".txt", x, delimiter=" ")
    #np.savetxt(out_path + "ydisplacement" + ".txt", y, delimiter=" ")
    #np.savetxt(out_path + "zdisplacement" + ".txt", z, delimiter=" ")


# not sure what to return