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

# this data can be inputed through a txt file
mesh_dir = "../mesh/0609G1/"
intput = ""
output = "../output_data/test"

# don't touch this yet, idk what it means.
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 5
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
set_log_level(20)

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

# don't need to touch this for now; high-level modularization of inputs will be done; this is the hard part.
def solver(u, mesh, domains, V0, V, B, T, f):
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
    
    #ef = MyFiber()
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

# defining contractile strength (f) as a class to give us the ability to spatially define it.

'''
class ff(UserExpression):
        
    def eval(self, value, x):
            
        # define contractile strenght spatially

        # start off with changing it with z position or something.

    def value_shape(self):
        return (1,)
        # not sure how to define the shape/how to access the shape of the stress field.
'''
def run():
    '''
    Define the solution, and external fields.
    Run solver to compute the solution and save it periodically.
    '''

    mesh, domains, V0, V = setup()
    u = Function(V)
    # time-dependent field of body force
    B = Expression(('0.','0.','t*0.'),t=0.,element=V.ufl_element())

    # EXTENDING THE FIBER CONTRACTION TO BE A FUNCTION OF POSITION
    # likely want to turn this into a class...not entirely sure how to do this.

    # changing the strength based on the z position of the coordinate.

    # having trouble accessing x outside of "solver"
    #yy = abs(x[2])

    # time-dependent field of active fiber contraction 
    #f = Expression(('t*5'),t=0.,element=V0.ufl_element())
    f = Expression(('t*5'),t=0.,element=V0.ufl_element())

    #f = ContractileStrength(t=0.,element=V0.ufl_element())
    # i think that this will be wrong. I never include the element=V0...... in it.


    # no traction boundary condition
    T = Constant((0.,0.,0.))

    # time increment
    step = 20
    dt = 1./step
    freq_checkout = 2 
    # file names to save the solutions
    vtkfile_x = File(output + '/solution_x.pvd')
    vtkfile_y = File(output + '/solution_y.pvd')
    vtkfile_z = File(output + '/solution_z.pvd')
    vtkfile_material = File(output + '/domains.pvd')
    vtkfile_material << domains    
    vtkfile_mx = File(output + '/deformed_fiber_x.pvd')
    vtkfile_my = File(output + '/deformed_fiber_y.pvd')
    vtkfile_mz = File(output + '/deformed_fiber_z.pvd')
    for n in range(step+1):
        print('n = %d'%(n))
        t = n*dt
        u, B, m = solver(u,mesh,domains,V0,V,B,T,f)

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

if __name__ == '__main__':
    run()
