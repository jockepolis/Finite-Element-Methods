from __future__ import print_function
from fenics import *

# Create mesh and define function space
mesh = Mesh("circle.xml.gz")

# Construct the finite element space
P1 = FiniteElement ("Lagrange", mesh.ufl_cell(), 1)
TH = P1 * P1 * P1
V = FunctionSpace(mesh, TH)

# Define parameters :
T = 500
dt = 0.5
numsteps = T/dt
delta1 = 1
delta2 = 1
delta3 = 1
alpha = 0.4
beta = 0.8
gamma = 0.8
zeta = 2
L_0 = 0.4
l = 0.6
m = 0.12

# Class representing the intial conditions
class InitialConditions(UserExpression):
    
    def eval(self, values, x):
        #values[0] = 0.1*((4/25)-2*pow(10,-7)*(x[0]-0.1*x[1]-225)*(x[0]-0.1*x[1]-675))
        values[0] = 0
        values[1] = (4/25)-2*pow(10,-7)*(x[0]-0.1*x[1]-225)*(x[0]-0.1*x[1]-675)
        values[2] = (22/45)-3*pow(10,-5)*(x[0]-450)-1.2*pow(10,-4)*(x[1]-150)
    
    def value_shape(self):
        return (3,)

# Define initial condition
indata = InitialConditions(degree=2)
u0 = Function(V)
u0 = interpolate(indata, V)

# Define test functions
v_1, v_2, v_3 = TestFunction(V)

# Define functions for velocity and concentrations
u = Function(V)
u_n = Function(V)

# Define the integrals
M0 = u[0] * dx
M1 = u[1] * dx
M2 = u[2] * dx

# Split system functions to access components
u_1 = u[0]
u_2 = u[1]
u_3 = u[2]
u_n1 = u_n[0]
u_n2 = u_n[1]
u_n3 = u_n[2]


# Non linear functions
def S_u(u_n1, u_n2):
    return (alpha* u_n1 ** 2)/(L_0 + l*u_n2)

def S_v(u_n1, u_n2, u_n3):
    return (beta * u_n2**2) + (u_n2*u_n3)/(alpha + u_n2 + m*u_n1)

def S_w(u_n1, u_n2, u_n3):
    return (zeta*u_n2*u_n3)/(alpha + u_n2 + m*u_n1)

# Define source terms
f_1 = Constant(0)
f_2 = Constant(0)
f_3 = Constant(0)

# Define expressions used in variational forms
k = Constant(dt)
delta1 = Constant(delta1)
delta2 = Constant(delta2)
delta3 = Constant(delta3)
alpha = Constant(alpha)
beta = Constant(beta)
gamma = Constant(gamma)

# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx + delta1*(inner(grad(v_1), grad(u_1)) + inner(grad(v_1), grad(u_n1))/2)*dx \
  - alpha*((inner(u_1, v_1) + inner(u_n1, v_1))/2)*dx + S_u(u_n1, u_n2)*inner(u_n1, v_1)*dx  \
  + ((u_2 - u_n2) / k)*v_2*dx + delta2*(inner(grad(v_2), grad(u_2)) + inner(grad(v_2), grad(u_n2))/2)*dx \
  - beta*((inner(u_2, v_2) + inner(u_n2, v_2))/2)*dx + S_v(u_n1, u_n2, u_n3)*inner(u_n2, v_2)*dx  \
  + ((u_3 - u_n3) / k)*v_3*dx + delta3*(inner(grad(v_3), grad(u_3)) + inner(grad(v_3), grad(u_n3))/2)*dx \
  + gamma*((inner(u_3, v_3) + inner(u_n3, v_3))/2)*dx - S_w(u_n1, u_n2, u_n3)*inner(u_n3, v_3)*dx  \
  - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx


# Create VTK files for visualization output
vtkfile_u = File('test/u.pvd')
vtkfile_v = File('test/v.pvd')
vtkfile_w = File('test/w.pvd')

# Compute the functional
population_u = assemble(M0)
population_v = assemble(M1)
population_w = assemble(M2)

# Time-stepping
t = 0
while t < T:

    # Update current time
    t += dt

    # Solve variational problem for time step
    solve(F == 0, u)

    if t in [100, 200, 300, 400, 500]:
        # Save solution to file (VTK)
        _u_1, _v_2, _w_3 = u.split()
        vtkfile_u << (_u_1, t)
        vtkfile_v << (_v_2, t)
        vtkfile_w << (_w_3, t)

    # Update previous solution
    u_n.assign(u)