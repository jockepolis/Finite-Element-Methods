from __future__ import print_function
from fenics import *
import csv

# Define parameters :
T = 500
dt = 0.5
delta1 = 1
delta2 = 1
delta3 = 1
alpha = 0.4
beta = 1
gamma = 0.8
zeta = 2
L_0 = 0.4
l = 0.4
m = 0.12

# Define expressions used in variational forms
k = Constant(dt)
delta1 = Constant(delta1)
delta2 = Constant(delta2)
delta3 = Constant(delta3)
alpha = Constant(alpha)
beta = Constant(beta)
gamma = Constant(gamma)
zeta = Constant(zeta)
L_0 = Constant(L_0)
l = Constant(l)
m = Constant(m)

# Create mesh and define function space
mesh = Mesh("circle.xml.gz")

# Construct the finite element space
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P1 * P1 * P1
V = FunctionSpace(mesh, TH)

# Define test functions
v = TestFunction(V)
v_1 = v[0]
v_2 = v[1]
v_3 = v[2]

# Define function
u = Function(V)

# Class representing the intial conditions
class InitialConditions(UserExpression):
    
    def eval(self, values, x):
        values[0] = 0
        values[1] = (4/15)-2*pow(10,-7)*(x[0]-0.1*x[1]-350)*(x[0]-0.1*x[1]-67)
        values[2] = (22/45)-3*pow(10,-5)*(x[0]-450)-1.2*pow(10,-4)*(x[1]-15)
    
    def value_shape(self):
        return (3,)

# Define initial condition
indata = InitialConditions(degree=2)
u0 = Function(V)
u0 = interpolate(indata, V)

# Split system functions to access components
# u_1, u_2, u_3 = split(u)
# u_n1, u_n2, u_n3 = split(u0)
u_1 = u[0]
u_2 = u[1]
u_3 = u[2]
u_n1 = u0[0]
u_n2 = u0[1]
u_n3 = u0[2]

# Define the integrals
M0 = u_1 * dx
M1 = u_2 * dx
M2 = u_3 * dx

# Define source terms
f_1 = Constant(0)
f_2 = Constant(0)
f_3 = Constant(0)

# Non linear functions
def S_u(u_n1, u_n2):
    return (alpha* u_n1 ** 2)/(L_0 + l*u_n2)

def S_v(u_n1, u_n2, u_n3):
    return (beta * u_n2**2) + (u_n2*u_n3)/(alpha + u_n2 + m*u_n1)

def S_w(u_n1, u_n2, u_n3):
    return (zeta*u_n2*u_n3)/(alpha + u_n2 + m*u_n1)


# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx + delta1*(dot(grad(v_1), grad(u_1)) + dot(grad(v_1), grad(u_n1))/2)*dx \
  - alpha*((u_1 + u_n1) / 2)*v_1*dx + S_u(u_n1, u_n2)*u_n1*v_1*dx  \
  + ((u_2 - u_n2) / k)*v_2*dx + delta2*(dot(grad(v_2), grad(u_2)) + dot(grad(v_2), grad(u_n2))/2)*dx \
  - beta*((u_2 + u_n2) / 2)*v_2*dx + S_v(u_n1, u_n2, u_n3)*u_n2*v_2*dx  \
  + ((u_3 - u_n3) / k)*v_3*dx + delta3*(dot(grad(v_3), grad(u_3)) + dot(grad(v_3), grad(u_n3))/2)*dx \
  + gamma*((u_3 + u_n3) / 2)*v_3*dx - S_w(u_n1, u_n2, u_n3)*u_n3*v_3*dx  \
  - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx

# Create VTK files for visualization output
vtkfile_u = File('LotkaVolterraTest/u.pvd')
vtkfile_u0 = File('LotkaVolterraTest/u0.pvd')

# Time-stepping
t = 0
vtkfile_u0 << (u0, t)
u.assign(u0)    
popu = []
popv = []
popw = []
while t < T:
    
    if t in [100, 200, 300, 400]:
    # Save solution to file (VTK)
        vtkfile_u << (u, t)

    population_u = assemble(M0)
    population_v = assemble(M1)
    population_w = assemble(M2)
    popu.append(population_u)
    popv.append(population_v)
    popw.append(population_w)
        
    # Update current time
    t += dt
    
    # Solve variational problem for time step
    solve(F == 0, u)
    
    # Update previous solution
    u0.assign(u)
    
print(u[0])
with open('populationstest.csv', 'w') as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow(popv)
    csv_writer.writerow(popw)