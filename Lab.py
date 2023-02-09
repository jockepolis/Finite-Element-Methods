from dolfin import *
import time

t = time.time()
# Create mesh and define function space
mesh = UnitSquareMesh(500, 500)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("8 * pow(pi, 2) * sin(2*pi*x[0]) * sin(2*pi*x[1])", degree=2)
g = Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Save solution in VTK format
file = File("poisson_new.pvd")
file << u

elapsed = time.time() - t
info ("My code took %f seconds " %elapsed)

# Plot solution
import matplotlib.pyplot as plt
plot(u)
plt.show()