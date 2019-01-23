from fenics import *
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

funnel = sys.argv[1]
ez = sys.argv[2]

mesh = Mesh('meshes/' + funnel + '_spine.xml')

coor = mesh.coordinates()
scale_factor = 1e-6
coor[:,:] *= scale_factor

min_z = coor[:,2].min()

class Absorb(SubDomain):
    def inside(self, x, on_boundary, min_z=min_z):
        return on_boundary and x[2] < min_z + DOLFIN_EPS

sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
sub_domains.set_all(0)
absorb = Absorb()
absorb.mark(sub_domains, 1)

P1 = FiniteElement('P', mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, P1)

bc = DirichletBC(V, Constant(0), sub_domains, 1)

E = Expression(("0", "0", ez), degree=4)
F = -E
R = 8.314   # Gas constant constant, J/(K*mol)
Fa = 9.648e4 # Faradays constant, C/mol
T = 300     # Temperature,  Kelvin
z = 2.
D = 4e-10
gamma = Constant(R*T/(D*Fa*z))

print(1/(100*R*T/Fa))

u = TrialFunction(V)
v = TestFunction(V)

form = (-D*inner(grad(u), grad(v)) - inner(F/gamma, grad(u))*v + 1*v)*dx

a, L = lhs(form), rhs(form)
u = Function(V)
solve(a==L, u, bc)
# f = plot(u)
# plt.colorbar(f)
# plt.show()
mfpt = u(Point(0,0,0))

with open('results/fenics/' + funnel + '.csv','a+') as fd:
    fd.write(ez + ',' + str(mfpt))
