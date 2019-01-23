from fenics import *
# import matplotlib.pyplot as plt
import numpy as np


mesh_list = ['funnel_spine',
             'sharp_spine']

e_list = np.linspace(-1.1e5, 3e4, 101)


mfpt_table = np.zeros([len(mesh_list), len(e_list)])
volume_table = np.zeros(len(mesh_list))

R = 8.314   # Rayleighs constant, J/(K*mol)
Fa = 9.648e4 # Faradays constant, C/mol
T = 300     # Temperature,  Kelvin
z = 2.
D = 4e-10
gamma = Constant(R*T/(D*Fa*z))

idx = 0
for mesh_idx, mesh_name in enumerate(mesh_list):
    mesh = Mesh('meshes/' + mesh_name + '.xml')
    mesh = refine(mesh)
    # mesh = UnitSquareMesh(10,10)
    scale_factor = 1e-6
    coor = mesh.coordinates()
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
    volume_table[mesh_idx] = assemble(project(1, V)*dx)

    for e_idx, e in enumerate(e_list):
        print idx
        idx += 1
        E = Expression(("0", "0", "e"), e=e, degree=4)
        F = -E

        u = TrialFunction(V)
        v = TestFunction(V)

        form = (-D*inner(grad(u), grad(v)) - inner(F/gamma, grad(u))*v + 1*v)*dx

        a, L = lhs(form), rhs(form)
        u = Function(V)
        solve(a==L, u, bc)
        file = File('results/' + mesh_name + '_sol_e' + str(e) + '.xml')

        file << u

        if 'neck' in mesh_name:
            z_p = coor[:,2].max()
        else:
            z_p = 0
        mfpt = u(Point(0, 0, z_p))
        mfpt_table[mesh_idx, e_idx] = mfpt

print(mfpt_table)
np.save('results/' + 'mfpt_table', mfpt_table)
np.save('results/' + 'volume_table', volume_table)
