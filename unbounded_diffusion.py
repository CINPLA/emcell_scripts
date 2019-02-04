import pymcell as m
import sys
import os

for mode in [1,2]:
    if mode==1:
        N = 2
        dt = 1e-4
    else:
        N = 101
        dt = 1e-6

    ex = 1e5
    ey = 3e5
    ez = -2e5

    def new_electric_field(x, y, z, ex=ex, ey=ey, ez=ez):
        # Efield should be volt per meter
        Ex = ex
        Ey = ey
        Ez = ez
        return (Ex, Ey, Ez)

    def create_vector3(x=0., y=0., z=0.):
        v = m.vector3()
        v.x = x
        v.y = y
        v.z = z
        return v

    m.initialize_electric_field()
    m.update_electric_field(new_electric_field)

    world = m.MCellSim(2)
    world.set_time_step(dt)
    world.set_iterations(iterations=N)

    spec = m.Species("Ca", 4e-6)
    spec.z = 2
    world.add_species(spec)
    shape = 'spherical'
    quantity = 10000
    diam = create_vector3()
    pos = create_vector3()
    world.create_release_site(spec, quantity, shape, pos, diam)

    f_n = 'results/output_' + str(mode) + '/'
    world.add_count(spec, folder_name=f_n)
    world.add_viz([spec], step=1, folder_name=f_n)
    for i in range(world._iterations):
        world.run_iteration()
    world.end_sim()
