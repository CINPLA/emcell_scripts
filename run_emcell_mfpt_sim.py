import pymcell as m
import sys
import os

funnel = sys.argv[1]
#reenter = bool(int(sys.argv[2]))
ez = float(sys.argv[2])
dt = float(sys.argv[3])

f_n = 'results/emcell_' + funnel + '_E' + sys.argv[2] + '_dt' + sys.argv[3] + '/'

if os.path.isdir(f_n):
    print("Simulation is already done!")
    sys.exit()

json_file_path = 'data_models/narrow_escape_' + funnel + '.json'

def new_electric_field(x, y, z, ez=ez):
    # Efield should be volt per meter
    Ex = 0
    Ey = 0
    Ez = ez
    return (Ex, Ey, Ez)

def create_spherical_release_sites_from_dm(data_model, world, species_dict):
    rel_site_dm_list = \
        data_model['mcell']['release_sites']['release_site_list']
    for rel_site_dm in rel_site_dm_list:
        if rel_site_dm['shape'] == 'SPHERICAL':
            pos = create_vector3(float(rel_site_dm['location_x']),
                                 float(rel_site_dm['location_y']),
                                 float(rel_site_dm['location_z']))
            diam = create_vector3(float(rel_site_dm['site_diameter']),
                                  float(rel_site_dm['site_diameter']),
                                  float(rel_site_dm['site_diameter']))
            shape = 'spherical'
            quantity = int(rel_site_dm['quantity'])
            vm = rel_site_dm['molecule']
            species = species_dict[vm]
            world.create_release_site(species, quantity, shape, pos, diam)


def create_vector3(x=0., y=0., z=0.):
    v = m.vector3()
    v.x = x
    v.y = y
    v.z = z
    return v

m.initialize_electric_field()
m.update_electric_field(new_electric_field)

world = m.MCellSim(1)
world.set_time_step(dt)
world.set_iterations(iterations=int(2/dt))
world.set_output_freq(1000)
data_model = m.read_json_data_model(json_file_path)

species_dict = m.create_species_from_dm(data_model)
species_dict['vm'].z = 2
for spec in species_dict:
    world.add_single_species(species_dict[spec])
meshobj_dict = m.create_meshobjs_from_dm(data_model)
world.add_geometry(meshobj_dict['spine'])

surface_classes = m.create_surface_classes_from_dm(data_model, world, species_dict)
absorb = m.SurfaceClass(m.SC(3), species_dict['vm'], 'abs')

create_spherical_release_sites_from_dm(data_model, world, species_dict)
world.assign_surf_class(absorb, meshobj_dict['spine'].regions[0])

world.add_count(species_dict['vm'], folder_name=f_n)
world.add_trigger(species_dict['vm'], mesh_obj=meshobj_dict['spine'], folder_name=f_n)

world.add_partitions('x', -0.5, 0.5, 0.02)
world.add_partitions('y', -0.5, 0.5, 0.02)
world.add_partitions('z', -1.5, 0.5, 0.02)

for i in range(world._iterations):
    world.run_iteration()
world.end_sim()
