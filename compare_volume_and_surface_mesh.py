import pymcell as m
from fenics import *
import numpy as np

json_file_path = 'data_models/narrow_escape_sharp.json'
data_model = m.read_json_data_model(json_file_path)

meshobj_dict = m.create_meshobjs_from_dm(data_model)
spine = meshobj_dict['spine']
vert_list = np.array(spine.vert_list)
print(vert_list[:,2].min())

mesh = Mesh('meshes/sharp_spine.xml')
coor = mesh.coordinates()
print(coor[:,2].min())
