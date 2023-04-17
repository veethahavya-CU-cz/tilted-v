#!/usr/bin/python
import numpy as np
import netCDF4 as nc
import pickle
import flopy as fp
import matplotlib.pyplot as plt


def flatten(xs):
    if isinstance(xs, (list, tuple)):
        return [x for sublist in xs for x in flatten(sublist)]
    else:
        return [xs]

def getDisV(mapNC_path='../../dfm/out/tilted-v_map.nc'):
    mapNC = nc.Dataset(mapNC_path)
    ugrid_dict = {}

    ugrid_dict['ne'] = mapNC.dimensions["mesh2d_nFaces"].size
    ugrid_dict['ex'] = np.array(mapNC["mesh2d_face_x"][:])
    ugrid_dict['ey'] = np.array(mapNC["mesh2d_face_y"][:])
    ugrid_dict['ez'] = np.array(mapNC["mesh2d_flowelem_bl"][:])
    ugrid_dict['ea'] = np.array(mapNC["mesh2d_flowelem_ba"][:])
    elems_ma_array = mapNC["mesh2d_face_nodes"][:]

    ugrid_dict['nn'] = mapNC.dimensions["mesh2d_nNodes"].size
    ugrid_dict['nx'] = np.array(mapNC["mesh2d_node_x"][:])
    ugrid_dict['ny'] = np.array(mapNC["mesh2d_node_y"][:])
    ugrid_dict['nz'] = np.array(mapNC["mesh2d_node_z"][:])

    mapNC.close()

    verteces = []
    for idx, (x,y) in enumerate(zip(ugrid_dict['nx'], ugrid_dict['ny'])):
        verteces.append([idx, x, y])
    
    elemental_nodes_ma = elems_ma_array - 1
    elemental_nodes = [en.compressed().tolist() for en in elemental_nodes_ma]

    cell2d = []
    for idx, (x,y) in enumerate(zip(ugrid_dict['ex'], ugrid_dict['ey'])):
        cell2d.append(flatten([idx, x, y, len(elemental_nodes[idx]), elemental_nodes[idx]]))

    return verteces, cell2d, ugrid_dict


name = 'tilted-v'
sim = fp.mf6.MFSimulation(sim_name=name, exe_name='/opt/mf6/bin/mf6.exe', version="mf6", sim_ws='../tilted-v')
tdis = fp.mf6.ModflowTdis(sim, time_units="seconds", nper=1, perioddata=[1, 1, 1.0]) # dummy
ims = fp.mf6.ModflowIms(sim)
gwf = fp.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

verteces, cell2d, ugrid_dict = getDisV()
bot = ugrid_dict['ez'] - 5
disv = fp.mf6.ModflowGwfdisv(gwf, length_units="meters", nlay=1, top=ugrid_dict['ez'], botm=bot, ncpl=ugrid_dict['ne'], nvert=ugrid_dict['nn'], vertices=verteces, cell2d=cell2d)
sim.write_simulation()

with open('../../common/data/ugrid.pkl', 'wb') as outfile:
    pickle.dump(ugrid_dict, outfile, pickle.HIGHEST_PROTOCOL)