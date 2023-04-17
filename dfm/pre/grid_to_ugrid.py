#!/usr/bin/python
import numpy as np
import netCDF4 as nc
import pickle
from shapely.geometry import Point, Polygon, LineString

from tqdm import tqdm


### Meshgrid creation and calculations for Ugrid-NC definition

## Define grid parameters
nx, ny = (110+1, 100+1)
lb_x, ub_x = (0, 110)
lb_y, ub_y = (0, 100)
slope_x, slope_y = (0.05, 0.02)
dx, dy = (ub_x-lb_x)/(nx-1), (ub_y-lb_y)/(ny-1)

# Generate node coordinates
x = np.linspace(lb_x, ub_x, nx)
y = np.linspace(lb_y, ub_y, ny)
xx, yy = np.meshgrid(x, y, indexing='ij')

# Generate zz coordinates
zz = np.zeros((nx, ny))
nch = 9
assert (nx-nch)%2 == 0, "adjust the number of cells for the channel so that the y-domain is equally divisible by the channel in the middle"
cs = np.concatenate([np.linspace((dx*slope_x*((nx-nch)//2)), 0, (nx-nch)//2), np.full(nch, 0), np.linspace(0, (dx*slope_x*((nx-nch)//2)), (nx-nch)//2)])
for j in range(ny):
    zz[:,j] = cs + slope_y*dy*j
datum = 100
zz = zz + datum
assert zz.min() > 10, "GWSWEX does not support -ve elevations"


## calculations for Ugrid-NC definition

# Generate node points
node_points = [Point(xx[i, j], yy[i, j], zz[i, j]) for i in range(nx) for j in range(ny)]

# Generate face polygons and node indices
face_polygons, face_nodes, face_coords_x, face_coords_y = [], [], [], []
print("Pass (1/4)")
for i in tqdm(range(nx-1)):
    for j in range(ny-1):
        idx1, idx2, idx3, idx4 = i*ny+j, i*ny+j+1, (i+1)*ny+j+1, (i+1)*ny+j
        poly_corners = [node_points[idx1], node_points[idx2], node_points[idx3], node_points[idx4]]
        face_nodes.append([idx1, idx2, idx3, idx4])
        face_polygons.append(Polygon(poly_corners))
        face_coords_x.append([p.x for p in poly_corners])
        face_coords_y.append([p.y for p in poly_corners])

# Generate edge points and maps
print("Pass (2/4)")
edge_points, edge_node_map, edge_poly_map = [], {}, {}
for face_idx, face in enumerate(tqdm(face_polygons)):
    for i in range(len(face.exterior.coords) - 1):
        start_point, end_point = Point(face.exterior.coords[i]), Point(face.exterior.coords[i+1])
        if (start_point, end_point) in edge_points or (end_point, start_point) in edge_points:
            continue
        edge_points.append((start_point, end_point))
        start_idx, end_idx = node_points.index(start_point), node_points.index(end_point)
        edge_node_map[(start_point, end_point)] = (start_idx, end_idx)
        edge_poly_map.setdefault(start_idx, []).append(face_idx)
        edge_poly_map.setdefault(end_idx, []).append(face_idx)

# Generate edge-face connectivity
print("Pass (3/4)")
edge_faces = [list(set(edge_poly_map[edge[0]]).intersection(set(edge_poly_map[edge[1]]))) for edge in tqdm(edge_node_map.values())]
edge_faces = [sublist + [-1] * (2 - len(sublist)) for sublist in edge_faces]

# Generate face centroids
print("Pass (4/4)")
face_centroids = [face.centroid for face in tqdm(face_polygons)]

# Calculate number of nodes, edges, and faces
nnodes, nedges, nfaces = nx * ny, len(edge_points), len(face_polygons)
nfaces_max = max(len(face) for face in face_nodes)



### start Ugrid-NC file creation
gridNC = nc.Dataset('../tilted-v_net.nc', 'w', format='NETCDF4')
gridDict = {}


## create dimensions
gridDict['D'] = {'Two': 2, 'mesh2d_nEdges': nedges, 'mesh2d_nNodes': nnodes, 'mesh2d_nFaces': nfaces, 'mesh2d_nMax_face_nodes': nfaces_max}

gridDict['dims'] = {}
for dname, dval in gridDict['D'].items():
    gridDict['dims'][dname] = gridNC.createDimension(dname, dval)


## create variables
gridDict['V'] = {}
gridDict['V']['projected_coordinate_system'] = ['i', ()]
gridDict['V']['mesh2d'] = ['i', ()]
gridDict['V']['mesh2d_node_x'] = ['f8', ('mesh2d_nNodes')]
gridDict['V']['mesh2d_node_y'] = ['f8', ('mesh2d_nNodes')]
gridDict['V']['mesh2d_node_z'] = ['f8', ('mesh2d_nNodes'), -999]
gridDict['V']['mesh2d_edge_x'] = ['f8', ('mesh2d_nEdges')]
gridDict['V']['mesh2d_edge_y'] = ['f8', ('mesh2d_nEdges')]
gridDict['V']['mesh2d_edge_nodes'] = ['i', ('mesh2d_nEdges','Two')]
gridDict['V']['mesh2d_face_nodes'] = ['i', ('mesh2d_nFaces', 'mesh2d_nMax_face_nodes'), -999]
gridDict['V']['mesh2d_edge_faces'] = ['i', ('mesh2d_nEdges', 'Two'), -999]
gridDict['V']['mesh2d_face_x'] = ['f8', ('mesh2d_nFaces')]
gridDict['V']['mesh2d_face_y'] = ['f8', ('mesh2d_nFaces')]
gridDict['V']['mesh2d_face_x_bnd'] = ['f8', ('mesh2d_nFaces', 'mesh2d_nMax_face_nodes'), -999]
gridDict['V']['mesh2d_face_y_bnd'] = ['f8', ('mesh2d_nFaces', 'mesh2d_nMax_face_nodes'), -999]

gridDict['vars'] = {}
for vname, attr in gridDict['V'].items():
    if len(attr) == 2:
        gridDict['vars'][vname] = gridNC.createVariable(vname, attr[0], attr[1])
    elif len(attr) == 3:
        gridDict['vars'][vname] = gridNC.createVariable(vname, attr[0], attr[1], fill_value=attr[2])


## create global attributes
gridNC.Conventions = 'CF-1.8 UGRID-1.0 Deltares-0.10'


## create variable attributes
gridDict['vars']['projected_coordinate_system']._name = "Unknown projected"
gridDict['vars']['projected_coordinate_system'].epsg = 0
gridDict['vars']['projected_coordinate_system'].grid_mapping_name = "Unknown projected"
gridDict['vars']['projected_coordinate_system'].longitude_of_prime_meridian = 0.0
gridDict['vars']['projected_coordinate_system'].EPSG_code = "EPSG:0"
gridDict['vars']['projected_coordinate_system'].value = "value is equal to EPSG code"
gridDict['vars']['projected_coordinate_system'].proj4_params = ""

gridDict['vars']['mesh2d'].cf_role = "mesh_topology"
gridDict['vars']['mesh2d'].long_name = "Topology data of 2D mesh"
gridDict['vars']['mesh2d'].topology_dimension = 2
gridDict['vars']['mesh2d'].node_coordinates = "mesh2d_node_x mesh2d_node_y"
gridDict['vars']['mesh2d'].node_dimension = "mesh2d_nNodes"
gridDict['vars']['mesh2d'].max_face_nodes_dimension = "mesh2d_nMax_face_nodes"
gridDict['vars']['mesh2d'].edge_node_connectivity = "mesh2d_edge_nodes"
gridDict['vars']['mesh2d'].edge_dimension = "mesh2d_nEdges"
gridDict['vars']['mesh2d'].edge_coordinates = "mesh2d_edge_x mesh2d_edge_y"
gridDict['vars']['mesh2d'].face_node_connectivity = "mesh2d_face_nodes"
gridDict['vars']['mesh2d'].face_dimension = "mesh2d_nFaces"
gridDict['vars']['mesh2d'].edge_face_connectivity = "mesh2d_edge_faces"
gridDict['vars']['mesh2d'].face_coordinates = "mesh2d_face_x mesh2d_face_y"

gridDict['vars']['mesh2d_node_x'].units = "m"
gridDict['vars']['mesh2d_node_x'].standard_name = "projection_x_coordinate"
gridDict['vars']['mesh2d_node_x'].long_name = "x-coordinate of mesh nodes"

gridDict['vars']['mesh2d_node_y'].units = "m"
gridDict['vars']['mesh2d_node_y'].standard_name = "projection_y_coordinate"
gridDict['vars']['mesh2d_node_y'].long_name = "y-coordinate of mesh nodes"

gridDict['vars']['mesh2d_node_z'].units = "m"
gridDict['vars']['mesh2d_node_z'].standard_name = "projection_z_coordinate"
gridDict['vars']['mesh2d_node_z'].long_name = "z-coordinate of mesh nodes"
gridDict['vars']['mesh2d_node_z'].mesh = "mesh2d"
gridDict['vars']['mesh2d_node_z'].location = "node"
gridDict['vars']['mesh2d_node_z'].coordinates = "mesh2d_node_x mesh2d_node_y"
gridDict['vars']['mesh2d_node_z'].grid_mapping = "projected_coordinate_system"

gridDict['vars']['mesh2d_edge_x'].units = "m"
gridDict['vars']['mesh2d_edge_x'].standard_name = "projection_x_coordinate"
gridDict['vars']['mesh2d_edge_x'].long_name = "characteristic x-coordinate of the mesh edge (e.g. midpoint)"

gridDict['vars']['mesh2d_edge_y'].units = "m"
gridDict['vars']['mesh2d_edge_y'].standard_name = "projection_y_coordinate"
gridDict['vars']['mesh2d_edge_y'].long_name = "characteristic y-coordinate of the mesh edge (e.g. midpoint)"

gridDict['vars']['mesh2d_edge_nodes'].long_name = "Start and end nodes of mesh edges"
gridDict['vars']['mesh2d_edge_nodes'].cf_role = "edge_node_connectivity"
gridDict['vars']['mesh2d_edge_nodes'].start_index = 1

gridDict['vars']['mesh2d_face_nodes'].long_name = "Vertex nodes of mesh faces (counterclockwise)"
gridDict['vars']['mesh2d_face_nodes'].cf_role = "face_node_connectivity"
gridDict['vars']['mesh2d_face_nodes'].start_index = 1

gridDict['vars']['mesh2d_edge_faces'].long_name = "Neighboring faces of mesh edges"
gridDict['vars']['mesh2d_edge_faces'].cf_role = "edge_face_connectivity"
gridDict['vars']['mesh2d_edge_faces'].start_index = 1

gridDict['vars']['mesh2d_face_x'].units = "m"
gridDict['vars']['mesh2d_face_x'].standard_name = "projection_x_coordinate"
gridDict['vars']['mesh2d_face_x'].long_name = "characteristic x-coordinate of the mesh face (e.g. centroid)"
gridDict['vars']['mesh2d_face_x'].bounds = "mesh2d_face_x_bnd"

gridDict['vars']['mesh2d_face_y'].units = "m"
gridDict['vars']['mesh2d_face_y'].standard_name = "projection_y_coordinate"
gridDict['vars']['mesh2d_face_y'].long_name = "characteristic y-coordinate of the mesh face (e.g. centroid)"
gridDict['vars']['mesh2d_face_y'].bounds = "mesh2d_face_y_bnd"

gridDict['vars']['mesh2d_face_x_bnd'].units = "m"
gridDict['vars']['mesh2d_face_x_bnd'].standard_name = "projection_x_coordinate"
gridDict['vars']['mesh2d_face_x_bnd'].long_name = "x-coordinate bounds of mesh faces (i.e. corner coordinates)"

gridDict['vars']['mesh2d_face_y_bnd'].units = "m"
gridDict['vars']['mesh2d_face_y_bnd'].standard_name = "projection_y_coordinate"
gridDict['vars']['mesh2d_face_y_bnd'].long_name = "y-coordinate bounds of mesh faces (i.e. corner coordinates)"


## write mesh data to NC
gridNC['mesh2d_node_x'][:] = [p.x for p in node_points]
gridNC['mesh2d_node_y'][:] = [p.y for p in node_points]
gridNC['mesh2d_node_z'][:] = [p.z for p in node_points]

gridNC['mesh2d_edge_nodes'][:] = np.array(list(edge_node_map.values())) + 1

gridNC['mesh2d_face_nodes'][:] = np.array(face_nodes) + 1

gridNC['mesh2d_face_x_bnd'][:] = np.array(face_coords_x)
gridNC['mesh2d_face_y_bnd'][:] = np.array(face_coords_y)

gridNC['mesh2d_face_x'][:] = np.array([p.x for p in face_centroids])
gridNC['mesh2d_face_y'][:] = np.array([p.y for p in face_centroids])

gridNC['mesh2d_edge_faces'][:] = np.array(edge_faces) + 1

gridNC['mesh2d_edge_x'][:] = np.array([p.x for p in [LineString([edge[0],edge[1]]).interpolate(0.5) for edge in edge_points]])
gridNC['mesh2d_edge_y'][:] = np.array([p.y for p in [LineString([edge[0],edge[1]]).interpolate(0.5) for edge in edge_points]])


## pickle the mesh data for outside use
pkl_dict = {'nx':nx, 'ny': ny, 'lb_x': lb_x, 'lb_y': lb_y, 'ub_x': ub_x, 'ub_y': ub_y, 'dx':dx, 'dy': dy, \
            'xx': xx, 'yy': yy, 'zz': zz, 'node_points': node_points, 'edge_points': edge_points, 'face_polygons':face_polygons, 'face_centroids': face_centroids, \
            'face_nodes': face_nodes, 'edge_faces': edge_faces, 'face_coords_x': face_coords_x, 'face_coords_y': face_coords_y}

with open('../../common/data/mesh.pkl', 'wb') as file:
    pickle.dump(pkl_dict, file, pickle.HIGHEST_PROTOCOL)
    
gridNC.close()