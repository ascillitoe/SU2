import pyvista as pv
import pandas as pd
import numpy as np
from tqdm import tqdm 

# This code converts data from a vtk file to a SU2 restart file, for use in a data-driven RANS version of SU2.
# It expects a (nnode,6) array named "deltas" in the vtk file. The vtk file containing the data, the target vtk file
# from the SU2 run, and the target SU2 csv file must be given.
# Author: Ashley Scillitoe (2020), ashleyscillitoe@googlemail.com.

# User inputs
tol = 1e-4 
data    = pv.read('2a_Y.vtk')  # Source vtk file
su2_vtk = pv.read('flow.vtu')  # Target vtk file (i.e. what SU2 writes out)
su2_pd  = pd.read_csv('solution_flow.csv') # Target csv file (also written out by SU2)

# Read files
ndata    = data.number_of_points
nsu2_vtk = su2_vtk.number_of_points
nsu2_pd  = len(su2_pd)

# Read data from ML vtk file
deltas_vtk = data.point_arrays['deltas']

# Get coordinates from files
data_coords = data.points[:,:2]
su2_vtk_coords = su2_vtk.points[:,:2]
su2_pd_coords = np.array([su2_pd['x'].to_numpy(),su2_pd['y'].to_numpy()]).T

# Interp from ML vtk file to su2 csv file
print('Interpolating from ML vtk file to su2 csv file...')
deltas = np.zeros([nsu2_pd,6])
count = 0
max_dist = 0.0
for row, point in enumerate(tqdm(su2_pd_coords)):
    dist = np.sqrt((data_coords[:,0]-point[0])**2 + (data_coords[:,1]-point[1])**2)
    idx = np.argmin(dist)
    deltas[row,:] = deltas_vtk[idx,:]

    max_dist = max(max_dist,dist[idx])
    if dist[idx] > tol: count += 1
print('Number of points above tolerance: ', count)
print('Max distance: ', max_dist)

# Interp from su2 csv file onto su2_vtk for debugging 
print('Interpolating su2 csv file to su2 vtk file for debugging...')
delta_vtk = np.zeros([nsu2_vtk,6])
count = 0
for row, point in enumerate(tqdm(su2_vtk_coords)):
    dist = np.sqrt((su2_pd_coords[:,0]-point[0])**2.0 + (su2_pd_coords[:,1]-point[1])**2.0)
    idx = np.argmin(dist)
    delta_vtk[row,:] = deltas[idx,:]

    if dist[idx] > tol: count += 1
print('Number of points above tolerance: ', count)
su2_vtk.point_arrays['delta'] = delta_vtk
su2_vtk.save('debug.vtk')

# Save new su2 restart file
print('Saving su2 restart file')
su2_pd.insert(8, 'delta xc',   deltas[:,0])
su2_pd.insert(9, 'delta yc',   deltas[:,1])
su2_pd.insert(10, 'delta h2',  deltas[:,2])
su2_pd.insert(11, 'delta h3',  deltas[:,3])
su2_pd.insert(12, 'delta h4',  deltas[:,4])
su2_pd.insert(13, 'delta logk',deltas[:,5])
su2_pd.to_csv('restart_flow.csv',index=False)

