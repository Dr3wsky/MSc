"""
This script performs data analysis on experimental and simulated data from OpenFOAM.
The OpenFOAM data is extracted from Paraview using plotOverLine method along the axis of a confined jet.
The python scrtipt herin loads data files, processes the information, and generates plots of the specified variable

Note: Ensure the necessary data files are available in the specified file paths before running this script.

Author: Andrew McFadden
Purpose: MSc. Data analysis
Date: 02/11/2023
Version: V1
"""

# Load external libraries and modules
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
# import data_handling
# from scipy.io import loadmat

# Assign file location and simulations
dir_home = 'C:\\Users\\drewm\\Documents\\2.0 MSc\\2.0 Simulations\\'
div = '//'
sim_locs = ['PIV_JetSim_M1','PIV_JetSim_M2', 'PIV_JetSim_M3', 'PIV_JetSim_M4']
sim_names = ['M1', 'M2', 'M3', 'M4']

# Load experimental Data
with open('exp_data.pkl', 'rb') as file:
    exp_data = pickle.load(file)
    file.close()

sim_data = {}
# Load simulation data
for i in range(len(sim_names)):
    sim_data.update({sim_names[i]: pd.read_csv(dir_home + sim_locs[i] + '\\axial_data.csv')})


    
    

# # Load k-w data
# loc_kw = f'C:/Users/drewm/Documents/2.0 MSc/2.0 Simulations/PIV_JetSim_M{mesh}/'
# type_kw = 'axial_data'
# data_kw = pd.read_csv(loc_kw + type_kw + '.csv')

# kw_tke_dimless = data_kw['kMean'] / (data_kw['UMean_0'] ** 2)

# # Load RSM data
# loc_rsm = f'C:/Users/drewm/Documents/2.0 MSc/2.0 Simulations/RSM_{model}/'
# type_rsm = 'axial_data'
# data_rsm = pd.read_csv(loc_rsm + type_rsm + '.csv')

# rsm_tke_dimless = data_rsm['kMean'] / (data_rsm['UMean_0'] ** 2)
# rsm_tke_dimless_unsteady = data_rsm['k'] / (data_rsm['U_0'] ** 2)

# # Other properties and data handling
# xd = (data_rsm['Points_0'] - 0.33533) / 0.0064
# start = np.where(xd > 0)[0][0]
# fin = np.where(xd > 40)[0][0]
# var = 'UMean_1'

# # Find starting point in experimental data
# exp_xd = exp_data['exp_x'][0] / 6.4
# exp_start = np.where(exp_data['exp_u_centreline'] == 0)[0][0]

# # Non-dimensionalize variables
# exp_u_centreline_nonDim = exp_data['exp_u_centreline'] / max(exp_data['exp_u_centreline'])
# kw_U_nonDim = data_kw[var] / np.mean(data_kw[var][:37])
# rsm_U_nonDim = data_rsm[var] / np.mean(data_rsm[var][:37])

# # Plotting
# colormap = ["#e74c3c", "#e67e22", "#2ecc71", "#8e44ad", "#3498db", "#f1c40f"]
# linetype = ["-", "--", ":", "-."]

# # Plot data
# plt.figure(1)
# plt.plot(xd, kw_tke_dimless, linetype[0], color=colormap[3], linewidth=1.25)
# plt.plot(xd, rsm_tke_dimless, linetype[1], color=colormap[4], linewidth=1.25)
# plt.plot(xd, rsm_tke_dimless_unsteady, linetype[2], color=colormap[1], linewidth=1.25)

# plt.grid(True)
# plt.grid(True, which='minor')
# plt.xlim([0, 100])
# plt.ylabel('$k/{u^2}_{jet}$', fontsize=14)
# plt.xlabel('$x/D_{jet}$', fontsize=14)
# plt.title('Turbulent Kinetic Energy Profile', fontsize=14)
# plt.legend(['k-w SST, M4', '7-equation RSM, LRR', '7-equation RSM-unsteady, LRR'], fontsize=14)

# plt.show()
