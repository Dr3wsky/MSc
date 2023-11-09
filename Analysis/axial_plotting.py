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
# import os
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

# Load and handle data
# ----------------------------------------------------------------------------   

# Assign file location and simulations
dir_home = 'C:\\Users\\drewm\\Documents\\2.0 MSc\\2.0 Simulations\\'
div = '//'
sim_locs = ['PIV_JetSim_M1','PIV_JetSim_M2', 'PIV_JetSim_M3', 'PIV_JetSim_M4']
sim_names = ['M1', 'M2', 'M3', 'M4']

# Load experimental Data
with open('exp_data.pkl', 'rb') as file:
    exp_data = pickle.load(file)
    file.close()
    
# Assign experimental data range and variable
exp_var = 'u_centreline'
exp_xd = exp_data['x'][0] / 6.4
exp_start = np.where(exp_data[f'{exp_var}'] > 0)[0][0]
exp_end = np.where(exp_data[f'{exp_var}'][exp_start:] < 100)[0][0]
    
sim_data = {}
# Load simulation data
for i in range(len(sim_names)):
    sim_data.update({sim_names[i]: pd.read_csv(dir_home + sim_locs[i] + '\\axial_data.csv')})

# Assign plotting range and variable of interest
sim_xd = (sim_data[f'{sim_names[0]}']['Points:0'] - 0.33533) / 0.0064
sim_end = np.where(sim_xd > 40 )[0][0]
sim_var = 'UMean:0'


# Plotting
# ----------------------------------------------------------------------------     
colormap = ["#e60049", "#0bb4ff", "#50e991", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"]
linetype = ["-", "--", ":", "-."]
legend = ['Experimental PIV']
styles = plt.style.available


# Plot data
plt.figure(1)
plt.plot(exp_xd[exp_start:exp_end], exp_data[f'{exp_var}'][exp_start:exp_end], color=colormap[0], linewidth=1)
i = 1
for name in sim_names:
    plt.plot(sim_xd[:sim_end], sim_data[f'{name}'][f'{sim_var}'][:sim_end], color=colormap[i], linewidth=0.8)
    legend += [f'CFD, k-omega SST {name}']
    i += 1        

# Plot settings
plt.grid(True)
plt.grid(True, which='minor')
plt.xlim([0, 40])
plt.ylabel('$U_x$', fontsize=14)
plt.xlabel('Axial Position $x/D_{jet}$', fontsize=14)
plt.title('Axial Velocity Profile', fontsize=14)
plt.legend(legend,fontsize=14)

for style in styles:
    plt.style.use(style)
    plt.show()

