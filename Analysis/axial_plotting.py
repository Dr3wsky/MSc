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
import seaborn as sns

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
exp_var = 'v_centreline'
exp_xd = exp_data['x'][0] / 6.4
# exp_start = np.where(exp_data[f'{exp_var}'] != 0)[0][0]
# exp_end = np.where(exp_data[f'{exp_var}'][exp_start:] < 100)[0][0]
exp_start = 10
exp_end = len(exp_xd)

    
sim_data = {}
# Load simulation data
for i in range(len(sim_names)):
    sim_data.update({sim_names[i]: pd.read_csv(dir_home + sim_locs[i] + '\\axial_data.csv')})

# Assign plotting range and variable of interest
sim_xd = (sim_data[f'{sim_names[0]}']['Points:0'] - 0.33533) / 0.0064
sim_end = np.where(sim_xd > 40 )[0][0]
sim_var = 'UMean:1'


# Plotting
# ----------------------------------------------------------------------------     
# legend = ['Experimental PIV']
legend = []
plt.clf()
# styles = plt.style.available
# print(styles)
bright = sns.color_palette("bright")
sns.set_palette(bright)
plt.style.use('seaborn-v0_8-bright')
plt.figure(dpi=1000)
fig, ax = plt.subplots()



# Plot data
# ax.plot(exp_xd[exp_start:exp_end-4], exp_data[f'{exp_var}'][0][exp_start:exp_end], linestyle='--', linewidth=0.8 )
i = 1
for name in sim_names:
    ax.plot(sim_xd[:sim_end], sim_data[f'{name}'][f'{sim_var}'][:sim_end], color=bright[i], linewidth=1.25)
    legend += [f'CFD, k-omega SST {name}']
    i += 1        

# Plot settings
ax.grid(True, which='minor')
ax.set_xlim([0, 25])
ax.set_ylabel('$U_y$\t\n [m/s]', fontsize=12, labelpad=10)
ax.yaxis.label.set(rotation='horizontal', ha='right');
ax.set_xlabel('Axial Position $x/D_{jet}$', fontsize=12, labelpad=10)
ax.set_title('Centreline Velocity Profile', fontsize=14, pad=10)
ax.legend(legend,fontsize=10)

plt.show()

fig.savefig('C:\\Users\\drewm\\Documents\\2.0 MSc\\2.0 Simulations\\Figures\\Python\\Grid Independence\\Uy_M1-4.png',
        dpi=1000 ,bbox_inches='tight', pad_inches=0.15)
