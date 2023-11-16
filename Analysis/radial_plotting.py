"""
This script interprets .csv data from a specified simulation, and plots the desired variable for a given domain.
This script will be expanded upon later to include importing and comparing experimental data, as well as calulations of further parameters. 

The OpenFOAM data is extracted from Paraview using an iterative plotOverLine method along specified axial locations of a confined jet.
The simulation data is loaded as a dictionary where the value is a Dataframe for each axial position corresponding key.
Plots for specified variables are generated with matplotlib.

Note: Ensure the necessary data files are available in the specified file paths before running this script.

Author: Andrew McFadden
Purpose: MSc. Data analysis
Date: 10/11/2023
Version: V1
"""

# Load external libraries and modules
import os
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

# # DEFINE WORKSPACE AND VARS
# # ---------------------------------------------------------------------------

dir_home = 'C:\\Users\\drewm\\Documents\\2.0 MSc\\2.0 Simulations\\'
sim_locs = ['PIV_JetSim_M4']
sim_names = ['kw']
files = os.listdir(dir_home + sim_locs[0])
all_pos = []

sim_data = {}
    
for num, name in enumerate(files):
    if name[:11] == 'radial_data':
        cur_pos = name[15:-4]
        all_pos.append(cur_pos)
        sim_data.update({ cur_pos : pd.read_csv(dir_home + sim_locs[0] + '\\' + name) })
        # # Add normalized data, currently on ly for velocity. Will calc excess velocity , RNorm and such in another script
        if 'UNorm:0' not in sim_data[cur_pos]:
            Ux_norm = sim_data[cur_pos]['UMean:0']/np.max(sim_data[cur_pos]['UMean:0'])
            Uy_norm = sim_data[cur_pos]['UMean:1']/np.max(sim_data[cur_pos]['UMean:0'])
            tke_norm = sim_data[cur_pos]['kMean']/(np.max(sim_data[cur_pos]['UMean:0'])**2)
            sim_data[cur_pos]['UNorm:0'] = Ux_norm
            sim_data[cur_pos]['UNorm:1'] = Uy_norm
            sim_data[cur_pos]['kNorm'] = tke_norm
            sim_data[cur_pos].to_csv(f'{dir_home}{sim_locs[0]}\\{name}',index=False)
## Add r and theta calcs for points and velocities, from axial calcs

# Select what positions and variable to plot        
xd_pos = np.linspace(0, 30, 13)
r_dimless = sim_data['0.0']['Points:1']/np.max(sim_data['0.0']['Points:1'])
sim_var = 'UMean:0'

# # PLOTTING
# # ---------------------------------------------------------------------------

# PLOT SETTINTGS
legend = []
plt.clf()
# styles = plt.style.available
# print(styles)
bright = sns.set_palette(sns.color_palette("bright"))
plt.style.use('seaborn-v0_8-bright')
plt.figure(dpi=1000)
mpl.rcParams['font.family'] = 'book antiqua'

# PLOT DATA
fig, ax = plt.subplots()
for pos in xd_pos:
    ax.plot(r_dimless, sim_data[pos][sim_var])
    legend.append(f'xd = {pos}')  

ax.grid(True, which='minor')
# ax.set_xlim([0, 25])
# ax.set_ylim([.25, 1])
ax.set_ylabel('Velocity: $u_i$', fontsize=14, labelpad=15)
# ax.yaxis.label.set(rotation='horizontal', ha='right');
ax.set_xlabel('Radial Position: $r/R_{tube}$', fontsize=12, labelpad=10)
ax.set_title('Radial Velocity Profile', fontsize=14, pad=12)
ax.legend(legend,fontsize=10)

fig.savefig('C:\\Users\\drewm\\Documents\\2.0 MSc\\2.0 Simulations\\Figures\\Python\\Radial Iso\\Ux_M4.png',
        dpi=1000 ,bbox_inches='tight', pad_inches=0.15)

