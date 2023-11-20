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

dir_home = r'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations'
sim_locs = ['PIV_JetSim_M4', 'kOmegaSSTdd_PIV_JetSim_M4']
sim_names = ['kw', 'kwSSTdd']
all_pos = []

sim_data = {}

# Read all radial_data.csv files, making dictionary with key of simulation, containing;
# a dictionary with key of each axial position, value of Dataframe for data

for idx, sim in enumerate(sim_names):
    sim_data.update({ sim : {} })
    files = os.listdir(fr'{dir_home}\{sim_locs[idx]}')
    # Loop through all radial position of each simulation folder    
    for filename in files:
        if filename[:11] == 'radial_data':
            cur_pos = filename[15:-4]
            if cur_pos not in all_pos:
                all_pos.append(cur_pos)
            sim_data[sim].update({ cur_pos : pd.read_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}') })
            # # Add normalized data for simulation outputs. Will calculate excess,half-width, etc in other scripts. 
            if 'UNorm:0' not in sim_data[sim][cur_pos]:
                u_jet = np.max(sim_data[sim][cur_pos]['UMean:0'])
                Ux_norm = sim_data[sim][cur_pos]['UMean:0'] / u_jet 
                Uy_norm = sim_data[sim][cur_pos]['UMean:1'] / u_jet
                tke_norm = sim_data[sim][cur_pos]['kMean']/ u_jet**2
                sim_data[sim][cur_pos]['UNorm:0'] = Ux_norm
                sim_data[sim][cur_pos]['UNorm:1'] = Uy_norm
                sim_data[sim][cur_pos]['kNorm'] = tke_norm
                sim_data[sim][cur_pos].to_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}',index=False)
                ## Add r and theta calcs for points and velocities, from axial calcs
                
# Select what positions and variable to plot        
xd_pos = np.linspace(12.5, 17.5, 3)
r_dimless = sim_data['kw']['0.0']['Points:1'] / np.max(sim_data['kw']['0.0']['Points:1'])
sim_var = 'UNorm:0'

# # PLOTTING
# # ---------------------------------------------------------------------------

# PLOT SETTINTGS
legend = []
plt.clf()
# styles = plt.style.available
# print(styles)
bright = sns.set_palette(sns.color_palette("bright"))
plt.style.use('seaborn-v0_8')
plt.figure(dpi=1000)
mpl.rcParams['font.family'] = 'book antiqua'

# PLOT DATA
fig, ax = plt.subplots()
# Loop through each axial position
for pos in xd_pos:
    hold_color=next(ax._get_lines.prop_cycler)['color']
    # Loop through simulations for that position
    for sim in sim_names:
        if sim == 'kw':
            ax.plot(r_dimless, sim_data[sim][f'{pos}'][sim_var], color=hold_color, linewidth=0.75)
        else:
            ax.plot(r_dimless, sim_data[sim][f'{pos}'][sim_var], color=hold_color, linestyle='--', linewidth=0.9)
        legend.append(fr'xd = {pos}     {sim}')  

ax.grid(True, which='minor')
# ax.set_xlim([0, 25])
# ax.set_ylim([.25, 1])
ax.set_ylabel('Axial Velocity: $\\frac{u_i}{u_jet}$', fontsize=14, labelpad=15)
# ax.yaxis.label.set(rotation='horizontal', ha='right');
ax.set_xlabel('Radial Position: $r/R_{tube}$', fontsize=12, labelpad=10)
ax.set_title('Radial Velocity Profiles', fontsize=14, pad=12)
ax.legend(legend,fontsize=10)

fig.savefig(r'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\Figures\Python\Radial Iso\test.png',
        dpi=1000 ,bbox_inches='tight', pad_inches=0.15)

