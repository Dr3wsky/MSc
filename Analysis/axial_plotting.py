"""
This script performs data analysis on experimental and simulated data from OpenFOAM.
The OpenFOAM data is extracted from Paraview using plotOverLine method along the centreline axis of a confined jet.
The simulation data is loaded as a dictionary where the value is a Dataframe for each Mesh as the corresponding key.
Due to the nature of experimental data, it is unpacked frfom a pickeld file and loaded as a dictionary where values are multi-dimentional arrays.
Plots for specified variables are generated with matplotlib.

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
import matplotlib as mpl

# # LOAD & HANDLE DATA
# # ---------------------------------------------------------------------------  

# Assign file location and simulations
dir_home = r'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations'
sim_locs = [ 'PIV_JetSim_M4', '15deg_PIVJetSim_M4']
sim_names = ['5deg wedge', '15deg wedge']

# LOAD EXP DATA
with open('exp_data.pkl', 'rb') as file:
    exp_data = pickle.load(file)
    file.close()
    
# Assign experimental data range and variable
exp_var = 'u_centreline_norm'
exp_xd = exp_data['x'][0] / 6.4
# # Indicies for x plotting
# exp_start = np.where(exp_data[f'{exp_var}'] > 0)[0][0]
# exp_end = len(exp_xd)
# # Indicies for y-plotting
exp_start = 10
exp_end = len(exp_xd)

    
sim_data = {}
# LOAD SIMULATION DATA
for i in range(len(sim_names)):
    sim_data.update({sim_names[i]: pd.read_csv(fr'{dir_home}\{sim_locs[i]}\axial_data.csv')})
    # Add normalized data
    if 'UNorm:0' not in sim_data[sim_names[i]]:
        Ux_norm = sim_data[sim_names[i]]['UMean:0'] / np.max(sim_data[f'{sim_names[i]}']['UMean:0'])
        Uy_norm = sim_data[sim_names[i]]['UMean:1'] / np.max(sim_data[f'{sim_names[i]}']['UMean:0'])
        Uz_norm = sim_data[sim_names[i]]['UMean:2'] / np.max(sim_data[f'{sim_names[i]}']['UMean:0'])
        tke_norm = sim_data[sim_names[i]]['kMean'] / (np.max(sim_data[f'{sim_names[i]}']['UMean:0'])**2)
        sim_data[sim_names[i]]['UNorm:0'] = Ux_norm
        sim_data[sim_names[i]]['UNorm:1'] = Uy_norm
        sim_data[sim_names[i]]['UNorm:2'] = Uz_norm
        sim_data[sim_names[i]]['kNorm'] = tke_norm
        sim_data[sim_names[i]].to_csv(fr'{dir_home}\{sim_locs[i]}\axial_data.csv',index=False)
        
    # Convert to polar coords    
    if 'Points:r' not in sim_data[sim_names[i]]:
        r = np.sqrt(sim_data[sim_names[i]]['Points:1']**2 + sim_data[sim_names[i]]['Points:2']**2)
        theta_rads = np.arctan(sim_data[sim_names[i]]['Points:2'] / sim_data[sim_names[i]]['Points:1'])
        sim_data[sim_names[i]]['Points:r'] = r
        sim_data[sim_names[i]]['Points:theta'] = theta_rads
        sim_data[sim_names[i]].to_csv(dir_home + sim_locs[i] + '\\axial_data.csv',index=False)
    if 'UMean:r' not in sim_data[sim_names[i]]:
        UMean_r = sim_data[sim_names[i]]['UMean:1'] * np.cos(sim_data[sim_names[i]]['Points:theta']) 
        UMean_theta = -sim_data[sim_names[i]]['UMean:1'] * np.sin(sim_data[sim_names[i]]['Points:theta']) 
        sim_data[sim_names[i]]['UMean:r'] = UMean_r
        sim_data[sim_names[i]]['UMean:theta'] = UMean_theta
        sim_data[sim_names[i]].to_csv(fr'{dir_home}\{sim_locs[i]}\axial_data.csv',index=False)

# # Re-save csv with normalized data and then remove the norm calcs. . . ? 
# df.to_csv('data.csv', index=False)
        

# Assign plotting range and variable of interest
sim_xd = (sim_data[f'{sim_names[0]}']['Points:0'] - 0.33533) / 0.0064
sim_end = np.where(sim_xd > 80 )[0][0]
sim_var = 'UMean:theta'


# # PLOTTING
# # ---------------------------------------------------------------------------   

# PLOT SETTINGS
# legend = ['Experimental PIV']
legend = []
plt.clf()
# styles = plt.style.available
# print(styles)
bright = sns.set_palette(sns.color_palette("bright"))
plt.style.use('seaborn-v0_8-bright')
plt.figure(dpi=1000)
mpl.rcParams['font.family'] = 'book antiqua'
fig, ax = plt.subplots()

# # PLOT DATA
# x data range experimental
# ax.plot(exp_xd[exp_start:exp_end], exp_data[f'{exp_var}'][exp_start:exp_end], linestyle='--', linewidth=1 )
# y data range experimental
# ax.plot(exp_xd[exp_start:exp_end-4], exp_data[f'{exp_var}'][0][exp_start:exp_end], linestyle='--', linewidth=1 )

# Simulation data
for name in sim_names:
    if name == '15deg wedge':
        ax.plot(sim_xd[:sim_end], sim_data[f'{name}'][f'{sim_var}'][:sim_end], linestyle='--', linewidth=0.9)
    else:
        ax.plot(sim_xd[:sim_end], sim_data[f'{name}'][f'{sim_var}'][:sim_end], linewidth=1.15)
    legend.append(f'{name}, CFD k-omega SST')   

# # PLOT SETTINGS
ax.grid(True, which='minor')
ax.set_xlim([0, 60])
# ax.set_ylim([.25, 1])
# ax.set_ylabel('$\\frac{u_{\\theta}}{u_{jet}}$', fontsize=20, labelpad=15)
ax.set_ylabel('$u_{\\theta}$\t \n[m/s]', fontsize=16, labelpad=15)
ax.yaxis.label.set(rotation='horizontal', ha='right');
ax.set_xlabel('Axial Position\t $x/D_{jet}$', fontsize=12, labelpad=10)
ax.set_title('Centreline Velocity Profile', fontsize=14, pad=12)
ax.legend(legend,fontsize=10)


plt.show()

fig.savefig(r'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\Figures\Python\Wedge Compare\UthetaNorm_5-15_deg_raw.png',
        dpi=1000 ,bbox_inches='tight', pad_inches=0.15)
