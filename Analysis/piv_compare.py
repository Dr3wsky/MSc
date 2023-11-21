"""
This script interprets .csv data from a specified simulation, and plots variables for each axial position specified.
It also import pickled experimental data and unpacks it as a  where the values are multi-dimaensional arrays. 
The arrays are parsed to make a dictionary object similar to the siulation data with a dataframe for each value at each axial position

The radial profile for various parameters of interest are then plotted as axial progressions and compared.

Author: Andrew McFadden
Purpose: MSc. Data analysis
Date: 10/16/2023
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

# # DEFINE WORKSPACE AND INITIALIZE VARS
# # ---------------------------------------------------------------------------

dir_home = r'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations'
sim_locs = ['PIV_JetSim_M4', 'kOmegaSSTdd_PIV_JetSim_M4']
sim_names = ['kwSST', 'kwSSTdd']
all_pos = []
r_tube = 0.025396
r_jet = 0.003175;
U_2 = {}
U_jet = {}
U_excess_centreline = {}
B = {}

sim_data = {}

# # READ SIMULATION DATA AND CLACULATE PARAMETERS
# # ---------------------------------------------------------------------------

# SIMULATION DATA
# Read all radial_data.csv files, making dictionary with key of simulation, containing;
# a dictionary with key of each axial position, value of Dataframe for data

for idx, sim in enumerate(sim_names):
    sim_data.update({ sim : {} })
    files = os.listdir(fr'{dir_home}\{sim_locs[idx]}')
    # Loop through all radial position of each simulation folder
    for filename in files:
        if filename[:11] == 'radial_data':
            cur_pos = filename[15:-4]
            all_pos.append(cur_pos)
            # Update DF to include data for current radial postion
            sim_data[sim].update({ cur_pos : pd.read_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}') })
            
            # Add normalized velocity data
            u_jet = np.max(sim_data[sim][cur_pos]['UMean:0'])
            U_jet.update({ cur_pos : u_jet })
            r_dimless = sim_data[sim][cur_pos]['Points:1'] / r_tube
            
            if 'UNorm:0' not in sim_data[sim][cur_pos]:
                Ux_norm = sim_data[sim][cur_pos]['UMean:0'] / u_jet 
                Uy_norm = sim_data[sim][cur_pos]['UMean:1']/ u_jet
                tke_norm = sim_data[sim][cur_pos]['kMean']/ u_jet**2
                # Save new vals to dataset
                sim_data[sim][cur_pos]['UNorm:0'] = Ux_norm
                sim_data[sim][cur_pos]['UNorm:1'] = Uy_norm
                sim_data[sim][cur_pos]['kNorm'] = tke_norm
                sim_data[sim][cur_pos].to_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}',index=False)
                
            # Calculate jet flow parameters
            if 'U_secondary' not in sim_data[sim][cur_pos] and cur_pos != '0.0':
                # Trim the wall effects and jet flow
                y_trim_start = np.where(r_dimless > 0.15)[0][0]
                y_trim_end = np.where(r_dimless < 0.9)[0][-1]
                
                # Use histogram to find common velocity for secondary stream
                hist, bins = np.histogram(sim_data[sim][cur_pos]['UMean:0'][y_trim_start:y_trim_end], bins=250) 
                bin_idx = np.where(hist == np.max(hist[:90]))
                U_secondary = bins[bin_idx[0]]
                # plt.hist(sim_data[sim][cur_pos]['UMean:0'][:y_trim_end], bins=250, edgecolor='black')
                # plt.show()
                U_2.update({ cur_pos : U_secondary[0] })
                sim_data[sim][cur_pos]['U_secondary'] = U_secondary[0]
                
                # Calculate excess velocities and normalize
                u_excess = sim_data[sim][cur_pos]['UMean:0'] - U_secondary[0]
                u_excess_centreline = np.max(u_excess)
                U_excess_centreline.update({ cur_pos : u_excess_centreline })
                f_zeta = u_excess / u_excess_centreline
                v_norm = sim_data[sim][cur_pos]['UMean:1'] / u_excess_centreline
                sim_data[sim][cur_pos]['u_excess'] = u_excess
                sim_data[sim][cur_pos]['u_excess_centreline'] = u_excess_centreline
                sim_data[sim][cur_pos]['f_zeta'] = f_zeta
                sim_data[sim][cur_pos]['v_norm'] = v_norm
                
                # Find half-width
                b_velocity = u_excess_centreline / 2
                b_idx = np.where(u_excess < b_velocity)[0][0]
                b = sim_data[sim][cur_pos]['Points:1'][b_idx]
                B.update({ cur_pos : b })
                zeta = sim_data[sim][cur_pos]['Points:1'] / b
                sim_data[sim][cur_pos]['b_idx'] = b_idx
                sim_data[sim][cur_pos]['b'] = b
                sim_data[sim][cur_pos]['zeta'] = zeta
                
                # Save calculated values
                sim_data[sim][cur_pos].to_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}',index=False)
                
# Select what positions and variable to plot        
xd_pos = np.linspace(15, 20, 3)
r_dimless = sim_data['kwSST']['0.0']['Points:1']/np.max(sim_data['kwSST']['0.0']['Points:1'])
sim_x_var = 'zeta'
sim_y_var = 'f_zeta'

# # LOAD EXPERIMENTAL DATA
# # ---------------------------------------------------------------------------
# Read data from pickled file
with open('exp_data.pkl', 'rb') as file:
    exp_data = pickle.load(file)
    file.close()  
    
# Use scale from raw exp to set and find xd locations
scale = exp_data['x'][0][1] - exp_data['x'][0][0]
exp_xd_ass = xd_pos
exp_x_loc= np.round( exp_xd_ass * (r_jet*2*1000) / scale )
exp_x_var = 'zeta'
exp_y_var = 'f_zeta'


              
# # PLOTTING
# # ---------------------------------------------------------------------------

# PLOT SETTINTGS
plt_xd =[17.5]
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
# Loop through radial poitions
for idx, pos in enumerate(plt_xd):
    # hold_color=next(ax._get_lines.prop_cycler)['color']
    cur_color = next(ax._get_lines.prop_cycler)['color']
    
    # plot experimental data
    ax.plot(exp_data[exp_x_var][:,int(exp_x_loc[idx])], exp_data[exp_y_var][:,int(exp_x_loc[idx])],
            marker='o', linestyle='', markeredgewidth=1,
            alpha=0.75, markerfacecolor='none', markeredgecolor=cur_color,  
            markersize=6)
    legend.append(fr'xd = {pos}     PIV Experiment') 
    # Loop through simulations alnd plot each
    for sim in sim_names:
        if sim == 'kwSST':
            ax.plot(sim_data[sim][f'{pos}'][sim_x_var], sim_data[sim][f'{pos}'][sim_y_var],
                    color='black', alpha=0.65, linewidth=1)
        else:
            ax.plot(sim_data[sim][f'{pos}'][sim_x_var], sim_data[sim][f'{pos}'][sim_y_var],
                    color=cur_color, linestyle='--', linewidth=1.25)
        legend.append(fr'xd = {pos}     {sim}') 

ax.grid(True, which='minor')
ax.set_xlim([0, 2])
# ax.set_ylim([.25, 1])
ax.set_ylabel('$\\frac{\\overline{U}}{\\overline{U}_m}$', fontsize=16, labelpad=15)
ax.yaxis.label.set(rotation='horizontal', ha='right');
ax.set_xlabel('$\\eta = \\frac{r}{b}$', fontsize=16, labelpad=10)
ax.set_title('Normalized Excess Velocity Profile', fontsize=14, pad=12)
ax.legend(legend,fontsize=10)

fig.savefig(r'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\Figures\Python\PIV Compare\test.png',
        dpi=1000 ,bbox_inches='tight', pad_inches=0.15)
