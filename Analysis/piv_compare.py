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
# r_jet = 0.003175;
r_jet = 0.0032
U2 = {}
U_jet = {}
U_excess_jet = {}
B = {}
Rs_max_norm = {}
sim_data = {}

R_stress_legend = ["Rxx", "Ryy", "Rzz", "Rxy", "Ryz", "Rxz"];
UPrime_2_Mean_legend = ["UPrime2Mean_xx", "UPrime2Mean_xy", "UPrime2Mean_xz", "UPrime2Mean_yy", "UPrime2Mean_yz", "UPrime2Mean_zz"];

# # READ SIMULATION DATA AND CLACULATE PARAMETERS
# # ---------------------------------------------------------------------------

# SIMULATION DATA
# Read all radial_data.csv files, making dictionary with key of simulation, containing;
# a dictionary with key of each axial position, value of Dataframe for data

for idx, sim in enumerate(sim_names):
    sim_data.update({ sim : {} })
    U_jet.update({ sim : {} })
    U2.update({ sim : {} })
    U_excess_jet.update({ sim : {} })
    B.update({ sim : {} })
    Rs_max_norm.update({ sim : {} })
    files = os.listdir(fr'{dir_home}\{sim_locs[idx]}')
    
    # Loop through all radial position of each simulation folder to load data and do calculations
    for filename in files:
        if filename[:11] == 'radial_data':
            cur_pos = filename[15:-4]
            all_pos.append(cur_pos)
            # Update DF to include data for current radial postion
            sim_data[sim].update({ cur_pos : pd.read_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}') })
            r_dimless = sim_data[sim][cur_pos]['Points:1'] / np.max(sim_data[sim][cur_pos]['Points:1']) 
            u_jet = np.max(sim_data[sim][cur_pos]['UMean:0'])
            
            # Make key dictionarys for quick debugging if vars already calculated
            if cur_pos != '0.0':
                U_jet[sim].update({ cur_pos : u_jet })
                U2[sim].update({ cur_pos : sim_data[sim][cur_pos]['U_secondary'][0] })
                U_excess_jet[sim].update({ cur_pos : sim_data[sim][cur_pos]['u_excess_centreline'][0] })
                B[sim].update({ cur_pos : sim_data[sim][cur_pos]['b'][0] })
                Rs_max_norm[sim].update({ cur_pos : sim_data[sim][cur_pos]['Rs_max_norm'][0] })
                
            # Add normalized velocity data
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
                sim_data[sim][cur_pos]['U_secondary'] = U_secondary[0]
                # Calculate excess velocities and normalize
                u_excess = sim_data[sim][cur_pos]['UMean:0'] - U_secondary[0]
                u_excess_centreline = np.max(u_excess)
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
                zeta = sim_data[sim][cur_pos]['Points:1'] / b
                sim_data[sim][cur_pos]['b_idx'] = b_idx
                sim_data[sim][cur_pos]['b'] = b
                sim_data[sim][cur_pos]['zeta'] = zeta
                # Save calculated values
                sim_data[sim][cur_pos].to_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}',index=False)
              
            # Normalize Reynolds Stresses with U_excess_jet^2     
            if 'RxxNorm' not in sim_data[sim][cur_pos] and cur_pos != '0.0':
                RxxNorm = sim_data[sim][cur_pos]['turbulenceProperties:R:0'] / U_excess_jet[sim][cur_pos]**2
                RyyNorm = sim_data[sim][cur_pos]['turbulenceProperties:R:1'] / U_excess_jet[sim][cur_pos]**2
                RxyNorm = sim_data[sim][cur_pos]['turbulenceProperties:R:3'] / U_excess_jet[sim][cur_pos]**2
                Rs_max = np.max(sim_data[sim][cur_pos]['turbulenceProperties:R:3'])
                Rs_max_norm= Rs_max / U_excess_jet[sim][cur_pos]**2
                sim_data[sim][cur_pos]['RxxNorm'] = RxxNorm
                sim_data[sim][cur_pos]['RyyNorm'] = RyyNorm
                sim_data[sim][cur_pos]['RxyNorm'] = RxyNorm
                sim_data[sim][cur_pos]['Rs_max_norm'] = Rs_max_norm
                # Save calculated values
                sim_data[sim][cur_pos].to_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}',index=False)
            
            if 'Mach_conv' not in sim_data[sim][cur_pos] and sim != 'kwSST' and cur_pos != '0.0':
                # Speed of sound
                a1_idx = np.where(sim_data[sim][cur_pos]['UMean:0'] == U_jet[sim][cur_pos])[0][0]
                a1 = U_jet[sim][cur_pos] / sim_data[sim][cur_pos]['Ma'][a1_idx]
                a2_idx = np.where(np.isclose(sim_data[sim][cur_pos]['UMean:0'], U2[sim][cur_pos], atol=0.25))[0][0]
                a2 = U2[sim][cur_pos] / sim_data[sim][cur_pos]['Ma'][a2_idx]
                
                # Convective terms
                U_conv = (a1*U_jet[sim][cur_pos] + a2*U2[sim][cur_pos]) / (a1 + a2)
                Mc_jet = (U_jet[sim][cur_pos]-U_conv) / a1
                Mc_secondary = (U_conv - U2[sim][cur_pos]) / a2
                M_conv = np.sqrt(Mc_jet*Mc_secondary)
                
                # Save calculated values
                sim_data[sim][cur_pos]['U_conv'] = U_conv
                sim_data[sim][cur_pos]['Ma_conv'] = M_conv
                sim_data[sim][cur_pos].to_csv(fr'{dir_home}\{sim_locs[idx]}\{filename}',index=False)                

# Calculate experimental fits
Mc_unity = np.linspace(0, 1, num=100)
phi_Mc_unity = 0.8 * np.exp(-3 * Mc_unity**2 ) + 0.2

# Re-order all_pos
all_pos = sorted(set([np.round(float(pos), 4) for pos in all_pos]))

# Select what positions and variable to plot        
xd_pos = np.linspace(10, 25, 7)
r_dimless = sim_data['kwSST']['0.0']['Points:1']/np.max(sim_data['kwSST']['0.0']['Points:1'])
sim_x_var = ''
sim_y_var = 'turbulenceProperties:R:3'


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
exp_x_var = 'r_norm'
exp_y_var = 'Rxy'


              
# # PLOTTING
# # ---------------------------------------------------------------------------

# PLOT SETTINTG
# styles = plt.style.available
# print(styles)
bright = sns.set_palette(sns.color_palette("bright"))
plt.style.use('seaborn-v0_8')
plt.figure(dpi=1000)
mpl.rcParams['font.family'] = 'book antiqua'
# cur_color = next(ax._get_lines.prop_cycler)['color']
# next_color = next(ax._get_lines.prop_cycler)['color']
# next_next_color = next(ax._get_lines.prop_cycler)['color']
# colors = [cur_color, next_color, next_next_color]

# PLOT DATA
# Loop through radial poitions
color = ['blue', 'green', 'magenta', 'red', 'purple', 'orange', 'teal']
mark = ['o', 's', '^', '<', '>', 'D', 'p']

for idx, pos in enumerate(xd_pos):
    plt.clf()
    legend = []    
    fig, ax = plt.subplots()
    # Plot experimental data
    ax.plot(exp_data[exp_x_var][:,int(exp_x_loc[idx])], exp_data[exp_y_var][:,int(exp_x_loc[idx])],
            marker=mark[idx], linestyle='', markeredgewidth=1,
            alpha=0.75, markerfacecolor='none', markeredgecolor=color[idx],  
            markersize=6)
    legend.append(fr'xd = {pos}     PIV Experiment') 
    
    # Loop through simulations alnd plot each
    for sim in sim_names:
        if sim == 'kwSST':
            ax.plot(r_dimless, sim_data[sim][f'{pos}'][sim_y_var],
                    color='black', alpha=0.75, linewidth=1, linestyle='--')
        else:
            ax.plot(r_dimless, sim_data[sim][f'{pos}'][sim_y_var],
                    color=color[idx], alpha=0.85, linewidth=1)
        legend.append(fr'xd = {pos}     {sim}') 

    # Plot rendering
    ax.grid(True, which='minor')
    ax.set_xlim([0, 1])
    ax.set_ylim([-25, 1250])
    ax.set_ylabel('$R_{xy}$     \n[$m^2/s^2$]', fontsize=16, labelpad=15)
    ax.yaxis.label.set(rotation='horizontal', ha='right');
    ax.set_xlabel('$r/R_{tube}$', fontsize=16, labelpad=10)
    ax.set_title('Reynolds Shear Stress Profile', fontsize=14, pad=12)
    ax.legend(legend,fontsize=10)

    plt.show()

    fig.savefig(fr'C:\Users\drewm\Documents\2.0 MSc\2.0 Simulations\Figures\Python\PIV Compare\Rxy_xd{pos:.3}_PIV.png',
                dpi=1000 ,bbox_inches='tight', pad_inches=0.15)
