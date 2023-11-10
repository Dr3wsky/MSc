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
sim_names = ['M4']
files = os.listdir(dir_home + sim_locs[0])
all_pos = []
sim_data = {}

# Read all radial_data.csv files, making dictionary with key of each axial position, value of Dataframe for data
for num, name in enumerate(files):
    if name[:11] == 'radial_data':
        all_pos.append( f'{name[15:-4]}')
        sim_data.update({f'{name[15:-4]}': pd.read_csv(dir_home + sim_locs[0] + '\\' + name)})
        
test = 'holder'




