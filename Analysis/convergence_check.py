import os
import sys
import pandas as pd
import data_handling

n_interval = 1000
n_runs = 10

axes_font_size = 15

# Data Location
# ----------------------------------------------------------------------------
# dir_home = os.getcwd()
dir_home = 'C:\\Users\\drewm\\Documents\\2.0 MSc\\2.0 Simulations\\'
# Enter simulation folder here
sim = 'PIV_JetSim_M4'
run_dir = dir_home + sim + '\\' + f'postProcessing-{n_runs}'
folds = os.listdir(f'{run_dir}')
all_times = os.listdir(f'{run_dir}/{folds[0]}')

# Check which simulation type
if len(all_times) == 1:
    mode = 'ss'
else:
    mode = 'trans'


datastore = {}


# 1) Domain Mass Flow Convergence    
# ----------------------------------------------------------------------------        
                           
ds_num = [0, 1, 8, 3]
patch = ['farFieldATMFlow', 'inletFlow', 'jetInletFlow', 'outletFlow']
# Extract data and convert dict to DataFrame
for idx in ds_num:
    datastore[folds[idx]] = data_handling.extract_postproc(run_dir + '\\' + folds[idx], folds[idx])
mass_frame = pd.DataFrame(datastore)

# Sum mass flow through domain boundaries 
mf_tot = 0
for column in mass_frame.columns:
    mf_tot += mass_frame.iloc[mass_frame.shape[0]-n_interval:][column].mean()
    if column == 'outletFlow':
        mf_out = mass_frame.iloc[mass_frame.shape[0]-n_interval:][column].mean()
        
mf_convergence = abs(mf_tot / mf_out)

if mf_convergence <= 0.005:
    print(f'Domain mass flow continuity converged: MF_err = {mf_convergence:.6f} < 0.005')
    print('----------------------------------------------------------------------')
else:
    print(f'Mass Flow continuity error = {mf_convergence:.6f}')
    print('Continue running simulation until convergence')
    sys.exit()

# Delete dataframe for rewtriting in next block
for num in ds_num:
    datastore.pop(folds[num])
    
    
# # 2) Residual Convergence    
# ----------------------------------------------------------------------------        

# Convert data to df                           
ds_num = 9
# Assign Column names
residuals = ['p', 'Ux', 'Uy', 'Uz', 'h', 'k', 'omega']
# Extract data and convert dict to DataFrame
for name in residuals:
    datastore[name] = data_handling.extract_postproc(run_dir + '\\' + folds[ds_num], name)
resid_frame = pd.DataFrame(datastore)

# Average last iters for each residual, and check
for column in resid_frame.columns:
    resid_mean = resid_frame.iloc[resid_frame.shape[0]-n_interval:][column].mean()
    if resid_mean <= 1e-5:
        print(f'{column} residuals converged:\t\t {resid_mean:1.3e} < 1e-5')
        continue
    else: 
        print(f'{column} > 1e-5, residuals NOT converged')
        print('Continue running simulation until convergence with lower residuals')
        sys.exit()

print('ALL RESIDUALS CONVERGED < 1e-5')    
print('----------------------------------------------------------------------')

# Delete dataframe for rewtriting in next block
for name in residuals:
    datastore.pop(name)

        
# 3) Solution Monitor Convergence    
# ----------------------------------------------------------------------------

monitors = ['T_max','rho_min', 'U_jet', 'MF_jet', 'entrainment']
conv_lim = 0.0005

# Extract data and convert to dataframe
for name in monitors:
    if name == 'T_max':
        ds_num = 6
        datastore[name] = data_handling.extract_postproc(run_dir + '\\' + folds[ds_num], name[0])
    elif name == 'rho_min':
        ds_num = 7
        datastore[name] = data_handling.extract_postproc(run_dir + '\\' + folds[ds_num], name[0:3])
    elif name =='U_jet':
        ds_num = 5
        datastore[name] = data_handling.extract_postproc(run_dir + '\\' + folds[ds_num], name[0])
    elif name == 'MF_jet':
        ds_num = 4
        datastore[name] = data_handling.extract_postproc(run_dir + '\\' + folds[ds_num], folds[ds_num])
    elif name == 'entrainment':
        ds_num = 10
        tubeInletFlow = data_handling.extract_postproc(run_dir + '\\' + folds[ds_num], folds[ds_num])
        datastore[name] = [a / b for a, b in zip(tubeInletFlow, datastore['MF_jet'])]
      
monit_frame = pd.DataFrame(datastore)

# Test convergence in solution monitors 
for column in monit_frame.columns:
    monit_mean = monit_frame.iloc[monit_frame.shape[0]-n_interval:][column].mean()
    monit_prev_mean =  monit_frame.iloc[monit_frame.shape[0]-2*n_interval:monit_frame.shape[0]-n_interval][column].mean()
    conv = abs(monit_mean - monit_prev_mean) / abs(monit_mean)
    if conv <= conv_lim:
        print(f'{column} converged:\t\t {conv:1.3e} < {conv_lim:.1e}')
        continue
    else: 
        print(f'{conv:1.3e} > {conv_lim:.1e}, {column} NOT converged')
        print('Continue running simulation . . . ')
        # sys.exit()

