import os
import sys
import numpy as np
import pandas as pd
import data_handling

n_interval = 1000
n_runs = 10

axes_font_size = 15

# Data Location
# ----------------------------------------------------------------------------
loc = os.getcwd()
folds = os.listdir(f'{loc}/postProcessing-{n_runs}/')
all_times = os.listdir(f'{loc}/postProcessing-{n_runs}/{folds[0]}')

# Check which simulation type
if len(all_times) == 1:
    mode = 'ss'
else:
    mode = 'trans'

dir_home = os.getcwd()
run_dir = dir_home + '\\' + f'postProcessing-{n_runs}'


# Domain Mass Flow Convergence    
# ----------------------------------------------------------------------------        

# Convert data to df                           
ds_num = [0, 1, 8, 3]
datastore = {}
for idx in ds_num:
    datastore[folds[idx]] = data_handling.extract(run_dir + '\\' + folds[idx], folds[idx])
mf = pd.DataFrame(datastore)

# Sum mass flow through domain boundaries      
mf_tot = mf.farfieldATMFlow.iloc[mf.shape[0]-n_interval:].mean() \
        + mf.inletFlow.iloc[mf.shape[0]-n_interval:].mean() \
        + mf.jetInletFlow.iloc[mf.shape[0]-n_interval:].mean() \
        + mf.outletFlow.iloc[mf.shape[0]-n_interval:].mean()
        
mf_out = mf.outletFlow.iloc[mf.shape[0]-n_interval:].mean()
mf_convergence = abs(mf_tot / mf_out)

if mf_convergence <= 0.005:
    print(f"Domain mass flow continuity converged. MF_err = {mf_convergence:.6f} < 0.005")
    print('----------------------------------------------------------------------')
else:
    print(f"Mass Flow continuity error = {mf_convergence:.6f}")
    print("Continue running simulation until convergence")
    sys.exit()

# Delete dataframe for rewtriting in next block
for num in ds_num:
    datastore.pop(folds[num])
    
    
# # Residual COnvergence    
# ----------------------------------------------------------------------------        

# Convert data to df                           
ds_num = 9
# Assign Column names
residuals = ['p', 'Ux', 'Uy', 'Uz', 'h', 'k', 'omega']
datastore = {}
for name in residuals:
    datastore[name] = data_handling.extract(run_dir + '\\' + folds[ds_num], name)

rf = pd.DataFrame(datastore)

# Average last iters for each residual, and check
for column in rf.columns:
    rf_mean = rf.iloc[rf.shape[0]-n_interval:][column].mean()
    if rf_mean <= 1e-5:
        print(f'{column} residuals converged < 1e-5')
        continue
    else: 
        print(f"{column} > 1e-5, residuals not converged")
        print("Continue running simulation until convergence with lower residuals")
        sys.exit()

print("ALL RESIDUALS CONVERGED < 1e-5")    
print('----------------------------------------------------------------------')
        
# Solution Monitor Convergence    
# ----------------------------------------------------------------------------

# # 3) Solution Monitor Convergence
# monit_count = 0
# monitors = ["T_max", "Rho_min", "Ujet", "MF_jet", 'Entrainment']
# ds_num = [0, 1, 6, 6, 9]
# conv_lim = 0.0005

# for i in ds_num:
#     monit_count += 1
#     data = datastores[i]

#     if folds[i] == 'Max':
#         lastIters = data.iloc[-n_interval:, 1]
#         prevIters = data.iloc[-(2 * n_interval):-n_interval, 1]
#         mean_maxT = np.mean(lastIters)

#     elif folds[i] == 'Min':
#         lastIters = data.iloc[-n_interval:, 2]
#         prevIters = data.iloc[-(2 * n_interval):-n_interval, 2]
#         mean_minRho = np.mean(lastIters)

#     elif folds[i] == 'jetOutlet_intB' and monitors[monit_count] == 'Ujet':
#         lastIters = (data.iloc[-n_interval:, 4].str.replace("(",
#                      "").str.replace(")", "", regex=False).astype(float)) ** 2
#         prevIters = (
#             data.iloc[-(2 * n_interval):-n_interval, 4]
#             .str.replace("(", "")
#             .str.replace(")", "", regex=False)
#             .astype(float)
#         ) ** 2
#         mean_Ujet = np.mean(lastIters)

#     elif folds[i] == 'jetOutlet_intB' and monitors[monit_count] == 'MF_jet':
#         lastIters = (
#             data.iloc[-n_interval:, 3]
#             * np.pi
#             * (0.0032 ** 2)
#             * np.sqrt(
#                 (data.iloc[-n_interval:, 4].str.replace("(",
#                  "").str.replace(")", "", regex=False).astype(float)) ** 2
#                 + (data.iloc[-n_interval:, 6]) ** 2
#                 + (data.iloc[-n_interval:, 7].str.replace("(",
#                    "").str.replace(")", "", regex=False).astype(float)) ** 2
#             )
#         )
#         prevIters = (
#             data.iloc[-(2 * n_interval):-n_interval, 3]
#             * np.pi
#             * (0.0032 ** 2)
#             * np.sqrt(
#                 (data.iloc[-(2 * n_interval):-n_interval, 4].str.replace("(",
#                  "").str.replace(")", "", regex=False).astype(float))
#                 ** 2
#                 + (data.iloc[-(2 * n_interval):-n_interval, 6]) ** 2
#                 + (data.iloc[-(2 * n_interval):-n_interval, 7].str.replace("(",
#                    "").str.replace(")", "", regex=False).astype(float))
#                 ** 2
#             )
#         )
#         mean_MFjet = np.mean(lastIters)
#         prev_mean_MFjet = np.mean(prevIters)

#     elif folds[i] == 'tubeInlet_intB' and monitors[monit_count] == 'Entrainment':
#         lastIters = (
#             data.iloc[-n_interval:, 3]
#             * np.pi
#             * ((0.025396 ** 2) - (0.0032 ** 2))
#             * np.sqrt(
#                 (data.iloc[-n_interval:, 4].str.replace("(",
#                  "").str.replace(")", "", regex=False).astype(float)) ** 2
#                 + (data.iloc[-n_interval:, 6]) ** 2
#                 + (data.iloc[-n_interval:, 7].str.replace("(",
#                    "").str.replace(")", "", regex=False).astype(float)) ** 2
#             )
#         )
#         prevIters = (
#             data.iloc[-(2 * n_interval):-n_interval, 3]
#             * np.pi
#             * ((0.025396 ** 2) - (0.0032 ** 2))
#             * np.sqrt(
#                 (data.iloc[-(2 * n_interval):-n_interval, 4].str.replace("(",
#                  "").str.replace(")", "", regex=False).astype(float))
#                 ** 2
#                 + (data.iloc[-(2 * n_interval):-n_interval, 6]) ** 2
#                 + (data.iloc[-(2 * n_interval):-n_interval, 7].str.replace("(",
#                    "").str.replace(")", "", regex=False).astype(float))
#                 ** 2
#             )
#         )
#         last_entrain = np.mean(lastIters) / mean_MFjet
#         prev_entrain = np.mean(prevIters) / prev_mean_MFjet

#     if monitors[monit_count] == 'Entrainment':
#         convergence = abs(last_entrain - prev_entrain) / last_entrain
#     else:
#         meanVal = np.mean(lastIters)
#         meanPrev = np.mean(prevIters)
#         convergence = abs(meanVal - meanPrev) / meanVal

#     if convergence < conv_lim:
#         print(
#             f"{monitors[monit_count]} Converged. Monitor variance = {convergence:.6f} over last 1k iterations")
#     else:
#         print(
#             f"{monitors[monit_count]} NOT converged. Continue simulation until monitor < 0.0005 variance")
