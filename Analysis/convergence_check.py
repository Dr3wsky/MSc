import os
import numpy as np
import pandas as pd

n_interval = 1000
n_runs = 10

axes_font_size = 15

# Data Read
all_folder = os.listdir(f'./postProcessing-{n_runs}/')
all_times = os.listdir(f'./postProcessing-{n_runs}/' + all_folder[2])

vars = []
times = []

datastores = []

for fold in range(0, len(all_folder)):
    if fold == 2:
        vars.append(str(all_folder[fold]))
    else:
        vars.append(vars[-1] + str(all_folder[fold]))

    for time in range(2, len(all_times)):
        if fold == 2 and time == 2:
            times.append(str(all_times[time]))
        else:
            times.append(times[-1] + str(all_times[time]))

        file = f'./postProcessing-{nRuns}/{all_folder[fold]}/{all_times[time]}'
        all_files = os.listdir(file)

        for n_files in range(2, len(all_files)):
            fileread = file + '/' + all_files[n_files]

            # Find number of header lines in file
            with open(fileread, 'r') as fileID:
                header = "true"
                numHead = 0
                while header == "true":
                    cur_line = fileID.readline()
                    if "#" in cur_line:
                        numHead += 1
                    else:
                        header = "false"

            ds = pd.read_csv(fileread, skiprows=numHead, delim_whitespace=True)
            datastores.append(ds)

# 1) Mass Flow Convergence
ds_num = [2, 3, 5, 7]
mf_tot = 0
mf_out = 0

for i in ds_num:
    data = datastores[i]
    mf_tot += data.iloc[-1, 1]
    exec(f'mf_{vars[i]} = data.iloc[-1, 1]')
    if vars[i] == 'outletFlow':
        mf_out += data.iloc[-1, 1]

mf_conv = abs(mf_tot / mf_out)

if mf_conv <= 0.005:
    print(f"Mass Flow error converged. {mf_conv:.6f} < 0.005")
else:
    print(f"Mass Flow error = {mf_conv:.6f}")
    print("Continue running simulation until convergence")

# 2) Residuals Limit Check
ds_num = [8]

for i in ds_num:
    data = datastores[i]
    p = data.iloc[-1, 28]
    Ux = data.iloc[-1, 12]
    Uy = data.iloc[-1, 15]
    Uz = data.iloc[-1, 18]
    h = data.iloc[-1, 2]
    k = data.iloc[-1, 23]
    omega = data.iloc[-1, 33]

residuals = pd.DataFrame({'p': p, 'Ux': Ux, 'Uy': Uy, 'Uz': Uz, 'h': h, 'k': k, 'omega': omega})

for column in residuals.columns:
    if residuals.iloc[-1][column] <= 1e-5:
        continue
    else:
        print(f"{column} > 1e-5, residuals not converged")
        print("Continue running simulation until convergence with lower residuals")

print("All residuals converged < 1e-5")

# 3) Solution Monitor Convergence
monit_count = 0
monitors = ["T_max", "Rho_min", "Ujet", "MF_jet", 'Entrainment']
ds_num = [0, 1, 6, 6, 9]
conv_lim = 0.0005

for i in ds_num:
    monit_count += 1
    data = datastores[i]

    if vars[i] == 'Max':
        lastIters = data.iloc[-n_interval:, 1]
        prevIters = data.iloc[-(2 * n_interval):-n_interval, 1]
        mean_maxT = np.mean(lastIters)

    elif vars[i] == 'Min':
        lastIters = data.iloc[-n_interval:, 2]
        prevIters = data.iloc[-(2 * n_interval):-n_interval, 2]
        mean_minRho = np.mean(lastIters)

    elif vars[i] == 'jetOutlet_intB' and monitors[monit_count] == 'Ujet':
        lastIters = (data.iloc[-n_interval:, 4].str.replace("(", "").str.replace(")", "", regex=False).astype(float)) ** 2
        prevIters = (
            data.iloc[-(2 * n_interval):-n_interval, 4]
            .str.replace("(", "")
            .str.replace(")", "", regex=False)
            .astype(float)
        ) ** 2
        mean_Ujet = np.mean(lastIters)

    elif vars[i] == 'jetOutlet_intB' and monitors[monit_count] == 'MF_jet':
        lastIters = (
            data.iloc[-n_interval:, 3]
            * np.pi
            * (0.0032 ** 2)
            * np.sqrt(
                (data.iloc[-n_interval:, 4].str.replace("(", "").str.replace(")", "", regex=False).astype(float)) ** 2
                + (data.iloc[-n_interval:, 6]) ** 2
                + (data.iloc[-n_interval:, 7].str.replace("(", "").str.replace(")", "", regex=False).astype(float)) ** 2
            )
        )
        prevIters = (
            data.iloc[-(2 * n_interval):-n_interval, 3]
            * np.pi
            * (0.0032 ** 2)
            * np.sqrt(
                (data.iloc[-(2 * n_interval):-n_interval, 4].str.replace("(", "").str.replace(")", "", regex=False).astype(float))
                ** 2
                + (data.iloc[-(2 * n_interval):-n_interval, 6]) ** 2
                + (data.iloc[-(2 * n_interval):-n_interval, 7].str.replace("(", "").str.replace(")", "", regex=False).astype(float))
                ** 2
            )
        )
        mean_MFjet = np.mean(lastIters)
        prev_mean_MFjet = np.mean(prevIters)

    elif vars[i] == 'tubeInlet_intB' and monitors[monit_count] == 'Entrainment':
        lastIters = (
            data.iloc[-n_interval:, 3]
            * np.pi
            * ((0.025396 ** 2) - (0.0032 ** 2))
            * np.sqrt(
                (data.iloc[-n_interval:, 4].str.replace("(", "").str.replace(")", "", regex=False).astype(float)) ** 2
                + (data.iloc[-n_interval:, 6]) ** 2
                + (data.iloc[-n_interval:, 7].str.replace("(", "").str.replace(")", "", regex=False).astype(float)) ** 2
            )
        )
        prevIters = (
            data.iloc[-(2 * n_interval):-n_interval, 3]
            * np.pi
            * ((0.025396 ** 2) - (0.0032 ** 2))
            * np.sqrt(
                (data.iloc[-(2 * n_interval):-n_interval, 4].str.replace("(", "").str.replace(")", "", regex=False).astype(float))
                ** 2
                + (data.iloc[-(2 * n_interval):-n_interval, 6]) ** 2
                + (data.iloc[-(2 * n_interval):-n_interval, 7].str.replace("(", "").str.replace(")", "", regex=False).astype(float))
                ** 2
            )
        )
        last_entrain = np.mean(lastIters) / mean_MFjet
        prev_entrain = np.mean(prevIters) / prev_mean_MFjet

    if monitors[monit_count] == 'Entrainment':
        convergence = abs(last_entrain - prev_entrain) / last_entrain
    else:
        meanVal = np.mean(lastIters)
        meanPrev = np.mean(prevIters)
        convergence = abs(meanVal - meanPrev) / meanVal

    if convergence < conv_lim:
        print(f"{monitors[monit_count]} Converged. Monitor variance = {convergence:.6f} over last 1k iterations")
    else:
        print(f"{monitors[monit_count]} NOT converged. Continue simulation until monitor < 0.0005 variance")
