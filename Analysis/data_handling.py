# Module extract_data.py
# Used to walk through file folders and extract data only when needed in data analysis to save memory and compute time
import os

def walk(home):
    for root, dirs, files, *extra in os.walk(home):
        for file in files:
            cur_file = os.path.join(root, file)
            # Open file and find number of header lines in file
            with open(cur_file, 'r') as full_file:
                return full_file.readlines(), file
            
def close_file(file):
    file.close() 

def extract_postproc(dir_loc, name):
    data, file = walk(dir_loc)
    for line in range(len(data)):
        if data[line][0] != '#':
            header_line = line - 1
            start = line
            break
    headers = data[header_line].split()

# Data extraction for convergence_check.py
# Takes specific data column(s) as raw data bassed on the underlying file structure
# ----------------------------------------------------------------------------
    # For mass flow patches
    if file == 'surfaceFieldValue.dat' and name[-4:] == 'Flow':
        raw_data = [float(row.split()[1]) for row in data[start:]]

            # For residuals file
    elif file == 'solverInfo.dat':
        idx = 0
        for col in headers:
            if col == f'{name}_initial':
                raw_data = [float(row.split()[idx-1]) for row in data[start:]]
            idx += 1

            # For other solution monitors
    elif file == 'volFieldValue.dat':
        idx = 0
        for col in headers:
            if col == f'max({name})' or col == f'min({name})':
                raw_data = [float(row.split()[idx-1]) for row in data[start:]]
            idx += 1

    elif file == 'surfaceFieldValue.dat' and name == 'U':
        idx = 0
        for col in headers:
            if col == f'weightedAreaAverage({name})':
                raw_data = [float(row.split()[idx-1][1:]) for row in data[start:]]
            idx += 1

    return raw_data
