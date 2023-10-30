# Module extract_data.py
# Used to walk through file folders and extract data only when needed in data analysis to save memory and compute time
import os

def extract(dir_loc, name):
    for root, dirs, files, *extra in os.walk(dir_loc):
        for file in files:
            cur_file = os.path.join(root, file)
            # Open file and find number of header lines in file
            with open(cur_file, 'r') as full_file:
                data = full_file.readlines()
                for line in range(len(data)):
                    if data[line][0] != '#':
                        header_line = line - 1
                        start = line
                        break
            headers = data[header_line].split()
            # Split data depending on underlying type, or folder name
            raw_data = [float(row.split()[1]) for row in data[start:]]
            full_file.close()

            return raw_data
