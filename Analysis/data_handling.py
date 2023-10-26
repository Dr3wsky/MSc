# Module extract_data.py
# Used to walk through file folders and extract data only when needed in data analysis to save memory and compute time
import os
import numpy as np
import pandas as pd


def extract(dir_loc):
    for root, dirs, files, *extra in os.walk(dir_loc):        
            for file in files:
                    cur_file = os.path.join(root, file)
                    # Open file and find number of header lines in file
                    with open(cur_file, 'r') as full_file:
                        data = full_file.readlines()
                        for line in range(len(data)):
                            if data[line][0] != '#':
                                start = line 
                                break
                            
                            
