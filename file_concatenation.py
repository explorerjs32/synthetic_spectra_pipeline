import numpy as np
import os

# Define the path to where the spectra is stored
spec_dir = '../sample_synth_spec/'

# Define the empty list to store the
spec_5000_5500 = []
spec_5500_6000 = []
spec_out = []

for file in os.listdir(spec_dir):
    # Define the conditions to store the file names to e=its respective list
    if file.endswith('5500_6000.spec'): spec_5500_6000.append(file)

    if file.endswith('5000_5500.spec'):
        # Split the file name
        name = file.split('5000_5500.spec')[0]
        sufix = '5000_6000.spec'

        spec_5000_5500.append(file)
        spec_out.append(name+sufix)

# create a table of the sorted file names
table_out = np.column_stack([spec_5000_5500.sort(), spec_5500_6000.sort(), spec_out.sort()])

# Svae the table as a text file
np.savetxt('./temp_dir/concatenation_spec.txt', table_out)
