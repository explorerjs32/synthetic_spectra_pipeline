import numpy as np
import os
from astropy import constants as const
# Define the path th the data
files_in_dir = '../smoothed_synthetic_spectra/'
files_out_dir = '../synthetic_spectra/'
temp_dir = './temp_dir/'


# Get the list of files
files_in = np.array([os.path.join(files_in_dir, f) for f in sorted(os.listdir(files_in_dir))])
files_out_txt = np.array([os.path.join(files_out_dir, f+'.txt') for f in sorted(os.listdir(files_in_dir))])
files_out_bin = np.array([os.path.join(files_out_dir, f) for f in sorted(os.listdir(files_in_dir))])

files_in = []
files_out_txt = []
files_out_bin = []

for f in sorted(os.listdir(files_in_dir)):

    if not os.path.isfile(files_out_dir+f):

        files_in.append(os.path.join(files_in_dir, f))
        files_out_txt.append(os.path.join(files_out_dir, f+'.txt'))
        files_out_bin.append(os.path.join(files_out_dir, f))

files_in = np.array(files_in)
files_out_txt = np.array(files_out_txt)
files_out_bin = np.array(files_out_bin)

# Calculate the instrumental broadening
inst_res = 120000.
dv_inst = (const.c.to('km/s').value/inst_res)*-1.

inst_val = np.zeros(files_out_txt.size) + dv_inst

# Save the data to a txt file
data_out = np.column_stack([files_in, files_out_txt, files_out_bin, inst_val])

np.savetxt(temp_dir+'inst_broad_files_fun.txt', data_out, fmt='%s')
