import numpy as np
import os
import sys

# Define the data directories
spec_dir = '../parameters_spectra/'
save_dir = './temp_dir/'

# Create the file list
spec_list = sorted(os.listdir(spec_dir))

# Reshape the file list by passing  in the arguments from SBATCH_chi2_minimization.sh
ncols = int(sys.argv[1])
nrows = int(sys.argv[2])

reshaped_spec_list = np.asarray(spec_list).reshape((nrows, ncols))

# Save this re-shped list to a text file in the temporary directory
np.savetxt(save_dir+'spec_to_fit.txt', reshaped_spec_list, fmt='%s')
