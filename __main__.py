import os
import numpy as np

print(f'############################################## \nWELCOME TO THE SYNTHETIC SPECTRAL PIPELINE \n\
############################################## \n')

# Decide whether to manually input the stellar parameters or use the MARCS default parameters
# for the spectra Generation
usr_input = str(input('Do you want to manualy input the stellar parameters? (Y/N): '))

if usr_input == 'Y':
    # Run the codes that creates the marcs models
    print('Create the atmospheric models')
    os.system('python3 usr_marcs_models.py')
    os.system('sbatch SBATCH_model_create.sh')
    print('Done!')

    # Run the code that creates the synthetic spectra and smoothing parameters
    print('Create the synthetic spectra')
    os.system('python3 usr_ts_spectra.py')
    os.system('sbatch SBATCH_spectra_create.sh')
    print('Done!')

    print('Create the smotthed spectra')
    os.system('sbatch SBATCH_spectra_smooth.sh')
    print('Done!')

elif usr_input == 'N': os.system('python3 black_box_parameters.py')
