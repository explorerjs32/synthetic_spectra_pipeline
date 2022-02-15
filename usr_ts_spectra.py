import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import copy

def abundance_range(element):
    # Define the bounds and steps for the abundance range
    low = element - .5
    up = element + .5
    step = ((up - low)/.25) + 1

    # Define the abundance abundance range
    abu_range = np.linspace(low, up, int(step))

    return abu_range

def ts_reformat(df):
    # Reformat the turbulent velocity
    df['Turb_vel'] = np.array(['{:.2f}'.format(float(v)) for v in df['Turb_vel']])

    # Reformat abundances
    df['Na'] = np.array(['{:.2f}'.format(float(na)) for na in df['Na']])
    df['Mg'] = np.array(['{:.2f}'.format(float(mg)) for mg in df['Mg']])

    # Compose the name out for the spectra
    names = []

    for i in range(df['Model'].size):
        # Extract the stellat parameters from the model filename
        T = re.split(r'[_gzvt]', df['Model'][i])[1]
        g = re.split(r'[_gzvt]', df['Model'][i])[2]
        z = re.split(r'[_gzvt]', df['Model'][i])[3]
        wmin = str(df['Min_wave'][i])
        wmax = str(df['Max_wave'][i])

        name_out = 'T'+T+'_g'+g+'_z'+z+'_vt'+df['Turb_vel'][i]+'_Na'+df['Na'][i]+'_Mg'+df['Mg'][i]+'_'+wmin+'_'+wmax+'.spec'
        names.append(name_out)

    df['Name_out'] = np.array(names)

    return df

def faltbon_format(files, vels):
    # Define the directory where the files are coming from and where they are going to be stored
    files_in_dir = '../Turbospectrum2019/COM-v19.1/syntspec/'
    files_out_dir = '../synthetic_spectra/'

    faltbon_params = []

    for i in range(len(files)):
        # Split the parameters from the file name
        T = re.split(r'_', files[i])[0]
        g = re.split(r'_', files[i])[1]
        z = re.split(r'_', files[i])[2]
        vt = re.split(r'_', files[i])[3]
        na = re.split(r'_', files[i])[4]
        mg = re.split(r'_', files[i])[5]
        minl = re.split(r'_', files[i])[6]
        maxl = re.split(r'_', files[i])[7]

        for v in vels:
            # Create the file name out of faltbon
            name_out = T+'_'+g+'_'+z+'_'+vt+'_'+'vr'+str(v)+'_'+na+'_'+mg+'_'+minl+'_'+maxl
            faltbon_params.append([files[i], v, name_out])

    # Create a dataframe with these parameters
    df = pd.DataFrame(data=faltbon_params, columns=['Faltbon_name_in', 'Rot_vels', 'Faltbon_name_out'])

    # Reformat the files in by adding the absolute directory path
    df['Faltbon_name_in'] = [files_in_dir+file_in for file_in in df['Faltbon_name_in']]
    df['Faltbon_name_out'] = [files_out_dir+file_out for file_out in df['Faltbon_name_out']]

    return df

# Set up the path and list of all the atmospheric models
model_dir = './test_models/'
model_list = sorted(os.listdir(model_dir))

# Enter the parameters you wnat to use in TurboSpectrum
# NOTE: The difference betweeen lam_ max and lam_min cannot be greater than 5000 A
lam_min = 5000
lam_max = 5500
delta_lam = 0.01
vt = 1.0

# Define the absolute abundances of the Sun
Na = 6.24
Mg = 7.60

# Calculate the abundances range using an arbitrary number
Na_abun = abundance_range(Na)
Mg_abun = abundance_range(Mg)

# Defone the rotational velocity range
rot_vel = np.linspace(5, 50, 10)*-1.

# Create the range of microturbulent volocities
vts = np.linspace(vt, vt + 1., int(1./.2)+1)

# Iterate over every model and extract the parameters required for TurboSpectrum along with all the
# possible combinationsof parameters and abundances
ts_parameters = []

for model in model_list:
    # Extract the metallicity value from the model
    z = re.split(r'[_gzvt]', model)[3]

    # Find all the possible combinations betweend the parameters and abundances
    params = [[model, lam_min, lam_max, delta_lam, z, v, na, mg] \
              for v in vts for na in Na_abun for mg in Mg_abun]

    ts_parameters.append(params)

# Reshape the array
spec_params_arr = np.array(ts_parameters)
spec_params = spec_params_arr.reshape(spec_params_arr.shape[0]*spec_params_arr.shape[1], spec_params_arr.shape[2])

# Create a dataframe of these parameters
columns = ['Model', 'Min_wave', 'Max_wave', 'Delta_wave', 'Z', 'Turb_vel', 'Na', 'Mg']
df = pd.DataFrame(data=spec_params, columns=columns)

# Create a copy of the dataframe and modify the wavelength columns, then append it at the end of dataframe 1
df2 = copy.deepcopy(df)
df2['Min_wave'] = int(df['Max_wave'][0])
df2['Max_wave'] = int(df['Max_wave'][0]) + 500

df = df.append(df2, ignore_index=True)

# Re format the spec parameters
spec_df = ts_reformat(df)

# Create the variations for the rotational velocity
faltbon_df = faltbon_format(spec_df['Name_out'], rot_vel)

# Save the parameters to a text file
save_dir = './temp_dir/'
if not os.path.isdir(save_dir): os.mkdir(save_dir)

np.savetxt(save_dir+'ts_parameters.txt', spec_df.values, fmt='%s')
np.savetxt(save_dir+'faltbon_parameters.txt', faltbon_df.values, fmt='%s')
