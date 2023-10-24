import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from subprocess import Popen, PIPE
import time
from tqdm.auto import tqdm
from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler
from astropy import units as u
import os
import sys
import argparse

import warnings
warnings.filterwarnings("ignore")

def marcs_reformat(df):
    # Re-format the temperatures
    df['Tmin'] = np.array([int(float(t)) for t in df['Tmin']])
    df['Tmax'] = np.array([int(float(t)) for t in df['Tmax']])
    df['Temp'] = np.array([int(float(t)) for t in df['Temp']])

    # Re-format the gravities
    df['gmin'] = np.array(['+'+g for g in df['gmin']])
    df['gmax'] = np.array(['+'+g for g in df['gmax']])
    df['logg'] = np.array(['+'+str(round(float(g), 2)) for g in df['logg']])

    # Re-format the metallicities
    for i in range(len(df['z'])):
        if float(df['zmin'][i]) < 0.: df['zmin'][i] = '{:.2f}'.format(float(df['zmin'][i]))
        else: df['zmin'][i] = '+{:.2f}'.format(float(df['zmin'][i]))
        if float(df['zmax'][i]) < 0.: df['zmax'][i] = '{:.2f}'.format(float(df['zmax'][i]))
        else: df['zmax'][i] = '+{:.2f}'.format(float(df['zmax'][i]))
        if float(df['z'][i]) < 0.: df['z'][i] = '{:.2f}'.format(float(df['z'][i]))
        else: df['z'][i] = '+{:.2f}'.format(float(df['z'][i]))

    # Re-format the turbulent velocity
    df['vt'] = np.array(['0'+v for v in df['vt']])

    return df

def generate_marcs_model(star, parameters, marcs_bash, marcs_script):

    # Define the grid of pararameters for the MARCS models
    Ts = np.linspace(4500, 6500, 9)
    gs = np.linspace(4.0, 5.0, 3)
    zs = np.linspace(-1.0, 1.0, 9)

    # Get the parameters for the MARCS models
    teff, logg, z = parameters['Teff'], parameters['logg'], parameters['Fe/H']

    # Find the closest interpolation parameters within the grid space
    min_teff, max_teff = Ts[np.argsort(abs(Ts - teff))[0]], Ts[np.argsort(abs(Ts - teff))[1]]
    min_logg, max_logg = 4.0, 5.0
    min_z, max_z = -1.0, 1.0

    # Create the alpha abundance dictionary and define the rest of the parameters for MARCS
    alphas = {'1.00' : '+0.00',
              '0.75' : '+0.00',
              '0.50' : '+0.00',
              '0.25' : '+0.00',
              '0.00' : '+0.00',
              '-0.25' : '+0.10',
              '-0.50' : '+0.20',
              '-0.75' : '+0.30',
              '-1.00' : '+0.40',
              '-1.50' : '+0.40',
              '-2.00' : '+0.40',
              '-2.50' : '+0.40',
              '-3.00' : '+0.40',
              '-4.00' : '+0.40',
              '-5.00' : '+0.40'}
    min_a, max_a = alphas['{:.2f}'.format(min_z)], alphas['{:.2f}'.format(max_z)]
    turb_vel = 1
    mass = '0.0'
    geom = 'p'
    comp = 'st'

    # Create the data frame and re-format the marameters
    columns = ['Name','Tmin', 'Tmax', 'gmin', 'gmax', 'zmin', 'zmax', 'vt', 'amin', 'amax', 'mass', 'geom', 'comp', 'Temp', 'logg', 'z']
    data = np.column_stack([star, min_teff, max_teff, min_logg, max_logg, min_z, max_z, turb_vel, min_a, max_a, mass, geom, comp, teff, logg, z])

    model_df = pd.DataFrame(data=data, columns=columns)

    # Re-format each of the colums of the columns so it has the MARCS model format
    model_df = marcs_reformat(model_df)

    # Create the atmospheric model
    model_args = ['./%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s > /dev/null 2>&1' \
            %(str(marcs_bash), str(marcs_script), str(model_df['Tmin'].values[0]), str(model_df['Tmax'].values[0]), model_df['gmin'].values[0], model_df['gmax'].values[0], \
              model_df['zmin'].values[0], model_df['zmax'].values[0], model_df['vt'].values[0], model_df['amin'].values[0], model_df['amax'].values[0], \
              model_df['mass'].values[0], model_df['geom'].values[0], model_df['comp'].values[0], str(model_df['Temp'].values[0]), model_df['logg'].values[0], \
              model_df['z'].values[0], model_df['Name'].values[0])]

    model_process = Popen(model_args, shell=True).wait()

    # Define the model name
    model = model_df['Name'].values[0]+'_T'+str(model_df['Temp'].values[0])+'g'+model_df['logg'].values[0]+'z'+model_df['z'].values[0]+'vt01.TURBO'

    return model_df, model

# Define all the parameters to parse into the code
parser = argparse.ArgumentParser(description='Parameters to parse to create the MARCS models for a set of stellar parameters')

parser.add_argument('-s', '--s', '--star_name', required=True, type=str, help='The star name to calculate the atmospheric models for')
parser.add_argument('-pd', '--pd', '-parameters_dir', default='/blue/rezzeddine/share/fmendez/results/stellar_parameters/', type=str, help='Directory where the stellar parameters are located')
parser.add_argument('-mod', '--mod', '-model_out_dir', default='/blue/rezzeddine/share/fmendez/results/marcs_models/', type=str, help='Directory where the output MARCS models are saved')
parser.add_argument('-mbs', '--mbs', '--marcs_bash_script', default='model_create.sh', type=str, help='Bash script to run the MARCS models script')
parser.add_argument('-ms', '--ms', '--marcs_script', default='scratch', type=str, help='Script to generate the MARCS models')

args = parser.parse_args()


# Open the asplund abundances
asplund_solar_abund = pd.read_csv('./asplund_solar_abund.csv', sep='\s+', names=['Atomic_Number','Element','Abundance','Abundance_err'])

# Open the measured stellar parameters
samples = np.genfromtxt('%s%s_corner_values_v3.txt' %(args.pd, args.s))
parameters = [np.percentile(samples[:,i], [16, 50,84]) for i in range(samples.shape[1] - 2)]
errors = [np.diff(np.percentile(samples[:,i], [16, 50,84])) for i in range(samples.shape[1] -2)]

# Define the directory where the iron abundances are being stored and list all the files in it
fe_abund_dir = '/blue/rezzeddine/share/fmendez/results/abundances/fe_abundances/'
measured_fe_abunds_list = sorted(os.listdir(fe_abund_dir))

# Check if the Fe abundances have been measured already
if '%s_Fe_abundance.csv' %args.s in measured_fe_abunds_list:

    print('Fe abundance has been measured! Re-write Fe/H')

    # Update the metallicity parameter
    asplund_abund = asplund_solar_abund[asplund_solar_abund['Element'].isin(['Fe'])][['Abundance','Abundance_err']]   
    measured_abund = pd.read_csv('%s%s_Fe_abundance.csv' %(fe_abund_dir, args.s), sep='\s+')

    new_met = measured_abund['Abundance'].values - asplund_abund['Abundance'].values

    parameters[2][1] = new_met

else:
    pass

# Save the parameters and respective uncertainties
parameters_out = []
uncertainties_out = []

for i in range(len(parameters)):

    # Get the index of largest uncertainty
    argmax_idx = np.array(errors[i]).argmax()

    # Add a negative or positive sign to the high uncertainty depending if it's an upper or lower value
    if argmax_idx == 0:
        uncertainty = errors[i][argmax_idx]*-1.

    if argmax_idx == 1:
        uncertainty = errors[i][argmax_idx]

    # Save the uncertainties and parameters to a list
    parameters_out.append(parameters[i][1])
    uncertainties_out.append(uncertainty)

# Add the uncertaintis to the parameters keeping one of them constant
parameters_and_uncertainties = [np.array(parameters_out)]

for i in range(len(parameters_out)):

    # Create a copy of the stellar parameters]
    parameters_to_edit = np.array(parameters_out)

    # Update the parameter value adding the uncertainty
    parameters_to_edit[i] = parameters_and_uncertainties[0][i] + uncertainties_out[i]

    # Add the modified parameters to the list
    parameters_and_uncertainties.append(parameters_to_edit)

# Convert the modified parameters to a DataFrame
star_parameters = pd.DataFrame(np.array(parameters_and_uncertainties), \
                               columns=['Teff','logg','Fe/H','turb_vel','vsini'])

for idx in tqdm(star_parameters.index):

    # Generate the MARCS atmospheric model
    model_params, model_out = generate_marcs_model(args.s, star_parameters.iloc[idx], args.mbs, args.ms)

print(star_parameters.iloc[0])