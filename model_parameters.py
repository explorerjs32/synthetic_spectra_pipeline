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

parameters_results_dir = '/blue/rezzeddine/share/fmendez/results/stellar_parameters/'
literature_values_dir = '/blue/rezzeddine/share/fmendez/literature_tables/'
models_out_dir = '/blue/rezzeddine/share/fmendez/results/marcs_models/'

# Define all the input variables
star = sys.argv[1]
marcs_bash = 'model_create.sh'
marcs_script = 'scratch'
sigma = 2.

# Extract the stellar parameters based on measurements or literature values
parameters_source = 'measured'

if parameters_source == 'literature':

    # Open the stellar parameters from the literature
    literature_parameters_df = pd.read_csv(literature_values_dir+'benchmark_stars_literature_parameters.txt', sep='\s+', \
                                           names=['Star','Name','Teff','e_Teff','logg','e_logg','Fe/H','e_Fe/H','turb_vel','e_turb_vel','vsini','e_vsini'])
    star_parameters = literature_parameters_df[literature_parameters_df['Name'] == star]

if parameters_source == 'measured':

    # Open the measured stellar parameters
    samples = np.genfromtxt(parameters_results_dir+star+'_corner_values.txt')
    parameters = [np.percentile(samples[:,i], [16, 50,84]) for i in range(samples.shape[1] - 2)]
    errors = [np.diff(np.percentile(samples[:,i], [16, 50,84])) for i in range(samples.shape[1] -2)]

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
        uncertainties_out.append(uncertainty*sigma)

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
    model_params, model_out = generate_marcs_model(star, star_parameters.iloc[idx], marcs_bash, marcs_script)
