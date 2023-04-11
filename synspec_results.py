import numpy as np
import pandas as pd
from subprocess import Popen, PIPE
import time
import os
import shutil
from tqdm.auto import tqdm
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

def generate_marcs_model(star, parameters):

    # Define the grid of pararameters for the MARCS models
    Ts = np.linspace(4500, 6500, 9)
    gs = np.linspace(4.0, 5.0, 3)
    zs = np.linspace(-1.0, 1.0, 9)

    # Get the parameters for the MARCS models
    teff, logg, z = parameters[0][1], parameters[1][1], parameters[2][1]

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
    model_args = ['./model_create.sh scratch %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s > /dev/null 2>&1' \
            %(str(model_df['Tmin'].values[0]), str(model_df['Tmax'].values[0]), model_df['gmin'].values[0], model_df['gmax'].values[0], \
              model_df['zmin'].values[0], model_df['zmax'].values[0], model_df['vt'].values[0], model_df['amin'].values[0], model_df['amax'].values[0], \
              model_df['mass'].values[0], model_df['geom'].values[0], model_df['comp'].values[0], str(model_df['Temp'].values[0]), model_df['logg'].values[0], \
              model_df['z'].values[0], model_df['Name'].values[0])]

    model_process = Popen(model_args, shell=True).wait()

    # Define the model name
    model = model_df['Name'].values[0]+'_'+str(model_df['Temp'].values[0])+'g'+model_df['logg'].values[0]+'z'+model_df['z'].values[0]+'vt01.TURBO'

    return model_df, model

def generate_synthetic_spectra(model_params, stellar_params, abund, model, element, atomic_number, linelist_dir):

    # Define the necesary parameters for the synthetic spectra
    model_out = model
    min_wl, max_wl, dwl = '5000', '6000', '0.01'
    z, turvel = str(model_params['z'].values[0]), str(round(stellar_params[3][1], 2))
    abundance = abund
    synth_spec_out = 'T'+str(model_params['Temp'].values[0])+'g'+str(model_params['logg'].values[0])+'z'+str(model_params['z'].values[0])+'vt'+turvel+'_'+element+str(abundance)+'_'+min_wl+'-'+max_wl+'.spec'
    linelist = element+'_lines.list'

    # Generate the synthetic spectra
    synthspec_args = ['./abund_spec_create.sh ts_script_abu.com %s %s %s %s %s %s %s %s %s %s > /dev/null 2>&1' \
                      %(model_out, min_wl, max_wl, dwl, z, turvel, synth_spec_out, str(atomic_number), \
                      str(abundance), str(linelist_dir+linelist))]
    synthspec_process = Popen(synthspec_args, shell=True).wait()

    return synth_spec_out

def spec_broadening(stellar_params, spec_out, temp_dir):

    # Define the broadening values
    inst_broadening = -1.2
    rv_broadening = round(stellar_params[3][1], 2)*-1.

    # Get the name for the broadened spectra
    split_name = spec_out.split('.spec')[0]
    spec_out_inst = split_name+'_inst'+str(inst_broadening)+'.spec'
    spec_out_rv = split_name+'_rv'+str(rv_broadening)+'.spec'

    # Apply instrumental broadening
    inst_spec_args = ['./spectra_smooth.sh inst_broadening.f90 %s %s %s inst_broadening > /dev/null 2>&1' \
                       %(temp_dir+spec_out, temp_dir+spec_out_inst, str(inst_broadening))]
    inst_process1 = Popen(inst_spec_args, shell=True).wait()

    # Apply instrumental broadening
    vr_args = ['./spectra_smooth.sh rot_vel.f90 %s %s %s rot_vel > /dev/null 2>&1' \
                       %(temp_dir+spec_out_inst, temp_dir+spec_out_rv, str(rv_broadening))]
    vr_process1 = Popen(vr_args, shell=True).wait()

    return spec_out_rv


# Define the directory paths
temp_dir = '/blue/rezzeddine/share/fmendez/temp_dir/'
synspec_save_dir = '/blue/rezzeddine/share/fmendez/results/synthetic_spectra/'
model_save_dir = '/blue/rezzeddine/share/fmendez/results/marcs_models/'
parameters_results_dir = '/blue/rezzeddine/share/fmendez/results/stellar_parameters/'
linelist_dir = '/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/linelists/'
abundances_save_dir = '/blue/rezzeddine/share/fmendez/results/abundances/'
atomic_numbers = pd.read_csv('./atomic_numbers.csv')

# Define the stars and elements you want to generate the synthetic spectra for
stars = ['Sun']
elements = ['Mg']

for star in tqdm(stars, position=0, leave=True):

    print(star)

    # Get the measured stellar parameters
    samples = np.genfromtxt(parameters_results_dir+star+'_corner_values.txt')
    stellar_parameters = [np.percentile(samples[:,i], [16, 50,84]) for i in range(samples.shape[1])]

    # Generate the MARCS model
    model_params, model_out = generate_marcs_model(star, stellar_parameters)

    for element in elements:

        # Get the atomic number for each element
        atomic_num = atomic_numbers[atomic_numbers['Symbol'] == element]['AtomicNumber'].values[0]

        # Extract the abundance value
        abundances_results = pd.read_csv(abundances_save_dir+star+'_'+element+'_abundances.csv')
        abu = abundances_results.sort_values(by='Chi').reset_index(drop=True)['Abundance'].astype(float).values[0]

        # Generate the synthetic spectra
        synth_spec_out = generate_synthetic_spectra(model_params, stellar_parameters, round(abu, 2), model_out, element, atomic_num, linelist_dir)

        # Apply instrumental and rotation broadening
        final_synth_spec = spec_broadening(stellar_parameters, synth_spec_out, temp_dir)

        # Move the final spectrum to the results folder
        mv = Popen(['mv %s %s' %(temp_dir+final_synth_spec, synspec_save_dir)], shell=True).wait()
        rename = Popen(['mv %s %s' %(synspec_save_dir+final_synth_spec, synspec_save_dir+star+'_'+final_synth_spec)], shell=True).wait()

        # Remove all the extra files
        rm = Popen(['rm %s' %(temp_dir+'*')], shell=True).wait()
