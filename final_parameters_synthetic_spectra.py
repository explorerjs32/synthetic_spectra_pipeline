import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from subprocess import Popen, PIPE
import time
from tqdm import tqdm
from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler
from astropy import units as u
import os
import sys
import shutil
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

def generate_marcs_model(star, parameters):

    # Define the grid of pararameters for the MARCS models
    Ts = np.linspace(4500, 6500, 9)
    gs = np.linspace(4.0, 5.0, 3)
    zs = np.linspace(-1.0, 1.0, 9)

    # Get the parameters for the MARCS models
    teff, logg, z = parameters['Teff'].values[0], parameters['logg'].values[0], parameters['Fe/H'].values[0]

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

    # Define the model name
    model = model_df['Name'].values[0]+'_T'+str(model_df['Temp'].values[0])+'g'+model_df['logg'].values[0]+'z'+model_df['z'].values[0]+'vt01.TURBO'

    return model_df, model

def generate_synthetic_spectra(model_params, stellar_params, model, ts_bash, ts_script):

    # Define the necesary parameters for the synthetic spectra
    model_out = model
    min_wl, max_wl, dwl = '5000', '6000', '0.01'
    z, turvel, na_abund, mg_abund = str(model_params['z'].values[0]), str(round(stellar_params['turb_vel'].values[0], 2)), str(round(stellar_params['Na'].values[0], 2)), str(round(stellar_params['Mg'].values[0], 2))
    
    synth_spec_out = 'T'+str(model_params['Temp'].values[0])+'g'+str(model_params['logg'].values[0])+'z'+str(model_params['z'].values[0])+'vt'+str(turvel)+'_'+min_wl+'-'+max_wl+'.spec'
    
    # Generate the synthetic spectra
    synthspec_args = ['./%s %s %s %s %s %s %s %s %s %s %s > /dev/null 2>&1' \
                      %(str(ts_bash), str(ts_script), model_out, min_wl, max_wl, dwl, z, turvel, synth_spec_out, na_abund, mg_abund)]
    synthspec_process = Popen(synthspec_args, shell=True).wait()

    return synth_spec_out


def spec_broadening(stellar_params, spec_out, temp_dir, broadening_bash, inst_broad_script, rot_vel_broad_script):

    # Define the broadening values
    inst_broadening = -1.2
    rv_broadening = round(stellar_params['vsini'].values[0], 2)*-1.

    # Get the name for the broadened spectra
    split_name = spec_out.split('.spec')[0]
    spec_out_inst = split_name+'_inst'+str(inst_broadening)+'.spec'
    spec_out_rv = split_name+'_rv'+str(rv_broadening)+'.spec'

    # Apply instrumental broadening
    inst_spec_args = ['./%s %s.f90 %s %s %s %s > /dev/null 2>&1' \
                       %(str(broadening_bash), str(inst_broad_script), \
                         temp_dir+spec_out, temp_dir+spec_out_inst, str(inst_broadening), \
                         str(inst_broad_script))]
    inst_process1 = Popen(inst_spec_args, shell=True).wait()

    # Apply instrumental broadening
    vr_args = ['./%s %s.f90 %s %s %s %s > /dev/null 2>&1' \
                       %(str(broadening_bash), str(rot_vel_broad_script), \
                         temp_dir+spec_out_inst, temp_dir+spec_out_rv, str(rv_broadening), \
                         str(rot_vel_broad_script))]
    vr_process1 = Popen(vr_args, shell=True).wait()

    return spec_out_rv

# Define the arguments to be parsed
parser = argparse.ArgumentParser(description='Parameters to parse to generate synthetic spectra')

parser.add_argument('-s', '--s', '--star_name', required=True, type=str, help='The star to derive the chemical abundances for')
parser.add_argument('-spd', '--spd', '--stellar_parameters_dir', default='/blue/rezzeddine/share/fmendez/results/stellar_parameters/', type=str, help='Directory where the stellar parameters are located')
parser.add_argument('-lld', '--lld', '--line_list_dir', default='/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/linelists/', type=str, help='Directory where the line list to generate the synthetic spectrum is stored')
parser.add_argument('-tod', '--tod', '--temporary_output_dir', default='/blue/rezzeddine/share/fmendez/temp_dir/', type=str, help='Directory where the temporary output synthetic spectra are stored')
parser.add_argument('-tsbs', '--tsbs', '--turbo_spectrum_bash_script', default='spectra_create.sh', type=str, help='Bash script to run the TurboSpectrum script')
parser.add_argument('-tss', '--tss', '--turbo_spectrum_script', default='ts_script.com', type=str, help='TurboSpectrum script to generate synthetic spectra')
parser.add_argument('-bbs', '--bbs', '--broadening_bash_script', default='spectra_smooth.sh', type=str, help='Bash script to run the script to apply broadening by instrumentation')
parser.add_argument('-ibs', '--ibs', '--instrumental_broadening_script', default='inst_broadening', type=str, help='Script to apply instrumental broadening')
parser.add_argument('-rbs', '--rbs', '--rotational_broadening_script', default='rot_vel', help='Script to apply rotational broadening')
parser.add_argument('-ssd','--ssd','--synthetic_spectra_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/synthetic_spectra/stellar_parameters/', help='Directory where to save the synthetic spectra')

args = parser.parse_args()

# Define additional variables
wl_range = [4200., 6600.]

# Open the measured stellar parameters
samples = np.genfromtxt('%s%s_corner_values.txt' %(args.spd, args.s))
parameters = np.array([np.percentile(samples[:,i], [16, 50,84])[1] for i in range(samples.shape[1])])

# Convert the list of parameters to a dataframe
star_parameters = pd.DataFrame(parameters.reshape(1,-1), columns=['Teff','logg','Fe/H','turb_vel','vsini','Na','Mg'])

# Generate the MARCS atmospheric model
model_params, model_out = generate_marcs_model(args.s, star_parameters)

# Generate the synthetic spectra
synth_spec_out = generate_synthetic_spectra(model_params, star_parameters, model_out, args.tsbs, args.tss)

# Apply instrumental and rotation broadening
final_synth_spec = spec_broadening(star_parameters, synth_spec_out, args.tod, args.bbs, args.ibs, args.rbs)

# Copy the final synthetic spectra to the results directory
shutil.copy(args.tod+final_synth_spec, args.ssd+args.s+'_'+final_synth_spec)

# Remove all the extra files
remove_synthspec = Popen(['rm %s*' %(args.tod)], shell=True).wait()
