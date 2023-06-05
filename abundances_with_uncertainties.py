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

    # Define the model name
    model = model_df['Name'].values[0]+'_T'+str(model_df['Temp'].values[0])+'g'+model_df['logg'].values[0]+'z'+model_df['z'].values[0]+'vt01.TURBO'

    return model_df, model

def generate_synthetic_spectra(model_params, stellar_params, abund, model, element, atomic_number, linelist_dir, ts_bash, ts_script, temp_dir):

    # Define the necesary parameters for the synthetic spectra
    model_out = model
    min_wl, max_wl, dwl = '4200', '6600', '0.01'
    z, turvel = str(model_params['z'].values[0]), str(round(stellar_params['turb_vel'], 2))
    abundance = abund
    synth_spec_out = 'T'+str(model_params['Temp'].values[0])+'g'+str(model_params['logg'].values[0])+'z'+str(model_params['z'].values[0])+'vt'+turvel+'_'+element+str(abundance)+'_'+min_wl+'-'+max_wl+'.spec'
    linelist = element+'_lines.list'
    print(synth_spec_out, model_out)
    # Generate the synthetic spectra
    synthspec_args = ['./%s %s %s %s %s %s %s %s %s %s %s %s > /dev/null 2>&1' \
                      %(str(ts_bash), str(ts_script), model_out, min_wl, max_wl, dwl, z, turvel, synth_spec_out, str(atomic_num), \
                      str(abundance), str(linelist_dir+linelist))]
    synthspec_process = Popen(synthspec_args, shell=True).wait()

    return synth_spec_out

def spec_broadening(stellar_params, spec_out, temp_dir, broadening_bash, inst_broad_script, rot_vel_broad_script):

    # Define the broadening values
    inst_broadening = -1.2
    rv_broadening = round(stellar_params['vsini'], 2)*-1.

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

def spec_resampler(obs_spec, sample_spec):
    '''
    The observed spectrum is not expected to be sampled in the same grid space as the synthetic spectra.
    In order to perform an acurate pixel to pixel fit we want to re-sample the observed spectrum to the
    same spectral grid than our library of synthetic spetra.
    '''
    # Convert the spectrum dataframe into the SpecUtils format
    obs_spec = Spectrum1D(spectral_axis=obs_spec['Wave'].values*u.AA, flux=obs_spec['Flux'].values*u.Jy)
    mock_spec = Spectrum1D(spectral_axis=sample_spec['Wave'].values*u.AA, flux=sample_spec['Wave'].values*u.Jy)

    # Re-sample the observe spectrum using linear interpolation
    new_spec_grid = np.linspace(4200, 6600, sample_spec.Wave.size)
    linear = LinearInterpolatedResampler()
    obs_spec_resamp = linear(obs_spec, new_spec_grid*u.AA)

    # Return a datafram of the re-sampled spectrum
    resamp_spec_df = pd.DataFrame(data=np.column_stack([obs_spec_resamp.spectral_axis.value, obs_spec_resamp.flux.value]), columns=['Wave','Flux'])

    return resamp_spec_df

def line_masker(obs_spec_df, synth_spec_df, spec_mask_df):
    '''
    This function masks the spectral lines that are going to be considered for the LM fit.
    It takes in both observed and synthetic spectra as parameters, as well as a pre-selected
    list of lines to mask. The output is the mask observed and synthetic spectra
    '''

    owl, ofl = obs_spec_df['Wave'].values, obs_spec_df['Flux'].values
    swl, sfl = synth_spec_df['Wave'].values, synth_spec_df['Flux'].values

    obs_regions = []
    synth_regions = []

    for i in range(spec_mask_df['Wave'].size):
        # Define the mask for each line
        lower = spec_mask_df['Wave'][i] - spec_mask_df['Region'][i]
        upper = spec_mask_df['Wave'][i] + spec_mask_df['Region'][i]

        obs_line_mask = (owl >= lower) & (owl <= upper)
        synth_line_mask = (swl >= lower) & (swl <= upper)

        mowl, mofl = owl[obs_line_mask], ofl[obs_line_mask]
        mswl, msfl = swl[synth_line_mask], sfl[synth_line_mask]

        # Clip the spectra with the highest amout of pixels
        if owl[obs_line_mask].size > swl[synth_line_mask].size:
            dif = abs(owl[obs_line_mask].size - swl[synth_line_mask].size)
            mowl = owl[obs_line_mask][:-dif]
            mofl = ofl[obs_line_mask][:-dif]

        if owl[obs_line_mask].size < swl[synth_line_mask].size:
            dif = abs(owl[obs_line_mask].size - swl[synth_line_mask].size)
            mswl = swl[synth_line_mask][:-dif]
            msfl = sfl[synth_line_mask][:-dif]

        obs_region = np.column_stack([mowl, mofl])

        synth_region = np.column_stack([mswl, msfl])

        obs_regions.append(obs_region)
        synth_regions.append(synth_region)

    mobs_spec = np.concatenate(obs_regions)
    msynth_spec = np.concatenate(synth_regions)

    # Create a dataframe for the masked spectra
    mobs_spec_df = pd.DataFrame(data=mobs_spec, columns=['Wave','Flux'])
    msynth_spec_df = pd.DataFrame(data=msynth_spec, columns=['Wave','Flux'])

    return mobs_spec_df, msynth_spec_df

# Define the arguments to be parsed
parser = argparse.ArgumentParser(description='Parameters to parse to derive chemical abundances')

parser.add_argument('-s', '--s', '--star_name', required=True, type=str, help='The star to derive the chemical abundances for')
parser.add_argument('-e', '--e', '--element', required=True, type=str, help='The element to measure the chemical abundance for')
parser.add_argument('-a', '--a', '--abundance', required=True, type=float, help='The abundance value to generate the synthetic spectrum for')
parser.add_argument('-nsd', '--nsd', '--normalized_spectra_dir', default='/blue/rezzeddine/share/fmendez/normalized_stellar_spectra/', type=str, help='Directory where the normalized stellar spectra is stored')
parser.add_argument('-spd', '--spd', '--stellar_parameters_dir', default='/blue/rezzeddine/share/fmendez/results/stellar_parameters/', type=str, help='Directory where the stellar parameters are located')
parser.add_argument('-lld', '--lld', '--line_list_dir', default='/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/linelists/', type=str, help='Directory where the line list to generate the synthetic spectrum is stored')
parser.add_argument('-tod', '--tod', '--temporary_output_dir', default='/blue/rezzeddine/share/fmendez/temp_dir/', type=str, help='Directory where the temporary output synthetic spectra are stored')
parser.add_argument('-od', '--od', '--output_dir', default='/blue/rezzeddine/share/fmendez/results/abundances/temp_abund/', type=str, help='Output directory to save the measured chemical abundances')
parser.add_argument('-slm', '--slm', '--spectral_line_mask', default='/blue/rezzeddine/share/fmendez/atomic_line_mask/chemical_abundances_lines.csv', type=str, help='List of atomic lines used to chemical abundances')
parser.add_argument('-tsbs', '--tsbs', '--turbo_spectrum_bash_script', default='abund_spec_create.sh', type=str, help='Bash script to run the TurboSpectrum script')
parser.add_argument('-tss', '--tss', '--turbo_spectrum_script', default='ts_script_abu.com', type=str, help='TurboSpectrum script to generate synthetic spectra')
parser.add_argument('-bbs', '--bbs', '--broadening_bash_script', default='spectra_smooth.sh', type=str, help='Bash script to run the script to apply broadening by instrumentation')
parser.add_argument('-ibs', '--ibs', '--instrumental_broadening_script', default='inst_broadening', type=str, help='Script to apply instrumental broadening')
parser.add_argument('-rbs', '--rbs', '--rotational_broadening_script', default='rot_vel', help='Script to apply rotational broadening')
parser.add_argument('-ps', '--ps', '--parameter_source', default='measured', type=str, help='Source of stellar parameters')
parser.add_argument('-S', '--S', '--sigma', default=1., type=float, help='Sigma factor to add to the parameters uncertainty')

args = parser.parse_args()


# Define the different directory paths
# observed_spectrum_dir = '/blue/rezzeddine/share/fmendez/normalized_stellar_spectra/' D
# parameters_results_dir = '/blue/rezzeddine/share/fmendez/results/stellar_parameters/' D
# linelist_dir = '/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/linelists/' D
# temp_dir = '/blue/rezzeddine/share/fmendez/temp_dir/' D
# save_dir = '/blue/rezzeddine/share/fmendez/results/abundances/temp_abund/' D

# Define the input variables
# marcs_bash = sys.argv[1] D
# marcs_script = sys.argv[2] D
# ts_bash = sys.argv[3] D
# ts_script = sys.argv[4] D
# abundance = np.round(float(sys.argv[5]), 2) D
# broadening_bash = sys.argv[6] D
# inst_broad_script = sys.argv[7] D
# inst_broad_script_out = sys.argv[8] D
# rot_vel_broad_script = sys.argv[9] D
# rot_vel_broad_script_out = sys.argv[10] D
# star = sys.argv[11] D
# element = sys.argv[12] D

# Define additional variables
wl_range = [4200., 6600.]
# parameters_source = 'measured' D
# sigma = 1. D

# Open the data we are going to use
asplund_solar_abund = pd.read_csv('./asplund_solar_abund.csv', sep='\s+', names=['Atomic_Number','Element','Abundance','Abundance_err'])
atomic_numbers = pd.read_csv('./atomic_numbers.csv')
atomic_num = atomic_numbers[atomic_numbers['Symbol'] == args.e]['AtomicNumber'].values[0]

if args.ps == 'measured':

    # Open the measured stellar parameters
    samples = np.genfromtxt('%s%s_corner_values.txt' %(args.spd, args.s))
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
        uncertainties_out.append(uncertainty*args.S)

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

# Open the observed spectrum and mask it based on the wavelength range
obs_spec = pd.read_csv('%s%s_spectrum.dat' %(args.nsd, args.s), delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python').fillna(1.0)
spec_mask = ((obs_spec['Wave'] >= wl_range[0]) & (obs_spec['Wave'] <= wl_range[1]))

obs_spec = obs_spec[spec_mask].reset_index(drop=True)

# Open the abundance line list and mask it based on element and wavelength range
if args.ps == 'measured':
    literature_star_linelist = pd.read_csv(args.slm, sep='\s+', names=['Wave','Potential','loggf','Element','Region']).fillna(0.5)
    linelist_mask = ((literature_star_linelist['Wave'] >= wl_range[0]) & (literature_star_linelist['Wave'] <= wl_range[1]))

    linelist_out = literature_star_linelist[linelist_mask].reset_index(drop=True)

if args.e+'I' in linelist_out['Element'] and args.e+'II' in linelist_out['Element']:
    lines_out = linelist_out[((linelist_out['Element'] == args.e+'I') & (linelist_out['Element'] == args.e+'II'))].reset_index(drop=True)

else:
    lines_out = linelist_out[linelist_out['Element'] == args.e+'I'].reset_index(drop=True)

# Iterate over all the combinations of stellar parameters and perform the abundance measurements
synspec_files = []
elements = []
abundances = []
chis2 = []

for idx in tqdm(star_parameters.index):

    # Generate the MARCS atmospheric model
    model_params, model_out = generate_marcs_model(args.s, star_parameters.iloc[idx])

    # Generate the synthetic spectra
    synth_spec_out = generate_synthetic_spectra(model_params, star_parameters.iloc[idx], round(args.a, 2), model_out, args.e, atomic_num, args.lld, args.tsbs, args.tss, args.tod)

    # Apply instrumental and rotation broadening
    final_synth_spec = spec_broadening(star_parameters.iloc[idx], synth_spec_out, args.tod, args.bbs, args.ibs, args.rbs)

    # Open the synthetic spectrum and re-sample it
    synth_spec = pd.read_csv('%s%s' %(args.tod, final_synth_spec), sep='\s+', names=['Wave','Flux','Int'])
    resamp_spec = spec_resampler(obs_spec, synth_spec)

    # Mask the absorption lines to use for the fit
    masked_obs_spec, masked_synth_spec = line_masker(resamp_spec, synth_spec, lines_out)

    # Calculate the chi^2
    chi2 = np.sum(((masked_obs_spec['Flux'].values - masked_synth_spec['Flux'].values)**2.)/masked_synth_spec['Flux'].values)

    # Add all the calculations to the lists
    synspec_files.append(final_synth_spec)
    elements.append(args.e)
    abundances.append(round(args.a,2))
    chis2.append(chi2)

# Save the data into a DataFrame
results_out = pd.DataFrame()
results_out['Syntspec_File'] = synspec_files
results_out['Element'] = elements
results_out['Abundance'] = abundances
results_out['Chi2'] = chis2

# Save the results
# results_out.to_csv(save_dir+star+'_'+element+'_'+str(abundance)+'_v3_2.csv', index=False, sep=' ')
results_out.to_csv('%s%s_%s_%s.csv' %(args.od, args.s, args.e, str(round(args.a, 2))), index=False, sep=' ')