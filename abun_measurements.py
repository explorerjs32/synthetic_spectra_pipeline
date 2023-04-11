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
                      %(model_out, min_wl, max_wl, dwl, z, turvel, synth_spec_out, str(atomic_num), \
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

def spec_resampler(obs_spec, sample_spec):
    '''
    The observed spectrum is not expected to be sampled in the same grid space as the synthetic spectra.
    In order to perform an acurate pixel to pixel fit we want to re-sample the observed spectrum to the
    same spectral grid than our library of synthetic spetra.
    '''
    # Convert the spectrum dataframe into the SpecUtils format
    obs_spec = Spectrum1D(spectral_axis=obs_spec.Wave.values*u.AA, flux=obs_spec.Flux.values*u.Jy)
    mock_spec = Spectrum1D(spectral_axis=sample_spec.Wave.values*u.AA, flux=sample_spec.Flux.values*u.Jy)

    # Re-sample the observe spectrum using linear interpolation
    new_spec_grid = np.linspace(5000, 6000, sample_spec.Wave.size)
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

def residual(fit_params, obs_spec, stellar_params, model_params, model, element, atomic_number, lines_out, linelist_dir, temp_dir):

    # Define the abundance
    abund = round(fit_params, 2)

    # Generate the synthetic spectra
    synth_spec_out = generate_synthetic_spectra(model_params, stellar_params, abund, model, element, atomic_number, linelist_dir)

    # Apply instrumental and rotation broadening
    final_synth_spec = spec_broadening(stellar_params, synth_spec_out, temp_dir)

    # Open the synthetic spectrum
    synth_spec = pd.read_csv(temp_dir+final_synth_spec, delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')

    # Resample the spectrum
    resamp_spec = spec_resampler(obs_spec, synth_spec)

    # Mask the absorption lines to use for the fit
    masked_obs_spec, masked_synth_spec = line_masker(resamp_spec, synth_spec, lines_out)

    # Define the chi2
    chi2 = np.sum(((masked_obs_spec['Flux'].values - masked_synth_spec['Flux'].values)**2.)/masked_synth_spec['Flux'].values)

    return chi2


# Define the different directories
temp_dir = '/blue/rezzeddine/share/fmendez/temp_dir/'
parameters_results_dir = '/blue/rezzeddine/share/fmendez/results/stellar_parameters/'
models_out_dir = '/blue/rezzeddine/share/fmendez/results/marcs_models/'
synthetic_spectra_out_dir = '/blue/rezzeddine/share/fmendez/results/synthetic_spectra/'
linelist_dir = '/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/linelists/'
observed_spectrum_dir = '/blue/rezzeddine/share/fmendez/normalized_stellar_spectra/'
save_dir = '/blue/rezzeddine/share/fmendez/results/abundances/'

# Open the data we are going to use
atomic_numbers = pd.read_csv('./atomic_numbers.csv')
asplund_solar_abund = pd.read_csv('./asplund_solar_abund.csv', sep='\s+', names=['Atomic_Number','Element','Abundance','Abundance_err'])

# Define the star we want to measure the abundances for
star = 'HD128621'
element = 'Mg'
atomic_num = atomic_numbers[atomic_numbers['Symbol'] == element]['AtomicNumber'].values[0]

# Open the observed spectrum
obs_spec = pd.read_csv(observed_spectrum_dir+star+'_5000-6000.dat', delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')

# Extract the measured stellar parameters
samples = np.genfromtxt(parameters_results_dir+star+'_corner_values.txt')
parameters = [np.percentile(samples[:,i], [16, 50,84]) for i in range(samples.shape[1])]

# Open the abundance line list and mask it based on element
linelist_out = pd.read_csv('./abundance_linelist_v2.csv').fillna(0.5)

if element+'I' in linelist_out['Element'] and element+'II' in linelist_out['Element']:
    lines_out = linelist_out[((linelist_out['Element'] == element+'I') & (linelist_out['Element'] == element+'II'))].reset_index(drop=True)

else:
    lines_out = linelist_out[linelist_out['Element'] == element+'I'].reset_index(drop=True)

# Define the range for the abundance
abundances = np.linspace(asplund_solar_abund[asplund_solar_abund['Element'] == element]['Abundance'].values[0] - 1.,
                         asplund_solar_abund[asplund_solar_abund['Element'] == element]['Abundance'].values[0] + 1.,
                         int(((asplund_solar_abund[asplund_solar_abund['Element'] == element]['Abundance'].values[0] + 1.) - (asplund_solar_abund[asplund_solar_abund['Element'] == element]['Abundance'].values[0] - 1.))/0.1) + 1)

# Generate the MARCS model
model_params, model_out = generate_marcs_model(star, parameters)

chis = []

for fit_params in tqdm(abundances):

    # Perform the LM Minimization
    out = residual(fit_params, obs_spec, parameters, model_params, model_out, element, atomic_num, lines_out, linelist_dir, temp_dir)

    chis.append(out)

# Create a dataframe to save the data
results_out = pd.DataFrame()
results_out['Abundance'] = abundances
results_out['Chi'] = chis

results_out.to_csv(save_dir+star+'_'+element+'_abundances_v2.csv', index=False)

# Remove temporary files
for f in os.listdir(temp_dir): os.remove(os.path.join(temp_dir, f))
print('DONE!')
