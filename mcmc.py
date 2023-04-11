import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import os
from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler
from astropy import units as u
from astropy import constants as const
import emcee
import corner
from multiprocessing import Pool
from subprocess import Popen, PIPE
import time
import random


def spec_resampler(obs_spec, sample_spec):
    '''
    The observed spectrum is not expected to be sampled in the same grid space as the synthetic spectra.
    In order to perform an acurate pixel to pixel fit we want to re-sample the observed spectrum to the
    same spectral grid than our library of synthetic spetra.
    '''
    # Convert the spectrum dataframe into the SpecUtils format
    obs_spec = Spectrum1D(spectral_axis=obs_spec_df.Wave.values*u.AA, flux=obs_spec_df.Flux.values*u.Jy)
    mock_spec = Spectrum1D(spectral_axis=sample_spec.Wave.values*u.AA, flux=sample_spec.Flux.values*u.Jy)

    # Re-sample the observe spectrum using linear interpolation
    new_spec_grid = np.linspace(5000, 6000, sample_spec.Wave.size)
    linear = LinearInterpolatedResampler()
    obs_spec_resamp = linear(obs_spec, new_spec_grid*u.AA)

    # Return a datafram of the re-sampled spectrum
    resamp_spec_df = pd.DataFrame(data=np.column_stack([obs_spec_resamp.spectral_axis.value, obs_spec_resamp.flux.value]), columns=['Wave','Flux'])

    return resamp_spec_df

def spec_out(params, synth_spec_dir):

    # Re-format the parameters
    Tmin, Tmax = str(int(params['T'][0])), str(int(params['T'][1]))
    loggmin, loggmax = '+{:.1f}'.format(params['logg'][0]), '+{:.1f}'.format(params['logg'][1])
    vtmin, vtmax =  '{:.2f}'.format(params['vt'][0]), '{:.2f}'.format(params['vt'][1])
    vrmin, vrmax = '{:.1f}'.format(params['vr'][0]), '{:.1f}'.format(params['vr'][1])
    Namin, Namax = '{:.2f}'.format(params['Na'][0]), '{:.2f}'.format(params['Na'][1])
    Mgmin, Mgmax = '{:.2f}'.format(params['Mg'][0]), '{:.2f}'.format(params['Mg'][1])
    if params['Z'][0] < 0.: Zmin = '{:.2f}'.format(params['Z'][0])
    else: Zmin = '+{:.2f}'.format(params['Z'][0])
    if params['Z'][1] < 0.: Zmax = '{:.2f}'.format(params['Z'][1])
    else: Zmax = '+{:.2f}'.format(params['Z'][1])

    # Create the file name for the input files
    file_min_spec_5000 = 'T%s_g%s_z%s_vt%s_vr-%s_Na%s_Mg%s_5000_5500.spec.feather' %(Tmin, loggmin, Zmin, vtmin, vrmin, Namin, Mgmin)
    file_min_spec_5500 = 'T%s_g%s_z%s_vt%s_vr-%s_Na%s_Mg%s_5500_6000.spec.feather' %(Tmin, loggmin, Zmin, vtmin, vrmin, Namin, Mgmin)
    file_max_spec_5000 = 'T%s_g%s_z%s_vt%s_vr-%s_Na%s_Mg%s_5000_5500.spec.feather' %(Tmax, loggmax, Zmax, vtmax, vrmax, Namax, Mgmax)
    file_max_spec_5500 = 'T%s_g%s_z%s_vt%s_vr-%s_Na%s_Mg%s_5500_6000.spec.feather' %(Tmax, loggmax, Zmax, vtmax, vrmax, Namax, Mgmax)

    # Open the binary files and extract the spectra
    min_spec_5000_df = pd.read_feather(synth_spec_dir+file_min_spec_5000, columns=['Wave','Flux'])
    min_spec_5500_df = pd.read_feather(synth_spec_dir+file_min_spec_5500, columns=['Wave','Flux'])
    max_spec_5000_df = pd.read_feather(synth_spec_dir+file_max_spec_5000, columns=['Wave','Flux'])
    max_spec_5500_df = pd.read_feather(synth_spec_dir+file_max_spec_5500, columns=['Wave','Flux'])

    return min_spec_5000_df, min_spec_5500_df, max_spec_5000_df, max_spec_5500_df

def spec_interpolate(min_spec_df, max_spec_df):
    # Set the flux from both spectra into a single dataframe along with the wavelength values
    synth_spectra_df = pd.concat([min_spec_df['Flux'], max_spec_df['Flux']], axis=1)
    synth_spectra_df.columns = ['Flux_min','Flux_max']

    # Create a column full of nan values to do the interpolation between the two spectra
    nans = np.empty(synth_spectra_df['Flux_min'].size)
    nans[:] = np.NaN

    synth_spectra_df.insert(1, 'New_Flux', nans)
    synth_spectra_df = synth_spectra_df.reset_index(drop=True)

    # Interpolate between the two spectra
    synth_spectra_df = synth_spectra_df.interpolate(method='linear', axis=1)

    # Create a dataframe for the interpolated spectra
    spec_interpolated_df = pd.DataFrame()
    spec_interpolated_df['Wave'] = min_spec_df['Wave'].values
    spec_interpolated_df['Flux'] = synth_spectra_df['New_Flux'].values

    # Plot the interpolated spectrum along with the spectra used for interpolation
    '''
    plt.plot(min_spec_df['Wave'], synth_spectra_df['Flux_min'], 'y-')
    plt.plot(min_spec_df['Wave'], synth_spectra_df['Flux_max'], 'r-')
    plt.plot(min_spec_df['Wave'], synth_spectra_df['New_Flux'], 'b-')
    plt.show()
    '''

    return spec_interpolated_df

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

    for i in range(spec_mask_df.wave.size):
        # Define the mask for each line
        lower = spec_mask_df.wave[i] - spec_mask_df.region[i]
        upper = spec_mask_df.wave[i] + spec_mask_df.region[i]

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


# Define the paths to the diffeent directories and flles
synth_spec_dir = '../synthetic_spectra/'
star_spec_dir = '../normalized_stellar_spectra/'
mock_sun_spec = '../mock_spectra/Sun_mock_T5777_g+4.4_z0.00_vt1.00_vr-2.0_inst-2.5_Na6.24_Mg7.60_5000-6000.spec'
line_mask = './spec_mask_v2.txt'
results_dir = '../results/'

# Open some of the files we are going to use
star = 'Sun'
mock_spec_df = pd.read_csv(mock_sun_spec, delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')
obs_spec_df = pd.read_csv(star_spec_dir+star+'_5000-6000.dat', delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')
spec_mask_df = pd.read_csv(line_mask, delimiter='\s+', engine='python', names=['wave','element','potential','loggf','region']).fillna(0.5)

# Mask the lines by loggf
sm_mask = ((spec_mask_df['element'] == 'FeI') & (spec_mask_df['loggf'] >= -1.5))

idx = spec_mask_df[sm_mask].index.values
spec_mask_df = spec_mask_df.drop(index=idx, axis=0).reset_index(drop=True)

# Resample the spectrum and create the flux error array
resamp_spec_df = spec_resampler(obs_spec_df, mock_spec_df)
resamp_spec_df['Flux_err'] = np.zeros([resamp_spec_df['Flux'].size]) + 0.05


# Define the parameter grid space
tgrid = np.linspace(4500, 6500, 21)
ggrid = np.linspace(4.0, 5.0, 6)
zgrid = np.linspace(-1., 1., 11)
vtgrid = np.linspace(1., 2., 6)
vrgrid = np.linspace(1., 10., 10)
nagrid = np.linspace(5.74, 6.74, 5)
mggrid = np.linspace(7.10, 8.10, 5)

matrix = {'T':tgrid,
          'logg':ggrid,
          'Z':zgrid,
          'vt':vtgrid,
          'vr':vrgrid,
          'Na':nagrid,
          'Mg':mggrid}

def log_likelihood(theta, obs_spec, matrix, spec_mask):
    # Define the spectral parameters
    T, logg, z, vt, vr, na, mg = theta

    params = {'T':T,
              'logg':logg,
              'Z':z,
              'vt':vt,
              'vr':vr,
              'Na':na,
              'Mg':mg}

    # Based on these parameters find the closest upper and lower value at each parameter grid
    params_pair = {}

    for key in matrix.keys():

        def closest_param(param_list, param):
            min_index = np.argsort(abs(param_list - param))[0]
            max_index = np.argsort(abs(param_list - param))[1]

            return sorted([param_list[min_index], param_list[max_index]])

        new_params = closest_param(matrix[key], params[key])
        params_pair[key] = new_params

    # Extract the synthetic flux and append it together
    min_spec_5000_df, min_spec_5500_df, max_spec_5000_df, max_spec_5500_df = spec_out(params_pair, synth_spec_dir)

    min_spec_df = pd.concat([min_spec_5000_df, min_spec_5500_df]).reset_index(drop=True)
    max_spec_df = pd.concat([max_spec_5000_df, max_spec_5500_df]).reset_index(drop=True)

    # Interpolate the spectrum for the given parameters
    spec_interpolated_df = spec_interpolate(min_spec_df.astype(float), max_spec_df.astype(float))

    # Mask the spectral lines we are going to use for the minimization
    masked_obs_spec_df, masked_interpolated_spec_df = line_masker(obs_spec, spec_interpolated_df, spec_mask)

    # Define the flux and the model
    flux = masked_obs_spec_df['Flux']
    model = masked_interpolated_spec_df['Flux']
    flux_err = np.zeros(model.size) + 0.05

    sigma2 = model**2. + flux_err**2.

    return -0.5 * np.sum((flux - model) ** 2 / sigma2 )

def log_prior(theta):
    T, logg, Z, vt, vr, na, mg = theta

    # Define the gridspace for each of the walkers at each parameter
    if (4500. <= T <= 6500. and 4.0 <= logg <= 5.0 and -1.0 <= Z <= 1.0 and 1.0 <= vt <= 2.0 \
    and 1.0 <= vr <= 10. and 5.74 <= na <= 6.74 and 7.10 <= mg <= 8.10):
        return 0.0
    return -np.inf

def log_probability(theta, obs_spec, matrix, spec_mask):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, obs_spec, matrix, spec_mask)

# Select the initial parameters where the sampling is going to start
T0, logg0, Z0, vt0, vr0, na0, mg0 = 5500., 4.5, 0.0, 1.5, 5.0, 6.24, 7.60

# Create the step range for each parameter
tstep = np.linspace(-100., 100., 2)
gstep = np.linspace(-.2, .2, 2)
zstep = np.linspace(-.25, .25, 2)
vtstep = np.linspace(-.2, .2, 2)
vrstep = np.linspace(-1., 1., 2)
nastep = np.linspace(-.25, .25, 2)
mgstep = np.linspace(-.25, .25, 2)

# Set up the starting parameter fo each walker
ndim = 7
nwalkers = 250
pos = np.array([np.array([T0, logg0, Z0, vt0, vr0, na0, mg0]) + \
               np.array([random.choice(tstep),
                         random.choice(gstep),
                         random.choice(zstep),
                         random.choice(vtstep),
                         random.choice(vrstep),
                         random.choice(nastep),
                         random.choice(mgstep)]) * np.random.randn(ndim) for i in range(nwalkers)])

# Run the MCMC code by applying paralelization
with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool, args=(resamp_spec_df, matrix, spec_mask_df))
    sampler.run_mcmc(pos, 1000, progress=True)

# Get the results from the MCMC
samples_corner = sampler.flatchain
results = samples_corner[np.argmax(sampler.flatlnprobability)]

# Get the walkers information
samples_walkers = sampler.get_chain()
samples_walkers_rs = samples_walkers.reshape(samples_walkers.shape[0], -1)

np.savetxt('../results/stellar_parameters/'+star+'_corner_values.txt', samples_corner)
np.savetxt('../results/stellar_parameters/'+star+'_walkers.txt', samples_walkers_rs)
