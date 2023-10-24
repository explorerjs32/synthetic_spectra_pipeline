import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import argparse
import time
from tqdm import tqdm
from astropy import units as u
from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler
import fnmatch
from subprocess import Popen, PIPE

import warnings
warnings.filterwarnings("ignore")


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

def abundances_out(obs_spec, synth_spec_dir, synth_spec_lib, element, line_list, abundance_range):

    # Re-sample the observed spectrum to the same wavelength of the synthetic spectra
    template_spec = pd.read_csv('%s%s' %(synth_spec_dir, synth_spec_lib[0]), sep='\s+', names=['Wave','Flux','Intensity'])
    resamp_spec = spec_resampler(obs_spec, template_spec)

    # Define the lists to store the data
    final_chis = []
    final_abundances = []

    # Interpolate over the wavelenghts of each line to mask the line region
    for i, wave in enumerate(tqdm(line_list['Wave'].values.astype(float))):

        # Define the lower and upper regions of the atomic line
        lower = wave - line_list['Region'][i]
        upper = wave + line_list['Region'][i]

        # Store the chi^2 measurements
        chis = []

        # Interpolate over the library of synthetic spectra
        for synth_spec_file in synth_spec_lib:

            # Open the synthetic spectra
            synth_spec = pd.read_csv('%s%s' %(synth_spec_dir, synth_spec_file), sep='\s+', names=['Wave','Flux','Intensity'])

            # Define and apply the wavelength masks for both the observed and synthetic spectra
            obs_line_mask = ((resamp_spec['Wave'].values >= lower) & (resamp_spec['Wave'] <= upper))
            synth_line_mask = ((synth_spec['Wave'].values >= lower) & (synth_spec['Wave'] <= upper))

            masked_obs_spec, masked_synth_spec = resamp_spec[obs_line_mask], synth_spec[synth_line_mask]

            # Clipped the spectra with the highest amount of pixels
            if masked_obs_spec['Wave'].size > masked_synth_spec['Wave'].size:
                diff = abs(masked_obs_spec['Wave'].size - masked_synth_spec['Wave'].size)
                masked_obs_spec = masked_obs_spec[:-diff]

            elif masked_obs_spec['Wave'].size < masked_synth_spec['Wave'].size:
                diff = abs(masked_obs_spec['Wave'].size - masked_synth_spec['Wave'].size)
                masked_synth_spec = masked_synth_spec[:-diff]

            # Calculate the chi^2
            chi2 = np.sum(((masked_obs_spec['Flux'].values - masked_synth_spec['Flux'].values)**2.)/masked_synth_spec['Flux'].values)
            chis.append(chi2)

        # Find the minimum chi2 value and its index
        min_idx = np.argmin(np.array(chis))

        min_chi = chis[min_idx]
        final_abund = abundance_range[min_idx]

        # Save the abundance and chi value to a list
        final_chis.append(min_chi)
        final_abundances.append(final_abund)

    # Save the rusults to a dataframe
    abundances_out_df = pd.DataFrame()
    abundances_out_df['Wave'] = line_list['Wave'].values
    abundances_out_df['Element'] = line_list['Element'].values
    abundances_out_df['Abundance'] = final_abundances
    abundances_out_df['Chi2'] = final_chis

    return abundances_out_df


# Define the arguments to be parsed
parser = argparse.ArgumentParser(description='Parameters to parse to derive chemical abundances')

parser.add_argument('-s', '--s', '--star_name', required=True, type=str, help='The star to derive the chemical abundances for')
parser.add_argument('-e', '--e', '--element', required=True, type=str, help='The element to measure the chemical abundance for')
parser.add_argument('-nsd', '--nsd', '--normalized_spectra_dir', default='/blue/rezzeddine/share/fmendez/normalized_stellar_spectra/', type=str, help='Directory where the normalized stellar spectra is stored')
parser.add_argument('-spd', '--spd', '--stellar_parameters_dir', default='/blue/rezzeddine/share/fmendez/results/stellar_parameters/', type=str, help='Directory where the stellar parameters are located')
parser.add_argument('-ssd', '--ssd', '--synthetic_spectra_dir', default='/blue/rezzeddine/share/fmendez/temp_dir/', type=str, help='Directory where the synthetic spectra are stored')
parser.add_argument('-od', '--od', '--output_dir', default='/blue/rezzeddine/share/fmendez/results/abundances/fe_abundances/', type=str, help='Output directory to save the measured chemical abundances')
parser.add_argument('-slm', '--slm', '--spectral_line_mask', default='/blue/rezzeddine/share/fmendez/atomic_line_mask/chemical_abundances_lines_v2.csv', type=str, help='List of atomic lines used to chemical abundances')
parser.add_argument('-tod', '--tod', '--temporary_output_dir', default='/blue/rezzeddine/share/fmendez/temp_dir/', type=str, help='Directory where the temporary output synthetic spectra are stored')

args = parser.parse_args()

# Define additional variables
wl_range = [4200., 6600.]

# Open the measured stellar parameters
samples = np.genfromtxt('%s%s_corner_values_v3.txt' %(args.spd, args.s))
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
    uncertainties_out.append(uncertainty)

# Add the uncertaintis to the parameters keeping one of them constant
parameters_and_uncertainties = [np.array(parameters_out)]

# Convert the modified parameters to a DataFrame
star_parameters = pd.DataFrame(np.array(parameters_and_uncertainties), \
                            columns=['Teff','logg','Fe/H','turb_vel','vsini'])

# Open the observed spectrum and mask it based on the wavelength range
obs_spec = pd.read_csv('%s%s_spectrum.dat' %(args.nsd, args.s), delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python').fillna(1.0)
spec_mask = ((obs_spec['Wave'] >= wl_range[0]) & (obs_spec['Wave'] <= wl_range[1]))

obs_spec = obs_spec[spec_mask].reset_index(drop=True)

# Open the abundance line list and mask it based on element and wavelength range
literature_star_linelist = pd.read_csv(args.slm, sep='\s+', names=['Wave','Potential','loggf','Element','Region']).fillna(0.5)
linelist_mask = ((literature_star_linelist['Wave'] >= wl_range[0]) & (literature_star_linelist['Wave'] <= wl_range[1]))

linelist_out = literature_star_linelist[linelist_mask].reset_index(drop=True)

if args.e+'II' in linelist_out['Element'].values:
    print('Line list has ion 2')

    # Define the element mask
    element_mask = (linelist_out['Element'].isin([args.e+'I', args.e+'II']))
    lines_out = linelist_out[element_mask].reset_index(drop=True)

else:
    print('Line list only has ion 1')
    lines_out = linelist_out[linelist_out['Element'] == args.e+'I'].reset_index(drop=True)

# Get the library of synthetic spectra for the respective star and element
synth_spec_library = sorted([file for file in os.listdir(args.ssd) if fnmatch.fnmatch(file, '%s_*%s*.spec' %(args.s, args.e))])

# Open the Asplund 2020 abundance and define the abundance values for the given element
asplund_solar_abund = pd.read_csv('./asplund_solar_abund.csv', sep='\s+', names=['Atomic_Number','Element','Abundance','Abundance_err'])
element_abundance = asplund_solar_abund[asplund_solar_abund['Element'] == args.e]['Abundance'].values[0]

# Create the bundance range
abundance_range = np.linspace(element_abundance - 1., element_abundance + 1., 21)

# Calculate the abundances for each atomic line
line_abundances = abundances_out(obs_spec, args.ssd, synth_spec_library, args.e, lines_out, abundance_range)

# Apply 3-sigma clipping to the abundances in order to get rid of outliers
abundance_mean = line_abundances['Abundance'].mean()
abundance_sigma = line_abundances['Abundance'].std()

sigma_clip = ((line_abundances['Abundance'] >= abundance_mean - 3.*abundance_sigma) & (line_abundances['Abundance'] <= abundance_mean + 3.*abundance_sigma))
masked_line_abundances = line_abundances[sigma_clip].reset_index(drop=True)

# Get the final abundance and uncertainty
final_abundance = masked_line_abundances['Abundance'].mean()
final_abundance_err = masked_line_abundances['Abundance'].std()

# Save the Fe abundances to the respective directory
fe_abundance_out = pd.DataFrame()
fe_abundance_out['Abundance'] = [final_abundance]
fe_abundance_out['Abundance_err'] = [final_abundance_err]

fe_abundance_out.to_csv('%s%s_%s_abundance.csv' %(args.od, args.s, args.e), index=False, sep=' ')

# Remove the extra files
remove_synthspec = Popen(['rm %s*' %(args.tod)], shell=True).wait()

