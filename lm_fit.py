import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd
import os
import lmfit as lm
from specutils import Spectrum1D
from specutils.manipulation import LinearInterpolatedResampler
from astropy import units as u
from astropy import constants as const
from subprocess import Popen, PIPE
import time


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

def synth_spec_in_out(parameters, inst_broad, synth_spec_dir, tmp_dir):
    '''
    This function outputs the synthetic spectra files that are going to be used for the LM Fit.
    It takes in the stellar parameters, then re-formats them so they match the format of the
    files in the synthetic spectra library
    '''

    # Re-format the injected parameters
    Teff = str(int(parameters['Teff'].value))
    logg = '+'+'{:.1f}'.format(parameters['logg'].value)
    if parameters['Z'].value < 0.: Z = '{:.2f}'.format(parameters['Z'].value)
    else: Z = '+'+'{:.2f}'.format(parameters['Z'].value)
    vt = '{:.2f}'.format(parameters['vt'].value)
    vr = '{:.1f}'.format(parameters['vr'].value)
    Na =  '{:.2f}'.format(parameters['Na'].value)
    Mg =  '{:.2f}'.format(parameters['Mg'].value)
    inst = '{:.1f}'.format(inst_broad)

    # Create the file name for the input files
    file_in_5000_bin = 'T%s_g%s_z%s_vt%s_vr%s_Na%s_Mg%s_5000_5500.spec.feather' %(Teff, logg, Z, vt, vr, Na, Mg)
    file_in_5500_bin = 'T%s_g%s_z%s_vt%s_vr%s_Na%s_Mg%s_5500_6000.spec.feather' %(Teff, logg, Z, vt, vr, Na, Mg)

    file_in_5000 = tmp_dir+file_in_5000_bin.split('.feather')[0]
    file_in_5500 = tmp_dir+file_in_5500_bin.split('.feather')[0]

    # Open the binary files and save them in a temporary directory
    file_in_5000_df = pd.read_feather(synth_spec_dir+file_in_5000_bin)
    file_in_5500_df = pd.read_feather(synth_spec_dir+file_in_5500_bin)

    file_in_5000_df.to_csv(file_in_5000, sep=' ', header=False, index=False)
    file_in_5500_df.to_csv(file_in_5500, sep=' ', header=False, index=False)

    # Create the name of the output files
    file_out_5000 = tmp_dir+'T%s_g%s_z%s_vt%s_vr%s_inst%s_Na%s_Mg%s_5000_5500.spec' %(Teff, logg, Z, vt, vr, inst, Na, Mg)
    file_out_5500 = tmp_dir+'T%s_g%s_z%s_vt%s_vr%s_inst%s_Na%s_Mg%s_5500_6000.spec' %(Teff, logg, Z, vt, vr, inst, Na, Mg)

    file_out_5000_bin = tmp_dir+'T%s_g%s_z%s_vt%s_vr%s_inst%s_Na%s_Mg%s_5000_5500.spec.feather' %(Teff, logg, Z, vt, vr, inst, Na, Mg)
    file_out_5500_bin = tmp_dir+'T%s_g%s_z%s_vt%s_vr%s_inst%s_Na%s_Mg%s_5500_6000.spec.feather' %(Teff, logg, Z, vt, vr, inst, Na, Mg)

    return file_in_5000, file_in_5500, file_out_5000, file_out_5500, file_out_5000_bin, file_out_5500_bin

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
tmp_dir = './temp_dir/'
synth_spec_dir = '../smoothed_synthetic_spectra/'
star_spec_dir = '../normalized_stellar_spectra/'
mock_sun_spec = '../mock_spectra/Sun_mock_T5777_g+4.4_z0.00_vt1.00_vr-2.0_inst-2.5_Na6.24_Mg7.60_5000-6000.spec'
line_mask = './spec_mask.txt'

# Open some of the files we are going to use
mock_spec_df = pd.read_csv(mock_sun_spec, delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')
obs_spec_df = pd.read_csv(star_spec_dir+'Sun_5000-6000.dat', delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')
spec_mask_df = pd.read_csv(line_mask, delimiter='\s+', engine='python', names=['wave','element','region']).fillna(0.3)

# Resample the spectrum
resamp_spec_df = spec_resampler(obs_spec_df, mock_spec_df)
'''
# Display the fit of the synthetic spectra
fig, ax = plt.subplots()
ax.legend()

for i in range(spec_mask_df.wave.size):
    # Define the mask for each line
    lower = spec_mask_df.wave[i] - spec_mask_df.region[i]
    upper = spec_mask_df.wave[i] + spec_mask_df.region[i]

    # Define the lines based on the elements
    if spec_mask_df.element[i] == 'MgI':
        ax.add_patch(Rectangle((lower, 0), spec_mask_df.region[i]*2., 2.1, facecolor='green', alpha=0.5, label=spec_mask_df.element[i]))

    if spec_mask_df.element[i] == 'NaI':
        ax.add_patch(Rectangle((lower, 0), spec_mask_df.region[i]*2., 2.1, facecolor='Blue', alpha=0.5, label=spec_mask_df.element[i]))

    if spec_mask_df.element[i] == 'FeI':
        ax.add_patch(Rectangle((lower, 0), spec_mask_df.region[i]*2., 2.1, facecolor='red', alpha=0.5, label=spec_mask_df.element[i]))

plt.plot(resamp_spec_df.Wave, resamp_spec_df.Flux+1, 'k.', label='Observed Spectrum')
plt.plot(mock_spec_df.Wave, mock_spec_df.Flux+1, 'r-',  label='Synthetic Spectrum')
plt.plot(resamp_spec_df.Wave, (resamp_spec_df.Flux - mock_spec_df.Flux), 'b.', label='Residuals')
plt.axhline(0., c='green', linestyle='--')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.tight_layout()
plt.legend(by_label.values(), by_label.keys())
plt.show()
'''
# Define the initial parameters for the fit interpolation
np.random.seed(0)
fit_params = lm.Parameters()
fit_params.add('Teff', value=5800, vary=True, min=4500, max=6500, brute_step=100.)
fit_params.add('logg', value=4.4, vary=True, min=4.0, max=5.0, brute_step=0.2)
fit_params.add('Z', value=0.00, vary=True, min = -1.00, max=1.00, brute_step=0.2)
fit_params.add('vt', value=1.00, vary=True, min=1.00, max=2.00, brute_step=0.2)
fit_params.add('vr', value=-3.0, vary=True, min=-10.0, max=-1.0, brute_step=1.)
fit_params.add('Na', value=6.24, vary=True, min=5.74, max=6.74, brute_step=0.25)
fit_params.add('Mg', value=7.60, vary=True, min=7.10, max=8.10, brute_step=0.25)

# Implement Instrumental broadening on the synthetic spectra by inputing the instrument resolution
inst_res = 250000.
dv_inst = (const.c.to('km/s').value/inst_res)*-1.

def residual(pars, inst_broad, synth_spec_dir, tmp_dir, obs_spec_df, spec_mask_df):
    print(pars.valuesdict())

    filein_5000, filein_5500, fileout_5000, fileout_5500, fileout_5000_bin, fileout_5500_bin = synth_spec_in_out(pars, inst_broad, synth_spec_dir, tmp_dir)

    process_5000 = Popen('./spectra_smooth.sh %s %s %s %s %s' % (str('inst_broadening_5000.f90'), str(filein_5000), str(fileout_5000), str('{:.1f}'.format(dv_inst)), str('inst_broadening_5000')), shell=True)
    process_5500 = Popen('./spectra_smooth.sh %s %s %s %s %s' % (str('inst_broadening_5500.f90'), str(filein_5500), str(fileout_5500), str('{:.1f}'.format(dv_inst)), str('inst_broadening_5500')), shell=True)

    time.sleep(3)

    # Open the created spectra and stack it together
    synth_spec_5000_df = pd.read_csv(fileout_5000, delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')
    synth_spec_5500_df = pd.read_csv(fileout_5500, delimiter='\s+', names=['Wave','Flux','Intensity'], engine='python')

    synth_spec_df = pd.concat([synth_spec_5000_df, synth_spec_5500_df])

    # Mask the spectral lines used for the Fitting
    masked_obs_spec_df, masked_synth_spec_df = line_masker(obs_spec_df, synth_spec_df, spec_mask_df)

    return masked_synth_spec_df['Flux'] - masked_obs_spec_df['Flux']

# Perform the LM Minimization

out = lm.minimize(residual, fit_params, method='leastsq', args=(dv_inst, synth_spec_dir, tmp_dir, resamp_spec_df, spec_mask_df))
print(lm.fit_report(out))
