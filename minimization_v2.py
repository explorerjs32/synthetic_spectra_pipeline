import numpy as np
import matplotlib.pyplot as plt
import os
import copy
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region
from specutils.fitting import estimate_line_parameters
from astropy import units as u
from astropy.modeling import models
from scipy.stats import chisquare


def HARPS_fiber_fill(obs_wl, obs_fl, synth_wl, synth_fl):
    # Figure out the max blue and min red wave from the spectrum
    max_blue = max(obs_wl[obs_wl <= 5320.])
    min_red = min(obs_wl[obs_wl >= 5320.])

    # Fill in the gap with ones
    wl_gap = np.linspace(max_blue, min_red, int((min_red - max_blue)/0.01))
    fl_gap = np.ones(wl_gap.size)

    # Concatenate the gap filler to the spectrum
    obs_wl, obs_fl = np.concatenate([obs_wl, wl_gap]), np.concatenate([obs_fl, fl_gap])

    # MAsk this region on the synthetic spectram
    region = (np.where(synth_wl >= max_blue)[0][0], np.where(synth_wl <= min_red)[0][-1])
    synth_fl[region[0]:region[1]] = 1.

    # remove any extra values if the obs_spec and synth_spec have not the same size
    if obs_wl.size > synth_wl.size: obs_wl = obs_wl[:-1]; obs_fl = obs_fl[:-1]
    if obs_wl.size < synth_wl.size: synth_wl = synth_wl[:-1]; synth_fl = synth_fl[:-1]

    return obs_wl, obs_fl, synth_wl, synth_fl

def region_fit(line):
    # Create a dictionary containing the line regions
    #regions = {'Mg':[5164.25, 5187.25]}
    regions = {'all':[5000., 6000.],
               'Mg':[5160., 5190.],
               'Na':[5887.5, 5898.]}

    return regions.get(str(line))

def spectra_clipper(obs_wl, obs_fl, synth_wl, synth_fl, region, lines_dict):
    obs_wls, obs_fls, synth_wls, synth_fls = [], [], [], []

    for r in region:
        # Based on the region_fit fuction mask the spectra to isolate the selected region
        mask_intervals = region_fit(r)

        obs_mask = (obs_wl >= mask_intervals[0]) & (obs_wl <= mask_intervals[1]+.01)
        #synth_mask = (synth_wl >= mask_intervals[0]) & (synth_wl <= mask_intervals[1])

        # Mask the observed and synthetic spectra
        clip_obs_wl, clip_obs_fl = obs_wl[obs_mask], obs_fl[obs_mask]
        clip_synth_wl, clip_synth_fl = synth_wl[obs_mask], synth_fl[obs_mask]

        # Set up the spectra for specutils
        obs_spec = Spectrum1D(flux=clip_obs_fl*u.Jy, spectral_axis=clip_obs_wl*u.AA)
        synth_spec = Spectrum1D(flux=clip_synth_fl*u.Jy, spectral_axis=clip_synth_wl*u.AA)

        lines = lines_dict[r]

        for l in lines:
            # Create a copy of the observed and synthetic spectra to modify it
            copy_obs_spec = copy.deepcopy(obs_spec)
            copy_synth_spec = copy.deepcopy(synth_spec)

            # Define line mask
            sub_region = SpectralRegion((l - 2.)*u.AA, (l + 2.)*u.AA)
            clip_spec = extract_region(copy_obs_spec, sub_region)

            # Estimate the line Gauss parameters
            line_stddev = estimate_line_parameters(clip_spec, models.Gaussian1D()).stddev.value

            # Isolate the lines and wings based on 1 sigma rejection
            line_mask = (copy_obs_spec.spectral_axis.value >= l - 2*line_stddev) & (copy_obs_spec.spectral_axis.value <= l + 2*line_stddev)

            # Isolate the individual lines
            obs_line_wl, obs_line_fl = copy_obs_spec.spectral_axis.value[line_mask], copy_obs_spec.flux.value[line_mask]
            synth_line_wl, synth_line_fl = copy_synth_spec.spectral_axis.value[line_mask], copy_synth_spec.flux.value[line_mask]

            # Make the rest of the spectrum equal to one
            cont_region = np.where(line_mask == False)[0]
            copy_obs_spec.flux.value[cont_region] = 1.
            copy_synth_spec.flux.value[cont_region] = 1.

            obs_wls.append(copy_obs_spec.spectral_axis.value), obs_fls.append(copy_obs_spec.flux.value)
            synth_wls.append(copy_synth_spec.spectral_axis.value), synth_fls.append(copy_synth_spec.flux.value)

    out_obs_wl, out_obs_fl = np.concatenate(obs_wls), np.concatenate(obs_fls)
    out_synth_wl, out_synth_fl = np.concatenate(synth_wls), np.concatenate(synth_fls)

    '''
    plt.plot(obs_wl, obs_fl, 'r-')
    plt.plot(out_obs_wl, out_obs_fl, 'k-')
    plt.show()
    '''

    return out_obs_wl, out_obs_fl, out_synth_wl, out_synth_fl


# Create a dictioary with the spectral lines of a deffined region
spectral_lines = {'Mg':[5167.322, 5172.684, 5183.604], \
                  'Na':[5889.950, 5895.924]}

# What star do you want to analyze
star = 'HD10700'

# Define the data directories
obs_data_dir = '/home/fmendez/Documents/Research/Stellar_Spectra/Normalized/'
#synth_data_dir = '/media/fmendez/FA62A68862A64969/Research/Stellar_Parameters/synthetic_spectra/'
synth_data_dir = '/media/fmendez/2512-CDD2/synth_spec/'

# List the synthetic spectra
synth_data_list = sorted(os.listdir(synth_data_dir))

# Define the spectral region to mask
region = ['Mg', 'Na']

# Create a list to save all the chi 2 measurements
chis = []
clipped_spectra = []
files = []

for num, file in enumerate(synth_data_list):\
    # Open the data
    obs_data = np.genfromtxt(obs_data_dir+'Mock_5000g+4.0z--0.25vt1.20_5000-6000.spec')
    synth_data = np.genfromtxt(synth_data_dir+'/'+file)

    obs_wl, obs_fl = obs_data[:,0], obs_data[:,1]
    synth_wl, synth_fl = synth_data[:,0], synth_data[:,1]

    files.append(file)

    # remove any extra values if the obs_spec and synth_spec have not the same size
    if obs_wl.size > synth_wl.size: obs_wl = obs_wl[:-1]; obs_fl = obs_fl[:-1]
    if obs_wl.size < synth_wl.size: synth_wl = synth_wl[:-1]; synth_fl = synth_fl[:-1]

    # For HARPS spectrum we have to account for the gap that is on the 5000 - 6000 region due to
    # the fiber change, So, here we figure out where that gap is happenning and we mask it
    obs_wl, obs_fl, synth_wl, synth_fl = HARPS_fiber_fill(obs_wl, obs_fl, synth_wl, synth_fl)

    # Match the wavelength values between obseved and synthetic spectra
    min_obs_wl, min_synth_wl, max_obs_wl, max_synth_wl = min(obs_wl), min(synth_wl), max(obs_wl), max(synth_wl)

    if min_obs_wl > min_synth_wl and max_obs_wl > max_synth_wl:
        wl_mask = (synth_wl >= min_obs_wl) & (obs_wl <= max_synth_wl)
        obs_wl, obs_fl, synth_wl, synth_fl = obs_wl[wl_mask], obs_fl[wl_mask], synth_wl[wl_mask], synth_fl[wl_mask]

    # Isolate the spectral region we are interested in analyzing
    obs_wl, obs_fl, synth_wl, synth_fl = spectra_clipper(obs_wl, obs_fl, synth_wl, synth_fl, region, spectral_lines)
    '''
    plt.plot(obs_wl, obs_fl)
    plt.plot(synth_wl, synth_fl)
    plt.show()
    '''

    # Stack the clipped spectra into columns to fit later
    clipped_stacked = np.column_stack([obs_wl, obs_fl, synth_wl, synth_fl])
    clipped_spectra.append(clipped_stacked)

    # perform the Chi 2 fitting
    print(f'Fitting synthetic spectra {num+1} of {len(synth_data_list)}')
    chi, p = chisquare(obs_fl, synth_fl)

    chis.append(chi)

print(min(chis), chis.index(min(chis)))

# After performing the chi^2 fitting, calculate the fit residuals
chi_residuals = (clipped_spectra[chis.index(min(chis))][:,1] - clipped_spectra[chis.index(min(chis))][:,3])/np.sqrt(clipped_spectra[chis.index(min(chis))][:,3])

# Plot the best fit and the associated residuals

plt.figure(figsize=(12,8))

plt.subplot(211)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.title(files[chis.index(min(chis))], size=25)
plt.xlabel(r'Wavelength $\AA$',size=15)
plt.ylabel('Flux')
plt.ylim(0.,1.2)
plt.text(5891.65, 1.05, r'$\chi^2=$' + '{:.2f}'.format(min(chis)), size=15)
plt.scatter(clipped_spectra[chis.index(min(chis))][:,0], clipped_spectra[chis.index(min(chis))][:,1], \
         s=1., c='black', label='Observed Spectrum')
#plt.scatter(clipped_spectra[chis.index(min(chis))][:,2], clipped_spectra[chis.index(min(chis))][:,3], \
         #s=1., c='red', label='Synthetic Spectrum')
plt.plot(clipped_spectra[chis.index(min(chis))][:,2], clipped_spectra[chis.index(min(chis))][:,3], \
         'r-', label='Synthetic Spectrum')
plt.legend()

plt.subplot(212)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel('Residuals', size=15)
plt.xlabel(r'Wavelength $\AA$',size=15)
plt.plot(clipped_spectra[chis.index(min(chis))][:,0], chi_residuals, 'k-')

plt.tight_layout()
plt.show()
