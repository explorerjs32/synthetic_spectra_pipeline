import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

    # Plot the spectra

    plt.plot(max_spec_5000_df['Wave'], max_spec_5000_df['Flux'], 'k-', label=file_max_spec_5000)
    plt.plot(max_spec_5500_df['Wave'], max_spec_5500_df['Flux'], 'k-')
    plt.plot(min_spec_5000_df['Wave'], min_spec_5000_df['Flux'], 'b-', label=file_min_spec_5000)
    plt.plot(min_spec_5500_df['Wave'], min_spec_5500_df['Flux'], 'b-')
    plt.legend()
    plt.tight_layout()
    plt.show()

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

    plt.plot(min_spec_df['Wave'], synth_spectra_df['Flux_min'], 'y-', label='Lower Spectrum T=5500 logg=4.4 z=0.0 vt=1.0 vr=-4.0 Na=6.24 Mg=7.10')
    plt.plot(min_spec_df['Wave'], synth_spectra_df['Flux_max'], 'r-', label='Higher Spectrum T=5600 logg=4.6 z=0.2 vt=1.2 vr=-5.0 Na=6.49 Mg=7.35')
    plt.plot(min_spec_df['Wave'], synth_spectra_df['New_Flux'], 'b-', label='Interpolated Spectrum T=5527 logg=4.5 z=0.1 vt=1.1 vr=-4.5 Na=6.30 Mg=7.15')
    plt.legend()
    plt.show()


    return spec_interpolated_df


# Define the data directories
synth_spec_dir = '../smoothed_synthetic_spectra/'
results_dir = '../results/'

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

# Define the original parameters
T, logg, z, vt, vr, na, mg = 5527., 4.5, 0.1, 1.1, 4.5, 6.30, 7.15

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

min_spec_df = pd.concat([min_spec_5000_df, min_spec_5500_df])
max_spec_df = pd.concat([max_spec_5000_df, max_spec_5500_df])

# Interpolate the spectrum for the given parameters
spec_interpolated_df = spec_interpolate(min_spec_df, max_spec_df)
spec_interpolated_df.to_csv(results_dir+'interpol_T5527_g+4.5_z+0.10_vt1.10_vr-4.5_Na6.30_Mg7.15.spec')

plt.plot(spec_interpolated_df['Wave'], spec_interpolated_df['Flux'])
plt.show()
