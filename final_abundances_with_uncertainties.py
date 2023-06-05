import numpy as np
import pandas as pd
from tqdm import tqdm
import shutil
import sys
import argparse
from subprocess import Popen, PIPE


# Define the different arguments to parse into the code
parser = argparse.ArgumentParser(description='Parameters to parse to save the stellar parameters, abundances, and save the synthetic specta')

parser.add_argument('-s','--s','--star_name', required=True, type=str, help='The star to save the results for')
parser.add_argument('-pd','--pd','--parameters_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/stellar_parameters/', help='The directory where the stellar parameters are saved')
parser.add_argument('-ad','--ad','--abundances_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/abundances/temp_abund/', help='Directory where the chemical abundances are saved')
parser.add_argument('-td','--td','--temp_dir', type=str, default='/blue/rezzeddine/share/fmendez/temp_dir/', help='Directory where the temporary files are saved')
parser.add_argument('-ssd','--ssd','--synthetic_spectra_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/synthetic_spectra/abundances/', help='Directory where to save the synthetic spectra')
parser.add_argument('-rsd','--rsd','--results_save_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/parameters_and_abundances/', help='Directrory where the final parameters and abundances are saved')

args = parser.parse_args()

# Define the elements to get the abundances for
elements = ['Mg','Si','Ca','Sc','Ti','V','Cr','Mn','Co','Ni']

# Get the stellar parameters of the input star
samples = np.genfromtxt('%s%s_corner_values.txt' %(args.pd, args.s))
parameters = [np.percentile(samples[:,i], [16, 50,84]) for i in range(samples.shape[1] - 2)]
errors = [np.diff(np.percentile(samples[:,i], [16, 50,84])) for i in range(samples.shape[1] -2)]

# Save the parameters and associated uncertainties to a list
values = [parameters[0][1], parameters[1][1], parameters[2][1], parameters[3][1], parameters[4][1]]
values_err = [np.mean(errors[0]), np.mean(errors[1]), np.mean(errors[2]), np.mean(errors[3]), np.mean(errors[4])]

# Save the abundance departure from each parameter for each element
abunds_delta_t, abunds_delta_g, abunds_delta_z, abunds_delta_vt, abunds_delta_vsini = [], [], [], [], []

for element in tqdm(elements):

    # Create the abundance range based on the chemical element
    asplund_abundances = pd.read_csv('./asplund_solar_abund.csv', sep=' ', \
                                     names=['Atomic_Number','Element','Abundance','Abundance_Error'])
    solar_abundance = asplund_abundances[asplund_abundances['Element'] == element]['Abundance'].values[0]

    abundance_range = np.linspace(solar_abundance - 1, solar_abundance + 1, 21)

    # Interpolate over the abundances and extract the best derived abundances
    parameters_abunds = []
    delta_t_abunds = []
    delta_g_abunds = []
    delta_z_abunds = []
    delta_vt_abunds = []
    delta_vsini_abunds = []

    for abundance in abundance_range:

        # Open the derived abundances for each set of parameters
        derived_abundances_df = pd.read_csv('%s%s_%s_%s.csv' %(args.ad, args.s, element, str(round(abundance, 2))), sep=' ')

        # Store each chi2 value into the respective list
        parameters_abunds.append(derived_abundances_df.iloc[0]['Chi2'])
        delta_t_abunds.append(derived_abundances_df.iloc[1]['Chi2'])
        delta_g_abunds.append(derived_abundances_df.iloc[2]['Chi2'])
        delta_z_abunds.append(derived_abundances_df.iloc[3]['Chi2'])
        delta_vt_abunds.append(derived_abundances_df.iloc[4]['Chi2'])
        delta_vsini_abunds.append(derived_abundances_df.iloc[5]['Chi2'])

    # Get the index of the smallest chi2 value per parameters set
    parameters_idx = np.array(parameters_abunds).argmin()
    delta_t_idx = np.array(delta_t_abunds).argmin()
    delta_g_idx = np.array(delta_g_abunds).argmin()
    delta_z_idx = np.array(delta_z_abunds).argmin()
    delta_vt_idx = np.array(delta_vt_abunds).argmin()
    delta_vsini_idx = np.array(delta_vsini_abunds).argmin()

    # Get the measured abundances per parameters set
    abundances_out = np.array([abundance_range[parameters_idx], abundance_range[delta_t_idx], abundance_range[delta_g_idx],
                               abundance_range[delta_z_idx], abundance_range[delta_vt_idx], abundance_range[delta_vsini_idx]])

    # Extract the final abundance and uncertainty
    final_abundance, final_abundance_uncertainty = abundances_out[0], np.sqrt(np.mean((abundances_out - abundances_out[0])**2.))

    # Add the measured abundances and uncertainties to the lists
    values.append(final_abundance)
    values_err.append(final_abundance_uncertainty)

    # Add the abundance departures based on stellar parameters
    abunds_delta_t.append(abundances_out[1] - final_abundance)
    abunds_delta_g.append(abundances_out[2] - final_abundance)
    abunds_delta_z.append(abundances_out[3] - final_abundance)
    abunds_delta_vt.append(abundances_out[4] - final_abundance)
    abunds_delta_vsini.append(abundances_out[5] - final_abundance)

    # Save each of the best fit synthetic spectra
    sub_dirs = ['measured','delta_t','delta_g','delta_z','delta_vt','delta_vsini']

    for i, abund_out in enumerate(abundances_out):

        # Extract the synthetic spectra of the measured abundance
        measured_abundance_df = pd.read_csv('%s%s_%s_%s.csv' %(args.ad, args.s, element, str(round(abund_out, 2))), sep=' ')
        spec_filename = measured_abundance_df['Syntspec_File'].iloc[i]

        # Move the synthetic spectra file to the results directory
        shutil.copy(args.td+spec_filename, args.ssd+sub_dirs[i]+'/'+args.s+'_'+sub_dirs[i]+'_'+spec_filename)

# Save the values and uncertainties to a Dataframe
parameters_and_abundances_df = pd.DataFrame()
parameters_and_abundances_df['Parameter'] = ['Teff','logg','Fe/H','vt','vsini'] + elements
parameters_and_abundances_df['Values'] = values
parameters_and_abundances_df['Uncertainties'] = values_err
parameters_and_abundances_df['Delta_Teff'] = [np.nan, np.nan, np.nan, np.nan, np.nan] + abunds_delta_t
parameters_and_abundances_df['Delta_logg'] = [np.nan, np.nan, np.nan, np.nan, np.nan] + abunds_delta_g
parameters_and_abundances_df['Delta_Fe/H'] = [np.nan, np.nan, np.nan, np.nan, np.nan] + abunds_delta_z
parameters_and_abundances_df['Delta_vt'] = [np.nan, np.nan, np.nan, np.nan, np.nan] + abunds_delta_vt
parameters_and_abundances_df['Delta_vsini'] = [np.nan, np.nan, np.nan, np.nan, np.nan] + abunds_delta_vsini

parameters_and_abundances_df.to_csv(args.rsd+args.s+'_parameters_and_abundances.csv', index=False, sep=' ')

# Remove all the extra files
remove_synthspec = Popen(['rm %s*' %(args.td)], shell=True).wait()
remove_abunds = Popen(['rm %s*' %(args.ad)], shell=True).wait()