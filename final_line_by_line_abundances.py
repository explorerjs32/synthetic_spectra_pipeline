import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse
from subprocess import Popen, PIPE
import os


# Define the parameters to parse into the code
parser = argparse.ArgumentParser(description='Parameters to parse to derive chemical abundances')

parser.add_argument('-s', '--s', '--star_name', required=True, type=str, help='The star to derive the chemical abundances for')
parser.add_argument('-pd','--pd','--parameters_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/stellar_parameters/', help='The directory where the stellar parameters are saved')
parser.add_argument('-ad','--ad','--abundances_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/abundances/temp_abund/', help='Directory where the chemical abundances are saved')
parser.add_argument('-rsd','--rsd','--results_save_dir', type=str, default='/blue/rezzeddine/share/fmendez/results/parameters_and_abundances/', help='Directrory where the final parameters and abundances are saved')

args = parser.parse_args()

# Define the elements to get the abundances for
elements = ['Mg','Si','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Fe']

# Open the data we are going to use
asplund_solar_abund = pd.read_csv('./asplund_solar_abund.csv', sep='\s+', names=['Atomic_Number','Element','Abundance','Abundance_err'])

# Get the stellar parameters of the input star
samples = np.genfromtxt('%s%s_corner_values_v3.txt' %(args.pd, args.s))
parameters = [np.percentile(samples[:,i], [16, 50,84]) for i in range(samples.shape[1] - 2)]
errors = [np.diff(np.percentile(samples[:,i], [16, 50,84])) for i in range(samples.shape[1] -2)]

# Define the directory where the iron abundances are being stored and list all the files in it
fe_abund_dir = '/blue/rezzeddine/share/fmendez/results/abundances/fe_abundances/'
measured_fe_abunds_list = sorted(os.listdir(fe_abund_dir))

# Check if the Fe abundances have been measured already
if '%s_Fe_abundance.csv' %args.s in measured_fe_abunds_list:

    print('Fe abundance has been measured! Re-write Fe/H')

    # Update the metallicity parameter
    asplund_abund = asplund_solar_abund[asplund_solar_abund['Element'].isin(['Fe'])][['Abundance','Abundance_err']]   
    measured_abund = pd.read_csv('%s%s_Fe_abundance.csv' %(fe_abund_dir, args.s), sep='\s+')

    new_met = measured_abund['Abundance'].values - asplund_abund['Abundance'].values

    parameters[2][1] = new_met

else:
    pass

# Save the parameters and associated uncertainties to a list
values = [parameters[0][1], parameters[1][1], parameters[2][1], parameters[3][1], parameters[4][1]]
values_err = [np.mean(errors[0]), np.mean(errors[1]), np.mean(errors[2]), np.mean(errors[3]), np.mean(errors[4])]

# Interpolate over the elements list to get the abundances
for element in tqdm(elements):

    # Open the measured line abundances
    line_abundances = pd.read_csv('%s%s_%s_line_abundances_v3.csv' %(args.ad, args.s, element), sep='\s+')

    if element == 'Mg':
        line_abundances = line_abundances.drop([0,2,4])

    # Apply 3-sigma clipping to the abundances in order to get rid of outliers
    abundance_mean = line_abundances['Abundance'].mean()
    abundance_sigma = line_abundances['Abundance'].std()

    sigma_clip = ((line_abundances['Abundance'] >= abundance_mean - 3.*abundance_sigma) & (line_abundances['Abundance'] <= abundance_mean + 3.*abundance_sigma))
    masked_line_abundances = line_abundances[sigma_clip].reset_index(drop=True)

    # Get the final abundance and uncertainty
    final_abundance = masked_line_abundances['Abundance'].mean()
    final_abundance_err = masked_line_abundances['Abundance'].std()

    # Append the abundances and uncertainties to the values lists
    values.append(final_abundance)
    values_err.append(final_abundance_err)

# Save the values and uncertainties to a Dataframe
parameters_and_abundances_df = pd.DataFrame()
parameters_and_abundances_df['Parameter'] = ['Teff','logg','Fe/H','vt','vsini'] + elements
parameters_and_abundances_df['Values'] = values
parameters_and_abundances_df['Uncertainties'] = values_err

parameters_and_abundances_df.to_csv(args.rsd+args.s+'_parameters_and_abundances_v3.csv', index=False, sep=' ')

