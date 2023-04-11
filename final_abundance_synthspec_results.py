import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil
from subprocess import Popen, PIPE
import sys

# Define the different directories
abundances_dir = '/blue/rezzeddine/share/fmendez/results/abundances/temp_abund/'
temp_dir = '/blue/rezzeddine/share/fmendez/temp_dir/'
synth_spec_save_dir = '/blue/rezzeddine/share/fmendez/results/synthetic_spectra/abundances/'

# Define the star and the element
star = str(sys.argv[1])
element = str(sys.argv[2])

# Combine the individual abundances for a single star and element
combine_args = ['cat %s%s_%s* >> %s%s_%s_combined_abundances_v2.txt' %(abundances_dir, star, element, abundances_dir, star, element)]
combine_process = Popen(combine_args, shell=True).wait()

# Open the measured abundances and sort them by lowest chi2
abundances = pd.read_csv(abundances_dir+star+'_'+element+'_combined_abundances_v2.txt', sep='\s+', names=['File_Name','Element','Abundance','Chi2'])
sorted_by_chis = abundances.sort_values(by='Chi2').reset_index(drop=True)

# Get the file name of the best fit spectrum, re-name it and move it to the respective folder
spec_filename = sorted_by_chis['File_Name'].values[0]
shutil.move(temp_dir+spec_filename, synth_spec_save_dir+star+'_measured_'+spec_filename)
