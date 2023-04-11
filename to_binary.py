import pandas as pd
import os

# Define the files directory
files_dir = '../Turbospectrum2019/COM-v19.1/syntspec/'

# Create the files list
file_list = os.listdir(files_dir)

# Convert each file in the list to binary
for file in file_list:
    df = pd.read_csv(files_dir+file, delimiter='\s+', names=['Wave','FLux','Intensity'], engine='python')
    df.to_feather(files_dir+file+'.feather')

    print(f'{file} converted to binary')

    # Remove the extra file
    os.remove(files_dir+file)
