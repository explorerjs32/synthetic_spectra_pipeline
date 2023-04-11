import pandas as pd
import sys
import os

# Define the file in and file out
filein = str(sys.argv[1])
fileout = str(sys.argv[2])

# Read the File in as a Dataframe and save it as a feather file
filein_df = pd.read_csv(filein, delimiter='\s+', names=['Wave', 'Flux', 'Intensity'], engine='python')
filein_df.to_feather(fileout)

# Remove the input file
os.remove(filein)
