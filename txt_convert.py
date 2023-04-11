import numpy as np
import pandas as pd
import sys

# Define the file in and fie out names
filein = str(sys.argv[1])
fileout = str(sys.argv[2])

# Read the binary file
df = pd.read_feather(filein)

# Save the dataframe as a text file
np.savetxt(fileout, df.values, fmt='%s')
