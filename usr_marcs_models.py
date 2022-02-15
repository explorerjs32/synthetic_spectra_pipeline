import os
import numpy as np
import pandas as pd

# Create the alpha abundance dictionary
alphas = {'1.00' : '+0.00',
          '0.75' : '+0.00',
          '0.50' : '+0.00',
          '0.25' : '+0.00',
          '0.00' : '+0.00',
          '-0.25' : '+0.10',
          '-0.50' : '+0.20',
          '-0.75' : '+0.30',
          '-1.00' : '+0.40',
          '-1.50' : '+0.40',
          '-2.00' : '+0.40',
          '-2.50' : '+0.40',
          '-3.00' : '+0.40',
          '-4.00' : '+0.40',
          '-5.00' : '+0.40'}

# Select the parameters you want to use for the model creation
'''
Tmin, Tmax = int(input('Choose a Tmin value (T > 2500 K): ')), int(input('Choose a Tmax value (T < 8000 K): '))
gmin, gmax = float(input('Choose a gmin value (logg > -1.0 dex): ')), float(input('Choose a gmax value (logg < 5.5 dex): '))
zmin, zmax = float(input('Choose a zmin value (Z > -5.00 dex): ')), float(input('Choose a zmax value (Z < 1.00 dex): '))
'''
Tmin, Tmax = 5000., 6000.
gmin, gmax = 4.0, 5.0
zmin, zmax = -1.0, 0.0
amin, amax = alphas['{:.2f}'.format(zmin)], alphas['{:.2f}'.format(zmax)]
turb_vel = 1
mass = '0.0'
geom = 'p'
comp = 'st'

# First, define the steps for each of the parameters
T_step = 250.
g_step = 0.2
z_step = 0.25

# Now, create all the parameters that can come from the window and the steps
temps = np.linspace(Tmin, Tmax, int((Tmax - Tmin)/T_step) + 1).astype('int')
gs = np.linspace(gmin, gmax, int((gmax - gmin)/g_step) + 1)
zs = np.linspace(zmin, zmax, int((zmax - zmin)/z_step) + 1)

# Find all the possible combinations between parametes
param_combine = [[Tmin, Tmax, gmin, gmax, zmin, zmax, turb_vel, amin, amax, mass, geom, comp, t, g, z] \
                 for t in temps for g in gs for z in zs]

# Reshape the array
model_params = np.reshape(np.asarray(param_combine), (len(param_combine), len(param_combine[0])))

# Create the data frame
columns = ['Tmin', 'Tmax', 'gmin', 'gmax', 'zmin', 'zmax', 'vt', 'amin', 'amax', 'mass', 'geom', 'comp', 'Temp', 'logg', 'z']
model_df = pd.DataFrame(data=model_params, columns=columns)

# Re-format each of the colums of the columns so it has the MARCS model format
def marcs_reformat(df):
    # Re-format the temperatures
    df['Tmin'] = np.array([int(float(t)) for t in df['Tmin']])
    df['Tmax'] = np.array([int(float(t)) for t in df['Tmax']])

    # Re-format the gravities
    df['gmin'] = np.array(['+'+g for g in df['gmin']])
    df['gmax'] = np.array(['+'+g for g in df['gmax']])
    df['logg'] = np.array(['+'+g for g in df['logg']])

    # Re-format the metallicities
    for i in range(len(df['z'])):
        if float(df['zmin'][i]) < 0.: df['zmin'][i] = '{:.2f}'.format(float(df['zmin'][i]))
        if float(df['zmin'][i]) >= 0.: df['zmin'][i] = '+{:.2f}'.format(float(df['zmin'][i]))
        if float(df['zmax'][i]) < 0.: df['zmax'][i] = '{:.2f}'.format(float(df['zmax'][i]))
        if float(df['zmax'][i]) >= 0.: df['zmax'][i] = '+{:.2f}'.format(float(df['zmax'][i]))
        if float(df['z'][i]) < 0.: df['z'][i] = '{:.2f}'.format(float(df['z'][i]))
        if float(df['z'][i]) >= 0.: df['z'][i] = '+{:.2f}'.format(float(df['z'][i]))

    # Re-format the turbulent velocity
    df['vt'] = np.array(['0'+v for v in df['vt']])

    return df

model_df = marcs_reformat(model_df)

# Save the parameters to a text file
save_dir = './temp_dir/'
if not os.path.isdir(save_dir): os.mkdir(save_dir)

np.savetxt(save_dir+'model_parameters.txt', model_df.values, fmt='%s')
