import numpy as np
import os

'''
These are all the conditions that must be met by the grid depending on the parameters window:
- If the temperature range is 2500K - 4000K, then the steps are 100K.
- If the temperature range is 4000K - 8000K the steps are 250K
- IF log)g) range is -1.0 to +5.5 with steps of 0.5
- If log(g) range is -1.0 to +3.5, then the geometry is spherical (s)
- If log(g) range is +3.5 to +5.5, then the geometry is plane parallel (p)
- If the metallicity range is -1.00 to +1.00, the steps are 0.25
- If the metallicity range is -3.00 to -1.00, the steps are 0.50
- If the metallicity range is -5.00 to -3.00, the steps are 1.00

- This table represents the alpha-abundance values depending on the metallicity values
  -------------------
  [Fe/H]        [a/H]
  -------------------
  +1.00         0.00
  +0.75         0.00
  +0.50         0.00
  +0.25         0.00
   0.00         0.00
  -0.25        +0.10
  -0.50        +0.20
  -0.75        +0.30
  -1.00        +0.40
  -1.50        +0.40
  -2.00        +0.40
  -2.50        +0.40
  -3.00        +0.40
  -4.00        +0.40
  -5.00        +0.40

After the parameters window has been created, now we create the list of parameters
we want to interpolate for.

This code generates all the possible stellar parameters based on the previously
generated in par_win.py and a grid for each of the parameters. The grid we are
using for each parameter is the following:
- Temp: 250K
- log(g): 0.2
- Z: 0.2

The patern of creating the model parameters goes as follows:
- Constant Temp - Constant log(g) - Vary Z
- Constant Temp - Vary log(g) - Constant Z
- Vary Temp - Constant log(g) - Constant Z
'''
# Create the alpha abundance dictionary
alphas = {'+1.00' : '+0.00',
          '+0.75' : '+0.00',
          '+0.50' : '+0.00',
          '+0.25' : '+0.00',
          '+0.00' : '+0.00',
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

# Select the parameters we want to start with
Tmin, Tmax = 5000, 6000
gmin, gmax = 4.0, 5.0
zmin, zmax = -1.00, 0.00
turb_vel = 1
mass = '0.0'

# Convert the injected parameters to the MARCS format so we can find the models
T_low, T_hi = str(Tmin), str(Tmax)
g_low, g_hi = '+'+str(gmin) if gmin >= 0. else str(gmin), '+'+str(gmax) if gmax > 0. else str(gmax)
z_low = '+'+'{:.2f}'.format(zmin) if zmin >= 0. else '{:.2f}'.format(zmin)
z_hi = '+'+'{:.2f}'.format(zmax) if zmax >= 0. else '{:.2f}'.format(zmax)
a_low, a_hi = alphas[z_low], alphas[z_hi]
vt = '0'+str(turb_vel)
mass = '0.0'
geom = 's' if gmin >= 1.0 and gmax <= 3.5 else 'p'
comp = 'st'

# Create the model combinations based on these parameters
# From the MARCS script, it looks as the following:
model1 = geom+T_low+'_g'+g_low+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_low+'_a'+a_low+'_c+0.00_n+0.00_o'+a_low+'_r+0.00_s+0.00.mod'
model2 = geom+T_low+'_g'+g_low+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_hi+'_a'+a_hi+'_c+0.00_n+0.00_o'+a_hi+'_r+0.00_s+0.00.mod'
model3 = geom+T_low+'_g'+g_hi+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_low+'_a'+a_low+'_c+0.00_n+0.00_o'+a_low+'_r+0.00_s+0.00.mod'
model4 = geom+T_low+'_g'+g_hi+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_hi+'_a'+a_hi+'_c+0.00_n+0.00_o'+a_hi+'_r+0.00_s+0.00.mod'
model5 = geom+T_hi+'_g'+g_low+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_low+'_a'+a_low+'_c+0.00_n+0.00_o'+a_low+'_r+0.00_s+0.00.mod'
model6 = geom+T_hi+'_g'+g_low+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_hi+'_a'+a_hi+'_c+0.00_n+0.00_o'+a_hi+'_r+0.00_s+0.00.mod'
model7 = geom+T_hi+'_g'+g_hi+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_low+'_a'+a_low+'_c+0.00_n+0.00_o'+a_low+'_r+0.00_s+0.00.mod'
model8 = geom+T_low+'_g'+g_hi+'_m'+mass+'_t'+vt+'_'+comp+'_z'+z_hi+'_a'+a_hi+'_c+0.00_n+0.00_o'+a_hi+'_r+0.00_s+0.00.mod'

mod_list = [model1, model2, model3, model4, model5, model6, model7, model8]

# Open the models file list
models = open('model_list.txt', 'r').read()

# Check if the models exist for these parameters
if model1 in models and model2 in models and model3 in models and model4 in models \
and model5 in models and model6 in models and model7 in models and model8 in models:
    print('All the models for these parameters exist, move on to create the parameters grid')

    # If all the models exist, then create the parameters grid
    # First, define the steps for each of the parameters
    T_step = 250.
    g_step = 0.2
    z_step = 0.25

    # Now, create all the parameters that can come from the window and the steps
    temps = np.linspace(Tmin, Tmax, int((Tmax - Tmin)/T_step) + 1).astype('int')
    gs = np.linspace(gmin, gmax, int((gmax - gmin)/g_step) + 1)
    zs = np.linspace(zmin, zmax, int((zmax - zmin)/z_step) + 1)

    # Finally, create all tyhe possible combinations of parameters that can be made
    # by varying one parameter at a time
    parameters1 = np.array([str(t) + ' +' + str(g) + ' +' + '{:.2f}'.format(z) for t in temps for g in gs for z in zs if g >= 0. and z >= 0.])
    parameters2 = np.array([str(t) + ' +' + str(g) + ' ' + '{:.2f}'.format(z) for t in temps for g in gs for z in zs if g >= 0. and z < 0.])
    parameters3 = np.array([str(t) + ' ' + str(g) + ' +' + '{:.2f}'.format(z) for t in temps for g in gs for z in zs if g < 0. and z >= 0.])
    parameters4 = np.array([str(t) + ' ' + str(g) + ' ' + '{:print(save_parmod).2f}'.format(z) for t in temps for g in gs for z in zs if g < 0. and z < 0.])

    # Vary the microturbulet velocitis and Zs that will be used in TS
    metvel = np.array(['{:.2f}'.format(z) + ' ' + '{:.2f}'.format(turb_vel) for z in zs])

    # Save the parameters window and model parameters
    save_dir = './temp_dir/'
    if not os.path.isdir(save_dir): os.mkdir(save_dir)

    save_parwin = np.column_stack(np.array([T_low, T_hi, g_low, g_hi,
                                      z_low, z_hi, a_low, a_hi,
                                      vt, mass, geom, comp]))

    save_parmod = np.vstack((parameters1.reshape(parameters1.size, 1), parameters2.reshape(parameters2.size, 1),
                             parameters3.reshape(parameters3.size, 1), parameters4.reshape(parameters4.size, 1)))

    save_parts = np.vstack((metvel.reshape(metvel.size, 1)))

    np.savetxt(save_dir+'parwin.txt', save_parwin, fmt='%s')
    np.savetxt(save_dir+'parmod.txt', save_parmod, fmt='%s')
    np.savetxt(save_dir+'parts.txt', save_parts, fmt='%s')

    print('DONE!')

else:
    print("These models don't exist!")
    pass
