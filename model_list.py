import os
import numpy as np

# Create a path to all the MARCS models
marcs_dir = '/home/fmendez/MARCS/marcs/'
models = []

# Create the atmospheric list
atms = [a for a in os.listdir(marcs_dir) if a != '.svn']

for atm in atms:
    # Create the geometries list
    geos = [g for g in os.listdir(marcs_dir+atm+'/') if g != '.svn'
            and os.path.isdir(marcs_dir+atm+'/'+g) == True]

    for geo in geos:
        # Create the temperature list
        temps = [t for t in os.listdir(marcs_dir+atm+'/'+geo+'/') if t.endswith("K")]

        for temp in temps:
            files = [f for f in os.listdir(marcs_dir+atm+'/'+geo+'/'+temp+'/')]
            [models.append(str(mod)) for mod in files if mod.endswith('mod')]


data_out = np.array(sorted(models))
print(data_out.size)
np.savetxt('model_list.txt', data_out, fmt='%s')
