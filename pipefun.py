import numpy as np
import copy
import os

def marcs_models_para_list(marcs_temperatures, marcs_gravities, marcs_metallicities, alphas, geom, comp, masses, turbvels):
    # Split each of the variables to the MARCS format
    T_low, T_hi = marcs_temperatures[0], marcs_temperatures[-1]
    g_low, g_hi = marcs_gravities[0], marcs_gravities[-1]
    z_low, z_hi = marcs_metallicities[0], marcs_metallicities[-1]
    a_low, a_hi = alphas[z_low], alphas[z_hi]
    geom = geom
    comp = 'st'
    mass = masses[0]
    turbvel = turbvels[0]

    # Create a list for the MARCS Parameters
    marcs_parameters_list = [T_low, T_hi, g_low, g_hi, z_low, z_hi, a_low, a_hi, geom, comp, mass, turbvel]

    return marcs_parameters_list

def marcs_models_components(modeln, mp):
    # mp list order [T_low, T_hi, g_low, g_hi, z_low, z_hi, a_low, a_hi, geom, comp, mass, turbvel]
    # Define all the different components depending on the model
    # Model 1
    if modeln == 'model1': mc = [mp[8], mp[0], '_g', mp[2], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[4], '_a', mp[6], '_c+0.00_n+0.00_o', mp[6], '_r+0.00_s+0.00.mod']

    # Model 2
    if modeln == 'model2': mc = [mp[8], mp[0], '_g', mp[2], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[5], '_a', mp[7], '_c+0.00_n+0.00_o', mp[7], '_r+0.00_s+0.00.mod']

    # Model 3
    if modeln == 'model3': mc = [mp[8], mp[0], '_g', mp[3], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[4], '_a', mp[6], '_c+0.00_n+0.00_o', mp[6], '_r+0.00_s+0.00.mod']

    # Model 4
    if modeln == 'model4': mc = [mp[8], mp[0], '_g', mp[3], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[5], '_a', mp[7], '_c+0.00_n+0.00_o', mp[7], '_r+0.00_s+0.00.mod']

    # Model 5
    if modeln == 'model5': mc = [mp[8], mp[1], '_g', mp[2], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[4], '_a', mp[6], '_c+0.00_n+0.00_o', mp[6], '_r+0.00_s+0.00.mod']

    # Model 6
    if modeln == 'model6': mc = [mp[8], mp[1], '_g', mp[2], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[5], '_a', mp[7], '_c+0.00_n+0.00_o', mp[7], '_r+0.00_s+0.00.mod']

    # Model 7
    if modeln == 'model7': mc = [mp[8], mp[1], '_g', mp[3], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[4], '_a', mp[6], '_c+0.00_n+0.00_o', mp[6], '_r+0.00_s+0.00.mod']

    # Model 8
    if modeln == 'model8': mc = [mp[8], mp[1], '_g', mp[3], '_m', mp[10], '_t', mp[11], '_', mp[9],'_z', \
    mp[5], '_a', mp[7], '_c+0.00_n+0.00_o', mp[7], '_r+0.00_s+0.00.mod']

    return mc

def param_swap(unmatch_id, params, temperatures, gravities, metallicities, turb_vels, masses, alphas):
    # Define the parameters list
    parameters = params

    # Set up the conditionald for each of the parameters
    if unmatch_id == 0:
        unmatch_index = temperatures.index(parameters[unmatch_id])
        new_param = temperatures[unmatch_index+1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    if unmatch_id == 1:
        unmatch_index = temperatures.index(parameters[unmatch_id])
        new_param = temperatures[unmatch_index-1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    # Gravity Conditionals
    if unmatch_id == 2:
        unmatch_index = gravities.index(parameters[unmatch_id])
        new_param = gravities[unmatch_index+1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    if unmatch_id == 3:
        unmatch_index = gravities.index(parameters[unmatch_id])
        new_param = gravities[unmatch_index-1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    # Metallicity Conditionals
    if unmatch_id == 4:
        unmatch_index = metallicities.index(parameters[unmatch_id])
        new_param = metallicities[unmatch_index+1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    if unmatch_id == 5:
        unmatch_index = metallicities.index(parameters[unmatch_id])
        new_param = metallicities[unmatch_index-1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    # Alpha enhancement
    if unmatch_id == 6:
        new_param = alphas[parameters[4]]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    if unmatch_id == 7:
        new_param = alphas[parameters[5]]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    # Mass Conditional
    if unmatch_id == 10:
        unmatch_index = masses.index(parameters[unmatch_id])
        new_param = masses[unmatch_index+1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    # Turbulent Velocity conditional
    if unmatch_id == 11:
        unmatch_index = turb_vels.index(parameters[unmatch_id])
        new_param = turb_vels[unmatch_index+1]
        parameters[unmatch_id] = new_param
        parameters = copy.deepcopy(parameters)

    return parameters

def parameter_matching(marcs_parameters_list, modeln, models, marcs_temperatures, marcs_gravities, marcs_metallicities, turb_vels, masses, alphas):
    # Create the MARCS models combinations based on the parameters
    model_list = marcs_models_components(modeln, marcs_parameters_list)

    # Define the model variables to be matched
    model = model_list[0]

    # Interpolate over the parameters of each model and match them one by one
    # with the model list
    for n in range(len(model_list)):
        try:
            # If the models is matching then add the next parameter
            if model in models:
                print(f'\n{modeln}: {model} is matching! Adding {model_list[n+1]}')
                model += model_list[n+1]

            # While the model does not match, identify which parameter does not match and swap it
            while model not in models:
                print(f'{model_list[n+1]} not matched')

                # Identify which parameter did not match and remove it from the model
                print(f'Removing {model_list[n+1]} from {model}')
                unmatch_id = marcs_parameters_list.index(model_list[n+1])

                model = model.replace(model_list[n]+model_list[n+1], model_list[n]+'')

                # Swap the unmatching parameter for a new one
                new_marcs_parameters_list = param_swap(unmatch_id, marcs_parameters_list, marcs_temperatures, \
                 marcs_gravities, marcs_metallicities, turb_vels, masses, alphas)

                model_list_new = marcs_models_components(modeln, new_marcs_parameters_list)

                model += model_list_new[n+1]

                marcs_parameters_list = copy.deepcopy(new_marcs_parameters_list)
                model_list = copy.deepcopy(model_list_new)

                # Test if the new parameter works for matching
                if model in models: break

        except(IndexError): pass

    return marcs_parameters_list

def format_to_MARCS(Tmin, Tmax, gmin, gmax, zmin, zmax, turb_vel, mass, abundances):

    # Format the MARCS parameters
    T_low, T_hi = str(Tmin), str(Tmax)
    g_low, g_hi = '{:+.1f}'.format(gmin) if gmin >= 0. else '{:.1f}'.format(gmin), '{:+.1f}'.format(gmax) if gmax >= 0. else '{:.1f}'.format(gmax)
    z_low = '{+:.2f}'.format(zmin) if zmin >= 0. else '{:.2f}'.format(zmin)
    z_hi = '+'+'{:.2f}'.format(zmax) if zmax >= 0. else '{:.2f}'.format(zmax)
    a_low, a_hi = abundances[z_low], abundances[z_hi]
    vt = '0'+str(turb_vel)
    mass = str(mass)
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

    marcs_par_list = [T_low, T_hi, g_low, g_hi, z_low, z_hi, a_low, a_hi, vt, mass, geom, comp]

    return marcs_par_list

def parameter_selection(abundances):

    # Select the parameters you want to use for the model creation
    Tmin, Tmax = int(input('Choose a Tmin value (T > 2500 K): ')), int(input('Choose a Tmax value (T < 8000 K): '))
    gmin, gmax = float(input('Choose a gmin value (logg > -1.0 dex): ')), float(input('Choose a gmax value (logg < 5.5 dex): '))
    zmin, zmax = float(input('Choose a zmin value (Z > -5.00 dex): ')), float(input('Choose a zmax value (Z < 1.00 dex): '))
    turb_vel = int(input('Chose a vt value (0 < vt < 5): '))
    mass = str(input('Choose a mass value (0.0 < mass < 5.0): '))

    # Convert these parameters to the MARCS models format
    marcs_par_list = format_to_MARCS(Tmin, Tmax, gmin, gmax, zmin, zmax, turb_vel, mass, abundances)
    par_list = [Tmin, Tmax, gmin, gmax, zmin, zmax, turb_vel, mass]

    return par_list, marcs_par_list
