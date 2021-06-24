#!/bin/bash

# Start the script by typying the name of the Star
read -p 'Enter the Star name: ' star
read -p 'Enter the starting wavlength for the synthetic spectrum [A}: ' wave_low
read -p 'Enter the ending wavelength for the synthetic spectrum [A]: ' wave_hi
read -p 'Enter the wavelength step for the synthetic spectrum [A]: ' wave_step

# This space is reserved for the Python parameters script
#
#
#

# Import the text file containing the window parameters and read through each line
cat ./temp_dir/parwin.txt | while read parline
do
  echo These are the parameters the MARCS script will use for interpolation:
  echo $parline

  # Split the parameters in aech of the lines
  par=( $parline )

  # Call in the variable values from the MARCS script
  # We want to know which is the current variable being used in the script
  # These are the model windows for each stellar parameter
  T_low=$(awk -F' ' '/^set T_low/ {print $4}' ./scratch )
  T_hi=$(awk -F' ' '/^set T_hi/ {print $4}' ./scratch )
  g_low=$(awk -F' ' '/^set g_low/ {print $4}' ./scratch )
  g_hi=$(awk -F' ' '/^set g_hi/ {print $4}' ./scratch )
  z_low=$(awk -F' ' '/^set z_low/ {print $4}' ./scratch )
  z_hi=$(awk -F' ' '/^set z_hi/ {print $4}' ./scratch )
  a_low=$(awk -F' ' '/^set a_low/ {print $4}' ./scratch )
  a_hi=$(awk -F' ' '/^set a_hi/ {print $4}' ./scratch )
  vt=$(awk -F' ' '/^set vt/ {print $4}' ./scratch )
  mass=$(awk -F' ' '/^set mass/ {print $4}' ./scratch )
  geom=$(awk -F' ' '/^set geom/ {print $4}' ./scratch )
  comp=$(awk -F' ' '/^set comp/ {print $4}' ./scratch )

  # Replace the current variables in the MARCS script by the ones in parwin.text
  sed -i "s/^set T_low = $T_low/set T_low = ${par[0]}/" ./scratch
  sed -i "s/^set T_hi = $T_hi/set T_hi = ${par[1]}/" ./scratch
  sed -i "s/^set g_low = $g_low/set g_low = ${par[2]}/" ./scratch
  sed -i "s/^set g_hi = $g_hi/set g_hi = ${par[3]}/" ./scratch
  sed -i "s/^set z_low = $z_low/set z_low = ${par[4]}/" ./scratch
  sed -i "s/^set z_hi = $z_hi/set z_hi = ${par[5]}/" ./scratch
  sed -i "s/^set a_low = $a_low/set a_low = ${par[6]}/" ./scratch
  sed -i "s/^set a_hi = $a_hi/set a_hi = ${par[7]}/" ./scratch
  sed -i "s/^set vt = $vt/set vt = ${par[8]}/" ./scratch
  sed -i "s/^set mass = $vt/set mass = ${par[9]}/" ./scratch
  sed -i "s/^set geom = $vt/set geom = ${par[10]}/" ./scratch
  sed -i "s/^set comp = $vt/set comp = ${par[11]}/" ./scratch

  # Import the text file including the model parameters and read through each line
  echo These are the model parameters the MARCS script will use for model creation:
  cat ./temp_dir/parmod.txt | while read modline
  do
    echo
    echo Creating model for these parameters: $modline

    # Split the ,odel parameters in each of the lines
    mod=( $modline )

    # Call in the variables present in the MARCS script that are used for model creation
    Tref=$(awk -F' ' '/^foreach Tref/ {print $4}' ./scratch )
    loggref=$(awk -F' ' '/^foreach loggref/ {print $4}' ./scratch )
    zref=$(awk -F' ' '/^foreach zref/ {print $4}' ./scratch )
    vtref=$(awk -F' ' '/^foreach vtref/ {print $4}' ./scratch )
    name=$(awk -F' ' '/^foreach name/ {print $4}' ./scratch )

    # Replace the model variables in MARCS by the ones in parmod.txt
    sed -i "s/^foreach Tref ( $Tref )/foreach Tref ( ${mod[0]} )/" ./scratch
    sed -i "s/^foreach loggref ( $loggref )/foreach loggref ( ${mod[1]} )/" ./scratch
    sed -i "s/^foreach zref ( $zref )/foreach zref ( ${mod[2]} )/" ./scratch
    
    sed -i "66s/.*/foreach name ( \"$star\" )/" ./scratch
   
    # Run the MARCS script
    ./scratch > /dev/null 2>&1
    
    model_out=${star}_${mod[0]}g${mod[1]}z${mod[2]}vt${vt}.TURBO
    
    mv ./test_models/$model_out /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/models/

    echo DONE CREATING THE MODEL $model_out !

    # Create the Synthetic Spectra on Turbo Spectrum
    # Import the text file including the parameters tovary in the TS script
    cat ./temp_dir/parts.txt | while read partsline
    do
      parts=( $partsline )

      echo Creating the synthetic spectrum fo z = ${mod[2]} vt = ${parts[1]} ...

      # Call in the variables from the TS script that are currently being used
      MODEL=$(awk -F' ' '/^foreach MODEL/ {print $3}' /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com )
      lam_min=$(awk -F' ' '/^set lam_min/ {print $4}' /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com )
      lam_max=$(awk -F' ' '/^set lam_max/ {print $4}' /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com )
      deltalam=$(awk -F' ' '/^set deltalam/ {print $4}' /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com )
      METALLIC=$(awk -F' ' '/^set METALLIC/ {print $4}' /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com )
      TURBVEL=$(awk -F' ' '/^set TURBVEL/ {print $4}' /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com )

      # Replace the synthetic spectrum paramters in TS
      sed -i "s/^foreach MODEL $MODEL/foreach MODEL ($model_out)/" /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com
      sed -i "s/^set lam_min = $lam_min/set lam_min = \'$wave_low\'/" /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com
      sed -i "s/^set lam_max = $lam_max/set lam_max = \'$wave_hi\'/" /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com
      sed -i "s/^set deltalam = $deltalam/set deltalam = \'$wave_step\'/" /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com
      sed -i "s/^set METALLIC = $METALLIC/set METALLIC = \'${mod[2]}\'/" /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com
      sed -i "s/^set TURBVEL = $METALLIC/set TURBVEL = \'${parts[1]}\'/" /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com

      # Run the TS script to create the synthetic spectra
      /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/ts_script.com > /dev/null 2>&1
      
      echo DONE CREATING THE SYNTHETIC SPECTRUM!
      echo It took $SECONDS seconds to create this spectrum
      echo
    done
    
  done
done

echo It took $SECONDS seconds to finish the process!


