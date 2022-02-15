#!/bin/bash

# Define the arguments to be called to the script
SCRIPT=${1}

# Define the path to the Marcs code
path='/blue/rezzeddine/share/fmendez/pipeline/'
path_to_fb="$path$SCRIPT"

# Call in the variables to be modifies in Faltv=bon
inspec=$(awk -F' ' '/^  inspec/ {print $3}' $path_to_fb )
outfil=$(awk -F' ' '/^  outfil/ {print $3}' $path_to_fb )
FWHM=$(awk -F' ' '/^  FWHM/ {print $3}' $path_to_fb )

# Replace the current variable in the falbon code for the ones in the faltbon parameters text file
sed -i "s|^  inspec = $inspec|  inspec = \'${2}\'|" $path_to_fb
sed -i "s|^  outfil = $outfil|  outfil = \'${3}\' |" $path_to_fb
sed -i "s/^  FWHM = $FWHM/  FWHM = ${4}/" $path_to_fb

# Compile the fortran script and run it
SCRIPTOUT=${5}
gfortran $path_to_fb -o $SCRIPTOUT
/blue/rezzeddine/share/fmendez/pipeline/$SCRIPTOUT > /dev/null 2>&1
