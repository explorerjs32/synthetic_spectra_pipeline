#!/bin/bash

# Define the arguments to be called to the script
SCRIPT=${1}

# Define the path to the TS code
path='/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/'
path_to_ts="$path$SCRIPT"

# Call in the variables we want to modify from the TurboSpectrum Code
model=$(awk -F' ' '/^foreach MODEL/ {print $3}' $path_to_ts )
lam_min=$(awk -F' ' '/^set lam_min/ {print $4}' $path_to_ts )
lam_max=$(awk -F' ' '/^set lam_max/ {print $4}' $path_to_ts )
deltalam=$(awk -F' ' '/^set deltalam/ {print $4}' $path_to_ts )
METALLIC=$(awk -F' ' '/^set METALLIC/ {print $4}' $path_to_ts )
TURBVEL=$(awk -F' ' '/^set TURBVEL/ {print $4}' $path_to_ts )
SUFFIX=$(awk -F' ' '/^set SUFFIX/ {print $4}' $path_to_ts )

# Replace the current variables in the Turv=boSpectrum script for the ones in ts_parameters.txt
sed -i "s/^foreach MODEL $model/foreach MODEL (${2})/" $path_to_ts
sed -i "s/^set lam_min = $lam_min/set lam_min = \'${3}\'/" $path_to_ts
sed -i "s/^set lam_max = $lam_max/set lam_max = \'${4}\'/" $path_to_ts
sed -i "s/^set deltalam = $deltalam/set deltalam = \'${5}\'/" $path_to_ts
sed -i "s/^set METALLIC = $METALLIC/set METALLIC = \'${6}\'/" $path_to_ts
sed -i "s/^set TURBVEL = $TURBVEL/set TURBVEL = \'${7}\'/" $path_to_ts
sed -i "s/^set SUFFIX = $SUFFIX/set SUFFIX = ${8}/" $path_to_ts
sed -i "58s/.*/${9} ${10}/" $path_to_ts
sed -i "64s#.*#${11}#" $path_to_ts

# Run the TurboSpectrum script
chmod a+x /blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/${SCRIPT}
/blue/rezzeddine/share/fmendez/Turbospectrum2019-master/COM-v19.1/${SCRIPT} #> /dev/null 2>&1
