#!/bin/bash

# Define the arguments to be called to the script
SCRIPT=${1}

# Define the path to the Marcs code
path='/blue/rezzeddine/share/fmendez/pipeline/'
path_to_marcs="$path$SCRIPT"

# Call in the variable values from the MARCS script
# We want to know which is the current variable being used in the script
# These are the model windows for each stellar parameter
T_low=$(awk -F' ' '/^set T_low/ {print $4}' $path_to_marcs )
T_hi=$(awk -F' ' '/^set T_hi/ {print $4}' $path_to_marcs )
g_low=$(awk -F' ' '/^set g_low/ {print $4}' $path_to_marcs )
g_hi=$(awk -F' ' '/^set g_hi/ {print $4}' $path_to_marcs )
z_low=$(awk -F' ' '/^set z_low/ {print $4}' $path_to_marcs )
z_hi=$(awk -F' ' '/^set z_hi/ {print $4}' $path_to_marcs )
a_low=$(awk -F' ' '/^set a_low/ {print $4}' $path_to_marcs )
a_hi=$(awk -F' ' '/^set a_hi/ {print $4}' $path_to_marcs )
vt=$(awk -F' ' '/^set vt/ {print $4}' $path_to_marcs )
mass=$(awk -F' ' '/^set mass/ {print $4}' $path_to_marcs )
geom=$(awk -F' ' '/^set geom/ {print $4}' $path_to_marcs )
comp=$(awk -F' ' '/^set comp/ {print $4}' $path_to_marcs )
Tref=$(awk -F' ' '/^foreach Tref/ {print $4}' $path_to_marcs )
loggref=$(awk -F' ' '/^foreach loggref/ {print $4}' $path_to_marcs )
zref=$(awk -F' ' '/^foreach zref/ {print $4}' $path_to_marcs )
vtref=$(awk -F' ' '/^foreach vtref/ {print $4}' $path_to_marcs )
name=$(awk -F' ' '/^foreach name/ {print $4}' $path_to_marcs )

# Replace the current variables in the MARCS script by the ones in parwin.text
sed -i "s/^set T_low = $T_low/set T_low = ${2}/" $path_to_marcs
sed -i "s/^set T_hi = $T_hi/set T_hi = ${3}/" $path_to_marcs
sed -i "s/^set g_low = $g_low/set g_low = ${4}/" $path_to_marcs
sed -i "s/^set g_hi = $g_hi/set g_hi = ${5}/" $path_to_marcs
sed -i "s/^set z_low = $z_low/set z_low = ${6}/" $path_to_marcs
sed -i "s/^set z_hi = $z_hi/set z_hi = ${7}/" $path_to_marcs
sed -i "s/^set vt = $vt/set vt = ${8}/" $path_to_marcs
sed -i "s/^set a_low = $a_low/set a_low = ${9}/" $path_to_marcs
sed -i "s/^set a_hi = $a_hi/set a_hi = ${10}/" $path_to_marcs
sed -i "s/^set mass = $vt/set mass = ${11}/" $path_to_marcs
sed -i "s/^set geom = $vt/set geom = ${12}/" $path_to_marcs
sed -i "s/^set comp = $vt/set comp = ${13}/" $path_to_marcs
sed -i "s/^foreach Tref ( $Tref )/foreach Tref ( ${14} )/" $path_to_marcs
sed -i "s/^foreach loggref ( $loggref )/foreach loggref ( ${15} )/" $path_to_marcs
sed -i "s/^foreach zref ( $zref )/foreach zref ( ${16} )/" $path_to_marcs

# Run the MARCS script
chmod a+x /blue/rezzeddine/share/fmendez/pipeline/$SCRIPT
/blue/rezzeddine/share/fmendez/pipeline/$SCRIPT > /dev/null 2>&1
