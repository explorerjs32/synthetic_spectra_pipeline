#!/bin/bash

# Define the element you want to get abundances for
element="Mg"

# Get the solar abundance from Axplund 2020
sol_abundance=`grep -H -m1 $element /home/fmendez/Documents/Research/Stellar_Spectra/Data_Tables/asplund_solar_abund.csv | cut -d' ' -f3`

# Define the abundances range
step=1
min=$(awk '{print $1-$2}' <<<"${sol_abundance} ${step}")
max=$(awk '{print $1+$2}' <<<"${sol_abundance} ${step}")

# Create an array of abundances
abund_arr=( $(seq -f "%f" ${min} 0.1 ${max}) )

# Assign each of the abundance to the code and create the synthetic spectra
idx_dif=1
index=$(awk '{print $1-$2}' <<<"21 ${idx_dif}")


echo ${abund_arr[$index]}
