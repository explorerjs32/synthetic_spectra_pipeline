#!/bin/bash
date;hostname;pwd


# Define the number of rows and columns the output text file is going to have
# Each row is going to be the files that are going to be read by the chi2 script
NROWS=1
NCOLS=5

# Run the python script to reshape the list of files
python3 spectra_list.py $NCOLS $NROWS

# After we created the text file with the spectra to be fitted, then we can interpolate over each of the lines
READLINE=1
FILEOUT=results_$READLINE.txt
LINEOUT=$(sed -n "$READLINE"'p' ./temp_dir/spec_to_fit.txt)
echo $LINEOUT

# Run the chi2 fitting script
python3 chi_minimization_.py $FILEOUT $LINEOUT

date
