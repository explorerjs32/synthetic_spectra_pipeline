#!/bin/bash
#SBATCH --job-name=parallel_job_test                   # Job name
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                    # Where to send mail
#SBATCH --nodes=1                                      # Run all processes on a single node
#SBATCH --ntasks=1                                     # Run a single task
#SBATCH --time=240:00:00                                # Time limit hrs:min:sec
#SBATCH --cpus-per-task=4                              # Number of CPU cores per task
#SBATCH --mem=100gb                                    # total memory limit
#SBATCH --output=../log_files/spec_out_%A_%a.out     # Standard output log
#SBATCH --error ../log_files/spec_error_%A_%a.error   # Standard error log
#SBATCH --array=1-100                                  # Array range
date;hostname;pwd

# Create the different copies of the scripts to run for each array
scp /blue/rezzeddine/share/fmendez/pipeline/spectra_create.sh /blue/rezzeddine/share/fmendez/pipeline/spectra_create_${SLURM_ARRAY_TASK_ID}.sh

scp /blue/rezzeddine/share/fmendez/pipeline/ts_script.com /blue/rezzeddine/share/fmendez/pipeline/ts_script_${SLURM_ARRAY_TASK_ID}.com

# The job array size will be the number of lines in the text file divided by the number of selected lines bellow
START=$SLURM_ARRAY_TASK_ID
NUMLINES=450
STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

for (( N = $START; N <= $STOP; N++ ))
do
    # Extract the parameters from each row in the parameters table
    LINE=( $(sed -n "$N"p ./temp_dir/ts_parameters.txt) )

    SCRIPT=ts_script_${SLURM_ARRAY_TASK_ID}.com
    MODEL=${LINE[0]}
    LAMMIN=${LINE[1]}
    LAMMAX=${LINE[2]}
    DELTALAM=${LINE[3]}
    METALLICITY=${LINE[4]}
    TURBVEL=${LINE[5]}
    SUFIX=${LINE[8]}
    NA=${LINE[6]}
    MG=${LINE[7]}

    echo Creating the synthetic spectra: $SUFIX running $SCRIPT in spectra_create_${SLURM_ARRAY_TASK_ID}.sh

    # Run the script
    chmod a+x /blue/rezzeddine/share/fmendez/pipeline/spectra_create_${SLURM_ARRAY_TASK_ID}.sh
    /blue/rezzeddine/share/fmendez/pipeline/spectra_create_${SLURM_ARRAY_TASK_ID}.sh $SCRIPT $MODEL $LAMMIN $LAMMAX $DELTALAM $METALLICITY $TURBVEL $SUFIX $NA $MG

    # Make sure that the created spectrum is not an empty file
    FILESIZE=$(wc -c /blue/rezzeddine/share/fmendez/Turbospectrum2019/COM-v19.1/syntspec/$SUFIX | awk '{print $1}')
    echo file size is  $FILESIZE
    
    # While the spectrum file is empty, brak the for loop
    while [$FILESIZE -eq 0]
    do
        echo file size is $FILESIZE
        break $N

        # If the spectrum file is not empty, then continue with the loop
        if [$FILESIZE -gt 0]
        then
            echo $SUFIX created sucessfully!
            continue
        fi

    echo
    
    done
done

# Remove the scripts for the array
rm spectra_create_${SLURM_ARRAY_TASK_ID}.sh
rm ts_script_${SLURM_ARRAY_TASK_ID}.com
​
date
