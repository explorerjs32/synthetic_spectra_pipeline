#!/bin/bash
#SBATCH --job-name=parallel_job_test                   # Job name
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                    # Where to send mail
#SBATCH --nodes=1                                      # Run all processes on a single node
#SBATCH --ntasks=1                                     # Run a single task
#SBATCH --time=240:00:00                                # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1                              # Number of CPU cores per task
#SBATCH --output=../log_files/smooth_out_%A_%a.out     # Standard output log
#SBATCH --error ../log_files/smooth_error_%A_%a.error   # Standard error log
#SBATCH --array=1-500                                # Array range
date;hostname;pwd

# Create the different copies of the scripts to run for each array
scp /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth.sh /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
scp /blue/rezzeddine/share/fmendez/pipeline/rot_vel.f90 /blue/rezzeddine/share/fmendez/pipeline/rot_vel_${SLURM_ARRAY_TASK_ID}.f90


# The job array size will be the number of lines in the text file divided by the number of selected lines bellow
START=$SLURM_ARRAY_TASK_ID
NUMLINES=900
STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

for (( N = $START; N <= $STOP; N++ ))
do
    # Extract the parameters from each row in the parameters table
    LINE=( $(sed -n "$N"p ./temp_dir/faltbon_parameters.txt) )

    SCRIPTIN=rot_vel_${SLURM_ARRAY_TASK_ID}.f90
    FILEIN=${LINE[0]}
    ROTVEL=${LINE[1]}
    FILEOUT=${LINE[2]}
    SCRIPTOUT=rot_vel_${SLURM_ARRAY_TASK_ID}

    echo smoothing the spectrum: $FILEIN running $SCRIPTIN in spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
    echo

    # Run the script
    chmod a+x /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
    /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh $SCRIPTIN $FILEIN $FILEOUT $ROTVEL $SCRIPTOUT

    # Make sure that the created spectrum is not an empty file
    FILESIZE=$(wc -c $FILEOUT | awk '{print $1}')
    echo file size is  $FILESIZE

    # While the spectrum file is empty, break the for loop
    while [$FILESIZE -eq 0]
    do
        echo file size is $FILESIZE
        break $N

        # If the spectrum file is not empty, then continue with the loop
        if [$FILESIZE -gt 0]
        then
            echo $FILEOUT created sucessfully!
            continue
        fi

    echo

    done

done

# Remove the scripts for the array
rm spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
rm rot_vel_${SLURM_ARRAY_TASK_ID}.f90
rm rot_vel_${SLURM_ARRAY_TASK_ID}

# Create temporary text files with the file names of the created spectra
# Separate these text files by wavelength range
find ../synthetic_spectra/ -name '*5000_5500.spec' >> ./temp_dir/spec_5000_5500.txt 
find ../synthetic_spectra/ -name '*5500_6000.spec' >> ./temp_dir/spec_5500_6000.txt

# Sort the files at each new text file
sort ./temp_dir/spec_5000_5500.txt > sorted_spec_5000_5500.txt
sort ./temp_dir/spec_5500_6000.txt > sorted_spec_5500_6000.txt

# stack the contents of both text files as columns in another file
paste -d' ' ./temp_dir/sorted_spec_5000_5500.txt ./temp_dir/sorted_spec_5500_6000.txt > ./temp_dir/spec_files.txt



​
​
date
