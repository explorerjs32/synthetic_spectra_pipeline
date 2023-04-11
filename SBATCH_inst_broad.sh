#!/bin/bash
#SBATCH --job-name=parallel_job_test                   # Job name
#SBATCH --qos=rezzeddine-b                             # Request normal or burst cores
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                    # Where to send mail
#SBATCH --nodes=1                                      # Run all processes on a single node
#SBATCH --ntasks=1                                     # Run a single task
#SBATCH --time=50:00:00                                # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1                              # Number of CPU cores per task
#SBATCH --output=../log_files/smooth_out_%A_%a.out     # Standard output log
#SBATCH --error ../log_files/smooth_error_%A_%a.error   # Standard error log
#SBATCH --array=1-500                                # Array range
date;hostname;pwd

# Create the different copies of the scripts to run for each array
scp /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth.sh /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
scp /blue/rezzeddine/share/fmendez/pipeline/inst_broadening.f90 /blue/rezzeddine/share/fmendez/pipeline/inst_broadening_${SLURM_ARRAY_TASK_ID}.f90

# The job array size will be the number of lines in the text file divided by the number of selected lines bellow
START=$SLURM_ARRAY_TASK_ID
NUMLINES=8316
STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

for (( N = $START; N <= $STOP; N++ ))
do
    # Extract the parameters from each row in the parameters table
    LINE=( $(sed -n "$N"p ./temp_dir/inst_broad_files.txt) )

    SCRIPTIN=inst_broadening_${SLURM_ARRAY_TASK_ID}.f90
    FILEIN=${LINE[0]}
    FILEOUT_TXT=${LINE[1]}
    INST_FILEOUT=${FILEOUT_TXT}_inst.txt
    FILEOUT_BIN=${LINE[2]}
    INST_BROAD=${LINE[3]}
    SCRIPTOUT=inst_broadening_${SLURM_ARRAY_TASK_ID}

    # Convert the binary files to txt
    python3 txt_convert.py $FILEIN $FILEOUT_TXT
    sleep 2

    # Apply the instrumental broadening
    chmod a+x /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
    /blue/rezzeddine/share/fmendez/pipeline/spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh $SCRIPTIN $FILEOUT_TXT $INST_FILEOUT $INST_BROAD $SCRIPTOUT
    sleep 3

    # Convert the instrumental broadened spectra to binary
    python3 binary_convert.py $INST_FILEOUT $FILEOUT_BIN
    sleep 2

    rm $FILEOUT_TXT
done

rm ../synthetic_spectra/*txt
rm spectra_smooth_${SLURM_ARRAY_TASK_ID}.sh
rm inst_broadening_${SLURM_ARRAY_TASK_ID}.f90
rm inst_broadening_${SLURM_ARRAY_TASK_ID}
