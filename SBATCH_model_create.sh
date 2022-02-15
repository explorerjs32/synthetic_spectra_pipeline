#!/bin/bash
#SBATCH --job-name=parallel_job_test                   # Job name
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                    # Where to send mail
#SBATCH --nodes=1                                      # Run all processes on a single node
#SBATCH --ntasks=1                                     # Run a single task
#SBATCH --cpus-per-task=4                              # Number of CPU cores per task
#SBATCH --output=../log_files/models_out_%A_%a.out     # Standard output log
#SBATCH --error ../log_files/model_error_%A_%a.error   # Standard error log
#SBATCH --array=1-10                                   # Array range
date;hostname;pwd

# Create the different copies of the scripts to run for each array
scp /blue/rezzeddine/share/fmendez/pipeline/model_create.sh /blue/rezzeddine/share/fmendez/pipeline/model_create_${SLURM_ARRAY_TASK_ID}.sh

scp /blue/rezzeddine/share/fmendez/pipeline/scratch /blue/rezzeddine/share/fmendez/pipeline/scratch_${SLURM_ARRAY_TASK_ID}

# The job array size will be the number of lines in the text file divided by the number of selected lines bellow
START=$SLURM_ARRAY_TASK_ID
NUMLINES=15
STOP=$((SLURM_ARRAY_TASK_ID*NUMLINES))
START="$(($STOP - $(($NUMLINES - 1))))"

for (( N = $START; N <= $STOP; N++ ))
do
    # Extract the parameters from each row in the parameters table
    LINE=( $(sed -n "$N"p ./temp_dir/model_parameters.txt) )
    
    SCRIPT=scratch_${SLURM_ARRAY_TASK_ID}
    TLOW=${LINE[0]}
    THI=${LINE[1]}
    GLOW=${LINE[2]}
    GHI=${LINE[3]}
    ZLOW=${LINE[4]}
    ZHI=${LINE[5]}
    VT=${LINE[6]}
    ALOW=${LINE[7]}
    AHI=${LINE[8]}
    MASS=${LINE[9]}
    GEOM=${LINE[10]}
    COMP=${LINE[11]}
    TREF=${LINE[12]}
    GREF=${LINE[13]}
    ZREF=${LINE[14]}
    
    echo Creating the atmospheric models for: T:$TREF logg:$GREF z:$ZREF running $SCRIPT in model_create_${SLURM_ARRAY_TASK_ID}.sh
    echo

    # Run the script
    chmod a+x /blue/rezzeddine/share/fmendez/pipeline/model_create_${SLURM_ARRAY_TASK_ID}.sh
    /blue/rezzeddine/share/fmendez/pipeline/model_create_${SLURM_ARRAY_TASK_ID}.sh $SCRIPT $TLOW $THI $GLOW $GHI $ZLOW $ZHI $VT $ALOW $AHI $MASS $GEOM $COMP $TREF $GREF $ZREF
done

# Remove the scripts for the array
rm model_create_${SLURM_ARRAY_TASK_ID}.sh
rm scratch_${SLURM_ARRAY_TASK_ID}

date
