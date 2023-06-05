#!/bin/bash
#SBATCH --job-name=MCMC                                # Job name
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                    # Where to send mail
#SBATCH --account=rezzeddine                           # Account name
#SBATCH --qos=rezzeddine-b                             # Request normal or burst cores
#SBATCH --nodes=1                                      # Run all processes on a single node
#SBATCH --ntasks=1                                     # Run a single task
#SBATCH --cpus-per-task=80                             # Number of CPU cores per task
#SBATCH --mem=8gb                                      # total memory limit
#SBATCH --time=24:00:00                                # Time limit hrs:min:sec
#SBATCH --output=../log_files/mcmc_out_%A_%a.out       # Standard output log
#SBATCH --error ../log_files/mcmc_error_%A_%a.error    # Standard error log

date;hostname;pwd

# Define the name of the star to measure the parameters
STAR=$1

# Create a copy of the MCMC file
scp /blue/rezzeddine/share/fmendez/pipeline/mcmc.py /blue/rezzeddine/share/fmendez/pipeline/mcmc_${SLURM_JOB_ID}.py 
chmod a+x /blue/rezzeddine/share/fmendez/pipeline/mcmc.py /blue/rezzeddine/share/fmendez/pipeline/mcmc_${SLURM_JOB_ID}.py

# Run the MCMC code
python3 /blue/rezzeddine/share/fmendez/pipeline/mcmc_${SLURM_JOB_ID}.py -s ${STAR}

# Remove the copy code
rm /blue/rezzeddine/share/fmendez/pipeline/mcmc_${SLURM_JOB_ID}.py 

date
