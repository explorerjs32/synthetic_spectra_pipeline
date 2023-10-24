#!/bin/bash
#SBATCH --job-name=master_abundances                           # Job name
#SBATCH --mail-type=END,FAIL                                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=fmendez@ufl.edu                            # Where to send mail
#SBATCH --account=rezzeddine                                   # Account name
#SBATCH --qos=rezzeddine-b                                     # Request normal or burst cores
#SBATCH --nodes=1                                              # Run all processes on a single node
#SBATCH --ntasks=1                                             # Run a single task
#SBATCH --cpus-per-task=1                                      # Number of CPU cores per task
#SBATCH --mem-per-cpu=4gb                                      # total memory limit
#SBATCH --time=24:00:00                                        # Time limit hrs:min:sec
#SBATCH --output=../log_files/master_abund_out_%A_%a.out       # Standard output log
#SBATCH --error ../log_files/master_abund_error_%A_%a.error    # Standard error log
#SBATCH --array=0-10                                           # Array range


date;hostname

# Define the star we want to do the analysis for
star=$1

# Define the index to input in the second sbatch script
input_idx=${SLURM_ARRAY_TASK_ID}

# Copy the main sbatch script
scp /blue/rezzeddine/share/fmendez/pipeline/SBATCH_abundance_syntspec_library_create.sh /blue/rezzeddine/share/fmendez/pipeline/SBATCH_abundace_syntspec_library_create_${input_idx}.sh

# Define the elements
elements=("Fe" "Mg" "Si" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Co" "Ni")
element=${elements[${input_idx}]}

# Run the script
echo Runing SBATCH_abundace_syntspec_library_create_${input_idx}.sh for ${star}

sbatch /blue/rezzeddine/share/fmendez/pipeline/SBATCH_abundace_syntspec_library_create_${input_idx}.sh ${star} ${element} ${input_idx}  &

bg_pid=$!
wait $bg_pid
wait

# Remove the extra files
rm /blue/rezzeddine/share/fmendez/pipeline/SBATCH_abundace_syntspec_library_create_${input_idx}.sh

date