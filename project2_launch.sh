#!/bin/bash
# Submission script for NIC4
#SBATCH --job-name=TestRun
#SBATCH --time=00:10:00 # hh:mm:ss
#
#SBATCH --mem-per-cpu=8000 # megabytes
#SBATCH --partition=defq
#
#SBATCH --mail-user=tomcrasset@gmail.com
#SBATCH --mail-type=ALL
#
#SBATCH --comment=HPCProject-2
#


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_NTASKS


output_filename=$1
rm $output_filename # Delete old filename

serverPath=/home/ulg/info0939/tcrasset/Project2

# Output filename header
echo "Scheme,Process number,Number of processes,Number of threads,Time per process,DeltaX,DeltaY,DeltaT,s,r_threshold" >> $output_filename

# Cycles through a list of filenames in namesOfParameterFiles.txt
# and runs the computation with it

while IFS="" read -r parameter_file || [ -n "$parameter_file" ]
do
    mpirun -np $SLURM_NTASKS waves $serverPath/$parameter_file $serverPath/refraction.dat 0 0 0 >> $output_filename
done < namesOfParameterFiles.txt
