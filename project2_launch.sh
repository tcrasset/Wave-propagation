#!/bin/bash
# Submission script for NIC4
#SBATCH --job-name=TestRun
#SBATCH --time=00:45:00 # hh:mm:ss
#
#SBATCH --mem-per-cpu=32000 # megabytes
#SBATCH --partition=defq
#
#SBATCH --mail-user=tomcrasset@gmail.com
#SBATCH --mail-type=ALL
#
#SBATCH --comment=HPCProject-2
#


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_NTASKS

############# VARIABLES TO MODIFY ############
serverPath=/home/ulg/info0939/tcrasset/Project2
map=refraction
map_filename="$map.dat"
##############################################

# Output filename header
echo "Scheme,Process number,Number of processes,Number of threads,Time per process,DeltaX,DeltaY,DeltaT,s,r_threshold" >> $output_filename

# Cycles through a list of filenames in namesOfParameterFiles.txt
# and runs the computation with it
while IFS="" read -r parameter_file || [ -n "$parameter_file" ]
do
    output_filename="$serverPath/Results/statistics_nbproc_$SLURM_NTASKS-nbthreads_$SLURM_CPUS_PER_TASK-$map-$parameter_file"
    rm $output_filename # Delete old filename

    #Running $map with $parameter_file on $SLURM_NTASKS processes with $SLURM_CPUS_PER_TASK threads
    #Output in $output_filename
    mpirun -np $SLURM_NTASKS waves $serverPath/Parameters/$parameter_file $serverPath/Maps/$map_filename 0 0 0 >> $output_filename
done < namesOfParameterFiles.txt
