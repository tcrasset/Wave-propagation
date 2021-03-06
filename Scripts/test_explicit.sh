#!/bin/bash
# Submission script for NIC4
#SBATCH --job-name=TestRun
#SBATCH --time=00:30:00 # hh:mm:ss
#
#SBATCH --mem-per-cpu=500 # megabytes
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
parameter_file=sriLanka_explicit_result_small.txt
map_file=sriLanka.dat
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
output_filename="$serverPath/Results/statistics_explicit_result_$SLURM_NTASKS-$SLURM_CPUS_PER_TASK-on-$TIMESTAMP"
##############################################

for ii in $(seq 2 2 $SLURM_NNODES);do sleep 1;echo $ii;srun -n $SLURM_NNODES -N $ii true;done

# Output filename header
echo "Scheme,Process number,Number of processes,Number of threads,Time per process,DeltaX,DeltaY,DeltaT,s,r_threshold" > $output_filename
mpirun -np $SLURM_NTASKS $serverPath/waves $serverPath/Parameters/$parameter_file $serverPath/Maps/$map_file 0 0 0 >> $output_filename
