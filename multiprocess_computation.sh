#!/bin/bash
module load openmpi
mpicc -g project2_Delaunoy_Crasset*.c -std=c99 -lm -fopenmp -o waves

echo "Launching multiprocess_computation.sh"

for nb_processes in 2 #1 2 4 8
do
    for nb_threads in 1 #1 2 4 8
    do
        echo "Submitting job with $nb_processes processes and $nb_threads threads"
        sbatch --ntasks=$nb_processes --cpus-per-task=$nb_threads project2_launch.sh
    done
done
