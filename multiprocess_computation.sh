#!/bin/bash
module load openmpi
mpicc -g project2_Delaunoy_Crasset.c -std=c99 -lm -fopenmp -o waves

serverPath=/home/ulg/info0939/tcrasset/Project2/Results


for nb_processes in 2 #1 2 4 8
do
    for nb_threads in 1 #1 2 4 8
    do
        filename=$serverPath/project2_out_nbproc_$nb_processes-nbthreads_$nb_threads.txt
        echo "HPC_Project2 $nb_processes processes and $nb_threads threads, Filename output: $filename"
        sbatch --ntasks=$nb_processes --cpus-per-task=$nb_threads project2_launch.sh $filename
    done
done
