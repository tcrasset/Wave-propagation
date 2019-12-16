#!/bin/bash

printf "Have you uploaded the needed files ? \n\n"
printf ""
printf "Have you put your uploaded files into the respective directories ? \n"
printf ""
printf "For example, parameters in Parameters/ \n\n"
printf "Did you change the files that you need to run ? \n\n"
printf "Did you update 'namesOfParameterFiles' ? \n\n"

read -r -p "Did you do all of the above ? [y/N]" response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
    
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
else
    printf "\nExiting\n"
    exit 1
fi
