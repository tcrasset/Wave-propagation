#!/bin/bash

############# VARIABLES TO MODIFY ############
serverPath=/home/ulg/info0939/tcrasset/Project2
##############################################


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

<<<<<<< HEAD
    for nb_processes in 1 2 4 8 16 32 64
=======
    for nb_processes in 2 4 6 8
>>>>>>> 7d73510... Change scripts and add test_map
    do
        for nb_threads in 1 2 4 8
        do
            # let CPU_MEM="$MAX_MEM_PER_NODE/$nb_processes/$nb_threads"
            # CPU_MEM="${CPU_MEM}G"
            echo "Submitting job with $nb_processes processes and $nb_threads threads"
            sbatch --ntasks=$nb_processes --cpus-per-task=$nb_threads $serverPath/Scripts/test_computing_configuration.sh
        done
    done
else
    printf "\nExiting\n"
    exit 1
fi
