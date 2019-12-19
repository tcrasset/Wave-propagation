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
    
    cd $serverPath
    module load openmpi
    mpicc -g $serverPath/project2_Delaunoy_Crasset*.c -std=c99 -lm -fopenmp -o $serverPath/waves

    echo "Launching multiprocess_computation.sh"

    # for nb_threads in 1 2 4 8
    # do
    #     echo "Submitting job with 1 processes and $nb_threads threads"
    #     sbatch --ntasks=1 --cpus-per-task=$nb_threads $serverPath/Scripts/test_computing_configuration.sh
    # done

    echo "Submitting job with 2 processes and 1 threads"
    sbatch --ntasks=2 --cpus-per-task=1 $serverPath/Scripts/test_computing_configuration.sh
    # echo "Submitting job with 2 processes and 2 threads"
    # sbatch --ntasks=2 --cpus-per-task=2 $serverPath/Scripts/test_computing_configuration.sh
    # echo "Submitting job with 4 processes and 1 threads"
    # sbatch --ntasks=4 --cpus-per-task=1 $serverPath/Scripts/test_computing_configuration.sh
else
    printf "\nExiting\n"
    exit 1
fi
