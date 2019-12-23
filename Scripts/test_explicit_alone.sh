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

    echo "Submitting job with 16 processes and 1 threads"
    sbatch --ntasks=16 --cpus-per-task=1 $serverPath/Scripts/test_explicit.sh

else
    printf "\nExiting\n"
    exit 1
fi
