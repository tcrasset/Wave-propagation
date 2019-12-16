#!/bin/bash


############# VARIABLES TO MODIFY ############
serverPath=/home/ulg/info0939/tcrasset/Project2
nb_processes=2
nb_threads=1
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
    
    # Load modules and compile
    module load openmpi
    mpicc -g project2_Delaunoy_Crasset*.c -std=c99 -lm -fopenmp -o waves

    # Cycle through parameter filenames in namesOfParameterFiles.txt
    # run a computation with it
    while IFS="" read -r parameter_file || [ -n "$parameter_file" ]
    do
        echo "Job with parameter file $parameter_file, $nb_processes processes and $nb_threads threads"
        sbatch --ntasks=$nb_processes --cpus-per-task=$nb_threads $serverPath/Scripts/test_parameter_file.sh $parameter_file

    done < $serverPath/Scripts/namesOfParameterFiles.txt
else
    printf "\nExiting\n"
    exit 1
fi
