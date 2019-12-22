#!/bin/bash

#########Change this with your personal information#########
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2
serverPath=/home/ulg/info0939/tcrasset/Project2
###############################################################

ssh-add $rsaPath

# Code
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_MAIN.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_EXPLICIT.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_EXPLICIT.h $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IMPLICIT.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IMPLICIT.h $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_SPARSE.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_SPARSE.h $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IO.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IO.h $user@$server:$serverPath

# Scripts
scp -i $rsaPath $homePath/Scripts/test_computing_configuration.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_parameter_file.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_single_implicit_weakscaling_processes.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_weak_scaling.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_parameters_files.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_implicit_weakscaling_processes.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_implicit_strongscaling_threads.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_explicit_weakscaling_threads.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_single_explicit_weakscaling_threads.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_implicit_strongscaling_processes.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_computing_configurations.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_implicit_weakscaling_threads.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_single_implicit_strongscaling.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_weak_scaling.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_single_explicit_weakscaling_processes.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_single_implicit_weakscaling_threads.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_explicit_strongscaling_processes.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_explicit_strongscaling_threads.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_explicit_weakscaling_processes.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_single_explicit_strongscaling.sh $user@$server:$serverPath/Scripts


# Parameters
# scp -i $rsaPath $homePath/Parameters/base_params_refraction.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/base_params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/process_1-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/process_2-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/process_4-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/process_8-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/process_16-params_sriLanka.txt $user@$server:$serverPath/Parameters

# scp -i $rsaPath $homePath/Parameters/threads_1-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/threads_2-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/threads_4-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/threads_8-params_sriLanka.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/threads_16-params_sriLanka.txt $user@$server:$serverPath/Parameters

scp -i $rsaPath $homePath/Parameters/sriLanka_strongscaling.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Parameters/sriLanka_implicit_25000.txt $user@$server:$serverPath/Parameters


scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_1-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_2-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_4-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_8-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_10-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_12-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_14-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/threads/explicit_threads_16-params_sriLanka.txt $user@$server:$serverPath/Parameters

scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_1-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_2-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_4-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_8-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_10-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_12-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_14-params_sriLanka.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Report/Parameters/weakscaling/process/explicit_process_16-params_sriLanka.txt $user@$server:$serverPath/Parameters

# scp -i $rsaPath $homePath/Parameters/mod_time_0,01.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_space_2000.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_time_1.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_space_25000.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_time_0,001.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_time_0,1.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_time_10.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_space_1000.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_space_10000.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_time_0,05.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_space_500.txt $user@$server:$serverPath/Parameters
# scp -i $rsaPath $homePath/Parameters/mod_space_250.txt $user@$server:$serverPath/Parameters

