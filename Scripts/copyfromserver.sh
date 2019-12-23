#!/bin/bash

#Change this with your personal information
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2
serverPath=/home/ulg/info0939/tcrasset/Project2

ssh-add $rsaPath
# Copy whole Results/ folder from the Cluster to the Results folder on the personal PC
# scp -i $rsaPath -C -r $user@$server:$serverPath/Results/ $homePath/Results

# scp -i $rsaPath $user@$server:$serverPath/Results/matrices_of_base_params_sriLanka_2_1/test.zip $homePath/Results


# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_for_configuration* $homePath/Results

# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_for_weakscaling* $homePath/Results
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_for_weakscaling_threads* $homePath/Results
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_strongscaling_explicit_10000_it* $homePath/Report/stats/strongscaling
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_explicit_weakscaling_* $homePath/Report/stats/weakscaling
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_strongscaling_implicit* $homePath/Report/stats/strongscaling/implicit/NEW2
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_implicit_weakscaling* $homePath/Report/stats/weakscaling/implicit
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_implicit_weakscaling_correct_processes* $homePath/Report/stats/weakscaling/implicit/
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_implicit_weakscaling_correct_threads* $homePath/Report/stats/weakscaling/implicit
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_strongscaling_implicit_correct_1* $homePath/Report/stats/strongscaling/implicit/threads
# scp -i $rsaPath $user@$server:$serverPath/Results/statistics_strongscaling_implicit_correct* $homePath/Report/stats/strongscaling/implicit/process

# scp -i $rsaPath $user@$server:$serverPath/Results/implicit_matrices_results.zip $homePath/Results
# scp -i $rsaPath $user@$server:$serverPath/Results/modspace_matrices.zip $homePath/Results
# scp -i $rsaPath $user@$server:$serverPath/Results/modtime_matrices.zip $homePath/Results
# scp -i $rsaPath $user@$server:$serverPath/Results/modtime_bis_matrices.zip $homePath/Results
# scp -i $rsaPath $user@$server:$serverPath/Results/modspace_bis_matrices.zip $homePath/Results
scp -i $rsaPath $user@$server:$serverPath/Results/mod_space_matrices_bisbis.zip $homePath/Results


# Implicit vs explicit comparison
# scp -i $rsaPath $user@$server:$serverPath/Results/matrices_of_sriLanka_explicit_result_16_1/eta_20000.dat $homePath/Results/comparison

