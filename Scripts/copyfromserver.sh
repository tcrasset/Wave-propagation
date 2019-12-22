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
scp -i $rsaPath $user@$server:$serverPath/Results/statistics_strongscaling* $homePath/Report/stats/strongscaling

# scp -i $rsaPath $user@$server:$serverPath/Results/mod_matrices.zip $homePath/Results
