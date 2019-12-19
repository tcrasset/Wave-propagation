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
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_MAIN.h $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IO.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IO.h $user@$server:$serverPath

# Scripts
scp -i $rsaPath $homePath/Scripts/test_computing_configuration.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_computing_configurations.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_multiple_parameters_files.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/test_parameter_file.sh $user@$server:$serverPath/Scripts
scp -i $rsaPath $homePath/Scripts/namesOfParameterFiles.txt $user@$server:$serverPath/Scripts

# Parameters
scp -i $rsaPath $homePath/Parameters/base_params_refraction.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Parameters/base_params_sriLanka.txt $user@$server:$serverPath/Parameters

