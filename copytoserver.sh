#!/bin/bash

#Change this with your personal information
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2
serverPath=/home/ulg/info0939/tcrasset/Project2

ssh-add $rsaPath

scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_MAIN.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_MAIN.h $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IO.c $user@$server:$serverPath
scp -i $rsaPath $homePath/project2_Delaunoy_Crasset_IO.h $user@$server:$serverPath

scp -i $rsaPath $homePath/project2_launch.sh $user@$server:$serverPath
scp -i $rsaPath $homePath/multiprocess_computation.sh $user@$server:$serverPath

scp -i $rsaPath $homePath/namesOfParameterFiles.txt $user@$server:$serverPath

scp -i $rsaPath $homePath/Parameters/base_params_refraction.txt $user@$server:$serverPath/Parameters
scp -i $rsaPath $homePath/Parameters/base_params_sriLanka.txt $user@$server:$serverPath/Parameters

