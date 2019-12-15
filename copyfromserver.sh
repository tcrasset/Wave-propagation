#!/bin/bash

#Change this with your personal information
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2/Results
serverPath=/home/ulg/info0939/tcrasset/Project2/Results

ssh-add $rsaPath
# Copy whole Results/ folder from the Cluster to the Results folder on the personal PC
scp -i $rsaPath -r $user@$server:$serverPath/matrices_of_base_params_refraction/all_matrices.zip $homePath
