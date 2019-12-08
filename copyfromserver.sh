#!/bin/bash

#Change this with your personal information
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2/Results
serverPath=/home/ulg/info0939/tcrasset/Project2/Results

ssh-add $rsaPath

scp -i $rsaPath $user@$server:$serverPath/* $homePath
