#!/bin/bash

#Change this with your personal information
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2/Results
serverPath=/home/ulg/info0939/tcrasset/Project2/Results

ssh-add $rsaPath

# Cycles through a list of filenames in filesToCopyFromServer.txt
# and copies them here from the server
while IFS="" read -r filename || [ -n "$filename" ]
do
scp -i $rsaPath $user@$server:$serverPath/$filename $homePath
done < filesToCopyFromServer.txt
