#!/bin/bash

#Change this with your personal information
rsaPath=/home/tom/.ssh/id_rsa.ceci
user=tcrasset
server=nic4.segi.ulg.ac.be
homePath=/home/tom/Documents/Uliege/Master2/HPC/Project2
serverPath=/home/ulg/info0939/tcrasset/Project2

ssh-add $rsaPath

# Cycles through a list of filenames in filesToCopyToServer.txt
# and copies them to the server
while IFS="" read -r filename || [ -n "$filename" ]
do
scp -i $rsaPath $homePath/$filename $user@$server:$serverPath
done < filesToCopyToServer.txt
