#!/bin/bash
while getopts v:r:i:t:e: option
do
case "${option}"
in
v)   volume=${OPTARG};;
r)      SRR=${OPTARG};;
i)   Series=${OPTARG};;
t)    Start=${OPTARG};;
e)      End=${OPTARG};;
esac
done
for i in $(seq ${Start} ${End})  
do
echo "Downloading ${SRR}${Series}${i}: ftp://ftp.sra.ebi.ac.uk/${volume}/fastq/${SRR}/$(printf "%03d" ${i})/${SRR}${Series}${i}/${SRR}${Series}${i}.fastq.gz";
ascp -v -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/${volume}/fastq/${SRR}/$(printf "%03d" ${i})/${SRR}${Series}${i}/${SRR}${Series}${i}.fastq.gz ./
sudo gunzip ${SRR}${Series}${i}.fastq.gz
done
