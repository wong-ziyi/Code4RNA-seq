#!/bin/bash
level=1
number=$(find . -maxdepth ${level} -type f -name "*.fastq"| wc -l)
core=$(grep -c ^processor /proc/cpuinfo)
echo "This folder has ${number} fastq files";
for i in $(seq 1 ${number})  
do
  seq1="NR==${i}";
  name1=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq1})
  name=${name1%".fastq"*}
  echo "Mapping sample $((${i}+1)) with ${core} CPU cores: ${name1} and ${name2}; to ${name}";
  salmon quant -p ${core} -i /mnt/d/TranscriptIndex_ENSEMBL_Salmon_Index -l A -r ${name1} --validateMappings -o ${name}
done
