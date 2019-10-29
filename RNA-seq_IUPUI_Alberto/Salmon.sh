#!/bin/bash
level=3
number=$(find . -maxdepth ${level} -type f -name "*.fastq"| wc -l)
core=$(grep -c ^processor /proc/cpuinfo)
echo "This folder has ${number} fastq files";

for i in $(seq 1 2 $((${number}-1)))  
do
  seq1="NR==${i}";
  seq2="NR==$((${i}+1))";
  name1=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq1})
  name2=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq2})
  name=${name1%"_L0"*}
  echo "Mapping sample $((${i}/2+1)) with ${core} CPU cores: ${name1} and ${name2}; to ${name}";
done

y="s";
while [[ $y != ["yYnN"] ]]
do
read -e -n 1 -p "Please check the file list, then [Yy] for continous or [Nn] to exit: " y
done

if [[ $y == ["yY"] ]]; then
	for i in $(seq 1 2 $((${number}-1)))  
	do
	  seq1="NR==${i}";
	  seq2="NR==$((${i}+1))";
	  name1=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq1})
	  name2=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq2})
	  name=${name1%"_L0"*}
	  echo "Mapping sample $((${i}/2+1)) with ${core} CPU cores: ${name1} and ${name2}; to ${name}";
	  salmon quant --validateMappings --hardFilter --dumpEq --seqBias --gcBias -p ${core} -i /mnt/g/TranscriptIndex_ENSEMBL_Salmon_Index -l A -1 ${name1} -2 ${name2} -o ${name}
	done
elif [[ $y == ["nN"] ]]; then
echo -e "\nExited";
fi