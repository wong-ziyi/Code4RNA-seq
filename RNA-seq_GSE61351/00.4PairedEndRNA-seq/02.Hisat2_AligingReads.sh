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
  if [[ $name == *"d4"* ]]; then
  name=${name//d4/d04}
  elif [[ $name == *"d9"* ]]; then
  name=${name//d9/d09}
  fi
  echo "We will map sample ${i} with ${core} CPU cores (${name}): ${name1} and ${name2}";
done

y="s";
while [[ $y != ["yYnN"] ]]
do
read -e -n 1 -p "Please check the file list, then y for continous or n to exit:" y
done

if [[ $y == ["yY"] ]]; then
	for i in $(seq 1 2 $((${number}-1)))  
	do
	  seq1="NR==${i}";
	  seq2="NR==$((${i}+1))";
	  name1=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq1})
	  name2=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq2})
	  name=${name1%"_L0"*}
	  if [[ $name == *"d4"* ]]; then
	  name=${name//d4/d04}
	  elif [[ $name == *"d9"* ]]; then
	  name=${name//d9/d09}
	  fi
	  echo "Mapping sample $((${i}/2+1)) with ${core} CPU cores (${name}): ${name1} and ${name2}";
	  #Hisat2: Aliging reads
	  hisat2 -p ${core} -x ./Hisat_Genome_Index_With_ss_exon/genome_tran -1 ${name1} -2 ${name2} -S ${name}.sam
	  #Samtools sorting and converting sam to bam file
	  samtools sort -@ ${core} -o ${name}.bam ${name}.sam
	  #delete old unsorted sam file
	  rm ${name}.sam
	done
elif [[ $y == ["nN"] ]]; then
echo -e "\nExited";
fi