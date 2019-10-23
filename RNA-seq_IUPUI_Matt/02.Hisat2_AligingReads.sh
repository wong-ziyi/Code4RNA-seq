#!/bin/bash
level=3
number=$(find . -maxdepth ${level} -type f -name "*.fastq"| wc -l)
core=$(grep -c ^processor /proc/cpuinfo)
echo "This folder has ${number} fastq files";
for i in $(seq 1 ${number})
do
  seq1="NR==${i}";
  name1=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq1})
  name=${name1%".fastq"*};
  echo "Mapping sample ${i} with ${core} CPU cores (${name}): ${name1}";

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
	  name1=$(find . -maxdepth ${level} -type f -name "*.fastq" | awk ${seq1})
	  name=${name1%".fastq"*};
	  echo "Mapping sample $((${i}/2+1)) with ${core} CPU cores (${name}): ${name1} and ${name2}";
	  #Hisat2: Aliging reads
	  hisat2 -p ${core} -x ./Hisat_Genome_Index_With_ss_exon/genome_tran -U ${name1} -S ./01.BAMfile/${name}.sam
	  #Samtools sorting and converting sam to bam file
	  samtools sort -@ ${core} -o ./01.BAMfile/${name}.bam ./01.BAMfile/${name}.sam
	  #delete old unsorted sam file
	  rm ./01.BAMfile/${name}.sam
	done
elif [[ $y == ["nN"] ]]; then
echo -e "\nExited";
fi