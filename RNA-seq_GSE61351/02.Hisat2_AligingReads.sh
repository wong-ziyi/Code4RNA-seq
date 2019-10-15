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
  #Hisat2: Aliging reads
  hisat2 -p ${core} -x ./Hisat_Genome_Index_With_ss_exon/genome_tran -U ${name1} -S ${name}.sam
  #Samtools sorting and converting sam to bam file
  samtools sort -@ ${core} -o ${name}.bam ${name}.sam
  #delete old unsorted sam file
  rm ${name}.sam
done
