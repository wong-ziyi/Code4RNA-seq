#!/bin/bash
while getopts n: option
do
case "${option}"
in
n)   Out_name=${OPTARG};;
esac
done

ENSEMBL_RELEASE=84
ENSEMBL_GRCm38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna
ENSEMBL_GRCm38_GTF_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus
GTF_FILE=Mus_musculus.GRCm38.${ENSEMBL_RELEASE}.gtf

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

if [ ! -f ./genome.gtf ] ; then
	get ${ENSEMBL_GRCm38_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
	sudo gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
	mv ${GTF_FILE} genome.gtf
fi

level=3
number=$(find . -maxdepth ${level} -type f -name "*.bam"| wc -l)
core=$(grep -c ^processor /proc/cpuinfo)
echo "This folder has ${number} bam files";

if [ -f ./02.AssembledGTF/GTFls.txt ] ; then
	rm ./02.AssembledGTF/GTFls.txt
fi

for i in $(seq 1 ${number})
do
  seq1="NR==${i}";
  name=$(find . -maxdepth ${level} -type f -name "*.bam" | awk ${seq1})
  name=${name#"./01.BAMfile/"*}
  name=${name%".bam"*}
  echo "Assemble transcriptome from sample ${name}";
  stringtie -p ${core} -G ./genome.gtf -o ./02.AssembledGTF/${name}.gtf -l ${name} ./01.BAMfile/${name}.bam
  echo ./02.AssembledGTF/${name}.gtf >> ./02.AssembledGTF/GTFls.txt; 
done
echo "Mergeing assembled transcriptome from ${number} samples";
stringtie --merge -G ./genome.gtf -o ${Out_name}.gtf ./02.AssembledGTF/GTFls.txt

F=Mus_musculus.GRCm38.dna.toplevel.fa
if [ ! -f ./genome.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi