#!/bin/sh

# Downloads sequence for the GRCm38 release 98 version of M. musculus (mouse) from Ensembl.

ENSEMBL_RELEASE=98
ENSEMBL_GRCm38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna
ENSEMBL_GRCm38_BASE_cDNA=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/cdna
ENSEMBL_GRCm38_BASE_ncRNA=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/ncrna
ENSEMBL_GRCm38_GTF_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus
F=Mus_musculus.GRCm38.dna.primary_assembly.fa
cDNA=Mus_musculus.GRCm38.cdna.all.fa
ncRNA=Mus_musculus.GRCm38.ncrna.fa
GTF_FILE=Mus_musculus.GRCm38.${ENSEMBL_RELEASE}.gtf

core=$(grep -c ^processor /proc/cpuinfo)

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


if [ ! -f ./genome.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi

if [ ! -f ./genocDNA.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE_cDNA}/${cDNA}.gz || (echo "Error getting $F" && exit 1)
	gunzip ${cDNA}.gz || (echo "Error unzipping $F" && exit 1)
	mv ${cDNA} genocDNA.fa
fi

if [ ! -f ./genoncRNA.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE_ncRNA}/${ncRNA}.gz || (echo "Error getting $F" && exit 1)
	gunzip ${ncRNA}.gz || (echo "Error unzipping $F" && exit 1)
	mv ${ncRNA} genoncRNA.fa
fi

if [ ! -f ./genome.gtf ] ; then
       get ${ENSEMBL_GRCm38_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
       gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
	   mv $GTF_FILE genome.gtf
fi

cat genocDNA.fa genoncRNA.fa > genomeTx.fa

echo "Making decoy file for Salmon software";
~/generateDecoyTranscriptome.sh -j ${core} -g ./genome.fa -t ./genomeTx.fa -a ./genome.gtf -o TranscriptIndex_ENSEMBL_Salmon_Decoy
echo "Build index file for Salmon software";
salmon index -p ${core} -t ./TranscriptIndex_ENSEMBL_Salmon_Decoy/gentrome.fa -i TranscriptIndex_ENSEMBL_Salmon_Index -d ./TranscriptIndex_ENSEMBL_Salmon_Decoy/decoys.txt -k 31

level=1
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
  echo "Mapping sample $((${i}/2+1)) with ${core} CPU cores: ${name1} and ${name2}; to ${name}";
  salmon quant -p ${core} -i ./TranscriptIndex_ENSEMBL_Salmon_Index -l A -1 ${name1} -2 ${name2} --seqBias --gcBias --posBias --mimicStrictBT2 -o ${name}
done