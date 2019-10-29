#!/bin/bash

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
	sudo gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi

if [ ! -f ./genocDNA.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE_cDNA}/${cDNA}.gz || (echo "Error getting $F" && exit 1)
	sudo gunzip ${cDNA}.gz || (echo "Error unzipping $F" && exit 1)
	mv ${cDNA} genocDNA.fa
fi

if [ ! -f ./genoncRNA.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE_ncRNA}/${ncRNA}.gz || (echo "Error getting $F" && exit 1)
	sudo gunzip ${ncRNA}.gz || (echo "Error unzipping $F" && exit 1)
	mv ${ncRNA} genoncRNA.fa
fi

if [ ! -f ./genome.gtf ] ; then
       get ${ENSEMBL_GRCm38_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
       sudo gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
	   mv $GTF_FILE genome.gtf
fi

if [ ! -f ./genomeTx.fa ] ; then
       cat genocDNA.fa genoncRNA.fa > genomeTx.fa
fi

if { [ ! -d ./TranscriptIndex_ENSEMBL_Salmon_Index ] || [ ! -d ./TranscriptIndex_ENSEMBL_Salmon_Decoy ]; }; then
       echo "Making decoy file for Salmon software";
	   ~/generateDecoyTranscriptome.sh -j ${core} -g ./genome.fa -t ./genomeTx.fa -a ./genome.gtf -o TranscriptIndex_ENSEMBL_Salmon_Decoy
	   echo "Build index file for Salmon software";
	   salmon index -p ${core} -t ./TranscriptIndex_ENSEMBL_Salmon_Decoy/gentrome.fa -i TranscriptIndex_ENSEMBL_Salmon_Index -d ./TranscriptIndex_ENSEMBL_Salmon_Decoy/decoys.txt -k 31
fi