#!/bin/sh

# Downloads sequence for the GRCm38 release 98 version of M. musculus (mouse) from Ensembl.

ENSEMBL_RELEASE=98
ENSEMBL_GRCm38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna
ENSEMBL_GRCm38_GTF_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/mus_musculus
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

F=Mus_musculus.GRCm38.dna.primary_assembly.fa
if [ ! -f ./genome.fa ] ; then
	get ${ENSEMBL_GRCm38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi

if [ ! -f genome.gtf ] ; then
       get ${ENSEMBL_GRCm38_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
       gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
	   mv $GTF_FILE genome.gft
fi

Rscript --vanilla ./DisGenels.R genome.gtf DisGenels.txt


echo "Remove the genes for Mt_rna, pseudogenes, immunoglobulin rna, and ribosomal rna from .gtf file";
mapfile -t DisGenels < ./DisGenels.txt;
echo "There are ${#DisGenels[@]} groups";
for i in $(seq 1 ${#DisGenels[@]})
do
	echo "${i} of ${#DisGenels[@]}";
	grep -v -E "${DisGenels[$((${i}-1))]}" genome.gtf > temp.gtf;
	rm genome.gtf;
	mv temp.gtf genome.gtf;
done
mv genome.gtf genomeTx.gtf
echo "Retrieve transcriptome .fa file based on merged assembled transcriptome gtf file";
gffread -w ./genomeTx.fa -g ./genome.fa genome.gtf;

echo "Making decoy file for Salmon software";
~/generateDecoyTranscriptome.sh -j ${core} -g ./genome.fa -t ./genomeTx.fa -a ./genomeTx.gtf -o TranscriptIndex_ENSEMBL_Salmon_Decoy
echo "Build index file for Salmon software";
salmon index -p ${core} -t ./TranscriptIndex_ENSEMBL_Salmon_Decoy/gentrome.fa -i TranscriptIndex_ENSEMBL_Salmon_Index -d ./TranscriptIndex_ENSEMBL_Salmon_Decoy/decoys.txt -k 31