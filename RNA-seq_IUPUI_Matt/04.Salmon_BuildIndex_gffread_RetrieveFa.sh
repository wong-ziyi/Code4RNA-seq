#!/bin/bash
echo "Retrieve transcriptome fa file based on merged assembled transcriptome gtf file";
gffread -w ${Out_name}.fa -g ./genome.fa ${Out_name}.gtf
echo "Making decoy file for Salmon software";
~/generateDecoyTranscriptome.sh -j ${core} -g ./genome.fa -t ./${Out_name}.fa -a ./${Out_name}.gtf -o TranscriptIndex_ENSEMBL_Salmon_Decoy_Assemble
echo "Build index file for Salmon software";
salmon index -p ${core} -t ./TranscriptIndex_ENSEMBL_Salmon_Decoy_Assemble/gentrome.fa -i TranscriptIndex_ENSEMBL_Salmon_Index_Assemble -d ./TranscriptIndex_ENSEMBL_Salmon_Decoy_Assemble/decoys.txt -k 31