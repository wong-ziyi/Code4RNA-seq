#!/bin/bash
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d4_1_S45 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d4_1_S45_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d4_1_S45_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d4_2_S46 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d4_2_S46_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d4_2_S46_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d4_3_S47 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d4_3_S47_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d4_3_S47_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d9_1_S48 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d9_1_S48_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d9_1_S48_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d9_2_S49 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d9_2_S49_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d9_2_S49_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d9_3_S50 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d9_3_S50_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d9_3_S50_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d18_1_S51 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d18_1_S51_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d18_1_S51_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d18_2_S52 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d18_2_S52_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d18_2_S52_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d18_3_S53 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d18_3_S53_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d18_3_S53_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d28_1_S54 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d28_1_S54_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d28_1_S54_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d28_2_S55 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d28_2_S55_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d28_2_S55_L006_R2_001.fastq
STAR --runThreadN 8 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFileNamePrefix ./SW3_d28_3_S56 --genomeDir /mnt/d/GenomeIndex_ENSEMBL_STAR_ReadLength74 --readFilesIn /mnt/d/RNA-seq_IUPUI_Matt/SW3_d28_3_S56_L006_R1_001.fastq /mnt/d/RNA-seq_IUPUI_Matt/SW3_d28_3_S56_L006_R2_001.fastq
samtools view -H SW3_d4_1_S45Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d4_1_S45Aligned.sortedByCoord.out.bam > SW3_d4_1_S45Aligned.sortedByCoord.out.chr.bam
rm SW3_d4_1_S45Aligned.sortedByCoord.out.bam
samtools view -H SW3_d4_2_S46Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d4_2_S46Aligned.sortedByCoord.out.bam > SW3_d4_2_S46Aligned.sortedByCoord.out.chr.bam
rm SW3_d4_2_S46Aligned.sortedByCoord.out.bam
samtools view -H SW3_d4_3_S47Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d4_3_S47Aligned.sortedByCoord.out.bam > SW3_d4_3_S47Aligned.sortedByCoord.out.chr.bam
rm SW3_d4_3_S47Aligned.sortedByCoord.out.bam
samtools view -H SW3_d9_1_S48Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d9_1_S48Aligned.sortedByCoord.out.bam > SW3_d9_1_S48Aligned.sortedByCoord.out.chr.bam
rm SW3_d9_1_S48Aligned.sortedByCoord.out.bam
samtools view -H SW3_d9_2_S49Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d9_2_S49Aligned.sortedByCoord.out.bam > SW3_d9_2_S49Aligned.sortedByCoord.out.chr.bam
rm SW3_d9_2_S49Aligned.sortedByCoord.out.bam
samtools view -H SW3_d9_3_S50Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d9_3_S50Aligned.sortedByCoord.out.bam > SW3_d9_3_S50Aligned.sortedByCoord.out.chr.bam
rm SW3_d9_3_S50Aligned.sortedByCoord.out.bam
samtools view -H SW3_d18_1_S51Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d18_1_S51Aligned.sortedByCoord.out.bam > SW3_d18_1_S51Aligned.sortedByCoord.out.chr.bam
rm SW3_d18_1_S51Aligned.sortedByCoord.out.bam
samtools view -H SW3_d18_2_S52Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d18_2_S52Aligned.sortedByCoord.out.bam > SW3_d18_2_S52Aligned.sortedByCoord.out.chr.bam
rm SW3_d18_2_S52Aligned.sortedByCoord.out.bam
samtools view -H SW3_d18_3_S53Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d18_3_S53Aligned.sortedByCoord.out.bam > SW3_d18_3_S53Aligned.sortedByCoord.out.chr.bam
rm SW3_d18_3_S53Aligned.sortedByCoord.out.bam
samtools view -H SW3_d28_1_S54Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d28_1_S54Aligned.sortedByCoord.out.bam > SW3_d28_1_S54Aligned.sortedByCoord.out.chr.bam
rm SW3_d28_1_S54Aligned.sortedByCoord.out.bam
samtools view -H SW3_d28_2_S55Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d28_2_S55Aligned.sortedByCoord.out.bam > SW3_d28_2_S55Aligned.sortedByCoord.out.chr.bam
rm SW3_d28_2_S55Aligned.sortedByCoord.out.bam
samtools view -H SW3_d28_3_S56Aligned.sortedByCoord.out.bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - SW3_d28_3_S56Aligned.sortedByCoord.out.bam > SW3_d28_3_S56Aligned.sortedByCoord.out.chr.bam
rm SW3_d28_3_S56Aligned.sortedByCoord.out.bam
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d4_1_S45Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d4_1_S45.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d4_2_S46Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d4_2_S46.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d4_3_S47Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d4_3_S47.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d9_1_S48Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d9_1_S48.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d9_2_S49Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d9_2_S49.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d9_3_S50Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d9_3_S50.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d18_1_S51Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d18_1_S51.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d18_2_S52Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d18_2_S52.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d18_3_S53Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d18_3_S53.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d28_1_S54Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d28_1_S54.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d28_2_S55Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d28_2_S55.counts
featurecounts /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d28_3_S56Aligned.sortedByCoord.out.chr.bam -T 32 -s 2 -p -Q 10 -a /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf -o /mnt/d/RNA-Seq_IUPUI_Matt/SW3_d28_3_S56.counts