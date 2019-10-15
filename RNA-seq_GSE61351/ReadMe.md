# Content Index
- [Introducation](#introducation)
- [00.Download raw fastq files](#00)
 - [00-01.Create a directory to store this example and then navigate into that directory](#00-01)
 - [00-02.Run ***00.Download raw fastq files*** with specific parameters to download the raw fastq files](#00-02)
- [01.Build or download the Hisat2 index](#01)
- [02.Hisat2 alinging reads](#02)

# Introducation

<p align="justify">
  Here is an example for study and practice the processing and analysis of RNA-seq data. In this section, we will build up our own transcriptome anotation based on the data from GSE61351 by assembly their .bam file. Then we will do quantitation for them on transcription resolution by Salmon software, and identifing significant transcripts switching, alternative splicing, differentiate expressed genes or transcripts, and finally do functional enrichment analysis by <a href="https://yulab-smu.github.io/clusterProfiler-book/">clusterProfiler</a> R software.
</p>

<h1 id="00">00.Download raw fastq files</h1>
<h2 id="00-01">00-01.Create a directory to store this example and then navigate into that directory</h2> 

```bash
$ mkdir ./RNA-seq_GSE61351
$ cd ./RNA-seq_GSE61351
```
<h2 id="00-02">
 00-02.Run 
  <a href="https://github.com/wong-ziyi/Code4RNA-seq/blob/master/RNA-seq_GSE61351/00.RawFASTqDownload.sh">
    <i>00.Download raw fastq files</i>
  </a>
 with specific parameters to download the raw fastq files
</h2>

```bash
$ ./00.RawFASTqDownload.sh -v vol1 -r SRR157 -i 367 -t 0 -e 5
```
<h1 id="01">01.Build or download the Hisat2 index</h1> 
# 

<p align="justify">
If you have a computer with 200GB physical memory, you could change the line 14 "ENSEMBL_RELEASE=84" in <a href="https://github.com/wong-ziyi/Code4RNA-seq/blob/master/RNA-seq_GSE61351/01.make_grcm38_tran.sh"><b><i>01.make_grcm38_tran.sh</i></b></a> into the newest release version (e.g. "ENSEMBL_RELEASE=98"), then excute it to build a Hisat2 index with ss (splice site) and exon infromation based on the newest ensembl genome anotation. Otherwise, you could download the pre-build Hisat2 index in 2016 with ensembl genome anoation release 84 from <a href="https://cloud.biohpc.swmed.edu/index.php/s/grch37_tran/download">here</a> (See all available pre-build index in <a href="https://ccb.jhu.edu/software/hisat2/index.shtml">here</a>).
</p>

Create a directory to store the index
```bash
$ mkdir ./Hisat_Genome_Index_With_ss_exon
```
Copy *01.make_grcm38_tran.sh* into that directory, navigate into that directory, 
```bash
$ mv 01.make_grcm38_tran.sh ./Hisat_Genome_Index_With_ss_exon
$ cd ./Hisat_Genome_Index_With_ss_exon
```
and run this bash file (**option#1**)
```bash
$ ./01.make_grcm38_tran.sh
```
or download the pre-build Hisat2 index then un-zip it (**option#2**)
```bash
$ wget -P ./ "https://cloud.biohpc.swmed.edu/index.php/s/grcm38_tran/download"
$ tar -xvzf grch37.tar.gz ./
```
<h1 id="02">02.Hisat2 alinging reads</h1> 

```bash
$ conda activate py3
$ ./02.Hisat2_AligingReads.sh
```
