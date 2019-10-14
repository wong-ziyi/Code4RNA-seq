# 00.Download raw fastq files
## 00-01.Create a directory to store this example and then go to that directory
```bash
$ mkdir ./RNA-seq_GSE61351
$ cd ./RNA-seq_GSE61351
```
## 00-02.Run *00.RawFASTqDownload.sh* with specific parameters to download the raw fastq files
```bash
$ ./00.RawFASTqDownload.sh -v vol1 -r SRR157 -i 367 -t 0 -e 5
```
# 01-00.Build or download the Hisat2 index
If you have a computer with 200GB physical memory, you could change the line 14 "ENSEMBL_RELEASE=84" in *make_grcm38_tran.sh* into the newest release version (e.g. "ENSEMBL_RELEASE=98"), then excute it to build a Hisat2 index with ss (splice site) and exon infromation based on the newest ensembl genome anotation. Otherwise, you could download the pre-build Hisat2 index in 2016 with ensembl genome anoation release 84 from [here](https://cloud.biohpc.swmed.edu/index.php/s/grch37_tran/download) (See all available pre-build index in [here](https://ccb.jhu.edu/software/hisat2/index.shtml)).
