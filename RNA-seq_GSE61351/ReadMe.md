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
