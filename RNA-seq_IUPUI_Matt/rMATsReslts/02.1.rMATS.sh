#!/bin/bash
python /mnt/d/rMATS/rMATS-turbo-Linux-UCS4/rmats.py --b1 rMATs_D04.txt --b2 rMATs_D09.txt --gtf /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf --od IUPUI_Matt_Dft_D09vsD04 -t paired --readLength 75 --cstat 0.0001 --libType fr-firststrand --nthread 8
python /mnt/d/rMATS/rMATS-turbo-Linux-UCS4/rmats.py --b1 rMATs_D04.txt --b2 rMATs_D18.txt --gtf /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf --od IUPUI_Matt_Dft_D18vsD04 -t paired --readLength 75 --cstat 0.0001 --libType fr-firststrand --nthread 8
python /mnt/d/rMATS/rMATS-turbo-Linux-UCS4/rmats.py --b1 rMATs_D04.txt --b2 rMATs_D28.txt --gtf /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf --od IUPUI_Matt_Dft_D28vsD04 -t paired --readLength 75 --cstat 0.0001 --libType fr-firststrand --nthread 8
python /mnt/d/rMATS/rMATS-turbo-Linux-UCS4/rmats.py --b1 rMATs_D09.txt --b2 rMATs_D18.txt --gtf /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf --od IUPUI_Matt_Dft_D18vsD09 -t paired --readLength 75 --cstat 0.0001 --libType fr-firststrand --nthread 8
python /mnt/d/rMATS/rMATS-turbo-Linux-UCS4/rmats.py --b1 rMATs_D09.txt --b2 rMATs_D28.txt --gtf /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf --od IUPUI_Matt_Dft_D28vsD09 -t paired --readLength 75 --cstat 0.0001 --libType fr-firststrand --nthread 8
python /mnt/d/rMATS/rMATS-turbo-Linux-UCS4/rmats.py --b1 rMATs_D18.txt --b2 rMATs_D28.txt --gtf /mnt/d/mm10_GTF/Mus_musculus.GRCm38.97.chr.gtf --od IUPUI_Matt_Dft_D28vsD18 -t paired --readLength 75 --cstat 0.0001 --libType fr-firststrand --nthread 8