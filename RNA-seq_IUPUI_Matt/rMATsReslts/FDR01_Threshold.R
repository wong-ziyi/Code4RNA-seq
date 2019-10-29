PSI_Calculation<-function(x){
  rownames(x)<-x$ID
  v1<-c()
  v2<-c()
  v3<-c()
  Par_pb <- txtProgressBar(min = 0, max = nrow(x), style = 3)
  for (i in 1:nrow(x)) {
    li<-x$IncFormLen[i]
    ls<-x$SkipFormLen[i]
    Psi_1<-suppressWarnings(as.numeric(unlist(strsplit(as.character(x$IncLevel1[i]),","))))
    Psi_2<-suppressWarnings(as.numeric(unlist(strsplit(as.character(x$IncLevel2[i]),","))))
    Mean_PSI_1<-mean(
      c(
        li*Psi_1[1]/(li*Psi_1[1]+ls*(1-Psi_1[1])),
        li*Psi_1[2]/(li*Psi_1[2]+ls*(1-Psi_1[2])),
        li*Psi_1[3]/(li*Psi_1[3]+ls*(1-Psi_1[3]))
      ),
      na.rm = TRUE
    )
    v1<-c(v1,Mean_PSI_1)
    names(v1)[i]<-x$ID[i]
    Mean_PSI_2<-mean(
      c(
        li*Psi_2[1]/(li*Psi_2[1]+ls*(1-Psi_2[1])),
        li*Psi_2[2]/(li*Psi_2[2]+ls*(1-Psi_2[2])),
        li*Psi_2[3]/(li*Psi_2[3]+ls*(1-Psi_2[3]))
      ),
      na.rm = TRUE
    )
    v2<-c(v2,Mean_PSI_2)
    names(v2)[i]<-x$ID[i]
    PSI_Diff<- Mean_PSI_1-Mean_PSI_2
    v3<-c(v3,PSI_Diff)
    names(v3)[i]<-x$ID[i]
    setTxtProgressBar(Par_pb, i)
  }
  y<-cbind(x[,2:3],v1,v2,v3)
  colnames(y)<-c("GeneID", "geneSymbol", "PSI_1", "PSI_2", "PSI_Difference")
  return(y)
}
SE.MATS.JCEC<-read.csv("SE.MATS.JCEC.txt", header = TRUE, sep = "\t")
PSI.SE.MATS.JCEC<-PSI_Calculation(SE.MATS.JCEC)
SE.MATS.JCEC<-SE.MATS.JCEC[abs(PSI.SE.MATS.JCEC$PSI_Difference)>=0.2,]
SE.MATS.JCEC<-SE.MATS.JCEC[SE.MATS.JCEC$FDR<=0.01,]
write.table(head(SE.MATS.JCEC,0), "FDR01_SE.MATS.JCEC.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(SE.MATS.JCEC, "FDR01_SE.MATS.JCEC.txt", quote = c(2,3), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
rownames(SE.MATS.JCEC)<-SE.MATS.JCEC$ID
SE.MATS.JCEC<-merge(SE.MATS.JCEC, PSI.SE.MATS.JCEC[rownames(PSI.SE.MATS.JCEC) %in% SE.MATS.JCEC$ID,3:5],by='row.names',all=TRUE)
rownames(SE.MATS.JCEC)<-SE.MATS.JCEC$Row.names
SE.MATS.JCEC<-SE.MATS.JCEC[,-1]

RI.MATS.JCEC<-read.csv("RI.MATS.JCEC.txt", header = TRUE, sep = "\t")
PSI.RI.MATS.JCEC<-PSI_Calculation(RI.MATS.JCEC)
RI.MATS.JCEC<-RI.MATS.JCEC[abs(PSI.RI.MATS.JCEC$PSI_Difference)>=0.2,]
RI.MATS.JCEC<-RI.MATS.JCEC[RI.MATS.JCEC$FDR<=0.01,]
write.table(head(RI.MATS.JCEC,0), "FDR01_RI.MATS.JCEC.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(RI.MATS.JCEC, "FDR01_RI.MATS.JCEC.txt", quote = c(2,3), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
rownames(RI.MATS.JCEC)<-RI.MATS.JCEC$ID
RI.MATS.JCEC<-merge(RI.MATS.JCEC, PSI.RI.MATS.JCEC[rownames(PSI.RI.MATS.JCEC) %in% RI.MATS.JCEC$ID,3:5],by='row.names',all=TRUE)
rownames(RI.MATS.JCEC)<-RI.MATS.JCEC$Row.names
RI.MATS.JCEC<-RI.MATS.JCEC[,-1]

MXE.MATS.JCEC<-read.csv("MXE.MATS.JCEC.txt", header = TRUE, sep = "\t")
PSI.MXE.MATS.JCEC<-PSI_Calculation(MXE.MATS.JCEC)
MXE.MATS.JCEC<-MXE.MATS.JCEC[abs(PSI.MXE.MATS.JCEC$PSI_Difference)>=0.2,]
MXE.MATS.JCEC<-MXE.MATS.JCEC[MXE.MATS.JCEC$FDR<=0.01,]
write.table(head(MXE.MATS.JCEC,0), "FDR01_MXE.MATS.JCEC.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(MXE.MATS.JCEC, "FDR01_MXE.MATS.JCEC.txt", quote = c(2,3), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
rownames(MXE.MATS.JCEC)<-MXE.MATS.JCEC$ID
MXE.MATS.JCEC<-merge(MXE.MATS.JCEC, PSI.MXE.MATS.JCEC[rownames(PSI.MXE.MATS.JCEC) %in% MXE.MATS.JCEC$ID,3:5],by='row.names',all=TRUE)
rownames(MXE.MATS.JCEC)<-MXE.MATS.JCEC$Row.names
MXE.MATS.JCEC<-MXE.MATS.JCEC[,-1]

A3SS.MATS.JCEC<-read.csv("A3SS.MATS.JCEC.txt", header = TRUE, sep = "\t")
PSI.A3SS.MATS.JCEC<-PSI_Calculation(A3SS.MATS.JCEC)
A3SS.MATS.JCEC<-A3SS.MATS.JCEC[abs(PSI.A3SS.MATS.JCEC$PSI_Difference)>=0.2,]
A3SS.MATS.JCEC<-A3SS.MATS.JCEC[A3SS.MATS.JCEC$FDR<=0.01,]
write.table(head(A3SS.MATS.JCEC,0), "FDR01_A3SS.MATS.JCEC.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(A3SS.MATS.JCEC, "FDR01_A3SS.MATS.JCEC.txt", quote = c(2,3), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
rownames(A3SS.MATS.JCEC)<-A3SS.MATS.JCEC$ID
A3SS.MATS.JCEC<-merge(A3SS.MATS.JCEC, PSI.A3SS.MATS.JCEC[rownames(PSI.A3SS.MATS.JCEC) %in% A3SS.MATS.JCEC$ID,3:5],by='row.names',all=TRUE)
rownames(A3SS.MATS.JCEC)<-A3SS.MATS.JCEC$Row.names
A3SS.MATS.JCEC<-A3SS.MATS.JCEC[,-1]

A5SS.MATS.JCEC<-read.csv("A5SS.MATS.JCEC.txt", header = TRUE, sep = "\t")
PSI.A5SS.MATS.JCEC<-PSI_Calculation(A5SS.MATS.JCEC)
A5SS.MATS.JCEC<-A5SS.MATS.JCEC[abs(PSI.A5SS.MATS.JCEC$PSI_Difference)>=0.2,]
A5SS.MATS.JCEC<-A5SS.MATS.JCEC[A5SS.MATS.JCEC$FDR<=0.01,]
write.table(head(A5SS.MATS.JCEC,0), "FDR01_A5SS.MATS.JCEC.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(A5SS.MATS.JCEC, "FDR01_A5SS.MATS.JCEC.txt", quote = c(2,3), sep = "\t", col.names = FALSE, row.names = FALSE, append = TRUE)
rownames(A5SS.MATS.JCEC)<-A5SS.MATS.JCEC$ID
A5SS.MATS.JCEC<-merge(A5SS.MATS.JCEC, PSI.A5SS.MATS.JCEC[rownames(PSI.A5SS.MATS.JCEC) %in% A5SS.MATS.JCEC$ID,3:5],by='row.names',all=TRUE)
rownames(A5SS.MATS.JCEC)<-A5SS.MATS.JCEC$Row.names
A5SS.MATS.JCEC<-A5SS.MATS.JCEC[,-1]

summary<-rbind(
  SE.MATS.JCEC[,c("ID","GeneID","geneSymbol","chr","IncLevelDifference","PValue","FDR","PSI_1","PSI_2","PSI_Difference")],
  RI.MATS.JCEC[,c("ID","GeneID","geneSymbol","chr","IncLevelDifference","PValue","FDR","PSI_1","PSI_2","PSI_Difference")],
  MXE.MATS.JCEC[,c("ID","GeneID","geneSymbol","chr","IncLevelDifference","PValue","FDR","PSI_1","PSI_2","PSI_Difference")],
  A3SS.MATS.JCEC[,c("ID","GeneID","geneSymbol","chr","IncLevelDifference","PValue","FDR","PSI_1","PSI_2","PSI_Difference")],
  A5SS.MATS.JCEC[,c("ID","GeneID","geneSymbol","chr","IncLevelDifference","PValue","FDR","PSI_1","PSI_2","PSI_Difference")]
)
summary<-cbind(
  summary,
  Catalog=c(rep("SE", nrow(SE.MATS.JCEC)),
            rep("RI", nrow(RI.MATS.JCEC)),
            rep("MXE", nrow(MXE.MATS.JCEC)),
            rep("A3SS", nrow(A3SS.MATS.JCEC)),
            rep("A5SS", nrow(A5SS.MATS.JCEC)))
)
