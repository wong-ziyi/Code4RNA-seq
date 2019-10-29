#||||||||||||||||||||||AS analysis||||||||||||||||||||||####
load("00.Results//DESeq2_Results.Rdata")
load("00.Results//Parameters.Rdata")
load("00.Results//DESeq2_Processing.Rdata")
AS_DEGsLs<-unique(
  c(
    as.vector(DEGs_ls$D09vsD04$DEGs_raw$ensembl_gene_id[DEGs_ls$D09vsD04$DEGs_raw$CutOff==1]),
    as.vector(DEGs_ls$D18vsD04$DEGs_raw$ensembl_gene_id[DEGs_ls$D18vsD04$DEGs_raw$CutOff==1]),
    as.vector(DEGs_ls$D28vsD04$DEGs_raw$ensembl_gene_id[DEGs_ls$D28vsD04$DEGs_raw$CutOff==1]),
    as.vector(DEGs_ls$D18vsD09$DEGs_raw$ensembl_gene_id[DEGs_ls$D18vsD09$DEGs_raw$CutOff==1]),
    as.vector(DEGs_ls$D28vsD09$DEGs_raw$ensembl_gene_id[DEGs_ls$D28vsD09$DEGs_raw$CutOff==1]),
    as.vector(DEGs_ls$D28vsD18$DEGs_raw$ensembl_gene_id[DEGs_ls$D28vsD18$DEGs_raw$CutOff==1])
  )
)

for (i in 1:6) {
  if(i==1){
    DEGs_Summary<-DEGs_ls[[1]]$DEGs_raw[,c(1:11,13:15)]
    colnames(DEGs_Summary)[12:14]<-paste0(colnames(DEGs_Summary)[12:14],"_",i)
    next()
  }else if(i==6){
    temp<-DEGs_ls[[i]]$DEGs_raw[,c(13:18)]
    CutOff<-as.numeric(rownames(temp)%in%All_DEGs)
    names(CutOff)<-rownames(temp)
    temp<-cbind(temp, CutOff=CutOff)
  }else{
    temp<-DEGs_ls[[i]]$DEGs_raw[,c(13:15)]
  }
  colnames(temp)[1:3]<-paste0(colnames(temp)[1:3],"_",i)
  DEGs_Summary<-merge(DEGs_Summary, temp, by='row.names', all=TRUE)
  rownames(DEGs_Summary)<-DEGs_Summary[,1]
  DEGs_Summary<-DEGs_Summary[,-1]
}

AS_AllGeLs<-DEGs_NorCounts
AS_Allgene<-DEGs_Summary

load("IUPUI_Matt_Dft_D09vsD04//.Rdata")
summary<-summary[order(abs(summary$PSI_Difference),decreasing=TRUE),]
summary<-summary[which((summary$GeneID %in% AS_DEGsLs)),]
summary<-lookup.wzy(summary,AS_AllGeLs,"GeneID")
summary<-cbind(MeanRawCount=rowMeans(summary[, 12:23]),summary)
AS_Sum.D09vsD04<-summary[summary$MeanRawCount>=500,]
load("IUPUI_Matt_Dft_D18vsD04//.Rdata")
summary<-summary[order(abs(summary$PSI_Difference),decreasing=TRUE),]
summary<-summary[which((summary$GeneID %in% AS_DEGsLs)),]
summary<-lookup.wzy(summary,AS_AllGeLs,"GeneID")
summary<-cbind(MeanRawCount=rowMeans(summary[, 12:23]),summary)
AS_Sum.D18vsD04<-summary[summary$MeanRawCount>=500,]
load("IUPUI_Matt_Dft_D28vsD04//.Rdata")
summary<-summary[order(abs(summary$PSI_Difference),decreasing=TRUE),]
summary<-summary[which((summary$GeneID %in% AS_DEGsLs)),]
summary<-lookup.wzy(summary,AS_AllGeLs,"GeneID")
summary<-cbind(MeanRawCount=rowMeans(summary[, 12:23]),summary)
AS_Sum.D28vsD04<-summary[summary$MeanRawCount>=500,]
load("IUPUI_Matt_Dft_D18vsD09//.Rdata")
summary<-summary[order(abs(summary$PSI_Difference),decreasing=TRUE),]
summary<-summary[which((summary$GeneID %in% AS_DEGsLs)),]
summary<-lookup.wzy(summary,AS_AllGeLs,"GeneID")
summary<-cbind(MeanRawCount=rowMeans(summary[, 12:23]),summary)
AS_Sum.D18vsD09<-summary[summary$MeanRawCount>=500,]
load("IUPUI_Matt_Dft_D28vsD09//.Rdata")
summary<-summary[order(abs(summary$PSI_Difference),decreasing=TRUE),]
summary<-summary[which((summary$GeneID %in% AS_DEGsLs)),]
summary<-lookup.wzy(summary,AS_AllGeLs,"GeneID")
summary<-cbind(MeanRawCount=rowMeans(summary[, 12:23]),summary)
AS_Sum.D28vsD09<-summary[summary$MeanRawCount>=500,]
load("IUPUI_Matt_Dft_D28vsD18//.Rdata")
summary<-summary[order(abs(summary$PSI_Difference),decreasing=TRUE),]
summary<-summary[which((summary$GeneID %in% AS_DEGsLs)),]
summary<-lookup.wzy(summary,AS_AllGeLs,"GeneID")
summary<-cbind(MeanRawCount=rowMeans(summary[, 12:23]),summary)
AS_Sum.D28vsD18<-summary[summary$MeanRawCount>=500,]
rm(summary,MXE.MATS.JCEC,RI.MATS.JCEC,SE.MATS.JCEC,A5SS.MATS.JCEC,A3SS.MATS.JCEC,
   PSI.A3SS.MATS.JCEC,PSI.A5SS.MATS.JCEC,PSI.MXE.MATS.JCEC,PSI.RI.MATS.JCEC,PSI.SE.MATS.JCEC)
AS_Sum.D09vsD04.SE<-as.vector(unique(AS_Sum.D09vsD04$GeneID[AS_Sum.D09vsD04$Catalog=="SE"]))
AS_Sum.D18vsD04.SE<-as.vector(unique(AS_Sum.D18vsD04$GeneID[AS_Sum.D18vsD04$Catalog=="SE"]))
AS_Sum.D28vsD04.SE<-as.vector(unique(AS_Sum.D28vsD04$GeneID[AS_Sum.D28vsD04$Catalog=="SE"]))
AS_Sum.D18vsD09.SE<-as.vector(unique(AS_Sum.D18vsD09$GeneID[AS_Sum.D18vsD09$Catalog=="SE"]))
AS_Sum.D28vsD09.SE<-as.vector(unique(AS_Sum.D28vsD09$GeneID[AS_Sum.D28vsD09$Catalog=="SE"]))
AS_Sum.D28vsD18.SE<-as.vector(unique(AS_Sum.D28vsD18$GeneID[AS_Sum.D28vsD18$Catalog=="SE"]))

#statistics
AS_Number.PSIdiff<-data.frame(
  PSI=c(">=20", ">=30", ">=40", ">=50", ">=60", ">=70", ">=80", ">=90"),
  D09vsD04=c(
    nrow(AS_Sum.D09vsD04),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.3,]),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.4,]),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.5,]),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.6,]),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.7,]),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.8,]),
    nrow(AS_Sum.D09vsD04[abs(AS_Sum.D09vsD04$IncLevelDifference)>=0.9,])
  ),
  D18vsD04=c(
    nrow(AS_Sum.D18vsD04),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.3,]),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.4,]),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.5,]),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.6,]),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.7,]),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.8,]),
    nrow(AS_Sum.D18vsD04[abs(AS_Sum.D18vsD04$IncLevelDifference)>=0.9,])
  ),
  D28vsD04=c(
    nrow(AS_Sum.D28vsD04),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.3,]),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.4,]),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.5,]),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.6,]),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.7,]),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.8,]),
    nrow(AS_Sum.D28vsD04[abs(AS_Sum.D28vsD04$IncLevelDifference)>=0.9,])
  )
)

AS_Number.Ratio<-data.frame(
  AS_Type=c("SE", "A5SS", "A3SS", "MXE", "RI"),
  D09vsD04=c(
    100*nrow(AS_Sum.D09vsD04[AS_Sum.D09vsD04$Catalog=="SE",])/nrow(AS_Sum.D09vsD04),
    100*nrow(AS_Sum.D09vsD04[AS_Sum.D09vsD04$Catalog=="A5SS",])/nrow(AS_Sum.D09vsD04),
    100*nrow(AS_Sum.D09vsD04[AS_Sum.D09vsD04$Catalog=="A3SS",])/nrow(AS_Sum.D09vsD04),
    100*nrow(AS_Sum.D09vsD04[AS_Sum.D09vsD04$Catalog=="MXE",])/nrow(AS_Sum.D09vsD04),
    100*nrow(AS_Sum.D09vsD04[AS_Sum.D09vsD04$Catalog=="RI",])/nrow(AS_Sum.D09vsD04)
  ),
  D18vsD04=c(
    100*nrow(AS_Sum.D18vsD04[AS_Sum.D18vsD04$Catalog=="SE",])/nrow(AS_Sum.D18vsD04),
    100*nrow(AS_Sum.D18vsD04[AS_Sum.D18vsD04$Catalog=="A5SS",])/nrow(AS_Sum.D18vsD04),
    100*nrow(AS_Sum.D18vsD04[AS_Sum.D18vsD04$Catalog=="A3SS",])/nrow(AS_Sum.D18vsD04),
    100*nrow(AS_Sum.D18vsD04[AS_Sum.D18vsD04$Catalog=="MXE",])/nrow(AS_Sum.D18vsD04),
    100*nrow(AS_Sum.D18vsD04[AS_Sum.D18vsD04$Catalog=="RI",])/nrow(AS_Sum.D18vsD04)
  ),
  D28vsD04=c(
    100*nrow(AS_Sum.D28vsD04[AS_Sum.D28vsD04$Catalog=="SE",])/nrow(AS_Sum.D28vsD04),
    100*nrow(AS_Sum.D28vsD04[AS_Sum.D28vsD04$Catalog=="A5SS",])/nrow(AS_Sum.D28vsD04),
    100*nrow(AS_Sum.D28vsD04[AS_Sum.D28vsD04$Catalog=="A3SS",])/nrow(AS_Sum.D28vsD04),
    100*nrow(AS_Sum.D28vsD04[AS_Sum.D28vsD04$Catalog=="MXE",])/nrow(AS_Sum.D28vsD04),
    100*nrow(AS_Sum.D28vsD04[AS_Sum.D28vsD04$Catalog=="RI",])/nrow(AS_Sum.D28vsD04)
  )
)
AS_Number.Number<-data.frame()
AS_Number.Number<-rbind(AS_Number.Number, (AS_Number.Ratio[1,2:4]/100)*AS_Number.PSIdiff[1,2:4])
AS_Number.Number<-rbind(AS_Number.Number, (AS_Number.Ratio[2,2:4]/100)*AS_Number.PSIdiff[1,2:4])
AS_Number.Number<-rbind(AS_Number.Number, (AS_Number.Ratio[3,2:4]/100)*AS_Number.PSIdiff[1,2:4])
AS_Number.Number<-rbind(AS_Number.Number, (AS_Number.Ratio[4,2:4]/100)*AS_Number.PSIdiff[1,2:4])
AS_Number.Number<-rbind(AS_Number.Number, (AS_Number.Ratio[5,2:4]/100)*AS_Number.PSIdiff[1,2:4])
AS_Number.Number<-cbind(AS_Type=c("SE", "A5SS", "A3SS", "MXE", "RI"),AS_Number.Number)

AS_VennPlot<-plotVenn(sets = list(AS_Sum.D09vsD04.SE,AS_Sum.D18vsD04.SE,AS_Sum.D28vsD04.SE),
                      systemShow=TRUE)

AS_ComGeLs<-AS_Allgene[AS_Allgene$ensembl_gene_id %in% AS_VennPlot$reg[[length(AS_VennPlot$reg)]],]
write.table(unique(AS_Sum.D28vsD04$GeneID[AS_Sum.D28vsD04$Catalog=="SE"]),"AS_2804.txt",sep = "\t", quote=FALSE)
save(AS_AllGeLs,AS_Allgene,AS_ComGeLs,AS_DEGsLs,AS_Number.Number,AS_Number.PSIdiff,AS_Number.Ratio,AS_VennPlot,
     AS_Sum.D09vsD04,AS_Sum.D09vsD04.SE,AS_Sum.D18vsD04,AS_Sum.D18vsD04.SE,AS_Sum.D18vsD09,AS_Sum.D18vsD09.SE,
     AS_Sum.D28vsD04,AS_Sum.D28vsD04.SE,AS_Sum.D28vsD09,AS_Sum.D28vsD09.SE,AS_Sum.D28vsD18,AS_Sum.D28vsD18.SE,
  file = "00.Resource//rMATS_AS_Analysis.Rdata")
load("00.Resource//rMATS_AS_Analysis.Rdata")
rm(AS_AllGeLs,AS_Allgene,AS_ComGeLs,AS_DEGsLs,AS_Number.Number,AS_Number.PSIdiff,AS_Number.Ratio,AS_VennPlot,
   AS_Sum.D09vsD04,AS_Sum.D09vsD04.SE,AS_Sum.D18vsD04,AS_Sum.D18vsD04.SE,AS_Sum.D18vsD09,AS_Sum.D18vsD09.SE,
   AS_Sum.D28vsD04,AS_Sum.D28vsD04.SE,AS_Sum.D28vsD09,AS_Sum.D28vsD09.SE,AS_Sum.D28vsD18,AS_Sum.D28vsD18.SE)
#||||||||||||||||||||||END||||||||||||||||||||||####