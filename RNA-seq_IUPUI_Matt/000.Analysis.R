#Packages & functions Loading####
invisible(
  lapply(
    c("Functions.R",
      "Functions_Query.R",
      "Function00_ParameterSetting.R",
      "Function01_Salmon_Tximport_RUVSeq.R",
      "Function02_DESeq2_PreProcessing.R",
      "Function03_DESeq2_GSEAenrich_Heatmap_PCA_MultiComparsion.R"),
    source
  )
)
invisible(PackageCheckInstall())
options(url.method="libcurl")
wd<-getwd()
#||||||||||||||||||||||(Salmon) Identify DEs by DESeq2 Packages||||||||||||||||||||||####
#Par_ConLs_N<-cbind(Load_ExDesign,model.matrix(~ CellType*Loading*Drug, Load_ExDesign)) #Numerical model matrices
P0_ParameterSetting(
  #|----Set up parameters----|####
  Par_Test="none",      #"LRT", "Wald", or "none" 
  Par_PorFRR="padj",    #"pvalue" or "padj"
  Par_Rep=3,            #How many replicates for each condition
  Par_ReN=5,            #Cut-off value of read number
  Par_SaN=2,            #Cut-off value for expressed samples number
  Par_FaN=3,            #Factors column numner in design table
  Par_LFC=1,            #Cut-off value for Log2FC
  Par_FDR=0.01,         #Cut-off value for FDR
  Par_FDR4FEA=0.05,     #Cut-off (qValue) value for GO enrichment by GSEA
  Par_FDR4KEG=0.05,     #Cut-off (qValue) value for KEGG enrichment by GSEA
  Par_FDR4FEA_OR=0.01,  #Cut-off (qValue) value for GO enrichment by Over Representation Test
  Par_FDR4KEG_OR=0.05,  #Cut-off (qValue) value for KEGG enrichment by Over Representation Test
  #|----Set up Design Table----|####
  Load_ExDesign=read.csv("ExDesign.csv"),
  Par_ConLs_C=list(
    c("RUV_x","D09","D04"),
    c("RUV_x","D18","D04"),
    c("RUV_x","D28","D04"),
    c("RUV_x","D18","D09"),
    c("RUV_x","D28","D09"),
    c("RUV_x","D28","D18")
  ) #Combined factors list
)
P1_Salmon_Tximport_RUVSeq()
dev.off()
P2_DESeq2_PreProcessing()
P3_DESeq2_GSEA_MultiComparsion()
setwd(wd)
#||||||||||||||||||||||Correlation to metabolic results||||||||||||||||||||||####
load("00.Results//DESeq2_Results.Rdata")
load("00.Results//Parameters.Rdata")
load("00.Results//DESeq2_Processing.Rdata")
Metabolic_Me<-read.csv("Metabolic_Media.csv")
rownames(Metabolic_Me)<-Metabolic_Me[,1]
Metabolic_Me<-Metabolic_Me[,-1]
Metabolic_Me<-as.matrix(Metabolic_Me)

Metabolic_Ce<-read.csv("Metabolic_CellLysis.csv")
rownames(Metabolic_Ce)<-Metabolic_Ce[,1]
Metabolic_Ce<-Metabolic_Ce[,-1]
Metabolic_Ce<-as.matrix(Metabolic_Ce)

MeanTPM<-DEGs_log2TPM                  %>%
         apply(., 2, wzy.2power)       %>%
         apply(., 1, BinMean, every=3) %>%
         t()
colnames(MeanTPM)<-c("D04_MeanTPM","D09_MeanTPM","D18_MeanTPM","D28_MeanTPM")

x<-matrix(nrow = nrow(MeanTPM), ncol = nrow(Metabolic_Ce))
for (i in 1:nrow(Metabolic_Ce)) {
  y<-c()
  for (j in 1:nrow(MeanTPM)) {
    y<-c(y, cor(MeanTPM[j,], Metabolic_Ce[i,]))
  }
  x[,i]<-y
}
colnames(x)<-rownames(Metabolic_Ce)
rownames(x)<-rownames(MeanTPM)
Correlations<-x
rm(x)
x<-matrix(nrow = nrow(MeanTPM), ncol = nrow(Metabolic_Me))
for (i in 1:nrow(Metabolic_Me)) {
  y<-c()
  for (j in 1:nrow(MeanTPM)) {
    y<-c(y, cor(MeanTPM[j,], Metabolic_Me[i,]))
  }
  x[,i]<-y
}
colnames(x)<-rownames(Metabolic_Me)
rownames(x)<-rownames(MeanTPM)
Correlations<-merge(x=Correlations, y=x, by='row.names', all=TRUE)
rownames(Correlations)<-Correlations[,1]
Correlations<-Correlations[,-1]



All_DEGs<-unique(
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


DEGs_Summary<-merge(DEGs_Summary, MeanTPM, by='row.names', all=TRUE)
rownames(DEGs_Summary)<-DEGs_Summary[,1]
DEGs_Summary<-DEGs_Summary[,-1]
DEGs_Summary<-merge(DEGs_Summary, Correlations, by='row.names', all=TRUE)
rownames(DEGs_Summary)<-DEGs_Summary[,1]
DEGs_Summary<-DEGs_Summary[,-1]


xlsx<-createWorkbook()
addWorksheet(xlsx, "Transcript_Metabolic")
writeDataTable(xlsx, 1, DEGs_Summary, startCol = 1, startRow = 2)
writeData(xlsx, 1, "Basic Information", startCol=1, startRow=1)
writeData(xlsx, 1, "Threshold", startCol=3, startRow=1)
writeData(xlsx, 1, 0.99, startCol=4, startRow=1)
style<-createStyle(fgFill="#ABBEDE", halign="LEFT", border="bottom", textDecoration="Bold")
addStyle(xlsx, 1, style=style, cols=c(1:11), rows=1)
j<-1
name<-c("D09 vs D04", "D18 vs D04", "D28 vs D04", "D18 vs D09", "D28 vs D09", "D28 vs D18")
color<-c("#EDEDED","#DBDBDB","#C9C9C9","#FCE4D6","#F8CBAD","#F4B084")
for (i in seq(12,27,3)) {
  mergeCells(xlsx, 1, cols=i:(i+2), rows=1)
  writeData(xlsx, 1, name[j] , startCol=i, startRow=1)
  style<-createStyle(fgFill=color[j], halign="CENTER", border="bottom", textDecoration="Bold")
  addStyle(xlsx, 1, style=style, cols=i, rows=1)
  j<-j+1
}
#Likelihood Ratio Test
mergeCells(xlsx, 1, cols=30:33, rows=1)
writeData(xlsx, 1, "Likelihood Ratio Test", startCol=30, startRow=1)
style<-createStyle(fgFill="#C7D3E9", halign="CENTER", border="bottom", textDecoration="Bold")
addStyle(xlsx, 1, style=style, cols=30, rows=1)
#Mean TPM
mergeCells(xlsx, 1, cols=34:37, rows=1)
writeData(xlsx, 1, "Mean TPM", startCol=34, startRow=1)
style<-createStyle(fgFill="#DDEBF7", halign="CENTER", border="bottom", textDecoration="Bold")
addStyle(xlsx, 1, style=style, cols=34, rows=1)
#Pearson's r to Metabolomics from cell lysis
mergeCells(xlsx, 1, cols=38:57, rows=1)
writeData(xlsx, 1, "Pearson's r to Metabolomics from cell lysis", startCol=38, startRow=1)
style<-createStyle(fgFill="#E3E9F4", halign="CENTER", border="bottom", textDecoration="Bold")
addStyle(xlsx, 1, style=style, cols=38, rows=1)
#Pearson's r to Metabolomics from media
mergeCells(xlsx, 1, cols=58:90, rows=1)
writeData(xlsx, 1, "Pearson's r to Metabolomics from media", startCol=58, startRow=1)
style<-createStyle(fgFill="#E2EFDA", halign="CENTER", border="bottom", textDecoration="Bold")
addStyle(xlsx, 1, style=style, cols=58, rows=1)
#Conditional Formatting red for positive green for negative
posStyle<- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
negStyle<- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
conditionalFormatting(xlsx, 1, cols=38:90, row=3:(nrow(DEGs_Summary)+3), rule="<=-ABS($D$1)", style = negStyle)
conditionalFormatting(xlsx, 1, cols=38:90, row=3:(nrow(DEGs_Summary)+3), rule=">=ABS($D$1)", style = posStyle)
freezePane(xlsx, 1, firstActiveRow = 3, firstActiveCol = 5, firstRow = FALSE, firstCol = FALSE)
saveWorkbook(xlsx, "Summary_Gene_Metabolomics_Correlation.xlsx", overwrite = TRUE)

for (i in 1:6) {
  name1<-c("D09vsD04", "D18vsD04", "D28vsD04", "D18vsD09", "D28vsD09", "D28vsD18")
  name2<-c("BP","MF","CC","KEGG")
  xlsx<-createWorkbook()
  for (j in 1:4) {
    if(nrow(as.data.frame(DEGs_ls[[i]]$GSEA_Enrich[[j]]$Sig))>0){
      temp1<-DEGs_ls[[i]]$GSEA_Enrich[[j]]$Sig[,1:11]
      temp2<-DEGs_ls[[i]]$GSEA_Enrich[[j]]$Sig[,12:15]
      temp2<-2^temp2
      temp2<-as.matrix(temp2)
      colnames(temp2)<-c("D04_MeanTPM","D09_MeanTPM","D18_MeanTPM","D28_MeanTPM")
      x<-matrix(nrow = nrow(temp2), ncol = nrow(Metabolic_Ce))
      for (m in 1:nrow(Metabolic_Ce)) {
        y<-c()
        for (n in 1:nrow(temp2)) {
          y<-c(y, cor(temp2[n,], Metabolic_Ce[m,]))
        }
        x[,m]<-y
      }
      colnames(x)<-rownames(Metabolic_Ce)
      rownames(x)<-rownames(temp2)
      Correlations<-x
      rm(x)
      x<-matrix(nrow = nrow(temp2), ncol = nrow(Metabolic_Me))
      for (m in 1:nrow(Metabolic_Me)) {
        y<-c()
        for (n in 1:nrow(temp2)) {
          y<-c(y, cor(temp2[n,], Metabolic_Me[m,]))
        }
        x[,m]<-y
      }
      colnames(x)<-rownames(Metabolic_Me)
      rownames(x)<-rownames(temp2)
      Correlations<-merge(x=Correlations, y=x, by='row.names', all=TRUE)
      rownames(Correlations)<-Correlations[,1]
      Correlations<-Correlations[,-1]
      temp<-merge(temp1,temp2, by='row.names', all=TRUE)
      rownames(temp)<-temp[,1]
      temp<-temp[,-1]
      temp<-merge(temp, Correlations, by='row.names', all=TRUE)
      rownames(temp)<-temp[,1]
      temp<-temp[,-1]
      addWorksheet(xlsx, name2[j])
      writeDataTable(xlsx, name2[j], temp, startCol = 1, startRow = 2)
      setColWidths(xlsx, name2[j], cols=1, widths = 12, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
      setColWidths(xlsx, name2[j], cols=2, widths = 45, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
      writeData(xlsx, name2[j], "Threshold", startCol=3, startRow=1)
      writeData(xlsx, name2[j], 0.99, startCol=4, startRow=1)
      writeData(xlsx, name2[j], "GSEA Enrichment Analysis", startCol=1, startRow=1)
      style<-createStyle(fgFill="#ABBEDE", halign="LEFT", border="bottom", textDecoration="Bold")
      addStyle(xlsx, name2[j], style=style, cols=1:11, rows=1)
      mergeCells(xlsx, name2[j], cols=12:15, rows=1)
      writeData(xlsx, name2[j], "Total Mean of genes of each condition", startCol=12, startRow=1)
      style<-createStyle(fgFill="#C9C9C9", halign="CENTER", border="bottom", textDecoration="Bold")
      addStyle(xlsx, name2[j], style=style, cols=12, rows=1)
      mergeCells(xlsx, name2[j], cols=16:35, rows=1)
      writeData(xlsx, name2[j], "Metabolomics from cell lysis", startCol=16, startRow=1)
      style<-createStyle(fgFill="#E3E9F4", halign="CENTER", border="bottom", textDecoration="Bold")
      addStyle(xlsx, name2[j], style=style, cols=16, rows=1)
      mergeCells(xlsx, name2[j], cols=36:68, rows=1)
      writeData(xlsx, name2[j], "Metabolomics from media", startCol=36, startRow=1)
      style<-createStyle(fgFill="#E2EFDA", halign="CENTER", border="bottom", textDecoration="Bold")
      addStyle(xlsx, name2[j], style=style, cols=36, rows=1)
      freezePane(xlsx, name2[j], firstActiveRow = 3, firstActiveCol = 5, firstRow = FALSE, firstCol = FALSE)
      #Conditional Formatting red for positive green for negative
      posStyle<- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      negStyle<- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      conditionalFormatting(xlsx, name2[j], cols = 16:67, rows = 3:(2+nrow(temp)), rule = "<=-ABS($D$1)", style = posStyle)
      conditionalFormatting(xlsx, name2[j], cols = 16:67, rows = 3:(2+nrow(temp)), rule = ">=ABS($D$1)", style = negStyle)
      rm(temp,temp1,temp2,Correlations)
    }
  }
  saveWorkbook(xlsx, paste0(i,"_",name1[i],"GSEA_Enrichment_Metabolomics.xlsx"),overwrite = TRUE)
}

All_DEGs_log2TPM<-DEGs_log2TPM[rownames(DEGs_log2TPM)%in%All_DEGs,]

#|----PCA----|####
PCA_DEGs<-All_DEGs_log2TPM
Par_Label_Row<-DEGs_Out_TPM$external_gene_name
PCA_DEGres<-prcomp(PCA_DEGs, scale. = TRUE)
PCA_I_PV<-summary(PCA_DEGres)
PCA_index<-get_dist(PCA_DEGs, method = "pearson")       %>%
  hclust( method = "ward.D2", members = NULL) %>%
  cutree(k=2)
svg(filename="Summary_PCA_plot.svg",width = 16, height = 9, pointsize = 12)
gplot<-fviz_pca_ind(PCA_DEGres, geom.ind = "point", pointshape = 21, 
                    pointsize = 1, invisible="quali",
                    fill.ind = paste0("Cluster", PCA_index), 
                    col.ind = "black", ellipse.level=0.95,
                    palette = c("red", "blue"), 
                    addEllipses = TRUE,
                    label = "all", #or "none"
                    col.var = "black",
                    repel = TRUE,
                    legend.title = "") +
  ggtitle(paste0("2D PCA-plot from all ", nrow(PCA_DEGs), " DEGs with triplicates")) +
  xlab(paste0("PC1 (",round(PCA_I_PV$importance[2,1]*100,2), "%)"))+
  ylab(paste0("PC2 (",round(PCA_I_PV$importance[2,2]*100,2), "%)"))+
  theme(
    legend.key.size = unit(1, "line"),
    legend.justification=c(1,0), 
    legend.position=c(0.65, 0.65),  
    legend.text = element_text(hjust = 0.5, color = "black", size = 26, family = "sans"),
    legend.spacing.y=unit(6,"cm"),
    plot.title = element_text(hjust = 0.5, color = "black", size = 26, family = "sans"),
    axis.line = element_line(colour = "black", size = 1.4),
    axis.ticks = element_line(colour = "black", size = 1.4),
    axis.ticks.length = unit(5, "points"),
    axis.title = element_text(color = "black", size = 26, family = "sans"),
    axis.text = element_text(color = "black", size = 26, family = "sans"),
    legend.title = element_text(color = "black", size = 26, family = "sans")
  )
print(gplot)
dev.off()
rm(gplot)

#|----Heatmap----|####
HMP<-as.matrix(All_DEGs_log2TPM)
#Color for column side
Par_CC_P<-rep(c("black","gray","cyan","blue"),each=Par_Rep)
Par_CatLab<-(c("D04", "D09", "D18", "D28")) # Category labels
Par_ColKey<-c("black","gray","cyan","blue") # color key
#Color for Scale bar
Par_Colors=c("blue","black", "red")
Par_Colors=colorRampPalette(Par_Colors)(1000)
png("Summary_HeatMap_Plot.png",  width = 1250, height = 800)
heatmap.2(HMP,col = Par_Colors, ColSideColors = Par_CC_P, 
          scale = "row", trace = "none",
          distfun = function(x) get_dist(x,method = 'spearman'),
          hclustfun = function(x) hclust(x,method = 'complete'),
          labRow = "", #cexRow = 2.8, 
          labCol = "", #cexCol = 2.8,
          key.xlab = "Standardized Log2TPM",
          key.par = list(cex=1.6),
          density.info = "none", 
          margins = c(0, 10)
)
legend(x=0.75, y=1.04, # Position
       legend = Par_CatLab, # Category labels
       col = Par_ColKey, # color key
       lty = 1, # line style
       lwd = 10, # line width
       cex = 2.3,
       bty = "n",
       seg.len = 1,
       x.intersp=0.5,
       y.intersp=0.8,
       xpd = TRUE
)
dev.off()


#|----Example scripts for intergration----|
temp<-DEGs_Summary[,c(1, 3, 38:89)]
temp<-melt(temp, c(1:2))
write.table(temp, file = "Genes_Metabolomics.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(DEGs_Summary$external_gene_name[rownames(DEGs_Summary)%in%All_DEGs] , 
            file = "SigGenes.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
gene.cor.meata<-temp

query<-c("mmu04921","mmu04310","mmu04360","mmu04024")
temp1<-DEGs_ls[[2]]$GSEA_Enrich[[4]]$Sig %>%
  .[.$ID%in%query,]
Cor.FEAGeneMeta<-data.frame()
for (i in 1:length(query)) {
  temp<-temp1[i,'core_enrichment'] %>%
    str_split(.,"/") %>%
    unlist() %>%
    match(as.character(DEGs_Summary$entrezgene_id),.,nomatch = 0)>0
  temp<-DEGs_Summary[temp,c('ensembl_gene_id', 'external_gene_name', 'log2FoldChange_2')]
  temp<-temp[temp$ensembl_gene_id%in%All_DEGs,]
  temp<-data.frame(
    ensembl_gene_id=temp$ensembl_gene_id,
    external_gene_name=temp$external_gene_name,
    variable=paste0(temp1$ID[i],":",temp1$Description[i])%>%
      rep(.,nrow(temp)),
    value=temp$log2FoldChange_2
  )
  Cor.FEAGeneMeta<-rbind(Cor.FEAGeneMeta,temp)
}
temp<-DEGs_Summary[DEGs_Summary$ensembl_gene_id%in%unique(Cor.FEAGeneMeta$ensembl_gene_id),c(1, 3, 38:89)]%>%
  melt(., c(1:2))
temp<-temp[abs(temp$value)>=0.99,]
Cor.FEAGeneMeta<-rbind(temp,Cor.FEAGeneMeta)
write.table(Cor.FEAGeneMeta, file = "D18vsD04_Cor.FEAGeneMeta.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



GeneExpress<-DEGs_Summary[!is.na(DEGs_Summary$entrezgene_id),
                         c('entrezgene_id','log2FoldChange_1', 'log2FoldChange_2', 'log2FoldChange_3')]
Dup_IDs<-unique(GeneExpress$entrezgene_id[duplicated(GeneExpress$entrezgene_id)])
temp1<-GeneExpress[!(GeneExpress$entrezgene_id%in%Dup_IDs),]
temp2<-GeneExpress[GeneExpress$entrezgene_id%in%Dup_IDs,]
for (m in 1:length(Dup_IDs)) {
  temp3<-c(
    entrezgene_id=Dup_IDs[m],
    log(
      colMeans(
        wzy.2power(
          temp2[temp2$entrezgene_id==Dup_IDs[m],2:4]
        )
      ),
      2
    )
  )
  temp1<-rbind(
    temp1,
    temp3
  )
}
rownames(temp1)<-temp1$entrezgene_id
GeneExpress<-temp1[,-1]
pathview(gene.data=GeneExpress,
         pathway.id="mmu04024",
         species="mmu")




temp<-wzy.correlation.FEA.query(i=2,j=4,  
                                name1=c("D09vsD04", "D18vsD04", "D28vsD04", "D18vsD09", "D28vsD09", "D28vsD18"),
                                Metabolic_Ce=Metabolic_Ce,
                                DEGs_ls=DEGs_ls)
temp$ID<-paste0(temp$ID,":",temp$Description)
write.table(temp, file = "D18vsD04_GSEA_KEGG_Metabolomics.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$variable), file = "Metabolomics.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$ID[abs(temp$value)>=0.99]), file = "D18vsD04_GSEA_KEGG_Metabolomics_Source.99.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$variable[abs(temp$value)>=0.99]), file = "D18vsD04_GSEA_KEGG_Metabolomics_Target.99.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

temp<-wzy.correlation.FEA.query(i=2,j=1,  
                                name1=c("D09vsD04", "D18vsD04", "D28vsD04", "D18vsD09", "D28vsD09", "D28vsD18"),
                                Metabolic_Ce=Metabolic_Ce,
                                DEGs_ls=DEGs_ls)
temp$ID<-paste0(temp$ID,":",temp$Description)
write.table(temp, file = "D18vsD04_GSEA_GOBP_Metabolomics.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$ID[abs(temp$value)>=0.99]), file = "D18vsD04_GSEA_GOBP_Metabolomics_Source.99.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$variable[abs(temp$value)>=0.99]), file = "D18vsD04_GSEA_GOBP_Metabolomics_Target.99.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



temp<-wzy.correlation.FEA.query(i=3,j=1,  
                                name1=c("D09vsD04", "D18vsD04", "D28vsD04", "D18vsD09", "D28vsD09", "D28vsD18"),
                                Metabolic_Ce=Metabolic_Ce,
                                DEGs_ls=DEGs_ls)
temp$ID<-paste0(temp$ID,":",temp$Description)
write.table(temp, file = "D28vsD04_GSEA_GOBP_Metabolomics.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$ID[abs(temp$value)>=0.99]), file = "D28vsD04_GSEA_GOBP_Metabolomics_Source.99.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(unique(temp$variable[abs(temp$value)>=0.99]), file = "D28vsD04_GSEA_GOBP_Metabolomics_Target.99.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)






