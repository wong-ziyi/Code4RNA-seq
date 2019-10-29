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
options(url.method="libcurl")
invisible(PackageCheckInstall())
wd<-getwd()
#||||||||||||||||||||||(Salmon) Identify DEs by DESeq2 Packages||||||||||||||||||||||####
#Par_ConLs_N<-cbind(Load_ExDesign,model.matrix(~ CellType*Loading*Drug, Load_ExDesign)) #Numerical model matrices
P0_ParameterSetting(
  #|----Set up parameters----|####
  Par_dist="pearson",  #one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
  Par_hclust= "ward.D2", #one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
  Par_SigNu=5,         #minimal number of DEGs for PCA and heatmap plot
  Par_Test="LRT",      #"LRT", "Wald", or "none"  
  Par_PorFRR="padj",  #"pvalue" or "padj"
  Par_Rep=3,            #How many replicates for each condition
  Par_ReN=10,            #Cut-off value of read number
  Par_SaN=4,            #Cut-off value for expressed samples number
  Par_FaN=3:4,          #Factors column numner in design table
  Par_LFC=0.6,          #Cut-off value for Log2FC
  Par_FDR=0.05,         #Cut-off value for FDR
  Par_FDR4FEA=0.05,     #Cut-off value for GO enrichment
  Par_FDR4KEG=0.05,     #Cut-off value for KEGG enrichment
  Par_FDR4FEA_OR=0.05,  #Cut-off (qValue) value for GO enrichment by Over Representation Test
  Par_FDR4KEG_OR=0.05,  #Cut-off (qValue) value for KEGG enrichment by Over Representation Test
  #|----Set up Design Table----|####
  Load_ExDesign=read.csv("ExDesign.csv"),
  Par_ConLs_C=list(
    c("RUV_x","Mloaded_Ctrl","Ctrl_Ctrl"),
    c("RUV_x","Ctrl_BAIBA","Ctrl_Ctrl"),
    c("RUV_x","Mloaded_BAIBA","Ctrl_Ctrl"),
    c("RUV_x","Mloaded_BAIBA","Mloaded_Ctrl"),
    c("RUV_x","Ctrl_BAIBA","Mloaded_Ctrl"),
    c("RUV_x","Mloaded_BAIBA","Ctrl_BAIBA")
  ) #Combined factors list
)
P1_Salmon_Tximport_RUVSeq()
dev.off()
P2_DESeq2_PreProcessing()
P3_DESeq2_GSEA_MultiComparsion()
setwd(wd)
#|||||||||||||||||Further|||||||||||||####
load("00.Results//DESeq2_Results.Rdata")
load("00.Resource//Annotation.Rdata")
load("00.Results//Parameters.Rdata")
load("00.Results//DESeq2_Processing.Rdata")

Baiba<-DEGs_ls$Ctrl_BAIBAvsCtrl_Ctrl$GSEA_Enrich$KEGG$Sig%>%
  .[.$ID=="mmu04310", "core_enrichment"]%>%
  strsplit(.,"/")%>%
  unlist()
Mload<-DEGs_ls$Mloaded_CtrlvsCtrl_Ctrl$GSEA_Enrich$KEGG$Sig%>%
  .[.$ID=="mmu04310", "core_enrichment"]%>%
  strsplit(.,"/")%>%
  unlist()
Combine<-DEGs_ls$Mloaded_BAIBAvsMloaded_Ctrl$GSEA_Enrich$KEGG$Sig%>%
  .[.$ID=="mmu04310", "core_enrichment"]%>%
  strsplit(.,"/")%>%
  unlist()
BaibaFC<-DEGs_ls$Ctrl_BAIBAvsCtrl_Ctrl$DEGs_raw%>%
  .[.$entrezgene_id%in%Combine,c(2,14)]
MloadFC<-DEGs_ls$Mloaded_CtrlvsCtrl_Ctrl$DEGs_raw%>%
  .[.$entrezgene_id%in%Combine,c(2,14)]
Combine<-Reduce(intersect,list(BaibaFC$entrezgene_id,MloadFC$entrezgene_id))
BaibaFC<-BaibaFC[BaibaFC$entrezgene_id%in%Combine,]
BaibaFC<-BaibaFC[match(Combine,BaibaFC$entrezgene_id),]
MloadFC<-MloadFC[MloadFC$entrezgene_id%in%Combine,]
MloadFC<-MloadFC[match(Combine,MloadFC$entrezgene_id),]

FC<-c()
for (i in 1:nrow(MloadFC)) {
  if(BaibaFC$log2FoldChange[i]>0 & MloadFC$log2FoldChange[i]>0){
    temp<-0
  }else if(BaibaFC$log2FoldChange[i]<0 & MloadFC$log2FoldChange[i]<0){
    temp<-0
  }else if(BaibaFC$log2FoldChange[i]>0 & MloadFC$log2FoldChange[i]<0){
    temp<-2
  }else if(BaibaFC$log2FoldChange[i]<0 & MloadFC$log2FoldChange[i]>0){
    temp<--2
  }
  FC<-c(FC,temp)
}
names(FC)<-MloadFC$entrezgene_id

x<-MloadFC$log2FoldChange
names(x)<-MloadFC$entrezgene_id
y<-BaibaFC$log2FoldChange
names(y)<-BaibaFC$entrezgene_id

ID<-DEGs_ls$Ctrl_BAIBAvsCtrl_Ctrl$GSEA_Enrich$KEGG$Raw@geneSets[["mmu04310"]]


Baiba<-DEGs_ls$Ctrl_BAIBAvsCtrl_Ctrl$DEGs_raw%>%
  .[.$entrezgene_id%in%ID,c(2,14)]

Mload<-DEGs_ls$Mloaded_CtrlvsCtrl_Ctrl$DEGs_raw%>%
  .[.$entrezgene_id%in%ID,c(2,14)]%>%
  .[.$entrezgene_id%in%Baiba$entrezgene_id,]
  

z<-cbind(Baiba$log2FoldChange,Mload$log2FoldChange)
rownames(z)<-Baiba$entrezgene_id

a<-pathview(gene.data=z, cpd.data = z,
         pathway.id="mmu04310", 
         species="mmu")


Load_fl[match(Load_ExDesign[,1],Load_fl[,1]),]


load("00.Results//Parameters.Rdata")


Mload_Ctrl<-wzy.gseaKEG.query(a=1,DEGs_ls = DEGs_ls)
BAIBA_Ctrl<-wzy.gseaKEG.query(a=2,DEGs_ls = DEGs_ls)
Mload_BAIBA_Ctrl<-wzy.gseaKEG.query(a=3,DEGs_ls = DEGs_ls)
Mload_BAIBA_Mload<-wzy.gseaKEG.query(a=4,DEGs_ls = DEGs_ls)

KEGG_ls<-c("mmu04020",	"mmu04024",	"mmu04022",	"mmu05418",	"mmu04510",	"mmu04390",
           "mmu04910",	"mmu04010",	"mmu04150",	"mmu04151",	"mmu04015",	"mmu04014",	
           "mmu04810",	"mmu03013",	"mmu03040",	"mmu04310")

temp<-Mload_Ctrl@result
temp<-temp[temp$qvalues<=Par_FDR4KEG,]
svg(
  "Mload_Ctr_gseKEGG_Pathway.svg", 
  width = 9, height = 6, pointsize = 12
)
gplot<-temp                                                               %>% 
  mutate(GO_Term=apply(.[,1:2], MARGIN = 1, FUN = paste, collapse = "_")) %>% #Create new variable
  mutate(
    Tags=strsplit(.$leading_edge, ",")%>%
      lapply(., function(x){as.numeric(strsplit(x,"[=\\%]")[[1]][2])})%>%
      unlist()
  )%>%
  ggplot(aes(x=NES, y=GO_Term, colour=pvalue, size=Tags)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Normalized Enrichment Score (NES)", y="KEGG ID", 
       title="KEGG pathway enrichment by GSEA",
       colour="P value", size="Tags (%)")
print(gplot)
dev.off()

temp1<-DEGs_ls$Mloaded_CtrlvsCtrl_Ctrl$DEGs_raw
temp2<-DEGs_ls$Ctrl_BAIBAvsCtrl_Ctrl$DEGs_raw
temp3<-DEGs_ls$Mloaded_BAIBAvsCtrl_Ctrl$DEGs_raw
temp4<-DEGs_ls$Mloaded_BAIBAvsMloaded_Ctrl$DEGs_raw

DEGs_Summary<-merge(
  temp1[rownames(temp1)%in%rownames(temp2),c('entrezgene_id','log2FoldChange') ],
  temp2[rownames(temp2)%in%rownames(temp1),c('ensembl_gene_id','log2FoldChange') ],
  by='row.names', all=TRUE
)
rownames(DEGs_Summary)<-DEGs_Summary[,1]
DEGs_Summary<-DEGs_Summary[,-1]
DEGs_Summary<-merge(
  DEGs_Summary[rownames(DEGs_Summary)%in%rownames(temp3), ],
  temp3[rownames(temp3)%in%rownames(DEGs_Summary),c('ensembl_gene_id','log2FoldChange') ],
  by='row.names', all=TRUE
)
rownames(DEGs_Summary)<-DEGs_Summary[,1]
DEGs_Summary<-DEGs_Summary[,-1]
DEGs_Summary<-merge(
  DEGs_Summary[rownames(DEGs_Summary)%in%rownames(temp4), ],
  temp4[rownames(temp4)%in%rownames(DEGs_Summary),c('ensembl_gene_id','log2FoldChange') ],
  by='row.names', all=TRUE
)
rownames(DEGs_Summary)<-DEGs_Summary[,1]
DEGs_Summary<-DEGs_Summary[,-1]
DEGs_Summary<-DEGs_Summary[,-c(3,5,7)]
colnames(DEGs_Summary)<-c('entrezgene_id','log2FoldChange_1', 'log2FoldChange_2', 'log2FoldChange_3','log2FoldChange_4')


GeneExpress<-DEGs_Summary[!is.na(DEGs_Summary$entrezgene_id), ]
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
setwd(file.path("03.Query", "All_genes"))
for (i in 1:length(KEGG_ls)) {
  pathview(gene.data=GeneExpress,
           pathway.id=KEGG_ls[i],
           species="mmu")
}
setwd(wd)

Mload_Ctrl@result<-Mload_Ctrl@result[Mload_Ctrl@result$ID%in%KEGG_ls,]
nrow(Mload_Ctrl@result)
BAIBA_Ctrl@result<-BAIBA_Ctrl@result[BAIBA_Ctrl@result$ID%in%KEGG_ls,]
nrow(BAIBA_Ctrl@result)
Mload_BAIBA_Ctrl@result<-Mload_BAIBA_Ctrl@result[Mload_BAIBA_Ctrl@result$ID%in%KEGG_ls,]
nrow(Mload_BAIBA_Ctrl@result)
Mload_BAIBA_Mload@result<-Mload_BAIBA_Mload@result[Mload_BAIBA_Mload@result$ID%in%KEGG_ls,]
nrow(Mload_BAIBA_Mload@result)

NES_Out<-data.frame()
core<-c()
CytoScape<-data.frame()
for (i in 1:length(KEGG_ls)) {
  temp1<-Mload_Ctrl@result[Mload_Ctrl@result$ID==KEGG_ls[i],]
  temp2<-BAIBA_Ctrl@result[BAIBA_Ctrl@result$ID==KEGG_ls[i],]
  temp3<-Mload_BAIBA_Ctrl@result[Mload_BAIBA_Ctrl@result$ID==KEGG_ls[i],]
  temp4<-Mload_BAIBA_Mload@result[Mload_BAIBA_Mload@result$ID==KEGG_ls[i],]
  temp0<-data.frame(paste0(temp1$ID,":",temp1$Description),temp1$NES,temp2$NES,temp3$NES,temp4$NES)
  colnames(temp0)<-c("ID","Mload vs Ctr","BAIBA vs Ctr","Mload+BAIBA vs Ctr","Mload+BAIBA vs Mload")
  NES_Out<-rbind(NES_Out,temp0)
  colnames(NES_Out)<-c("ID","Mload vs Ctr","BAIBA vs Ctr","Mload+BAIBA vs Ctr","Mload+BAIBA vs Mload")
  
  svg(
    file.path("03.Query", make.names(paste0(NES_Out$ID[i],"Mload_Ctrl.svg"))), 
    width = 9, height = 6, pointsize = 12
  )
  gplot<-gseaplot(Mload_Ctrl, geneSetID=temp1$ID[1],title = NES_Out$ID[i])
  print(gplot)
  dev.off()
  
  svg(
    file.path("03.Query", make.names(paste0(NES_Out$ID[i],"BAIBA_Ctrl.svg"))), 
    width = 9, height = 6, pointsize = 12
  )
  gplot<-gseaplot(BAIBA_Ctrl, geneSetID=temp1$ID[1],title = NES_Out$ID[i])
  print(gplot)
  dev.off()
  
  svg(
    file.path("03.Query", make.names(paste0(NES_Out$ID[i],"Mload_BAIBA_Ctrl.svg"))), 
    width = 9, height = 6, pointsize = 12
  )
  gplot<-gseaplot(Mload_BAIBA_Ctrl, geneSetID=temp1$ID[1],title = NES_Out$ID[i])
  print(gplot)
  dev.off()
  
  svg(
    file.path("03.Query", make.names(paste0(NES_Out$ID[i],"Mload_BAIBA_Mload.svg"))), 
    width = 9, height = 6, pointsize = 12
  )
  gplot<-gseaplot(Mload_BAIBA_Mload, geneSetID=temp1$ID[1],title = NES_Out$ID[i])
  print(gplot)
  dev.off()
  
  temp1.core<-temp1$core_enrichment %>%
    str_split(.,"/") %>%
    unlist()
  temp2.core<-temp2$core_enrichment %>%
    str_split(.,"/") %>%
    unlist()
  temp3.core<-temp3$core_enrichment %>%
    str_split(.,"/") %>%
    unlist()
  temp4.core<-temp4$core_enrichment %>%
    str_split(.,"/") %>%
    unlist()
  
  core.temp<-temp1.core[temp1.core%in%temp2.core & temp1.core%in%temp3.core & temp1.core%in%temp4.core]
  core<-c(core,core.temp)
  
  CytoScape.temp<-GeneExpress[rownames(GeneExpress)%in%core.temp,]
  
  setwd(file.path("03.Query", "Core_Shared_genes"))
  pathview(gene.data=CytoScape.temp,
           pathway.id=KEGG_ls[i],
           species="mmu")
  setwd(wd)
  
  CytoScape.temp<-data.frame(
    Source=rownames(CytoScape.temp),
    Targe=rep(x=as.character(NES_Out$ID[i]), times=nrow(CytoScape.temp)),
    Value=rep(1, times=nrow(CytoScape.temp)),
    log2FoldChange_1=CytoScape.temp$log2FoldChange_1,
    log2FoldChange_2=CytoScape.temp$log2FoldChange_2,
    log2FoldChange_3=CytoScape.temp$log2FoldChange_3,
    log2FoldChange_4=CytoScape.temp$log2FoldChange_4
  )
  CytoScape<-rbind(CytoScape,CytoScape.temp)
}
CytoScape<-merge(CytoScape, FEA_GeneLs[FEA_GeneLs$entrezgene_id%in%CytoScape$Source,2:3], 
                 by.x='Source',by.y='entrezgene_id',
                 all=TRUE)
write.table(CytoScape, "CytoScape_Alberto.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(CytoScape$external_gene_name[duplicated(CytoScape$external_gene_name)], 
            "CytoScape_Alberto_shared.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)




CytoScape.2<-data.frame()
for (i in 1:length(KEGG_ls)) {
  print(i)
  setwd(file.path("03.Query", "Core_Genes"))
  a<-pathview(gene.data=GeneExpress[rownames(GeneExpress)%in%unique(CytoScape$Source),],
              pathway.id=KEGG_ls[i],
              species="mmu")
  setwd(wd)
  
  
  setwd(file.path("03.Query", "Shared_genes"))
  pathview(gene.data=GeneExpress[rownames(GeneExpress)%in%unique(CytoScape$Source[duplicated(CytoScape$Source)]),],
           pathway.id=KEGG_ls[i],expand.node=TRUE,
           species="mmu")
  setwd(wd)
  
  b<-unlist(a$plot.data.gene$kegg.names[!is.na(a$plot.data.gene$log2FoldChange_1)])
  
  CytoScape.temp<-GeneExpress[rownames(GeneExpress)%in%b,]
  CytoScape.temp<-data.frame(
    Source=rownames(CytoScape.temp),
    Targe=rep(x=paste0(KEGG_ls[i],":", FEA_KEGG@result$Description[FEA_KEGG@result$ID==KEGG_ls[i]]), 
              times=nrow(CytoScape.temp)),
    Value=rep(1, times=nrow(CytoScape.temp)),
    log2FoldChange_1=CytoScape.temp$log2FoldChange_1,
    log2FoldChange_2=CytoScape.temp$log2FoldChange_2,
    log2FoldChange_3=CytoScape.temp$log2FoldChange_3,
    log2FoldChange_4=CytoScape.temp$log2FoldChange_4
  )
  CytoScape.2<-rbind(CytoScape.2,CytoScape.temp)
  
}

CytoScape.2<-merge(CytoScape.2, FEA_GeneLs[FEA_GeneLs$entrezgene_id%in%CytoScape.2$Source,2:3], 
                 by.x='Source',by.y='entrezgene_id',
                 all=TRUE)
write.table(CytoScape.2, "CytoScape.2_Alberto.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(CytoScape.2$external_gene_name[duplicated(CytoScape.2$external_gene_name)], 
            "CytoScape.2_Alberto_shared.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)


GeneExpress.all<-merge(DEGs_Summary[rownames(DEGs_Summary)%in%rownames(DEGs_log2TPM),2:1],
                       DEGs_log2TPM[rownames(DEGs_log2TPM)%in%rownames(DEGs_Summary),], by='row.names', all=TRUE)
GeneExpress.all<-GeneExpress.all[,-2]
GeneExpress.all<-GeneExpress.all[!is.na(GeneExpress.all$entrezgene_id), ]
Dup_IDs<-unique(GeneExpress.all$entrezgene_id[duplicated(GeneExpress.all$entrezgene_id)])
temp1<-GeneExpress.all[!(GeneExpress.all$entrezgene_id%in%Dup_IDs),]
temp2<-GeneExpress.all[GeneExpress.all$entrezgene_id%in%Dup_IDs,]
for (m in 1:length(Dup_IDs)) {
  temp3<-c(
    Row.names=temp2$Row.names[temp2$entrezgene_id==Dup_IDs[m]],
    entrezgene_id=Dup_IDs[m],
    log(
      colMeans(
        wzy.2power(
          temp2[temp2$entrezgene_id==Dup_IDs[m],3:14]
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
GeneExpress.all<-temp1[,-c(1:2)]


#|----PCA----|####
PCA_DEGs<-as.matrix(GeneExpress.all[rownames(GeneExpress.all)%in%unique(CytoScape$Source),])
#PCA_DEGs<-as.matrix(GeneExpress.all[rownames(GeneExpress.all)%in%CytoScape$Source[duplicated(CytoScape$Source)],])
PCA_DEGs<-apply(X = PCA_DEGs,MARGIN = 2, FUN = as.numeric)
Par_Label_Row<-DEGs_Out_TPM$external_gene_name
PCA_DEGres<-prcomp(PCA_DEGs, scale. = TRUE)
PCA_I_PV<-summary(PCA_DEGres)
PCA_index<-get_dist(PCA_DEGs, method = "euclidean")       %>%
  hclust( method = "ward.D", members = NULL) %>%
  cutree(k=3)
svg(filename="Summary_PCA_plot.svg",width = 16, height = 9, pointsize = 12)
gplot<-fviz_pca_ind(PCA_DEGres, geom.ind = "point", pointshape = 21, 
                    pointsize = 3, invisible="quali",
                    col.ind = "black",
                    fill.ind = paste0("Cluster", PCA_index), #
                    llipse.level=0.95, #
                    palette = c("red", "green","yellow"),# 
                    addEllipses = TRUE, #
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
HMP<-PCA_DEGs
#Color for column side
Par_CC_P<-rep(c("black","gray","cyan","blue"),each=Par_Rep)
Par_CatLab<-c("Ctr","Mload","BAIBA","Mload+BAIBA") # Category labels
Par_ColKey<-c("black","gray","cyan","blue") # color key
#Color for Scale bar
Par_Colors=c("blue","black", "red")
Par_Colors=colorRampPalette(Par_Colors)(1000)
png("Summary_HeatMap_Plot.png",  width = 1250, height = 800)
heatmap.2(HMP,col = Par_Colors, ColSideColors = Par_CC_P, 
          #Colv=FALSE,
          #dendrogram="row",
          scale = "row", trace = "none",
          distfun = function(x) get_dist(x,method = 'euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D'),
          labRow = "", #cexRow = 2.8, 
          labCol = "", #cexCol = 2.8,
          key.xlab = "Standardized Log2TPM",
          key.par = list(cex=1.6),
          density.info = "none", 
          margins = c(0, 10)
)
legend(x=0.8, y=1.1, # Position
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
