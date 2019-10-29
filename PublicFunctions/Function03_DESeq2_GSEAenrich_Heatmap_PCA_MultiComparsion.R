P3_DESeq2_GSEA_MultiComparsion<-function(){
  if(!file.exists("00.Resource//Annotation.Rdata")){
    stop("no such file: Annotation.Rdata")
  }
  if(!file.exists("00.Results//Parameters.Rdata")){
    stop("no such file: Parameters.Rdata")
  }
  if(!file.exists("00.Results//SalmonRawData.Rdata")){
    stop("no such file: almonRawData.Rdata")
  }
  if(!file.exists("00.Results//DESeq2_Processing.Rdata")){
    stop("no such file: DESeq2_Processing.Rdata")
  }
  #|----DESeq2 Differenate Expressed Genes Visulazing Heatmap Enrichment Analysis----|####
  load("00.Resource//Annotation.Rdata")
  load("00.Results//Parameters.Rdata")
  load("00.Results//DESeq2_Processing.Rdata")
  DEGs_ls<-list()
  for (i in 1:length(Par_ConLs_C)){
    if(Par_Test=="none"){
      Sheet1<-"LRT_Test"
      DEGs_res<-results(DEGs_dds, contrast = Par_ConLs_C[[i]])
    }else if(Par_Test=="Wald"){
      Sheet1<-"Wald_Test"
      DEGs_res<-results(DEGs_dds, test="Wald",contrast = Par_ConLs_C[[i]])
    }else if(Par_Test=="LRT"){
      Sheet1<-"LRT_Test"
      load("00.Results//SalmonRawData.Rdata")
      index<-c(
        which(RUV_x%in%Par_ConLs_C[[i]][2]),
        which(RUV_x%in%Par_ConLs_C[[i]][3])
      )
      temp1<-counts(RUV_set)[,index]
      temp1<-temp1[apply(temp1, 1, function(x) length(x[x>Par_ReN])>=Par_SaN),]
      temp2<-pData(RUV_set)[index,]
      temp2$RUV_x<-factor(temp2$RUV_x,levels = unique(temp2$RUV_x))
      DEGs_dds <- DESeqDataSetFromMatrix(
        countData = temp1, 
        colData = temp2,         #cbind(pData(RUV_set),Load_ExDesign[,Par_FaN])
        design = ~ W_1 + RUV_x           #W_1 + CellType*Loading*Drug
      ) #attr(DEGs_dds, "modelMatrixType") View(attr(DEGs_dds, "modelMatrix"))
      DEGs_dds <- DESeq(
        DEGs_dds, test = "LRT", reduced = as.formula("~ W_1")
        #minReplicatesForReplace=Par_Rep
      )
      DEGs_res<-results(DEGs_dds,contrast = Par_ConLs_C[[i]])
      rm(RUV_set,Txi,temp1,temp2)
    }
    Par_name<-paste(Par_ConLs_C[[i]][2], "vs", Par_ConLs_C[[i]][3])
    dir.create("01.PlotOutput", showWarnings = FALSE)
    #
    #|----PlotMA----|####
    svg(filename = file.path("01.PlotOutput", 
                             paste0("Salmon_DEGs_PlotMA_",sprintf("%02d", i),"_",gsub(" ", "", Par_name),".svg")), 
        width = 10, height = 10, pointsize = 12)
    par(cex=1.5, cex.lab=1.5, cex.axis=1.5)
    options(scipen=5)
    plotMA(DEGs_res,ylim=c(-round(quantile(abs(DEGs_res$log2FoldChange),0.99)),
                           round(quantile(abs(DEGs_res$log2FoldChange),0.99))), alpha=Par_FDR,
           main=Par_name, ylab="Log2 Fold Change", xlab="Mean of Normalized Counts")
    abline(h=c(-Par_LFC,Par_LFC), col="dodgerblue", lwd=2)
    dev.off()
    #
    #|----Annotation----|####
    temp<-merge(Anno,Anno_Summary[,2:8],by="row.names", all=TRUE)
    rownames(temp)<-temp[,1]
    temp<-temp[,-1]
    DEGs_res_out<-as.data.frame(DEGs_res)
    DEGs_res_out<-cbind(baseMean=DEGs_res_out[,1],FoldChange=2^DEGs_res_out$log2FoldChange,DEGs_res_out[,2:6])
    DEGs_res_out<-merge(temp[rownames(temp)%in%rownames(DEGs_res_out),], 
                        DEGs_res_out, by="row.names", all=TRUE)
    rownames(DEGs_res_out)<-DEGs_res_out[,1]
    DEGs_res_out<-DEGs_res_out[,-1]
    DEGs_CutOff<-as.numeric(abs(DEGs_res_out$log2FoldChange)>=Par_LFC & 
                              DEGs_res_out[,Par_PorFRR]<=Par_FDR)#DEGs_CutOff for indentify significance
    names(DEGs_CutOff)<-rownames(DEGs_res_out)
    DEGs_res_out<-cbind(DEGs_res_out, CutOff=DEGs_CutOff)#Log2FC and P value
    #
    #|----Construction----|####
    DEGs_Div<-list(
      DEGs_raw=DEGs_res_out,
      GSEA_Enrich=list(GO_BP=list(),GO_MF=list(),GO_CC=list(),KEGG=list()),
      Over_Represent=list(
        All=list(GO_BP=list(),GO_MF=list(),GO_CC=list(),KEGG=list()),
        Up_regulated=list(GO_BP=list(),GO_MF=list(),GO_CC=list(),KEGG=list()),
        Down_regulated=list(GO_BP=list(),GO_MF=list(),GO_CC=list(),KEGG=list())
      )
    )
    #
    if(sum(DEGs_CutOff,na.rm = TRUE)>=Par_SigNu){
      #|----PCA----|####
      PCA_DEGs<-cbind(
        DEGs_Out_TPM[
          rownames(DEGs_Out_TPM) %in% rownames(DEGs_Div$DEGs_raw[DEGs_Div$DEGs_raw$CutOff==1,]),
          grep(
            paste0("log2TPM_", Par_ConLs_C[[i]][3]),
            colnames(DEGs_Out_TPM), value = TRUE
          )
          ],
        DEGs_Out_TPM[
          rownames(DEGs_Out_TPM) %in% rownames(DEGs_Div$DEGs_raw[DEGs_Div$DEGs_raw$CutOff==1,]),
          grep(
            paste0("log2TPM_", Par_ConLs_C[[i]][2]),
            colnames(DEGs_Out_TPM), value = TRUE
          )
          ]
      )
      Par_Label_Row<-DEGs_Out_TPM$external_gene_name
      PCA_DEGres<-prcomp(PCA_DEGs, scale. = TRUE)
      PCA_I_PV<-summary(PCA_DEGres)
      PCA_index<-get_dist(PCA_DEGs, method = Par_dist)       %>%
        hclust( method = Par_hclust, members = NULL) %>%
        cutree(k=2)
      svg(filename=file.path("01.PlotOutput", 
                             paste0("Clst_PCA_AllDEGs_",sprintf("%02d", i),"_",gsub(" ", "", Par_name),".svg")), 
          width = 16, height = 9, pointsize = 12)
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
      #
      #|----Heatmap----|####
      HMP<-as.matrix(PCA_DEGs)
      #Color for column side
      Par_CC_P<-rep(c("black","gray"),each=Par_Rep)
      Par_CatLab<-c(Par_ConLs_C[[i]][3],Par_ConLs_C[[i]][2]) # Category labels
      Par_ColKey<-c("black","gray") # color key
      #Color for Scale bar
      Par_Colors=c("blue","black", "red")
      Par_Colors=colorRampPalette(Par_Colors)(1000)
      png(
        file.path("01.PlotOutput", 
                  paste0("Clst_Heatmap_AllDEGs_",sprintf("%02d", i),"_",gsub(" ", "", Par_name),".png")),  
        width = 1250, height = 800
      )
      heatmap.2(HMP,col = Par_Colors, ColSideColors = Par_CC_P, 
                scale = "row", trace = "none",
                distfun = function(x) get_dist(x,method = Par_dist),
                hclustfun = function(x) hclust(x,method = Par_hclust),
                labRow = "", #cexRow = 2.8, 
                labCol = "", #cexCol = 2.8,
                key.xlab = "Standardized Log2TPM",
                key.par = list(cex=1.6),
                density.info = "none", 
                margins = c(0, 10)
      )
      legend(x=0.72, y=1.04, # Position
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
      #
    }
    #|----Functional Enrichment----|####
    #|----Prepare folder----|####
    FolderName<-paste0("FEA_Plot_",sprintf("%02d", i),"_",gsub(" ", "", Par_name))
    dir.create(file.path("01.PlotOutput",FolderName), showWarnings = FALSE)
    if(file.exists(file.path("01.PlotOutput",FolderName,"GSEA_Enrich"))){
      unlink(file.path("01.PlotOutput",FolderName,"GSEA_Enrich"), recursive=TRUE)
    }
    dir.create(file.path("01.PlotOutput",FolderName,"GSEA_Enrich"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"GSEA_Enrich","GO_BP"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"GSEA_Enrich","GO_MF"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"GSEA_Enrich","GO_CC"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"GSEA_Enrich","KEGG"),showWarnings = FALSE)
    if(file.exists(file.path("01.PlotOutput",FolderName,"Over_Represent"))){
      unlink(file.path("01.PlotOutput",FolderName,"Over_Represent"), recursive=TRUE)
    }
    dir.create(file.path("01.PlotOutput",FolderName,"Over_Represent"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"Over_Represent","All_KEGG"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"Over_Represent","Up_regulated_KEGG"),showWarnings = FALSE)
    dir.create(file.path("01.PlotOutput",FolderName,"Over_Represent","Down_regulated_KEGG"),showWarnings = FALSE)
    #
    #|----Prepare data for GSEA enrichment----|####
    FEA_GeneLs<-DEGs_Div$DEGs_raw[!is.na(DEGs_Div$DEGs_raw$entrezgene_id),]
    FEA_GeneLsFC<-FEA_GeneLs$log2FoldChange[!is.na(FEA_GeneLs$entrezgene_id)]
    names(FEA_GeneLsFC)<-FEA_GeneLs$entrezgene_id
    #|----get log2 mean fold change for duplicated IDs
    Dup_IDs<-unique(names(FEA_GeneLsFC)[duplicated(names(FEA_GeneLsFC))])
    temp1<-FEA_GeneLsFC[!(names(FEA_GeneLsFC)%in%Dup_IDs)]
    temp2<-FEA_GeneLsFC[names(FEA_GeneLsFC)%in%Dup_IDs]
    for (m in 1:length(Dup_IDs)) {
      temp3<-log(mean(wzy.2power(temp2[names(temp2)==Dup_IDs[m]])),2)
      names(temp2)<-Dup_IDs[m]
      temp1<-c(
        temp1,
        temp3
      )
    }
    #FEA_GeneLsFC_temp<-Anno2$entrezgene_id[!(Anno2$entrezgene_id%in%names(temp1))]
    #FEA_GeneLsFC_temp2<-rep(0,length(FEA_GeneLsFC_temp))
    #names(FEA_GeneLsFC_temp2)<-FEA_GeneLsFC_temp
    #FEA_GeneLsFC_temp2<-FEA_GeneLsFC_temp2[!duplicated(names(FEA_GeneLsFC_temp2))]
    #FEA_GeneLsFC<-sort(c(temp1,FEA_GeneLsFC_temp2), decreasing = TRUE)
    FEA_GeneLsFC<-sort(temp1, decreasing = TRUE)
    #
    #|----Gene Ontology (GO) enrichment (GSEA enrichment)----|####
    for (j in c("BP","MF","CC")) {
      FEA_GO <- gseGO(geneList     = FEA_GeneLsFC,
                      OrgDb        = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont          = j,
                      nPerm        = 1000,
                      minGSSize    = 20,
                      maxGSSize    = 2000,
                      pvalueCutoff = 0.1,
                      verbose      = FALSE)
      if(is.null(FEA_GO)){next()}
      temp1<-as.data.frame(FEA_GO@result)
      if(sum(is.na(temp1$qvalues))>0){
        FEA_GO_MostDescend<-temp1$ID[temp1$p.adjust<Par_FDR4FEA]
      }else{
        FEA_GO_MostDescend<-temp1$ID[temp1$qvalues<Par_FDR4FEA]
      }
      #
      #|----Get most descend significant term----|####
      if(length(FEA_GO_MostDescend)>0){
        FEA_GO_MostDescend<-GOMostDescend(x=FEA_GO_MostDescend)
        FEA_GO_MostDescend<-GOMostDescend2(x=FEA_GO_MostDescend)
        FEA_GO_MostDescend<-temp1[temp1$ID %in% FEA_GO_MostDescend, ]
        rownames(FEA_GO_MostDescend)<-FEA_GO_MostDescend$ID
        DEGs_Div$GSEA_Enrich[[match(j,c("BP","MF","CC"))]]<-list(Raw=FEA_GO,Sig=FEA_GO_MostDescend)
      }else{
        DEGs_Div$GSEA_Enrich[[match(j,c("BP","MF","CC"))]]<-list(Raw=list(),Sig=list())
      }
      #
      #|----Plot Out----|####
      if(nrow(as.data.frame(FEA_GO_MostDescend))>0){
        svg(
          file.path("01.PlotOutput",FolderName,"GSEA_Enrich",paste0("gseGO_",j,".svg")), 
          width = 9, height = 6, pointsize = 12
        )
        gplot<-FEA_GO_MostDescend                                                 %>% 
          top_n(30, wt=abs(.$NES))                                                %>%
          mutate(GO_Term=apply(.[,1:2], MARGIN = 1, FUN = paste, collapse = "_")) %>% #Create new variable
          mutate(
            Tags=strsplit(.$leading_edge, ",")%>%
              lapply(., function(x){as.numeric(strsplit(x,"[=\\%]")[[1]][2])})%>%
              unlist()
          )%>%
          ggplot(aes(x=NES, y=GO_Term, colour=pvalue, size=Tags)) +
          geom_point() +
          expand_limits(x=0) +
          labs(x="Normalized Enrichment Score (NES)", y=paste0("GO term (",j,")"), 
               title="The most descend significant GO term",
               colour="P value", size="Tags (%)")
        print(gplot)
        dev.off()
        temp<-top_n(FEA_GO_MostDescend,30, wt=abs(FEA_GO_MostDescend$NES))[,1:2]
        for (k in 1:nrow(temp)) {
          svg(
            filename = file.path(
              "01.PlotOutput", 
              FolderName,"GSEA_Enrich",
              paste0("GO_",j),
              paste0(substr(make.names(paste0("GSEA_",paste0(temp[k,1:2],collapse = "_"))),0,65),".svg")
            ), 
            width = 9, height = 6, pointsize = 12
          )
          gplot<-gseaplot(FEA_GO, geneSetID=temp[k,1],title = paste0(temp[k,1:2],collapse = " "))
          print(gplot)
          dev.off()
        }
        rm(gplot,temp)
      }
    }
    temp<-as.data.frame(DEGs_Div$GSEA_Enrich$GO_BP$Sig)
    temp<-temp[grepl("[Pp]athway",temp$Description),]
    if(nrow(as.data.frame(temp))>0){
      svg(
        file.path("01.PlotOutput",FolderName,"GSEA_Enrich",paste0("gseGO_goPathwayTerm.svg")), 
        width = 9, height = 6, pointsize = 12
      )
      gplot<-temp                                                               %>% 
        mutate(GO_Term=apply(.[,1:2], MARGIN = 1, FUN = paste, collapse = "_")) %>% #Create new variable 
        ggplot(aes(x=NES, y=GO_Term, colour=pvalue, size=setSize)) +
        geom_point() +
        expand_limits(x=0) +
        labs(x="Normalized Enrichment Score (NES)", y=paste0("GO term (",j,")"), 
             title='GO term with "pathway" keyword',
             colour="P value", size="Tags (%)")
      print(gplot)
      dev.off()
    }
    #|----KEGG pathway enrichment (GSEA enrichment)----|####
    FEA_KEGG <- gseKEGG(geneList     = FEA_GeneLsFC,
                        organism     = 'mmu',
                        nPerm        = 1000,
                        minGSSize    = 120,
                        pvalueCutoff = 0.1,
                        verbose      = TRUE)
    if(is.null(FEA_KEGG)){
      DEGs_Div$GSEA_Enrich$KEGG<-list(Raw=list(),Sig=list())
    }else{
      temp1<-as.data.frame(FEA_GO@result)
      if(sum(is.na(temp1$qvalues))>0){
        FEA_GO_MostDescend<-temp1$ID[temp1$p.adjust<Par_FDR4FEA]
      }else{
        FEA_GO_MostDescend<-temp1$ID[temp1$qvalues<Par_FDR4FEA]
      }
      temp<-as.data.frame(FEA_KEGG@result)
      temp<-temp[temp$qvalues<=Par_FDR4KEG,]
      if(nrow(temp)>0){
        wd<-getwd()
        setwd(file.path("01.PlotOutput",FolderName,"GSEA_Enrich","KEGG"))
        if(length(temp$ID)!=0){
          for (m in 1:length(temp$ID)) {
            tryCatch({
              pathview(gene.data=FEA_GeneLsFC,
                       pathway.id=temp$ID[m],#[!(temp[43:331]%in%c("mmu05206","mmu00511"))], 
                       species="mmu")
            }, error=function(e){})
          }
        }
        setwd(wd)
        for (k in 1:nrow(temp)) {
          svg(
            file.path(
              "01.PlotOutput",FolderName,"GSEA_Enrich","KEGG",
              make.names(paste0(paste0(temp$ID[k],collapse = "_"),".svg"))
            ), 
            width = 9, height = 6, pointsize = 12
          )
          gplot<-gseaplot(FEA_KEGG, geneSetID=temp[k,1],title = paste0(temp[k,1:2],collapse = " "))
          print(gplot)
          dev.off()
        }
        DEGs_Div$GSEA_Enrich$KEGG<-list(Raw=FEA_KEGG,Sig=temp)
        svg(
          file.path("01.PlotOutput",FolderName,"GSEA_Enrich",paste0("gseKEGG_Pathway.svg")), 
          width = 9, height = 6, pointsize = 12
        )
        gplot<-temp                                                               %>% 
          top_n(30, wt=abs(.$NES))                                                %>%
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
      }else{
        DEGs_Div$GSEA_Enrich$KEGG<-list(Raw=list(),Sig=list())
      }
    }
    #|----Over Representation Test----|####
    for (l in c("All","Up_regulated","Down_regulated")) {
      if(l=="All"){
        gene<-FEA_GeneLs$entrezgene_id[FEA_GeneLs$CutOff==1]
        gene<-unique(gene)
      }else if(l=="Up_regulated"){
        gene<-FEA_GeneLs$entrezgene_id[FEA_GeneLs$CutOff==1 & FEA_GeneLs$log2FoldChange>0]
        gene<-unique(gene)
      }else if(l=="Down_regulated"){
        gene<-FEA_GeneLs$entrezgene_id[FEA_GeneLs$CutOff==1 & FEA_GeneLs$log2FoldChange<0]
        gene<-unique(gene)
      }
      #|----Gene Ontology (GO) enrichment (Over Representation Test)----|####
      for (j in c("BP","MF","CC")) {
        FEA_GO_OR <-enrichGO(gene          = gene,
                             OrgDb         = org.Mm.eg.db,
                             keyType       = "ENTREZID",
                             ont           = j,
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
        if(is.null(FEA_GO_OR)){next()}
        temp2<-as.data.frame(FEA_GO_OR@result)
        if(sum(is.na(temp2$qvalue))>0){
          FEA_GO_OR_MostDescend<-temp2$ID[temp2$p.adjust<Par_FDR4FEA_OR]
        }else{
          FEA_GO_OR_MostDescend<-temp2$ID[temp2$qvalue<Par_FDR4FEA_OR]
        }
        #
        #|----Get most descend significant term----|####
        if(length(FEA_GO_OR_MostDescend)>0){
          FEA_GO_OR_MostDescend<-GOMostDescend(x=FEA_GO_OR_MostDescend)
          FEA_GO_OR_MostDescend<-GOMostDescend2(x=FEA_GO_OR_MostDescend)
          FEA_GO_OR_MostDescend<-temp2[temp2$ID %in% FEA_GO_OR_MostDescend, ]
          rownames(FEA_GO_OR_MostDescend)<-FEA_GO_OR_MostDescend$ID
          DEGs_Div$Over_Represent[[l]][[match(j,c("BP","MF","CC"))]]<-list(Raw=FEA_GO_OR,Sig=FEA_GO_OR_MostDescend)
        }else{
          DEGs_Div$Over_Represent[[l]][[match(j,c("BP","MF","CC"))]]<-list(Raw=list(),Sig=list())
        }
        #
        #|----Plot Out----|####
        if(nrow(as.data.frame(FEA_GO_OR_MostDescend))>0){
          svg(
            file.path("01.PlotOutput",FolderName,"Over_Represent",paste0(l,"_GO_",j,".svg")), 
            width = 9, height = 6, pointsize = 12
          )
          temp<-FEA_GO_OR_MostDescend                         %>% 
            top_n(30, wt=-FEA_GO_OR_MostDescend$pvalue)
          x<-c()
          for (k in 1:nrow(temp)) {
            x<-c(x,100*temp$Count[k]/length(FEA_GO_OR@geneSets[[temp$ID[k]]]))
          }
          temp<-cbind(temp,Hits=x)
          gplot<-temp                                                               %>%
            mutate(GO_Term=apply(.[,1:2], MARGIN = 1, FUN = paste, collapse = "_")) %>% #Create new variable
            ggplot(aes(x=Hits, y=GO_Term, colour=pvalue, size=Count)) +
            geom_point() +
            expand_limits(x=0) +
            labs(x="Hits (%)", y=paste0("GO term (",j,")"), 
                 title="The most descend significant GO term",
                 colour="P value", size="Gene Number")
          print(gplot)
          dev.off()
        }
      }
      temp<-as.data.frame(DEGs_Div$Over_Represent[[l]]$GO_BP$Sig)
      temp<-temp[grepl("[Pp]athway",temp$Description),]
      if(nrow(as.data.frame(temp))>0){
        svg(
          file.path("01.PlotOutput",FolderName,"Over_Represent",paste0(l,"_GO_goPathwayTerm.svg")), 
          width = 9, height = 6, pointsize = 12
        )
        x<-c()
        for (k in 1:nrow(temp)) {
          x<-c(x,100*temp$Count[k]/length(DEGs_Div$Over_Represent[[l]]$GO_BP$Raw@geneSets[[temp$ID[k]]]))
        }
        temp<-cbind(temp,Hits=x)
        gplot<-temp                                                               %>%
          mutate(GO_Term=apply(.[,1:2], MARGIN = 1, FUN = paste, collapse = "_")) %>% #Create new variable
          ggplot(aes(x=Hits, y=GO_Term, colour=pvalue, size=Count)) +
          geom_point() +
          expand_limits(x=0) +
          labs(x="Hits (%)", y=paste0("GO term (",j,")"), 
               title='GO term with "pathway" keyword',
               colour="P value", size="Gene Number")
        print(gplot)
        dev.off()
      }
      #|----KEGG pathway enrichment (Over Representation Test)----|####
      FEA_KEGG_OR <- enrichKEGG(gene      = gene,
                                organism  = 'mmu',
                                pvalueCutoff = Par_FDR4KEG_OR)
      if(is.null(FEA_KEGG_OR)){next()}
      temp<-as.data.frame(FEA_KEGG_OR@result)
      if(sum(is.na(temp$qvalue))>0){
        temp<-temp[temp$p.adjust<Par_FDR4KEG_OR,]
      }else{
        temp<-temp[temp$qvalue<Par_FDR4KEG_OR,]
      }
      if(nrow(temp)>0){
        wd<-getwd()
        setwd(file.path("01.PlotOutput",FolderName,"Over_Represent",paste0(l,"_KEGG")))
        if(length(rownames(temp))!=0){
          for (m in 1:length(temp$ID)) {
            tryCatch({
              pathview(gene.data=FEA_GeneLsFC,
                       pathway.id=temp$ID[m],#[!(temp[43:331]%in%c("mmu05206","mmu00511"))], 
                       species="mmu")
            }, error=function(e){})
          }
        }
        setwd(wd)
        DEGs_Div$Over_Represent[[l]]$KEGG<-list(Raw=FEA_KEGG_OR,Sig=temp)
        svg(
          file.path("01.PlotOutput",FolderName,"Over_Represent",paste0(l,"_KEGG_Pathway_.svg")), 
          width = 9, height = 6, pointsize = 12
        )
        temp<-top_n(temp,30, wt=-temp$pvalue)
        x<-c()
        for (k in 1:nrow(temp)) {
          x<-c(x,100*temp$Count[k]/length(FEA_KEGG_OR@geneSets[[temp$ID[k]]]))
        }
        temp<-cbind(temp,Hits=x)
        gplot<-temp                                                               %>%
          mutate(GO_Term=apply(.[,1:2], MARGIN = 1, FUN = paste, collapse = "_")) %>% #Create new variable
          ggplot(aes(x=Hits, y=GO_Term, colour=pvalue, size=Count)) +
          geom_point() +
          expand_limits(x=0) +
          labs(x="Hits (%)", y="KEGG ID", 
               title="Over represented KEGG pathway",
               colour="P value", size="Gene Number")
        print(gplot)
        dev.off()
      }else{
        DEGs_Div$Over_Represent$KEGG<-list(Raw=list(),Sig=list())
      }
    }
    #|----Make Excel .xlsx format output----|####
    #|----All DEGs output----|####
    xlsx<-createWorkbook()
    addWorksheet(xlsx,Sheet1)
    writeDataTable(xlsx,Sheet1,DEGs_res_out,startCol = 1,startRow = 2)
    setColWidths(xlsx,Sheet1,cols=1,widths=20,hidden=rep(FALSE,length(cols)),ignoreMergedCells=FALSE)
    setColWidths(xlsx,Sheet1,cols=4,widths=14,hidden=rep(FALSE,length(cols)),ignoreMergedCells=FALSE)
    mergeCells(xlsx,Sheet1, cols=1:11, rows=1)
    mergeCells(xlsx,Sheet1, cols=12:19, rows=1)
    writeData(xlsx, Sheet1, "Basic Information", startCol=1, startRow=1)
    style<-createStyle(fgFill="#ABBEDE", halign="LEFT", border="bottom", textDecoration="Bold")
    addStyle(xlsx, Sheet1, style=style, cols=1, rows=1)
    writeData(xlsx, Sheet1, "DESeq2 Analysis", startCol=12, startRow=1)
    style<-createStyle(fgFill="#C7D3E9", halign="CENTER", border="bottom", textDecoration="Bold")
    addStyle(xlsx, Sheet1, style=style, cols=12, rows=1)
    freezePane(xlsx, Sheet1, firstActiveRow = 3, firstActiveCol = 5, firstRow = FALSE, firstCol = FALSE)
    saveWorkbook(
      xlsx, 
      file = file.path(
        "02.TableOutput", 
        paste0("Salmon_DEGs_ResOut_",sprintf("%02d", i),"_",gsub(" ", "", Par_name),".xlsx")
      ),overwrite = TRUE
    )
    rm(xlsx)
    #|----All enrichment resuts output (GSEA)----|####
    xlsx<-createWorkbook()
    for (j in 1:4) {
      if(nrow(as.data.frame(DEGs_Div$GSEA_Enrich[[j]]$Sig))>0){
        addWorksheet(xlsx,c("GO_BP","GO_MF","GO_CC","KEGG")[j])
        temp<-DEGs_Div$GSEA_Enrich[[j]]$Raw@geneSets
        temp2<-DEGs_Div$GSEA_Enrich[[j]]$Sig
        temp3<-data.frame()
        for (k in 1:nrow(temp2)) {
          ID<-temp[[rownames(temp2)[k]]]
          temp4<-MeanRetrieve(ID,Anno2,DEGs_log2TPM,RUV_x,Par_Rep)
          temp3<-rbind(temp3,temp4)
          if(k==1){colnames(temp3)<-names(temp4)}
        }
        DEGs_Div$GSEA_Enrich[[j]]$Sig<-cbind(temp2,temp3)
        writeDataTable(xlsx,c("GO_BP","GO_MF","GO_CC","KEGG")[j],DEGs_Div$GSEA_Enrich[[j]]$Sig,startCol = 1,startRow = 2)
      }
    }
    if(length(names(xlsx))>0){
      for (k in 1:length(names(xlsx))) {
        setColWidths(xlsx, sheet=k, cols=1, widths = 12, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
        setColWidths(xlsx, sheet=k, cols=2, widths = 45, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
        freezePane(xlsx, sheet=k, firstActiveRow = 3, firstActiveCol = 3, firstRow = FALSE, firstCol = FALSE)
        mergeCells(xlsx, sheet=k, cols=1:11, rows=1)
        writeData(xlsx, sheet=k, "GSEA Enrichment Analysis", startCol=1, startRow=1)
        style<-createStyle(fgFill="#ABBEDE", halign="LEFT", border="bottom", textDecoration="Bold")
        addStyle(xlsx, sheet=k, style=style, cols=1, rows=1)
        mergeCells(xlsx, sheet=k, cols=12:(11+length(levels(RUV_x))), rows=1)
        writeData(xlsx, sheet=k, "Total Mean of genes of each condition", startCol=12, startRow=1)
        style<-createStyle(fgFill="#C7D3E9", halign="CENTER", border="bottom", textDecoration="Bold")
        addStyle(xlsx, sheet=k, style=style, cols=12, rows=1)
        mergeCells(xlsx, sheet=k, cols=(12+length(levels(RUV_x))):(11+length(levels(RUV_x))+length(RUV_x)), rows=1)
        writeData(xlsx, sheet=k, "Total Mean of genes of each sample", startCol=(12+length(levels(RUV_x))), startRow=1)
        style<-createStyle(fgFill="#E3E9F4", halign="CENTER", border="bottom", textDecoration="Bold")
        addStyle(xlsx, sheet=k, style=style, cols=12+length(levels(RUV_x)), rows=1)
      }
    }else{
      addWorksheet(xlsx,"NULL")
    }
    saveWorkbook(
      xlsx, 
      file.path(
        "02.TableOutput",
        paste0("GSEA_EnrichmentResults_",sprintf("%02d", i),"_",gsub(" ", "", Par_name),".xlsx")
      ),overwrite = TRUE
    )
    rm(xlsx)
    #
    #|----All enrichment resuts output (Over Representation Test)----|####
    for (l in c("All","Up_regulated","Down_regulated")) {
      xlsx<-createWorkbook()
      for (j in 1:4) {
        if(nrow(as.data.frame(DEGs_Div$Over_Represent[[l]][[j]]$Sig))>0){
          addWorksheet(xlsx,c("GO_BP","GO_MF","GO_CC","KEGG")[j])
          temp<-DEGs_Div$Over_Represent[[l]][[j]]$Raw@geneSets
          temp2<-DEGs_Div$Over_Represent[[l]][[j]]$Sig
          temp3<-data.frame()
          for (k in 1:nrow(temp2)) {
            ID<-temp[[rownames(temp2)[k]]]
            temp4<-MeanRetrieve(ID,Anno2,DEGs_log2TPM,RUV_x,Par_Rep)
            temp3<-rbind(temp3,temp4)
            if(k==1){colnames(temp3)<-names(temp4)}
          }
          DEGs_Div$Over_Represent[[l]][[j]]$Sig<-cbind(temp2,temp3)
          writeDataTable(xlsx,c("GO_BP","GO_MF","GO_CC","KEGG")[j],
                         DEGs_Div$Over_Represent[[l]][[j]]$Sig,startCol = 1,startRow = 2)
        }
      }
      if(length(names(xlsx))>0){
        for (k in 1:length(names(xlsx))) {
          setColWidths(xlsx, sheet=k, cols=1, widths = 12, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
          setColWidths(xlsx, sheet=k, cols=2, widths = 45, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
          freezePane(xlsx, sheet=k, firstActiveRow = 3, firstActiveCol = 3, firstRow = FALSE, firstCol = FALSE)
          mergeCells(xlsx, sheet=k, cols=1:9, rows=1)
          writeData(xlsx, sheet=k, "GSEA Enrichment Analysis", startCol=1, startRow=1)
          style<-createStyle(fgFill="#ABBEDE", halign="LEFT", border="bottom", textDecoration="Bold")
          addStyle(xlsx, sheet=k, style=style, cols=1, rows=1)
          mergeCells(xlsx, sheet=k, cols=10:(9+length(levels(RUV_x))), rows=1)
          writeData(xlsx, sheet=k, "Total Mean of genes of each condition", startCol=10, startRow=1)
          style<-createStyle(fgFill="#C7D3E9", halign="CENTER", border="bottom", textDecoration="Bold")
          addStyle(xlsx, sheet=k, style=style, cols=10, rows=1)
          mergeCells(xlsx, sheet=k, cols=(10+length(levels(RUV_x))):(9+length(levels(RUV_x))+length(RUV_x)), rows=1)
          writeData(xlsx, sheet=k, "Total Mean of genes of each sample", startCol=(10+length(levels(RUV_x))), startRow=1)
          style<-createStyle(fgFill="#E3E9F4", halign="CENTER", border="bottom", textDecoration="Bold")
          addStyle(xlsx, sheet=k, style=style, cols=10+length(levels(RUV_x)), rows=1)
        }
      }else{
        addWorksheet(xlsx,"NULL")
      }
      saveWorkbook(
        xlsx, 
        file.path(
          "02.TableOutput",
          paste0("OverRepresentTest_",sprintf("%02d", i),"_",gsub(" ", "", Par_name),"_",l,".xlsx")
        ),overwrite = TRUE
      )
      rm(xlsx)
    }
    #
    #|----Output----|####
    DEGs_ls[[i]]<-DEGs_Div
    names(DEGs_ls)[i]<-gsub(" ", "", Par_name)
  }
  rm(i)
  save(DEGs_ls, file = "00.Results//DESeq2_Results.Rdata")
}