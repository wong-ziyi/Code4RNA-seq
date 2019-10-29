P2_DESeq2_PreProcessing<-function(){
  if(!file.exists("00.Resource//Annotation.Rdata")){
    stop("no such file: Annotation.Rdata")
  }
  if(!file.exists("00.Results//Parameters.Rdata")){
    stop("no such file: Parameters.Rdata")
  }
  if(!file.exists("00.Results//SalmonRawData.Rdata")){
    stop("no such file: almonRawData.Rdata")
  }
  #|----DESeq2 linear model Pre-possesing----|####
  #|----DESeq2 linear model----|####
  load("00.Results//Parameters.Rdata")
  load("00.Results//SalmonRawData.Rdata")
  load("00.Resource//Annotation.Rdata")
  DEGs_dds <- DESeqDataSetFromMatrix(
    countData = counts(RUV_set), 
    colData = pData(RUV_set),         #cbind(pData(RUV_set),Load_ExDesign[,Par_FaN])
    design = ~ W_1 + RUV_x           #W_1 + CellType*Loading*Drug
  ) #attr(DEGs_dds, "modelMatrixType") View(attr(DEGs_dds, "modelMatrix"))
  if(Par_Test=="none"){
    DEGs_dds <- DESeq(
      DEGs_dds, test = "LRT", reduced = as.formula("~ W_1"), betaPrior=FALSE
      #minReplicatesForReplace=Par_Rep
    )
  }else if(Par_Test=="Wald"){
    DEGs_dds <- DESeq(
      DEGs_dds, test = "LRT", reduced = as.formula("~ W_1"), betaPrior=FALSE
      #minReplicatesForReplace=Par_Rep
    )
  }else if(Par_Test=="LRT"){
    DEGs_dds <- DESeq(
      DEGs_dds, test = "LRT", reduced = as.formula("~ W_1"), betaPrior=FALSE
      #minReplicatesForReplace=Par_Rep
    )
  }
  #
  #|----Getting Normalized Read Counts (NorCounts)----|####
  temp<-as.data.frame(counts(DEGs_dds, normalize = TRUE))
  DEGs_NorCounts<-merge(
    temp,
    Txi$length[rownames(Txi$length)%in%rownames(temp),],
    by='row.names', all=TRUE
  )
  rownames(DEGs_NorCounts)<-DEGs_NorCounts[,1]
  DEGs_NorCounts<-DEGs_NorCounts[,-1]
  colnames(DEGs_NorCounts)<-c(paste("NorCounts", RUV_x, c(1:Par_Rep), sep="_"), 
                              paste("EffLength",RUV_x, c(1:Par_Rep), sep="_"))
  #
  #|----Transcript Million Mapped reads (log2TPM)----|####
  DEGs_log2TPM<-DEGs_NorCounts[,1:(ncol(DEGs_NorCounts)/2)]/
    DEGs_NorCounts[,(ncol(DEGs_NorCounts)/2+1):ncol(DEGs_NorCounts)]
  Par_a<-(colSums(DEGs_log2TPM)/1000000)
  Par_cl<-makeCluster(detectCores())
  registerDoSNOW(Par_cl)
  Par_pb <- tkProgressBar("Parallel task", min=1, max=length(rownames(DEGs_log2TPM)))
  progress <- function(n) setTkProgressBar(Par_pb, n)
  Par_opts<-list(progress = progress)
  temp<-foreach(i=1:nrow(DEGs_log2TPM), .combine=rbind, 
                .options.snow=Par_opts, .packages = "tcltk") %dopar% {
                  log2(DEGs_log2TPM[i,1:(ncol(DEGs_NorCounts)/2)]/Par_a)
                }
  close(Par_pb)
  stopCluster(Par_cl)
  DEGs_log2TPM<-temp
  rm(temp, Par_a, Par_cl, Par_opts, Par_pb)
  #|____File with log2TPM only____|####
  DEGs_log2TPM[DEGs_log2TPM=="-Inf"]<-0
  colnames(DEGs_log2TPM)<-paste("log2TPM", RUV_x, c(1:Par_Rep), sep="_")
  DEGs_Out_TPM<-merge(DEGs_log2TPM[rownames(DEGs_log2TPM)%in%rownames(DEGs_NorCounts),], 
                      DEGs_NorCounts, by="row.names", all=TRUE)
  rownames(DEGs_Out_TPM)<-DEGs_Out_TPM[,1]
  DEGs_Out_TPM<-DEGs_Out_TPM[,-1]
  DEGs_Out_TPM<-merge(Anno[rownames(Anno)%in%rownames(DEGs_Out_TPM),], 
                      DEGs_Out_TPM, by="row.names", all=TRUE)
  rownames(DEGs_Out_TPM)<-DEGs_Out_TPM[,1]
  DEGs_Out_TPM<-DEGs_Out_TPM[,-1]
  #
  #|----Xlsx output----|####
  n1<-length(RUV_x)
  n2<-length(levels(RUV_x))
  xlsx<-createWorkbook()
  addWorksheet(xlsx,"ReadsCount")
  writeDataTable(xlsx,"ReadsCount",DEGs_Out_TPM[,1:4],startCol = 1,startRow = 2)
  writeData(xlsx,"ReadsCount",DEGs_Out_TPM[,-(1:4)],startCol = 5,startRow = 2)
  setColWidths(xlsx, "ReadsCount", cols=1, widths = 20, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
  setColWidths(xlsx, "ReadsCount", cols=4, widths = 14, hidden = rep(FALSE,length(cols)), ignoreMergedCells = FALSE)
  mergeCells(xlsx,"ReadsCount", cols=1:4, rows=1)
  writeData(xlsx, "ReadsCount", "Basic Information", startCol=1, startRow=1)
  style<-createStyle(fgFill="#ABBEDE", halign="LEFT", border="bottom", textDecoration="Bold")
  addStyle(xlsx, "ReadsCount", style=style, cols=1, rows=1)
  mergeCells(xlsx,"ReadsCount", cols=5:(4+n1), rows=1)
  writeData(xlsx, "ReadsCount", "Log2TPM", startCol=5, startRow=1)
  style<-createStyle(fgFill="#C7D3E9", halign="CENTER", border="bottom", textDecoration="Bold")
  addStyle(xlsx, "ReadsCount", style=style, cols=5, rows=1)
  mergeCells(xlsx,"ReadsCount", cols=(5+n1):(4+2*n1), rows=1)
  writeData(xlsx, "ReadsCount", "Normalized Raw Reads Count", startCol=(5+n1), startRow=1)
  style<-createStyle(fgFill="#E3E9F4", halign="CENTER", border="bottom", textDecoration="Bold")
  addStyle(xlsx, "ReadsCount", style=style, cols=(5+n1), rows=1)
  mergeCells(xlsx,"ReadsCount", cols=(5+2*n1):(4+3*n1), rows=1)
  writeData(xlsx, "ReadsCount", "Effective gene length from Salmon software", startCol=(5+2*n1), startRow=1)
  style<-createStyle(fgFill="#E9DDC7", halign="CENTER", border="bottom", textDecoration="Bold")
  addStyle(xlsx, "ReadsCount", style=style, cols=(5+2*n1), rows=1)
  freezePane(xlsx, "ReadsCount", firstActiveRow = 3, firstActiveCol = 5, firstRow = FALSE, firstCol = FALSE)
  saveWorkbook(
    xlsx, 
    file = file.path(
      "02.TableOutput", 
      "00.DEGs_Out_TPM.xlsx"
    ),overwrite = TRUE
  )
  #|----Finishing----|####
  save(DEGs_dds,DEGs_log2TPM,DEGs_NorCounts,DEGs_Out_TPM,file = "00.Results\\DESeq2_Processing.Rdata")
}