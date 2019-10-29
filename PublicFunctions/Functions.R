PackageCheckInstall<-function(){
  list.of.packages <- c("locfit","hwriter","coda","Hmisc","jsonlite","askpass","pracma","stringi",
                        "plotrix","pkgconfig","backports","memoise","cli","ashr","numDeriv",
                        "gplots", "plyr","scales","zeallot","RMariaDB","ggplot2","nVennR",
                        "RColorBrewer","data.table", "RSpectra","MatrixCorrelation","resample", 
                        "processx","pheatmap","xml2", "devtools","resample", "foreach",
                        "iterators","tcltk","parallel","doParallel", "stringr","doSNOW",
                        "dplyr","readr","factoextra","cluster","pvclust","openxlsx")
  list.of.packages.Bio <- c("DOSE","clusterProfiler","geneplotter","genefilter","SummarizedExperiment","Rsamtools",
                            "Biostrings","AnnotationDbi","GenomicRanges","GenomeInfoDb","S4Vectors","BiocGenerics",
                            "doSNOW","mygene","apeglm","RUVSeq","EnhancedVolcano","GO.db", "goseq", 
                            "biomaRt", "tximport","rtracklayer","geneLenDataBase", "ensembldb","gage",
                            "gageData","pathview","org.Mm.eg.db","DESeq2")
  #Packages that does not install yet
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  new.packages.Bio <- list.of.packages.Bio[!(list.of.packages.Bio %in% installed.packages()[, "Package"])]
  #install required packages
  if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
  if(length(new.packages.Bio)) {BiocManager::install(new.packages.Bio)}
  if(length(new.packages)) {install.packages(new.packages)}
  devtools::install_github("grimbough/biomaRt")
  #Loading all packages
  lapply(c(list.of.packages,list.of.packages.Bio), library, character.only = TRUE)
}
Rtrend<-function(x){
  if(x>0){
    return(1)
  }else if(x<0){
    return(-1)
  }else if(x==0){
    return(0)
  }
}
pattern<-function(x){
  for (l in 1:length(x)) {
    x[l]<-Rtrend(x[l])
  }
  return(x)
}
GOMostDescend2<-function(x){
  Par_GOCC <- as.list(GOCCOFFSPRING)
  Par_GOBP <- as.list(GOBPOFFSPRING)
  Par_GOMF <- as.list(GOMFOFFSPRING)
  GOlist<-c(Par_GOCC, Par_GOBP, Par_GOMF)
  for (i in length(x):1) {
    if(sum(x %in% GOlist[[x[i]]])!=0){
      x<-x[-i]
    }
  }
  return(x)
}
GOMostDescend<-function(x){
  Par_GOCC <- as.list(GOCCANCESTOR)
  Par_GOBP <- as.list(GOBPANCESTOR)
  Par_GOMF <- as.list(GOMFANCESTOR)
  GOlist<-c(Par_GOCC, Par_GOBP, Par_GOMF)
  for (i in 1:length(x)) {
    x<-x[!(x %in% GOlist[[x[i]]])]
  }
  return(x)
}
GOMostAscend<-function(x, GOlist){
  while (i<=length(x)) {
    x<-x[!(x %in% GOlist[[x[i]]])]
    i<-i+1
  }
  return(x)
}
lookup.wzy<-function(a,b, cola){
  temp<-c()
  for (i in 1:nrow(a)) {
    temp<-rbind(temp, b[rownames(b)==as.character(a[i, cola]),])
  }
  a<-cbind(a,temp)
  return(a)
}
Df4Plot<-function(Gene,Type="log2TPM",RUV_x=RUV_x,DEGs_Out_TPM=DEGs_Out_TPM,colname="external_gene_name"){
  temp<-as.numeric(DEGs_Out_TPM[DEGs_Out_TPM[,colname]==Gene,
                           c(grep(Type, colnames(DEGs_Out_TPM), value=TRUE))])
  temp<-data.frame(as.character(RUV_x), as.numeric(temp))
  colnames(temp)<-c("x",Type)
  temp<-as.data.frame(temp)
  temp$x<-factor(temp$x,unique(temp$x))
  return(temp)
}
BinMean <- function (vec, every, na.rm = FALSE) {
  n <- length(vec)
  x <- .colMeans(vec, every, n %/% every, na.rm)
  r <- n %% every
  if (r) x <- c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
  x
}
FUN<-function(x){
  if(is.null(x)){
    return(NA)
  }else{
    substr(paste(apply(x, MARGIN = 1, FUN = paste, collapse = "_"),collapse = " "),0,32700)
  }
}
wzy.2power<-function(x){
  if(is.null(dim(x))){
    x[x==0]<--Inf
  }else{
    x<-apply(x, 2, function(x){x[x==0]<--Inf;return(x)})
  }
  x<-2^x
  return(x)
}
ping <- function(x, stderr = FALSE, stdout = FALSE, ...){
  pingvec <- system2("ping", x,
                     stderr = FALSE,
                     stdout = FALSE,...)
  if (pingvec == 0) TRUE else FALSE
}
MeanRetrieve<-function(ID,Anno2,DEGs_log2TPM,RUV_x,Par_Rep){
  ID<-rownames(Anno2)[Anno2$entrezgene_id%in%ID]
  temp<-DEGs_log2TPM        %>%
    apply(., 2, wzy.2power) %>%
    .[rownames(.)%in%ID,]   %>%
    colMeans(.)             %>%
    c(BinMean(.,Par_Rep),.)    %>%
    log(.,2)
  names(temp)<-ifelse(names(temp)=="",paste0(levels(RUV_x),"_Log2MeanTPM"),names(temp))
  return(temp)
}
