P0_ParameterSetting<-function(
  Ann_0mart=Ann_0mart,
  #|----Set up parameters----|####
  Par_dist,  #one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
  Par_hclust, #one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
  Par_SigNu,         #minimal number of DEGs for PCA and heatmap plot 
  Par_Test,       #"Wald" or "none" 
  Par_PorFRR,     #"pvalue" or "padj"
  Par_Rep,        #How many replicates for each condition
  Par_ReN,        #Cut-off value of read number
  Par_SaN,        #Cut-off value for expressed samples number
  Par_FaN,        #Factors column numner in design table
  Par_LFC,        #Cut-off value for Log2FC
  Par_FDR,        #Cut-off value for FDR
  Par_FDR4FEA,    #Cut-off (qValue) value for GO enrichment
  Par_FDR4KEG,    #Cut-off (qValue) value for KEGG enrichment
  Par_FDR4FEA_OR, #Cut-off (qValue) value for GO enrichment by Over Representation Test
  Par_FDR4KEG_OR, #Cut-off (qValue) value for KEGG enrichment by Over Representation Test
  #|----Set up Design Table----|####
  Load_ExDesign=read.csv("ExDesign.csv"),
  Par_ConLs_C
){
  if(!file.exists("00.Results")){
    dir.create("00.Results")
  }
  if(!file.exists("01.PlotOutput")){
    dir.create("01.PlotOutput")
  }
  if(!file.exists("02.TableOutput")){
    dir.create("02.TableOutput")
  }
  #|----Set up Biomart parameters----|####
  View(listMarts()) #Check database list
  Ann_0mart <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
  View(listDatasets(mart=Ann_0mart)) #Check speices list
  Ann_0mart <- useDataset("mmusculus_gene_ensembl", Ann_0mart)
  #|========Parameter Setting========|####
  #Loading experiment design table
  Load_ExDesign[]<-lapply(Load_ExDesign,as.character) #Loading experiment design table
  RUV_x<-c()
  for (i in 1:nrow(Load_ExDesign)) {
    RUV_x<-c(RUV_x,paste(Load_ExDesign[i,Par_FaN], collapse = "_"))
  }
  RUV_x<-factor(RUV_x, levels = unique(RUV_x))
  for (i in 1:ncol(Load_ExDesign)) {
    Load_ExDesign[,i]<-factor(Load_ExDesign[,i], levels = unique(Load_ExDesign[,i]))
  }
  rm(i)
  #|----Finishing----|####
  save(Par_ConLs_C,Par_FaN,Par_FDR,Par_FDR4KEG,Par_FDR4FEA,Par_LFC,Ann_0mart,Par_Test,Par_dist,Par_hclust,Par_SigNu,
       Par_PorFRR,Par_ReN,Par_Rep,Par_SaN,Load_ExDesign,RUV_x,Par_FDR4FEA_OR,Par_FDR4KEG_OR,
       file = "00.Results//Parameters.Rdata")
}