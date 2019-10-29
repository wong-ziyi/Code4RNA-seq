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
