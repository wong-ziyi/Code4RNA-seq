#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "DisGenels.txt"
}
library(biomaRt)
library(refGenome)
library(DESeq2)
#|----Parapare Annotation----|####
View(listMarts()) #Check database list
Ann_0mart <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
View(listDatasets(mart=Ann_0mart)) #Check speices list
Ann_0mart <- useDataset("mmusculus_gene_ensembl", Ann_0mart)
#|----Load gtf file
Anno_gtf<-ensemblGenome()
read.gtf(Anno_gtf, args[1]) #"SW3_D04091828_merged.gtf"
Anno_gtf<-Anno_gtf@ev$gtf
genels<-unique(Anno_gtf$ref_gene_id)
#Extract unwanted IG, pseudogene, Mt RNA, 
Anno_biotype<-getBM(
  mart=Ann_0mart,
  attributes=c("ensembl_gene_id","gene_biotype","entrezgene_id"),
  filter="ensembl_gene_id",
  values=genels,
  uniqueRows = TRUE,useCache = TRUE)
table(Anno_biotype$gene_biotype)
DisGenels<-Anno_biotype[
  Anno_biotype$gene_biotype%in%c(
    "IG_C_gene","IG_C_pseudogene","IG_D_gene","IG_D_pseudogene","IG_J_gene","IG_LV_gene","IG_pseudogene","IG_V_gene",
    "IG_V_pseudogene","Mt_rRNA","Mt_tRNA","polymorphic_pseudogene","processed_pseudogene","pseudogene","rRNA",
    "TR_C_gene","TEC","TR_C_gene","TR_D_gene","TR_J_gene","TR_J_pseudogene","TR_V_gene","TR_V_pseudogene",
    "transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene",
    "translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene"
  ),
]
table(DisGenels$gene_biotype)
a<-c()
n<-nrow(DisGenels)%/%500
for (i in 1:n) {
  if(i < n){
    a<-c(
      a,
      paste0(DisGenels$ensembl_gene_id[(((i-1)*500)+1):(i*500)],collapse = "|")
    )
  }else{
    a<-c(
      a,
      paste0(DisGenels$ensembl_gene_id[(((i-1)*500)+1):nrow(DisGenels)],collapse = "|")
    )
  }
}
write.table(x = a, file = args[2], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")