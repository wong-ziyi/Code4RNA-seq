#|----Query Pathview----|####
pathview.query.wzy<-function(Path_name=Path_name,Path_ID=Path_ID,DEGs_ls=DEGs_ls,kegg.sets.mm=kegg.sets.mm){
  if(!file.exists("04.PathwayQuery")){
    dir.create("04.PathwayQuery")
  }
  wd<-getwd()
  if(!file.exists(
    file.path("04.PathwayQuery", Path_name)
  )){
    dir.create(file.path("04.PathwayQuery", Path_name))
  }
  if(!file.exists(
    file.path("04.PathwayQuery", Path_name,"All")
  )){
    dir.create(file.path("04.PathwayQuery", Path_name,"All"))
  }
  if(!file.exists(
    file.path("04.PathwayQuery", Path_name,"Sig")
  )){
    dir.create(file.path("04.PathwayQuery", Path_name,"Sig"))
  }
  setwd(file.path("04.PathwayQuery", Path_name))
  temp<-DEGs_ls[[Path_name]][["DEGs_raw"]]
  if(length(grep(Path_ID, names(kegg.sets.mm), value = TRUE))==0){
    print("No such ID in current local database")
    stop()
  }
  write.csv(
    temp[as.character(temp$entrezgene_id)%in%kegg.sets.mm[grep(Path_ID,names(kegg.sets.mm),value = TRUE)][[1]],],
    file = paste0(Path_ID,".csv"),
    row.names = FALSE
  )
  setwd(wd)
  setwd(file.path("04.PathwayQuery", Path_name,"All"))
  temp<-DEGs_ls[[Path_name]][["DEGs_raw"]]
  Path_GeneLsFC<-temp$log2FoldChange
  names(Path_GeneLsFC)<-temp$entrezgene_id
  pathview(gene.data=Path_GeneLsFC, #kegg.dir = "Data", 
           pathway.id=Path_ID,
           species="mmu")
  setwd(wd)
  setwd(file.path("04.PathwayQuery", Path_name,"Sig"))
  temp<-DEGs_ls[[Path_name]][["All_DEGs"]][[1]]
  Path_GeneLsFC<-temp$log2FoldChange
  names(Path_GeneLsFC)<-temp$entrezgene_id
  pathview(gene.data=Path_GeneLsFC, #kegg.dir = "Data", 
           pathway.id=Path_ID, 
           species="mmu")
  setwd(wd)
  return()
}
wzy.gseaKEG.query<-function(a=1,DEGs_ls=DEGs_ls){
  #|----Prepare data for GSEA enrichment----|####
  FEA_GeneLs<-DEGs_ls[[a]]$DEGs_raw[!is.na(DEGs_ls[[a]]$DEGs_raw$entrezgene_id),]
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
  FEA_GeneLsFC<-sort(temp1, decreasing = TRUE)
  FEA_KEGG <- gseKEGG(geneList     = FEA_GeneLsFC,
                      organism     = 'mmu',
                      nPerm        = 1000,
                      minGSSize    = 120,
                      pvalueCutoff = 1,
                      verbose      = TRUE)
}