P1_Salmon_Tximport_RUVSeq<-function(){
  if(!file.exists("00.Results//Parameters.Rdata")){
    stop("no such file: Parameters.Rdata")
  }
  #|========Load file & Remove unwanted varations within groups========|####
  load("00.Results//Parameters.Rdata")
  #|========Load .sf Salmon file========|####
  Load_fl<-list.dirs(path = ".", full.names = TRUE, recursive = TRUE) %>%
    list.files(pattern = "\\.sf$",full.names = TRUE)
  Load_fl<-strsplit(Load_fl, "/")                                     %>%
    lapply(FUN = nth, n=-2)                                    %>%
    as.data.frame()                                            %>%
    t()                                                        %>%
    cbind(Load_fl)
  rownames(Load_fl)<-Load_fl[,1]
  colnames(Load_fl)<-c("Name", "Path")
  sum(!(Load_ExDesign[,1]%in%Load_fl[,1]))==0
  Load_fl<-as.data.frame(Load_fl[match(Load_ExDesign[,1],Load_fl[,1]),])
  #
  #|========Anno_Make Tx2Gene file========|####
  if(!file.exists("00.Resource//Annotation.Rdata")){
    #|----Load first .sf to get gene list----|####
    Load_Tx2Gene<-Load_fl$Path[1]              %>%
      as.character()                           %>%
      read.table(header = TRUE, row.names = 1) %>%
      rownames()                               %>%
      as.character()
    #|----Build up Tx2Gene file----|####
    #|----Annotation retrieve setting----|####
    Anno_tx2gene <- getBM(
      mart=Ann_0mart,
      attributes=c("ensembl_transcript_id_version","ensembl_gene_id","entrezgene_id","gene_biotype",
                   "external_gene_name","transcript_biotype","description"),
      filter="ensembl_transcript_id_version",
      values=Load_Tx2Gene,
      uniqueRows = FALSE)
    #
    #|----make annotation file----|####
    Anno<-Anno_tx2gene[!duplicated(Anno_tx2gene$ensembl_gene_id),]
    rownames(Anno)<-Anno$ensembl_gene_id
    Anno<-data.frame(
      ensembl_gene_id=Anno$ensembl_gene_id,
      entrezgene_id=Anno$entrezgene_id,
      external_gene_name=Anno$external_gene_name,
      gene_biotype=Anno$gene_biotype,
      row.names = rownames(Anno)
    )
    Anno2<-Anno[!is.na(Anno$entrezgene_id),] #For EntrezGeneID
    Anno2<-Anno2[!duplicated(Anno2[,'ensembl_gene_id']),] #For EntrezGeneID
    rownames(Anno2)<-Anno2$ensembl_gene_id
    #|----get summary----|####
    Anno_Summary<-getGenes(geneid = rownames(Anno) , fields = c("symbol","summary","generif"))
    temp<-lapply(Anno_Summary$generif,FUN = FUN)
    temp<-unlist(temp)
    Anno_Summary<-data.frame(
      ensembl_gene_id=Anno_Summary@listData$query,
      summary=Anno_Summary@listData$summary,
      generif=temp
    )
    Anno_Summary<-Anno_Summary[!duplicated(Anno_Summary$ensembl_gene_id),]
    rownames(Anno_Summary)<-Anno_Summary$ensembl_gene_id
    temp2<-data.frame()
    for (i in 1:length(Anno_Summary$ensembl_gene_id)) {
      temp2<-suppressWarnings(
        rbind(
          temp2,
          cbind(
            ensembl_gene_id=Anno_Summary$ensembl_gene_id[i],
            fread(
              paste0(
                "https://www.uniprot.org/uniprot/?query=",
                Anno_Summary$ensembl_gene_id[i],
                "+organism:10090&format=tab&limit=1&columns=id,comment(FUNCTION),go,database(KEGG),comment(PATHWAY)"
              ),verbose=FALSE,showProgress=FALSE
            )
          ),
          fill=TRUE
        )
      )
      cat(i, paste0("of ",length(Anno_Summary$ensembl_gene_id),"\r"))
      flush.console()
    }
    rm(i)
    temp2<-as.data.frame(temp2)
    rownames(temp2)<-temp2$ensembl_gene_id
    Anno_Summary<-merge(Anno_Summary[rownames(Anno_Summary)%in%rownames(temp2),],temp2[,2:6],by="row.names",all=TRUE)
    Anno_Summary<-data.frame(
      ensembl_gene_id=Anno_Summary$ensembl_gene_id,
      UniProtKBID=Anno_Summary$Entry,
      UniProtFunction=Anno_Summary$`Function [CC]`,
      UniProtPathway=Anno_Summary$Pathway,
      RefSeqSummary=Anno_Summary$summary,
      KEGG=Anno_Summary$`Cross-reference (KEGG)`,
      GO=Anno_Summary$`Gene ontology (GO)`,
      GeneRif=Anno_Summary$generif,
      row.names = Anno_Summary$Row.names
    )
    rm(temp,temp2)
    #
    #|----Finishing----|####
    pattern<-'(\\\\r)|(\\\\n)|(\\\\a)|(\\\\f)|(\\\\v)|(\\\\b)|(\\\\s)|(\\\\t)|(\\\\)'
    Anno_Summary$GeneRif<-gsub(pattern,'', Anno_Summary$GeneRif)
    Anno_Summary$UniProtFunction<-gsub(pattern,'', Anno_Summary$UniProtFunction)
    Anno_Summary$UniProtPathway<-gsub(pattern,'', Anno_Summary$UniProtPathway)
    Anno_Summary$RefSeqSummary<-gsub(pattern,'', Anno_Summary$RefSeqSummary)
    Anno_Summary$GO<-gsub(pattern,'', Anno_Summary$GO)
    Anno_Summary$KEGG<-gsub(pattern,'', Anno_Summary$KEGG)
    Anno_Summary$GeneRif<-substr(gsub('"',"'", Anno_Summary$GeneRif),0,32700)
    dir.create("00.Resource", showWarnings = FALSE)
    save(Anno_tx2gene, Anno, Anno2, Anno_Summary, file = "00.Resource//Annotation.Rdata")
    load("00.Resource//Annotation.Rdata")
    write.xlsx(
      list(
        "Annotation"=Anno_Summary,"ConvertID"=Anno_tx2gene
      ),file = "Annotation.xlsx",asTable=FALSE,keepNA =FALSE
    )
    rm(Anno_tx2gene, Anno, Anno2, Anno_Summary)
  }
  #
  #|========load Salmon quant.sf file========|####
  load("00.Resource//Annotation.Rdata")
  Txi<-as.character(Load_fl$Path) %>%
    tximport(type = "salmon", tx2gene = Anno_tx2gene[,1:2])
  mode(Txi$counts)<-"integer"
  colnames(Txi$counts)<-Load_fl$Name
  head(Txi$counts,1)
  #
  #|========Remove unwanted variation by RUVSeq method (DOI:10.1038/nbt.2931)========|####
  Par_colors <- brewer.pal(ncol(Txi$counts)/Par_Rep, "Dark2")
  #|----Filter out non-expressed genes, by requiring more than 5 reads in at least two samples for each gene----|####
  RUV_filter <- apply(Txi$counts, 1, function(x) length(x[x>Par_ReN])>=Par_SaN)
  RUV_filtered <- Txi$counts[RUV_filter,]
  RUV_genes <- rownames(RUV_filtered)
  #|----Construct a matrix specifying the replicates----|####
  RUV_difference <- makeGroups(RUV_x)
  RUV_difference
  #|----Store the data in an object of S4 class----|####
  RUV_set <- newSeqExpressionSet(as.matrix(RUV_filtered),
                                 phenoData = data.frame(RUV_x, row.names = colnames(RUV_filtered)))
  RUV_set
  colMeans(counts(RUV_set))
  sd(colMeans(counts(RUV_set)))
  plotRLE(RUV_set, outline=FALSE, ylim=c(-4, 4), col=Par_colors[RUV_x])
  plotPCA(RUV_set, col=Par_colors[RUV_x], cex=1.2)
  #|----Use upper-quartile (UQ) normalization----|####
  RUV_set <- betweenLaneNormalization(RUV_set, which = "upper")
  plotRLE(RUV_set, outline=FALSE, ylim=c(-4, 4), col=Par_colors[RUV_x])
  plotPCA(RUV_set, col=Par_colors[RUV_x], cex=1.2)
  #|----Use all the genes as control genes for the estimation of the factors of unwanted variation----|####
  RUV_set <- RUVs(RUV_set, RUV_genes, k=1, RUV_difference)
  pData(RUV_set)
  colMeans(counts(RUV_set)*pData(RUV_set)$W_1)
  sd(colMeans(counts(RUV_set)*pData(RUV_set)$W_1))
  plotRLE(RUV_set, outline=FALSE, ylim=c(-4, 4), col=Par_colors[RUV_x])
  plotPCA(RUV_set, col=Par_colors[RUV_x], cex=1.2)
  #|----Finishing----|####
  dir.create("00.Results", showWarnings = FALSE)
  save(Txi,RUV_set,file = "00.Results//SalmonRawData.Rdata")
}