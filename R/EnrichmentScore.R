
calculateESx_upGene<-function(calc_r,calc_g){
  total<-inner_join(calc_r,calc_g,by = "Repeats")
  names(total)<-c("Repeats","r","g")
  sumR<-sum(total$r)
  sumG<-sum(total$g)
  total$R<-c(sumR)
  total$G<-c(sumG)
  total$ESx<-c((total$r/total$R)/(total$g/total$G))
  p<-apply(as.matrix(total[,2:5]), 1, function(x)    # to calculate pvalue with using fisher exact test
    fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)
  total$pValue<-signif(p,digits = 6)
  pAdjust<-p.adjust(p, method = "BH", n = length(p))     # to calculate p adjust value with using binom test
  total$pAdjust<-signif(pAdjust, digits = 6)
  return(total)
}

calculateESx_upCluster<-function(calc_r,calc_g){
  total<-inner_join(calc_r,calc_g,by = "Repeats")
  names(total)<-c("Repeats","Modules","r","g")
  u<-total$Modules   
  u<-as.data.frame(table(u)) #get unique fields of modules and their counts
  colnames(u)<-c("Modules","R")
  total<-merge(total,u,by="Modules")
  sumG<-sum(total$g)
  total$G<-c(sumG)
  total$ESx<-c((total$r/total$R)/(total$g/total$G))
  p<-apply(as.matrix(total[,3:6]), 1, function(x)         # to calculate pvalue with using fisher exact test
    fisher.test(matrix(round(x), ncol=2), workspace=1e9)$p.value)
  total$pValue<-signif(p,digits = 6)
  pAdjust<-p.adjust(p, method = "BH", n = length(p))
  total$pAdjust<-signif(pAdjust, digits = 6)
  return(total)
}

calculateShuffleEsx<-function(calcCountShfl_T2,calcCountIntrsct_T2,calcCountGenome_T2,nshuffle){
  orderList<-inner_join(calcCountShfl_T2,calcCountIntrsct_T2,by="Repeats")
  orderList<-inner_join(orderList,calcCountGenome_T2,by="Repeats")
  p<-(orderList$Shuffle)/((orderList$Genome)*nshuffle)
  bt <- function(x, n, p ) {binom.test(x, n, p, alternative="greater", conf.level = 0.95)$p.value}    # to calculate p adjust value with using binom test
  orderList$pValue <- mapply(bt, orderList$Shuffle, orderList$Genome ,p)
  pAdjust<-p.adjust(orderList$pValue, method = "BH", n = length(orderList$pValue))   
  orderList$pAdjust<-signif( -log10(pAdjust), digits = 6)
  orderList$pValue <-signif( -log10(orderList$pValue), digits = 6)
  return(orderList)
}

enrichmentAnalysis<-function(x){
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  eg <- bitr(x, fromType="SYMBOL", toType= c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  
  #GO over-representation test
  
  eg_ensembl <- eg$ENSEMBL
  ego <- enrichGO(gene         = eg_ensembl,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable = TRUE,
                  keyType = "ENSEMBL")
  
  return(ego)
  
}

makeGeneModules<-function(dtego){
  
  
  dtego$geneID<-strsplit(dtego$geneID,split="/")
  dflast<-data.frame()
  
  for (index in 1:nrow(dtego)) {
    if(index==1){
      df1<-data.frame("Genes"=dtego$geneID[[index]], "Clusters"=dtego$Description[index])
      df2<-data.frame("Genes"=dtego$geneID[[index+1]], "Clusters"=dtego$Description[index+1])
      dflast<-rbind(df1,df2)
    }else if(index>2){
      dfnew<-data.frame("Genes"=dtego$geneID[[index]], "Clusters"=dtego$Description[index])
      dflast<-rbind(dflast,dfnew)
    }
  }
  
  return(dflast)
}


makeResultTable<-function(dtego,RepeatsAndGenes){
  
  dt<-data.frame()
  uGenes<-unique(RepeatsAndGenes$genes)
  dtego$geneID<-strsplit(dtego$geneID,split="/")
  
  for (i in 1:nrow(dtego)){
    
    geneList<-dtego$geneID[[i]]
    vector<-character()
    for (j in 1:length(geneList)){
      
      index<-geneList[j]==RepeatsAndGenes$genes
      v<-RepeatsAndGenes$info[index]
      if(j==1){
        vector<-as.vector(v)
      }else{
        vector<-c(vector, as.vector(v))
      }
      
    }
    
    countOfRepeats<-length(vector)
    repeatsList<-unique(vector)
    if(i==1){
      dt<-data.frame("Clusters"=dtego$Description[i],"CountOfGenes"=length(geneList), "CountOfRepeats"= countOfRepeats, "Genes"=I(list(geneList)), "Repeats"=I(list(repeatsList)))
    }else{
      dtnew<-data.frame("Clusters"=dtego$Description[i],"CountOfGenes"=length(geneList), "CountOfRepeats"= countOfRepeats, "Genes"=I(list(geneList)), "Repeats"=I(list(repeatsList)))
      dt<-rbind(dt,dtnew)
    }
    
  }
  
  return(dt)
}