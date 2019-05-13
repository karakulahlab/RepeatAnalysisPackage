library(DT)
library(dplyr)
library(devtools)



source("R/FileOperations.R")
source("R/GenomicRanges.R")
source("R/EnrichmentScore.R")
source("R/MakePlotDiagrams.R")
source("R/Methylation.R")
source("R/biomartr.R")


#to filter input genes from known gene file
filterSubset<-function(input, assembly, type){
  filtered_df_kn <- data.frame()
  if(!is.null(assembly)){
    if(type=="geneNames"){
      df<-scan(input,character())
      filtered_df_kn <- filterGenes(assembly,df)
    }else if(type=="geneOntologyID"){
      filtered_df_kn <- filterGOID(assembly,input)
    }else if(type=="geneModules"){
      df<-read.csv(input,header=FALSE)
      filtered_df_kn <- filterGenes(assembly,df$V1)}
    head(filtered_df_kn)
    return(filtered_df_kn)
  }else{return(NULL)}}


data<-function(assembly){
  information<-readRepeatMasker(assembly)
  return(information)
}

#
# #to filter input genes from known gene file
# filterSubset<-function(input, genome, type, species){
#   filtered_df_kn <- data.frame()
#   if(!is.null(genome) & species=="Human_hg19"){
#     if(type=="geneNames"){
#       df<-scan(input,character())
#       filtered_df_kn <- filterAsGene(genome,df)
#       }else if(type=="geneOntologyID"){
#       filtered_df_kn <- filterGOID(species,input)
#       }else if(type=="geneModules"){
#         df<-read.csv(input,header=TRUE)
#         filtered_df_kn <- filterAsGene(genome,df$genes)}
#     head(filtered_df_kn)
#     return(filtered_df_kn)
#   }else{return(NULL)}}


#to get interval region as regions of input genes
getRegion<-function(input,length,region,isGenebody){
  if(region=="upstream"){
    md <- getUpstream(input, length, isGenebody)
    return(md)
  }else if(region=="downstream"){
    md <- getDownstream(input, length, isGenebody)
    return(md)
  }else if(region=="upstream and downstream"){
    md <- getDownAndUpStream(input, length, length)
    return(md)
  }else{
    print("please control the region")
    return(NULL)
  }
}


#to return intersect between genes which are given as input and repeat annotations
getOverlap<-function(genes,repeats,type,species){
  gr_genes<-makeGRangeObj(genes)
  gr_repeats<-makeGRangeObj(repeats)
  if(type=="geneNames"){
    overlaps<-returnOverlapUpGene(gr_repeats,gr_genes)
    suppressMessages(library(ggplot2))
    suppressMessages(library(interactiveDisplay))
    suppressMessages(library(Biobase))
    return(overlaps)
  }else if(type=="geneModules"){
    overlaps<-returnOverlapUpCluster(gr_repeats,gr_genes)
    suppressMessages(library(ggplot2))
    suppressMessages(library(interactiveDisplay))
    suppressMessages(library(Biobase))
    return(overlaps)
  }
}









