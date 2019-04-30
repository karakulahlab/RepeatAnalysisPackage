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
filterSubset<-function(input,genome, type, specie){
  filtered_df_kn <- data.frame()
  if(!is.null(genome) & specie=="Human_hg19"){
    if(type=="geneNames"){
      df<-scan(input,character())
      filtered_df_kn <- filterAsGene(genome,df)
      }else if(type=="geneOntologyID"){
      filtered_df_kn <- filterGOID(specie,input)
      }else if(type=="geneModules"){
        df<-read.csv(input,header=TRUE)
        filtered_df_kn <- filterAsGene(genome,df$genes)}
    head(filtered_df_kn)
    return(filtered_df_kn)
  }else{return(NULL)}}


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












