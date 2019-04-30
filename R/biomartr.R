
library(biomaRt)
#for humans
#listDatasets(ensembl)


filterGOID<-function(specie,goID){
  if(length(goID)!=0 && specie=="Human_hg19"){
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host = "https://grch37.ensembl.org", GRCh = 37, verbose = TRUE)
    #version 78
    # goID <-c("GO:0003700")
    data <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                 "chromosome_name", "start_position", "end_position",
                                 "strand"),

                  filters = c("go"), #listFilters(ensembl)
                  values = goID,
                  mart=ensembl,
                  verbose = TRUE)

    names(data)<-c("geneID","geneName","chr","start","end","strand")

    hit1<-data$strand=="1"
    data$strand[hit1]<-"+"

    hit2<-data$strand=="-1"
    data$strand[hit2]<-"-"


    return(data)
  }else return(NULL)


}




