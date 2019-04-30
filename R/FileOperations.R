
library(dplyr)
require(GenomicRanges)

readfile <- function(filepath){    #to read file and return data frame
  df <- read.table(filepath,header = FALSE)
  return(df)
}

writefile <- function(file,filepath){  #to save file to given filepath
  write.table(file, filepath, sep="\t",row.names=FALSE,col.names = FALSE,quote = FALSE)
}


dropcols <- function(df,vectColms){     #to delete columns which you can see as a vector you want delete
  df =df[,!(names(df) %in% vectColms)]
  return(df)
}

setHeadrs <- function(df,vectHeaders){ #to set header to columns
  names(df) <- vectHeaders[1:length(names(df))]
  return(df)
}

filterAsRepeats <- function(df,vRepeats){    #to select repeats

  header <- c('chr','start','end','strand','repeats')
  df<-setHeadrs(df,header)

  filtered_df <- df %>%
    dplyr::select(chr, start, end, strand, repeats) %>%
    filter(repeats %in% vRepeats)

  return(filtered_df)
}

filterAsGene <- function(df,vGenes){  #to select genes
  filtered_df <- df %>%
    dplyr::select(chr, start, end, strand, genes) %>%
    filter(genes %in% vGenes)
  return(filtered_df)
}

filterGenes <-function(genes){ #to select genes from knowngenes for human
  df_known <- readKnownGenes()
  filtered_df_kn <- filterAsGene(df_known,genes)
  return(filtered_df_kn)
}



prepKnownGenes <- function(){  # to prepeare grange object of known genes for human on same strand
  df_known <- readKnownGenes()
  k <- makeGRangeObj(df_known)
  return(k)
}

prepKnownGenes2 <- function(){  # to prepeare grange object of known genes for human on strandness
  df_known <- readKnownGenes()
  k <- makeGrObj_Unstrand(df_known)
  return(k)
}

readKnownGenes <- function(specie){  #to read known genes for human
  if(specie=="Human_hg19"){
    fknown <- read.table("~/R_codes/genomeArithmetic/RepeatAnalysis/Data/known.canonical.bed", header=FALSE)
    header <- c('chr','start','end','strand','genes')
    df_known <- setHeadrs(fknown,header)
    return(df_known)
  }else{return(NULL)}
}


saveCheckFile<-function(df_overlap,df_genome){   #to save file to check later
  df<-merge(df_overlap,df_genome,by="genes")
  df$GeneDescription<-""
  df$RepeatDescription<-""
  df<-df[c(1,12,8,9,10,11,6,13,2,3,4,5)]
  colnames(df)<-c('GeneID','GeneDescription','Chr','Start','End','Strand','RepeatID','RepeatDescription','Chr','Start','End','Strand')
  write.table(df, file="~/R_codes/genomeArithmetic/RepeatAnalysis/Data/overlapped_gene_repeat_information.bed", quote=F, sep="\t", row.names=F, col.names=T)
}

rGrAsDataTable<-function(gr){    #to convert from grange object to data table
  df1 <- data.frame(seqnames=seqnames(gr),
                    starts=start(gr)-1,
                    ends=end(gr),
                    strands=strand(gr))
  df2<-elementMetadata(gr)
  df<-cbind(df1,df2)
  return(df)
}

rGrAsDataTableForGenome<-function(gr){    #to convert from grange object to data table for Genome Annotation
  df1 <- data.frame(seqnames=seqnames(gr),
                    starts=start(gr)-1,
                    ends=end(gr),
                    strands=strand(gr))
  genes<-elementMetadata(gr)@listData[["id"]]
  df<-cbind(df1,genes)
  return(df)
}


toMakeShuffle<-function(genes,nshuffle){   #to get shuffle genome interval with using bedr shuffle function
  genome <- trbl_genome(
    ~chrom, ~size,
    # "chr1", 1e6,
    # "chr2", 2e6
    "chr1", 1e6
  )
  new_gr<-bed_shuffle(genes, genome, seed = nshuffle)
  new_gr<-dropcols(new_gr,c("name","score"))
  colnames(new_gr)<-c("chr","start","end","strand")
  return(new_gr)
}
