library(dplyr)
require(GenomicRanges)


getUpstream <- function(df,length, isWithGeneBody){ # to get genome interval of up stream of given genome
  
  dftemp <- df
  dfnew <- df
  
  indexP <- dftemp$strand == "+"
  indexN <- dftemp$strand == "-"
  
  dfnew$start[indexP] <- c(dftemp$start[indexP] - length)
  dfnew$end[indexN] <- c(dftemp$end[indexN] + length)
  
  if(isWithGeneBody == FALSE)
  {
    dfnew$end[indexP] <- c(dftemp$start[indexP])
    dfnew$start[indexN] <- c(dftemp$end[indexN] )
  }
  
  return(dfnew)
}

getDownstream <- function(df,length,isWithGeneBody){ # to get genome interval of down stream of given genome
  dftemp <- df
  dfnew <- df
  indexP <- dftemp$strand == "+"
  indexN <- dftemp$strand == "-"
  dfnew$start[indexN] <- dftemp$start[indexN] - length
  dfnew$end[indexP] <- dftemp$end[indexP] + length
  if(isWithGeneBody == FALSE)
  {
    dfnew$end[indexN] <- dftemp$start[indexN]
    dfnew$start[indexP] <- dftemp$end[indexP]
  }
  return(dfnew)
}

getDownAndUpStream <- function(df,len_up,len_down){ # to get genome interval of up and down stream of given genome with its
  dftemp <- df
  dfnew <- df
  
  indexP <- dftemp$strand == "+"
  indexN <- dftemp$strand == "-"
  
  dfnew$start[indexP] <- dftemp$start[indexP]-len_up
  dfnew$end[indexP] <- dftemp$end[indexP]+len_down
  
  dfnew$start[indexN] <- dftemp$start[indexN]-len_down
  dfnew$end[indexN] <- dftemp$end[indexN]+len_up
  
  return(dfnew)
}


makeGRangeObj <- function(df){  # to make grange object to analysis genome arithmetics well
  library(GenomicRanges)
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = strand))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = strand, id=genes))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = strand, id=geneName))
  }
  return(gr)
}

makeGrObj_Unstrand <- function(df){  #within strandness
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = "*"))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), strand = "*", id=genes))
  } 
  return(gr)
}

returnOverlapUpGene<-function(ref,test){  # to return overlap for uploaded genes 
  hits1 <- findOverlaps(ref, test)
  overlaps1 <- pintersect(ref[queryHits(hits1)], test[subjectHits(hits1)])
  hits2 <- findOverlaps(test, ref)
  overlaps2 <- pintersect(test[queryHits(hits2)], ref[subjectHits(hits2)])
  overlaps <- overlaps1
  id <- na.omit(match(overlaps1, overlaps2))
  GenomicRanges::values(overlaps) <- cbind(GenomicRanges::values(overlaps1), elementMetadata(overlaps2)[id, ])
  f <-elementMetadata(overlaps)
  f <-f[,c(1,3,4)]
  header<-c('repeats','genes','hit')
  f<-setHeadrs(f,header)
  elementMetadata(overlaps)<-f
  return(overlaps)
}

returnOverlapUpCluster<-function(ref,test){ # to return overlap for uploaded genes with their modules (clustered genes) 
  hits1 <- findOverlaps(ref, test)
  overlaps1 <- pintersect(ref[queryHits(hits1)], test[subjectHits(hits1)])
  hits2 <- findOverlaps(test, ref)
  overlaps2 <- pintersect(test[queryHits(hits2)], ref[subjectHits(hits2)])
  overlaps <- overlaps1
  id <- na.omit(match(overlaps1, overlaps2))
  GenomicRanges::values(overlaps) <- cbind(GenomicRanges::values(overlaps1), elementMetadata(overlaps2)[id, ])
  f <-elementMetadata(overlaps)
  f <-f[,c(1,3,4)]
  header<-c('repeats','genes','hit')
  f<-setHeadrs(f,header)
  elementMetadata(overlaps)<-f
  return(overlaps)
}

toFindOverlaps<-function(gr_repeats,gr_genome){  #to get overlap in general function
  hits<-findOverlaps(gr_repeats, gr_genome)
  overlaps<-pintersect(gr_repeats[queryHits(hits)], gr_genome[subjectHits(hits)])
  return(overlaps)
}