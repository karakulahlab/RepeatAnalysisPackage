library(dplyr)

readfile <- function(filepath){    #to read file and return data frame
  df <- read.table(filepath,header = FALSE)
  return(df)
}

setHeadrs <- function(df,vectHeaders){ #to set header to columns
  names(df) <- vectHeaders[1:length(names(df))]
  return(df)
}


filterAsGene <- function(df,vGenes){  #to select genes
  filtered_df <- df %>%
    dplyr::select(chr, start, end, strand, genes) %>%
    filter(genes %in% vGenes)
  return(filtered_df)
}




