library(dplyr)

readfile <- function(filepath){    #to read file and return data frame
  df <- read.table(filepath,header = FALSE)
  return(df)
}

setHeadrs <- function(df,vectHeaders){ #to set header to columns
  names(df) <- vectHeaders[1:length(names(df))]
  return(df)
}

#
# filterAsGene <- function(df,vGenes){  #to select genes
#   filtered_df <- df %>%
#     dplyr::select(chr, start, end, strand, genes) %>%
#     filter(genes %in% vGenes)
#   return(filtered_df)
# }
# library(biomartr)
# #https://drive.google.com/open?id=14-CbLlCgnYRr_GX6jiu-Ixn-5MziNcnQ
#
readRepeatMasker<-function(assembly){
library(googledrive)
 dt<-drive_download(
  as_id("14-CbLlCgnYRr_GX6jiu-Ixn-5MziNcnQ"), type = ".csv", overwrite = TRUE)
 dt<-read.csv(dt$local_path,sep = "\t")
 return(dt)
}
