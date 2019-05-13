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
# readRepeatMasker<-function(assembly){
#   datatable(readr::read_csv("https://docs.google.com/spreadsheets/d/1f98kFdj9mxVDc1dv4xTZdx8iWgUiDYO-qiFJINvmTZs/export?format=csv&gid=2046985454"),
#             filter = 'top',
#             options = list(scrollX = TRUE, keys = TRUE, pageLength = 40),
#             rownames = FALSE)
#
# }
