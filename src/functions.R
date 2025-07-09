read_narrowPeak <- function(path, id, group, file){
    #' @description Reads narrowPeak file.
    #' Adds metadata columns of id (unique sample identifer) 
    #' and group (a shared grouping variable to organize distinct replicates).
  
  Require('dplyr')
  Require('magrittr')

  peaks = read.table(path,  
                      header = F, 
                      sep = '\t') %>%
    dplyr::mutate(id = id,
                  group = group)

  ## set standard macs2 narrowPeak column names 
  colnames(peaks) = c('seqnames', 
                        'start', 
                        'end', 
                        'name', 
                        'score', 
                        'strand', 
                        'signalValue', 
                        'pValue', 
                        'qValue', 
                        'peak',
                        'id',
                        'group')
  return(peaks)
}