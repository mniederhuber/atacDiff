

message("Reading peaks...")
## map each unique sample (id) to read_narrow_peak()
## returns a list of granges objects - one for each unique sample described by id column
peaks <- purrr::map(sample_table$id, function(x){
  peak_df = read_narrowPeak(sample_table[sample_table$id == x,]$peak_file, x, sample$condition)
  gr = GRanges(peak_df) # convert the peak dataframe to GRanges object
  return(gr)
}) %>%
   set_names(sample_table$id) %>%
  GRangesList() 

message("Filtering peaks...")
## peaks.rowBound is a dataframe of all peaks from all samples
## the column `qc_keep` is added to mark which peaks are above the bottom 10% qValues for each sample
## and that peaks are on chromosomes of interest, which is all autosomes, X, and Y
## this is subsequently used for plotting cumsum of qvalues and bar plot for qc filtering
peaks.rowBound <- names(peaks) %>%
  purrr::map(function(x){
    peaks[[x]] %>%
      data.frame() %>%
      dplyr::mutate(qc_keep = ifelse(qValue > quantile(.$qValue, probs = params$qvalue_percentile)[[1]] & seqnames %in% params$chromosomes, T, F),
                    id = x) 
}) %>%
  GRangesList() %>%
  unlist()

message("Splitting peaks by condition...")
## peaks.condition takes the rowBound peaks and applies the filter, then splits out by group(condition) 
peaks.condition <- peaks.rowBound[peaks.rowBound$qc_keep,] %>%
  split(mcols(.)$group)

message("Finding high-confidence peaks...")
## peaks.hc is a further annotation of our peaks to just those that are reproducible within conditions
## we map makeUnion to each of the peak conditions 
## makeUnion splits each condition into its constituent samples (reps)
## and then marks those peaks that are reproduced (overlap) in at least 2 replicates
peaks.hc <- names(peaks.condition) %>%
  purrr::map(function(x){
    makeUnion(peaks.condition[[x]], x)
  }) %>%
  set_names(names(peaks.condition))

## we then apply the filter in the above step to keep only the final set of high confidence (hc) peaks in each condition
## This is useful for finding overlaps and peak categories
main <- names(peaks.hc) %>%
  purrr::map(function(x){
    peaks.hc[[x]] %>%
      dplyr::filter(keep)
  }) %>%
  set_names(names(peaks.hc))

### union for deseq, note that we don't use the high confidence set here because we don't want to bias our sampling at all
# we use a more liberal union set of all peaks that pass our initial qc filtering 
union <- peaks.condition %>%
  unlist() %>%
  GenomicRanges::reduce()
 
# convert the union set to a SAF format, which has specific column names
# saf format is required by rsubread
union.saf <- union %>%
   data.frame() %>%
   dplyr::mutate(GeneID = 1:nrow(.)) %>%
   dplyr::rename('Chr' = 'seqnames', 'Start' = 'start', 'End' = 'end', 'Strand' = 'strand') %>%
   dplyr::select(GeneID,Chr,Start,End,Strand)

message('Peaks loaded and processed!')