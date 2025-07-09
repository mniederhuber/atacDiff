###
# setting column data 
###

## column_data should look like:
##          group <Additional metadata columns, ie timepoint, batch, sex>
## WT_rep1  WT
## WT_rep2  WT
## KO_rep1  KO
##   etc...
column_data <- sample_table %>% # sample table is loaded in loadPeaks.R
  dplyr::select(id, group) %>% ##! select all metadata necessary for deseq model
  tibble::column_to_rownames('id') %>%
  dplyr::mutate(group = factor(group),
 ##! If there are other metadata features - like batch - that should be included in deseq design model, 
 ##! they should be added here and set as factors if they are categorical variables.
               # batch = factor(batch)
                )
## make sure that the control group set in params file is present in column data group
if(!(params$control_group %in% column_data$group){
    stop("Control group must be a grouping variable in column data -- make sure control group is set correctly in params file.")
}


### 
# building count table
### 
message('Building counts table...')
# rsubread will take all of the bam files and count reads within our features, in this case the union peak set
# rsubread is slow so only run if necessary - setting more threads in params will increase speed 
##! but make sure resource request on HPC matches params threads
counts <- Rsubread::featureCounts(sample_table$bam,
                                  annot.ext = union.saf, # union.saf is made in loadPeaks.R
                                  isPairedEnd = T, 
                                  nthreads = params$nthreads)

## set the column names for the counts matrix to match the sample ids in sample_table
colnames(counts$counts) <- sample_table$id
head(counts$counts)


##! It's critical that the order of rows in column_data matches the order of columns in the counts matrix
##  ie. for the example column_data above the counts matrix should have columns ordered: WT_rep1, WT_rep2, KO_rep1 ...
## this checks that everything is ordered correctly, but it's advisable to double check manually
if(!all(colnames(counts$counts) == rownames(column_data))){
  stop("!Something's wrong with coldata or counts columns for DESeq!")
}

## set design as needed
design <- ~ group 

dds <- DESeq2::DESeqDataSetFromMatrix(counts$counts, column_data, design = design)
dds <- DESeq2::DESeq(dds)
## vst normalization is used only for pca and any plotting of relative counts, not for actual DESeq differential contrast
vsd <- DESeq2::vst(dds, blind = F) # blind indicates if the design should be considered for vst normalization

###
# Run contrasts 
### 

message("Running differential contrasts...")
## test groups are all groups that are not defined as control 
test_groups <- setdiff(column_data$group, params$control_group)

## Create all contrasts
## loop through test_groups , returns a list of lists, each list element is a different contrast result with contrast name
contrast_results <- lapply(test_groups, function(test) {
  contrast_name <- paste0(test, '-', control_group) ## used for naming purposes
  message("Running contrast: ", contrast_name)
  ##! if other metadata columns have been added and a more complex design is being used then this may need to be adjusted
  result <- DESeq2::results(dds, 
                            contrast = c('group', test, control_group)) %>%
                    data.frame()
  ## add suffix to describe contrast name for each results column
  ## this is optional, but useful with multiple contrasts, especially if you want to combine results into a single dataframe later
  colnames(result) <- paste0(colnames(result), '.', contrast_name)
  return(list(name = contrast_name, result = result))
})

message('Saving results...')
purrr::map(contrast_results, function(x){
    result = x$result %>%
        dplyr::arrange(desc(log2FoldChange)) 

    fp = paste0(gsub('/$','',params$outdir),'/',x$name,'.csv')
    write.table(result, file = fp, quote = F, row.names = F, sep = ',')
})

message('Done!')