library(magrittr)
library(purrr)
library(dplyr)
library(tibble)
library(GenomicRanges)
library(DESeq2)
library(Rsubread)
library(rtracklayer)
#library(biomaRt)
#library(universalmotif)
#library(BSgenome)
#library(memes)
#library(Biostrings)
source('src/functions.R')

## Set custom biomaRt cache directory
#if (!dir.exists(".cache/biomaRt")) {
#  dir.create(".cache/biomaRt", recursive = TRUE)
#}
#options(biomart.cache.dir = file.path(getwd(), ".cache/biomaRt"))
#Sys.setenv(BIOMART_CACHE = file.path(getwd(), ".cache/biomaRt"))

# Check if the parameters script is provided as an argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  params_file <- args[1]
  # Source the parameters file
  params <- yaml::read_yaml(params_file)
} else if (file.exists('config/R.yaml')){
  params <- yaml::read_yaml('config/R.yaml')
} else {
  stop('No params file!')
}

set.seed(params$seed)

## input - sample table with paths to peak files and bam files, with condition and sample metadata
## columns: id, condition, rep, peak_file, bam_file
sampleTable <- read.csv(params$sample_table)
print(sample_table)

source('src/loadPeaks.R')
source('src/deseq.R')