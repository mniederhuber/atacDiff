# atacDiff

A simple R pipeline for differential analysis of ATAC-seq peak data using DESeq2.

## Dependencies

- R packages:
  - magrittr
  - purrr
  - dplyr
  - tibble
  - GenomicRanges
  - DESeq2
  - Rsubread
  - rtracklayer
  - yaml

## Usage

```
Rscript src/main.R config/R.yaml
```

## Input

### Sample Table

A CSV file with columns:
- `id`: Unique sample identifier
- `group`: Experimental condition/group
- `rep`: Replicate number
- `peak_file`: Path to MACS2 narrowPeak file
- `bam_file`: Path to aligned BAM file

```
id,group,rep,peak_file,bam_file
WT_rep1,untreated,1,data/peaks/WT_rep1_peaks.narrowPeak,data/bams/WT_rep1.bam
WT_rep2,untreated,2,data/peaks/WT_rep2_peaks.narrowPeak,data/bams/WT_rep2.bam
WT_rep3,untreated,3,data/peaks/WT_rep3_peaks.narrowPeak,data/bams/WT_rep3.bam
KO_rep1,knockout,1,data/peaks/KO_rep1_peaks.narrowPeak,data/bams/KO_rep1.bam
KO_rep2,knockout,2,data/peaks/KO_rep2_peaks.narrowPeak,data/bams/KO_rep2.bam
KO_rep3,knockout,3,data/peaks/KO_rep3_peaks.narrowPeak,data/bams/KO_rep3.bam 
```

### Configuration

Parameters in `config/R.yaml`:

```yaml
sample_table: 'sample_table.csv'
outdir: 'rOut'
chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, X, Y]
qvalue_percentile: 0.1
nthreads: 12
control_group: 'untreated'
seed: 1234
```

`chromosomes` -- sets which chromomes to *keep* -- note that there is not automated check to ensure seqname equivalency between params file and peak ranges!! 
`qvalue_percentile` -- Macs2 qValue percentile, peaks with qValue below threshold will be filtered out, filter is on per sample basis
`control_group` -- defines a grouping variable to use as control for `DESeq2`
`nthreads` -- set the number of threads to use for featureCounts ** make sure to also request enough threads/cpus if running script on shared resources **

## Output

Results are saved to the directory specified by `outdir` parameter. For each contrast, a CSV file is generated with differential analysis results. 