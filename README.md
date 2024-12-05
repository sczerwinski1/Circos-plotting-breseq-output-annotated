# Circos-plotting-breseq-output-annotated
 Plotting with circos of single nucleotide polymorphisms determinued using breseq with additional SNP annotations

## Dependencies
This script requires the following R packages:
- VariantAnnotation
- Rsamtools
- circlize
- GenomicRanges
- rtracklayer
- Biostrings
- RColorBrewer
- readr
- rvest
- dplyr
- tidyr

Please ensure these are installed before running the script. Most can be installed via CRAN or Bioconductor:
```bash
install.packages(c("circlize", "RColorBrewer", "readr", "rvest", "dplyr", "tidyr"))

## For Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation","Rsamtools","GenomicRanges","rtracklayer","Biostrings"))
