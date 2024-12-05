# Libs
library(VariantAnnotation)
library(Rsamtools)
library(circlize)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(RColorBrewer)
library(readr)
library(rvest)
library(dplyr)
library(tidyr)

# Define file paths to breseq outputs
vcf_files <- c(
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5265_1_WTnPAO1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5265_2_PAt10a.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5265_3_PAt10b.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5265_4_PAt10c.2.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5265_5_PAt10d.1.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5265_6_PAt10e.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_11_PA+SAt10a.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_12_PA+SAt10b.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_13_PA+SAt10c.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_14_PA+SAt10d.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_15_PA+SAt10e.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_1_L19PAt10a.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_2_L19PAt10b.3/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_3_L19PAt10c.3/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_4_L19PAt10d.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_5_L19PAt10e.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_6_L19PA+SAt10a.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_7_L19PA+SAt10b.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_8_L19PA+SAt10c.1/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_9_L19PA+SAt10d.2/data/output.vcf",
  "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_10_L19PA+SAt10e.3/data/output.vcf"
)

# function to read VCF files and handle errors and duplicates
read_vcf_with_header_check_and_handle_duplicates <- function(vcf_file) {
  tryCatch({
    cat("Reading VCF file:", vcf_file, "\n")
    vcf <- readVcf(vcf_file)
    if (!inherits(vcf, "VCF")) {
      stop("The object read from the file is not a valid VCF object.")
    }
    vcf_ranges <- rowRanges(vcf)
    
    # Detect and handle duplicate row names
    df <- as.data.frame(vcf_ranges, stringsAsFactors = FALSE)
    df$seqnames <- as.character(df$seqnames)
    df$unique_id <- paste(df$seqnames, df$start, df$end, sep = ":")
    
    duplicated_rows <- duplicated(df$unique_id)
    if (any(duplicated_rows)) {
      cat("Duplicates found in VCF file:", vcf_file, "\n")
      df <- df[!duplicated_rows, ]
      vcf_ranges <- GRanges(seqnames = df$seqnames, ranges = IRanges(start = df$start, end = df$end))
    }
    vcf
  }, error = function(e) {
    cat("Error reading VCF file:", vcf_file, "\n", e$message, "\n")
    NULL
  })
}

# load VCF files
vcf_list <- lapply(vcf_files, read_vcf_with_header_check_and_handle_duplicates)

# remove any NULL entries from the list
vcf_list <- Filter(Negate(is.null), vcf_list)

# ensure all VCF files are read correctly
if (length(vcf_list) != length(vcf_files)) {
  cat("Warning: Not all VCF files were loaded correctly.\n")
}

# additional Debug Information
cat("Number of valid VCF objects read:", length(vcf_list), "\n")
for (i in seq_along(vcf_list)) {
  cat("VCF file", i, "has", length(rowRanges(vcf_list[[i]])), "SNPs\n")
}

# define file paths for the reference genome and annotation file
reference_genome <- "~/../../Volumes/Extreme_SSD/seq/PAO1/Pseudomonas_aeruginosa_PAO1_107.fna"
reference_gff3 <- "~/../../Volumes/Extreme_SSD/seq/BreseqOut/BsqOutDirs/5238_1_L19PAt10a.1/data/reference.gff3"

# load the reference genome using Biostrings
ref_genome <- readDNAStringSet(reference_genome, format = "fasta")

# ensure reference genome is in a data frame format suitable for circos.initialize
genome_df <- data.frame(chr = "NC_002516.2", start = 1, end = width(ref_genome))

# load the GFF3 file
ref_gff <- import(reference_gff3, format = "gff3")

# read the tabular annotation file using tidyverse
annotation_file_path <- "~/../../Volumes/Extreme_SSD/seq/BreseqOut/CoevoComparisons2/AnnotationsFileFIN.txt"
annotation_data <- read_tsv(annotation_file_path)

# ensure gene_name column is character and fix any special character issues
annotation_data <- annotation_data %>%
  mutate(gene_name = as.character(gene_name))

# adjust annotation positions to avoid overlap
adjust_text_positions <- function(annotation_df, buffer = 0.1, proximity_threshold = 5000, additional_spacing = 0.05) {
  # ensure start_position and end_position are numeric and handle non-numeric values
  annotation_df$start_position <- as.numeric(gsub(",", "", annotation_df$start_position))
  annotation_df$end_position <- as.numeric(gsub(",", "", annotation_df$end_position))
  
  # remove rows with missing or non-numeric start_position values
  annotation_df <- annotation_df[!is.na(annotation_df$start_position) & !is.na(annotation_df$end_position), ]
  
  # deduplicate based on gene_name, arbitrarily choosing the first occurrence
  annotation_df <- annotation_df[!duplicated(annotation_df$gene_name), ]
  
  # order the dataframe by start_position
  annotation_df <- annotation_df[order(annotation_df$start_position), ]
  
  y_positions <- rep(buffer, nrow(annotation_df))  # initialize y-positions with a buffer for extra spacing
  
  for (i in 2:nrow(annotation_df)) {
    for (j in 1:(i-1)) {
      if (abs(annotation_df$start_position[i] - annotation_df$start_position[j]) < proximity_threshold) {
        if (y_positions[i] == y_positions[j]) {
          y_positions[i] <- y_positions[j] + buffer
        } else {
          y_positions[i] <- y_positions[j] + additional_spacing
        }
      }
    }
  }
  
  # ensure only one of the overlapping annotations is moved outward
  for (i in 2:nrow(annotation_df)) {
    if (y_positions[i] == y_positions[i - 1]) {
      y_positions[i] <- y_positions[i] + buffer
    }
  }
  
  annotation_df$y <- y_positions
  return(annotation_df)
}

# adjust the positions for the annotation data with a buffer, proximity threshold, and additional spacing
buffer <- 0.1
proximity_threshold <- 5000
additional_spacing <- 0.05

adjusted_annotation_track <- adjust_text_positions(annotation_data, buffer = buffer, proximity_threshold = proximity_threshold, additional_spacing = additional_spacing)

# save the plot to a file with high resolution and white background
save_circos_plot <- function(filename, format = c("pdf", "png"), width = 20, height = 17, res = 200) {
  format <- match.arg(format)
  if (format == "pdf") {
    pdf(file = filename, width = width, height = height, bg = "white") 
  } else if (format == "png") {
    png(file = filename, width = width, height = height, units = "in", res = res, bg = "white") 
  }
}

# call this function with the desired format before plotting
save_circos_plot("CircosPlot8-15-24_WhtBckgrnd_1.4.2.png", format = "png")

# set the font
par(family = "Arial")

# initialize the circos plot
circos.clear()
circos.par("track.height" = 0.05, "cell.padding" = c(0.005, 0.01, 0.005, 0.01), "track.margin" = c(0.01, 0.01))
circos.initialize(factors = genome_df$chr, xlim = genome_df[, c("start", "end")])

# plot the annotations first
circos.genomicTrackPlotRegion(
  adjusted_annotation_track,
  ylim = c(0, max(adjusted_annotation_track$y) + 0.1), 
  bg.border = NA, 
  panel.fun = function(region, value, ...) {
    circos.text(
      x = region$start_position, 
      y = value$y,                     
      labels = value$gene_name,        
      facing = "clockwise", 
      niceFacing = TRUE, 
      adj = c(0, 0.5), 
      cex = 0.8,
      col = "black" 
    )
  }, 
  track.height = 0.4 
)

# create an empty track for the axis
circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.border = NA, panel.fun = function(region, value, ...) {
  circos.axis(h = "top", major.at = seq(0, max(genome_df$end), by = 250000), labels = format(seq(0, max(genome_df$end), by = 250000), scientific = FALSE), labels.cex = 0.7, lwd = 0.7, col = "black", labels.col = "black") # Changed axis and labels to black
}, track.height = 0.03)


# group samples by condition
conditions <- list(
  condition1 = 1,
  condition2 = 2:6,
  condition3 = 7:11,
  condition4 = 12:16,
  condition5 = 17:21
)

# function to get SNPs appearing in all samples
get_common_snps <- function(vcf_list) {
  snp_counts <- list()
  
  for (vcf in vcf_list) {
    vcf_ranges <- rowRanges(vcf)
    df <- as.data.frame(vcf_ranges, stringsAsFactors = FALSE)
    df$unique_id <- paste(df$seqnames, df$start, df$end, sep = ":")
    snp_counts[[length(snp_counts) + 1]] <- df$unique_id
  }
  
  common_snps <- Reduce(intersect, snp_counts)
  return(common_snps)
}

# get SNPs common to all samples
common_snps <- get_common_snps(vcf_list)

# function to plot SNPs by condition with varying line heights and filtering
plot_snp_condition <- function(condition_indices, col, exclude_snps, track_height = 0.1) {
  if (any(condition_indices > length(vcf_list))) {
    stop("Condition index out of range. Ensure the condition indices are within the range of the VCF list.")
  }
  
  all_snp_data <- list()
  
  for (i in condition_indices) {
    vcf <- vcf_list[[i]]
    vcf_ranges <- rowRanges(vcf)
    
    # exclude SNPs found in exclude_snps
    df <- as.data.frame(vcf_ranges, stringsAsFactors = FALSE)
    df$unique_id <- paste(df$seqnames, df$start, df$end, sep = ":")
    df <- df[!(df$unique_id %in% exclude_snps), ]
    
    # check if vcf_ranges is empty
    if (nrow(df) == 0) next
    
    snp_data <- data.frame(
      chrom = df$seqnames,
      start = df$start,
      end = df$end,
      value = 1,
      stringsAsFactors = FALSE
    )
    
    all_snp_data[[length(all_snp_data) + 1]] <- snp_data
  }
  
  merged_snp_data <- do.call(rbind, all_snp_data)
  
  # explicitly remove row names after merging
  row.names(merged_snp_data) <- NULL
  
  # aggregate SNPs by position and sum their values
  snp_count <- aggregate(value ~ chrom + start + end, data = merged_snp_data, FUN = sum)
  
  circos.genomicTrackPlotRegion(snp_count, ylim = c(0, max(snp_count$value)), bg.border = NA, panel.fun = function(region, value, ...) {
    circos.genomicRect(region, ybottom = 0, ytop = value$value, col = adjustcolor(col, alpha.f = 1), border = col, lwd =2, ...)
  }, track.height = track_height)
}

# manually specify colors using RColorBrewer
# colors <- c("#f58231", "#74B435", "#6E2AAD", "#46f0f0", "#f032e6")
# darker colors
colors <- c("#F06A0B", "#68CB06", "#7A11DC", "#00D3D3", "#E011D5")

# plot SNPs for each condition with specified track heights
plot_snp_condition(conditions$condition5, colors[1], common_snps, track_height = 0.07) 
plot_snp_condition(conditions$condition4, colors[2], common_snps, track_height = 0.07)
plot_snp_condition(conditions$condition3, colors[3], common_snps, track_height = 0.07)
plot_snp_condition(conditions$condition2, colors[4], common_snps, track_height = 0.07)
plot_snp_condition(conditions$condition1, colors[5], common_snps, track_height = 0.07)


# add a legend to identify the colors
legend(x = 0.75, y = 1, legend = c("+Luz19 +S.aureus", "+Luz19", "+S.aureus", "PAO1 t10", "wt Ancestral"), fill = colors, title = "Conditions", cex = 1.2, bty = "n", text.col = "black") # Changed legend text to black

# finalize the plot
dev.off()

circos.clear()
