# Script for coverting bigWig + bed files to counts

# Partially adapted from https://lcolladotor.github.io/protocols/bigwig_DEanalysis/

# set up
rm(list = ls())
library(rtracklayer)
library(GenomicRanges)

# control panel
raw_data_folder <- '/Users/sobahytm/BESE_394A/class_3/mycode/new_dataset/big_wig_bed_files'
read_length <- 36
chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))
reference_file <- '/Users/sobahytm/BESE_394A/class_3/mycode/new_dataset/by_conditions/bed_files/reference.bed'

# bw files
bw_files <- file.path(raw_data_folder, dir(raw_data_folder, pattern = '*.bw'))
bw_files <- bw_files[1:6]

# loading reference bed
peaks <- import(reference_file)

# count matrix
count_matrix <- matrix(0, length(peaks), length(bw_files))
rownames(count_matrix) <- paste0(seqnames(peaks), '_', start(peaks), '_', end(peaks))
colnames(count_matrix) <- letters[1:length(bw_files)]

# looping over files
for(i in 1:length(bw_files)){
  
  # current files
  print(paste0('sample ', i, ' out of ', length(bw_files)))
  bw_file <- bw_files[i]
  
  # sample name
  sample_name <- gsub(raw_data_folder, '', bw_file, fixed = TRUE)
  sample_name <- gsub('.bw', '', sample_name, fixed = TRUE)
  sample_name <- gsub('/', '', sample_name, fixed = TRUE)
  sample_name <- strsplit(sample_name, '_')[[1]][2]
  if(grepl('HL60', sample_name)){
    sample_name <- paste0('T0h-', sample_name)
  }else{
    sample_name <- paste0('T', sample_name) 
  }
  
  # loadind and downsizing the bigwigfile
  bw_file_list <- BigWigFileList(bw_file)
  coverage <- import(bw_file_list[[1]], as = 'RleList')
  # head(coverage)
  coverage <- coverage[names(coverage) %in% chromosomes]
  
  # split the peaks across chromosomes
  peaks_list <- split(peaks, seqnames(peaks))
  #to_keep <- which(sapply(peaks_list, length) > 0)
  #peaks_list <- peaks_list[to_keep]
  
  # ensuring peaks list and coverage have the same chromosomes
   common_chromosomes <- intersect(names(coverage), 
                                   as.character(unique(seqnames(peaks_list))))
   peaks_list <- peaks_list[names(peaks_list) %in% common_chromosomes]
   coverage <- coverage[names(coverage) %in% common_chromosomes]
  
  # coverage per peak
  coverage <- coverage[names(peaks_list)]
  peaks_coverage <- Views(coverage, ranges(peaks_list))
  # peaks_coverage$chr1
  # peaks_coverage$chr1[[1]]
  
  # count values
  counts <- sapply(peaks_coverage, sum)
  # head(counts$chr1) 
  # sum(sapply(counts, length))
  # length(peaks)
  
  # ensuring to have the right peak information
  chrs <- rep(names(peaks_coverage), sapply(peaks_coverage, length))
  starts <- sapply(peaks_coverage, start)
  ends <- sapply(peaks_coverage, end)

  # converting to vector
  counts <- unlist(counts)
  names(counts) <- paste0(chrs, '_', unlist(starts), '_', unlist(ends))
  # head(counts)
  
  # rounding up
  counts <- round(counts / read_length)
  # head(counts)
  
  # fractions of reads in peaks
  sum(counts)
  all_counts <- sapply(coverage, sum)
  all_counts <- sum(all_counts)
  all_counts <- round(all_counts / read_length)
  all_counts
  sum(counts) / all_counts
  
  # count as data frame
  count_matrix[names(counts), i] <- counts
  colnames(count_matrix)[i] <- sample_name
  # head(count_matrix)
  
}

# writing
count_matrix <- as.data.frame(count_matrix)
count_matrix <- cbind(peak = rownames(count_matrix), count_matrix)
head(count_matrix)
write.csv(count_matrix, row.names = FALSE,
          file = 'count_matrix.csv')
