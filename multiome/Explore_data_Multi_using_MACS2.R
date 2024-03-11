#load library...
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
#load data...
counts <- Read10X("/Users/sobahytm/BESE_394A/class_5/data/filtered_feature_bc_matrix")
fragpath <- '/Users/sobahytm/BESE_394A/class_5/data/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz'
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
# create a Seurat object containing the RNA adata
mouse_brain <- CreateSeuratObject(
  counts = counts[['Gene Expression']],
  assay = "RNA"
)
# create ATAC assay and add it to the object
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
atac_counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

mouse_brain[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
################################################################################
# Calculate metrics for QC
################################################################################
DefaultAssay(mouse_brain) <- "ATAC"
mouse_brain <- NucleosomeSignal(mouse_brain)
mouse_brain <- TSSEnrichment(mouse_brain)

VlnPlot(
  object = mouse_brain,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
#save...

total_fragments <- CountFragments(fragments = fragpath)
rownames(total_fragments) <- total_fragments$CB

mouse_brain$fragments <- total_fragments[colnames(mouse_brain), "frequency_count"]
mouse_brain <- FRiP(object = mouse_brain, assay = 'ATAC', total.fragments = 'fragments')
mouse_brain$nucleosome_group <- ifelse(mouse_brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

mouse_brain$high.tss <- ifelse(mouse_brain$TSS.enrichment > 3, 'High', 'Low')
################################################################################
# Filter cells based on QC metrics
################################################################################
mouse_brain <- subset(
  x = mouse_brain,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 & 
    FRiP > 0.2
)
################################################################################
# Call peaks using MACS2
################################################################################

peaks_macs2 <- CallPeaks(mouse_brain)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks_macs2 <- keepStandardChromosomes(peaks_macs2, pruning.mode = "coarse")
peaks_macs2 <- subsetByOverlaps(x = peaks_macs2, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak called by MACS2
macs2_counts <- FeatureMatrix(
  fragments = Fragments(mouse_brain),
  features = peaks_macs2,
  cells = colnames(mouse_brain)
)

################################################################################
# Call peaks using FindTopFeatures
################################################################################
DefaultAssay(mouse_brain) <- "ATAC"
mouse_brain_Signac <- FindTopFeatures(mouse_brain, min.cutoff = 5)
mouse_brain_Signac <- RunTFIDF(mouse_brain_Signac)
mouse_brain_Signac <- RunSVD(mouse_brain_Signac)
# Extract the top features (peaks)
top_features <- VariableFeatures(object = mouse_brain_Signac)
################################################################################
# Compare peaks by peaks numbers
################################################################################
# Convert peak coordinates to GRanges objects
#macs peaks
peaks_macs2_granges <- GRanges(
  seqnames = Rle(seqnames(peaks_macs2)),
  ranges = IRanges(start = start(peaks_macs2), end = end(peaks_macs2)),
  strand = strand(peaks_macs2)
)
#FindTopFeatures paeks
# Split the strings in top_features to extract chromosome, start, and end
split_features <- strsplit(top_features, "-")
# Extract chromosome, start, and end from split_features
chr <- sapply(split_features, `[`, 1)
start <- as.numeric(sapply(split_features, `[`, 2))
end <- as.numeric(sapply(split_features, `[`, 3))
# Create a strand vector (assuming all strands are positive)
strand <- rep("+", length(chr))
# Create GRanges object for top_features
top_features_granges <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end), strand = strand)
# Compare peaks called by MACS2 and top features
overlapping_peaks <- findOverlaps(peaks_macs2_granges, top_features_granges)
unique_peaks_macs2 <- setdiff(1:length(peaks_macs2), subjectHits(overlapping_peaks))
unique_top_features <- setdiff(1:length(top_features), queryHits(overlapping_peaks))
# Print the numbers
cat("Number of peaks called by MACS2:", length(peaks_macs2), "\n")
cat("Number of top features identified:", length(top_features), "\n")
cat("Number of overlapping peaks:", length(overlapping_peaks), "\n")
cat("Number of unique peaks called by MACS2:", length(unique_peaks_macs2), "\n")
cat("Number of unique top features:", length(unique_top_features), "\n")
#Plot the numbers ...
# Create a data frame with the counts
data <- data.frame(
  Type = c("MACS2 Peaks", "Top Features", "Overlapping Peaks", "Unique MACS2 Peaks", "Unique Top Features"),
  Count = c(length(peaks_macs2), length(top_features), length(overlapping_peaks), length(unique_peaks_macs2), length(unique_top_features))
)

# Create the bar plot
bar_plot <- ggplot(data, aes(x = Type, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Comparison of Peak Counts", y = "Count", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show the plot
print(bar_plot)
################################################################################
# Compare peaks by peaks distribution
################################################################################
# Calculate sample size
sample_size <- min(length(peaks_macs2_granges), length(top_features_granges), 100)  # Sample only 100 peaks for visualization

# Subsample data
peaks_macs2_granges_sample <- sample(peaks_macs2_granges, sample_size)
top_features_granges_sample <- sample(top_features_granges, sample_size)
# Check if the subsampled peak data is not empty or NULL
if (!is.null(peaks_macs2_granges_sample) && length(peaks_macs2_granges_sample) > 0 &&
    !is.null(top_features_granges_sample) && length(top_features_granges_sample) > 0) {
  
  # Plot genomic locations of subsampled peaks identified by MACS2
  peaks_macs2_plot <- ggplot(as.data.frame(peaks_macs2_granges_sample)) +
    geom_segment(aes(x = start, xend = end, y = seqnames, yend = seqnames), size = 2, color = "blue") +
    geom_point(aes(x = (start + end) / 2, y = seqnames), color = "blue", size = 2) +  # Add points
    labs(title = "MACS2 Peaks Distribution", x = "Genomic Position", y = "Chromosome") +
    theme_minimal()
  
  # Plot genomic locations of subsampled peaks identified by FindTopFeatures
  top_features_plot <- ggplot(as.data.frame(top_features_granges_sample)) +
    geom_segment(aes(x = start, xend = end, y = seqnames, yend = seqnames), size = 2, color = "red") +
    geom_point(aes(x = (start + end) / 2, y = seqnames), color = "red", size = 2) +  # Add points
    labs(title = "FindTopFeatures Peaks Distribution", x = "Genomic Position", y = "Chromosome") +
    theme_minimal()
  
  # Combine both plots
  grid.arrange(peaks_macs2_plot, top_features_plot, ncol = 2)
  
} else {
  cat("Subsampled peak data is empty or NULL. Please check the data.")
}



