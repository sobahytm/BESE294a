#####################################################
#This is part TWO of homework#2.

#Here, we use Limma-trend with ORA & GSEA.

#This would allow us to compare :
#limma-trend + ORA V. voom + ORA
#limma-trend + GSEA + voom + GSEA
#####################################################
# Read data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Read Meta data
library(GEOquery)
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Group <- Meta_GSE198256[,c("disease state:ch1")]

dim(GSE198256_count)
Group
#######################
#2 DATA ANALYSIS-LIMMA
#######################
# set DGE class
require(limma)
require(edgeR)
dge <- DGEList(counts=GSE198256_count)
# Make sure on the metadata
rownames(Meta_GSE198256)==colnames(GSE198256_count)
Group[Group=="Covid19: Acute infection"] <- "Covid19AI"
Group[Group=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
Group[Group=="Covid19: Recovery 6Mo"] <- "Covid196Mo"
design <- model.matrix(~ Group )

# Filter
keep <- filterByExpr(dge, design=design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Normalization
dge <- calcNormFactors(dge)
#limma-trend
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(magrittr)

logCPM <- cpm(dge, log=TRUE, prior.count=3)

colnames(design) <- c("Intercept","Covid196Mo","Covid19AI","Healthy")
fit <- lmFit(logCPM, design)

contrast.matrix <- makeContrasts(Covid19AI-Healthy, Healthy, 
                                 Covid196Mo-Healthy,    
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_intercept_top <- topTable(fit2) 
res_Covid196Mo_top <- topTable(fit2,coef=1) 
res_Covid19AI_top <- topTable(fit2,coef=2) 
res_Healthy_top <- topTable(fit2,coef=3) 

ENSEMBL_vector <- mapIds(
  org.Hs.eg.db,
  keys = rownames(GSE198256_count),
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)

gene_key_df <- data.frame(
  ensembl_id = ENSEMBL_vector,
  entrez_id = names(ENSEMBL_vector),
  stringsAsFactors = FALSE
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(ensembl_id))

######
#GSEA
######
# Step 1: determine genes of interest.
diff_table_all <- topTable(fit2,coef=1,p.value=1,number=nrow(logCPM)) 
# Step 2: determine background.

# Step 3: Determine gene sets.

msigdbr_species()
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
head(hs_msigdb_df)

hs_kegg_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", 
    gs_subcat == "CP:KEGG" 
  )

# Step 4: conduct GSEA

list_ordered <- diff_table_all[,"B"]
names(list_ordered) <- rownames(diff_table_all)

gsea_results <- GSEA(
  geneList = list_ordered, 
  minGSSize = 25, 
  maxGSSize = 500, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hs_kegg_df,
    gs_name,
    human_entrez_gene
  )
)
# Step 5: Visualize / explore

head(gsea_results@result)
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "KEGG_RIBOSOME",
  title = "KEGG_RIBOSOME",
  color.line = "#0d76ff"
)
most_positive_nes_plot

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "KEGG_RIBOSOME",
  title = "KEGG_RIBOSOME",
  color.line = "#0d76ff"
)
most_negative_nes_plot

#ORA
# Step 1: determine genes of interest.

diff_table <- topTable(fit2,coef=1,p.value=0.01,number=10000) 

genes_dif<- rownames(diff_table )
# Step 2: determine background.

background_set <- unique(rownames(logCPM))
# Step 3: Determine gene sets.

msigdbr_species()
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
head(hs_msigdb_df)

hs_kegg_df <- hs_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat == "CP:KEGG"# This is because we only want KEGG pathways 
  )
# Step 4: conduct ORA.
kegg_ora_results <- enricher(
  gene = genes_dif, # A vector of your genes of interest
  pvalueCutoff = 0.1, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = background_set, # A vector containing your background set genes
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_kegg_df,
    gs_name,
    human_entrez_gene
  )
)
enrich_plot <- enrichplot::dotplot(kegg_ora_results)
enrich_plot

upset_plot <- enrichplot::upsetplot(kegg_ora_results)
upset_plot