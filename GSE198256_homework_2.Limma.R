#####################################################
#This is part TWO of homework#2.
#Here, we will apply learning concepts from the class.

#Here, we use DESLimma-trend and vroom.
#Also, we do impllement ORA & GSEA to gain 
#mechanistical interpretation for DEG.
#####################################################
#######################
#1 DATA LOADING-RNA-Seq
#######################
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
#using the defalut values [p-value  = 0.05, LFC = 0, top 10 genes]...
res_limma_trend <- topTable(fit, coef=ncol(design), sort.by="B" ,number=100)
#voom
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
#using the defalut values [p-value  = 0.05, LFC = 0, top 10 genes]...
res_voom <- topTable(fit, coef=ncol(design), , sort.by="B" ,number=100)

```

## ACTIVITY 1:

-   How would you compare the results between voom and trend?
-   Is it required to run more than analysis?
-   What exactly are we asking with this differential expression?

```{r ACTIVITY 1}




```
#how many shared genes between limma & voom in our dataset?
#using the top 100 genes when the results are sorted by p-value or log-odds of differential expression...
#Limma
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
#set parameters for sorting and genes selection...
res_limma_trend <- topTable(fit, coef=ncol(design), sort.by="B" ,number=100)
#voom
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
#set parameters for sorting and genes selection...
res_voom <- topTable(fit, coef=ncol(design), , sort.by="B" ,number=100)
genes_res_limma_trend <- rownames(res_limma_trend)
genes_res_limma_trend
genes_res_voom <- rownames(res_voom)
genes_res_voom

match <- 0
 
for (gene in genes_res_limma_trend) { 
    if (gene %in% genes_res_voom) { 
        match <- match +1 
        } 
    }
match_percentage <- match / length(genes_res_limma_trend) *100
print(paste("% of matches is", match_percentage))
############################
#3 BIOLOGICAL INTERPRETATION
############################
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(magrittr)
# Add more contrasts

v <- voom(dge, design, plot=TRUE)
colnames(design) <- c("Intercept","Covid196Mo","Covid19AI","Healthy")
fit <- lmFit(v, design)

contrast.matrix <- makeContrasts(Covid19AI-Healthy, Healthy, 
                                 Covid196Mo-Healthy,    
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_intercept_top <- topTable(fit2) 
res_Covid196Mo_top <- topTable(fit2,coef=1) 
res_Covid19AI_top <- topTable(fit2,coef=2) 
res_Healthy_top <- topTable(fit2,coef=3) 

```

## ORA and Gene Set Enrichment analysis.

-   What do we need to do the analysis?
-   What are the tools required?
-   

```{r Prepare ORA and GSEA}

# If we want to shift annotations:
ENSEMBL_vector <- mapIds(
  org.Hs.eg.db,
  keys = rownames(GSE198256_count),
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)


# We would like a data frame we can join to the differential expression stats
gene_key_df <- data.frame(
  ensembl_id = ENSEMBL_vector,
  entrez_id = names(ENSEMBL_vector),
  stringsAsFactors = FALSE
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(ensembl_id))

######
#ORA
######
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
# Step 5: Visualize / explore

enrich_plot <- enrichplot::dotplot(kegg_ora_results)
enrich_plot

upset_plot <- enrichplot::upsetplot(kegg_ora_results)
upset_plot

# Step 6: EXERCISE: alternatives to KEGG?


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
  geneSetID = "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  title = "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
  color.line = "#0d76ff"
)
most_positive_nes_plot

gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "KEGG_SPLICEOSOME",
  title = "KEGG_SPLICEOSOME",
  color.line = "#0d76ff"
)
most_negative_nes_plot

# Step 7: EXERCISE: compare GSEA vs ORA?
#by associated pathways...
#1 get ORA & GSEA associated pathways...
ora_enriched_pathways <- kegg_ora_results[, c("ID", "Description", "pvalue", "p.adjust")]
gsea_enriched_gene_sets <- gsea_results[, c("ID", "Description", "pvalue", "p.adjust")]
#2 how many are shared ?

for (i in 1:nrow(ora_enriched_pathways)) {
  if (ora_enriched_pathways$pvalue[i] < 0.05) {
    # Check if pathway is also enriched in the longer list with p.adjust < 0.01
    matching_pathway <- gsea_enriched_gene_sets$Description[gsea_enriched_gene_sets$ID == ora_enriched_pathways$ID[i] & gsea_enriched_gene_sets$pvalue < 0.05]
    if (length(matching_pathway) > 0) {
      print(paste("Pathway found in both lists:", ora_enriched_pathways$Description[i]))
    }
  }
}

###############
#MORE DETAILED CODE FOR THE COMPARISON IS IN ....
###############
#delete this...
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

#change the code below
logCPM <- cpm(dge, log=TRUE, prior.count=3) #v <- voom(dge, design, plot=TRUE)
colnames(design) <- c("Intercept","Covid196Mo","Covid19AI","Healthy")
fit <- lmFit(logCPM, design) #fit <- lmFit(v, design)

contrast.matrix <- makeContrasts(Covid19AI-Healthy, Healthy, 
                                 Covid196Mo-Healthy,    
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res_intercept_top <- topTable(fit2) 
res_Covid196Mo_top <- topTable(fit2,coef=1) 
res_Covid19AI_top <- topTable(fit2,coef=2) 
res_Healthy_top <- topTable(fit2,coef=3) 
#V
ENSEMBL_vector <- mapIds(
  org.Hs.eg.db,
  keys = rownames(GSE198256_count),
  keytype = "ENTREZID",
  column = "ENSEMBL",
  multiVals = "first"
)


#STOP