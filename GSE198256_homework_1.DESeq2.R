#####################################################
#This is part ONE of homework#2.
#Here, we revise the class activity & use DESeq2
#####################################################
#######################
#1 DATA LOADING-RNA-Seq
#######################
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1);
#######################
#1 DATA LOADING-META
#######################
library(GEOquery)
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]
Factors_GSE198256 <- Meta_GSE198256[,c("disease state:ch1")]
################################
#1 DATA LOADING-GENE ANNOTATIONS
################################
annotgene <- read.csv("/Users/sobahytm/BESE_394A/class_1/mart_export.txt",sep="\t",header = T)
annotgene <- annotgene[annotgene$Chromosome %in% c(as.character(1:22) ,"X","Y"),]
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)
#any duplications!
rownames(annotgene) <- annotgene$Entrezgene
annotgene_filt <- annotgene[!duplicated(annotgene$Entrezgene),]

sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)
sum(annotgene_filt$Entrezgene %in% rownames(GSE198256_count))

rownames(annotgene_filt) <- as.character(annotgene_filt$Entrezgene)
sum(as.character(rownames(annotgene_filt)) %in% rownames(GSE198256_count))
#only genes with full annotations!
GSE198256_count_filt <- GSE198256_count[rownames(GSE198256_count) %in% rownames(annotgene_filt),]
annotgene_ord <- annotgene_filt[rownames(GSE198256_count_filt ),]
sum(rownames(annotgene_ord)==rownames(GSE198256_count_filt))
#########################
#2 DATA INSPECTION-NOISEQ
#########################
library(NOISeq)

Factors_GSE198256 <- data.frame(Meta_GSE198256 [ colnames(GSE198256_count_filt),c("disease state:ch1")])
colnames(Factors_GSE198256)[1]<- "Group"
#get the data in the proper format for NOISEQ readData function...
lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
names(lengthuse) <- rownames(annotgene_ord)
gc <- annotgene_ord$GC
names(gc) <- rownames(annotgene_ord)
biotype <-annotgene_ord$type
names(biotype) <- rownames(annotgene_ord)
chromosome <- annotgene_ord[,c("Chromosome","start","end")]
#process the data with the package...
data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)
#check for data bias prioir normaliztion!
#2.1 information about the abunces of the different RNA molecules types..
myexplodata <- dat(data_NOISEQ, type = "biodetection")#specify the type...
explo.plot(myexplodata, plottype = "persample")
#2.2 samples comparisions...
#all counts/expersions (CPM)....
mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = "protein_coding", samples = 1:34, plottype = "barplot")
#2.3 features gain relative to sequencing depth...
mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:34)
#2.4testing for length bais... 
mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")
#2.5testing for GC% bais... 
myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")
#2.6testing for RNA compistion bais... 
myRNAc = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(myRNAc,samples = 1:12)
#testing for the study cohorts segragation based on non-normalized counts...
myPCA = dat(data_NOISEQ, type = "PCA")
explo.plot(myPCA, factor = "Group")
############################
#3 DATA NORMALIZATION-NOISEQ
############################
library(Biobase)
myRPKM_matrix <- rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myRPKM <- ExpressionSet(assayData = myRPKM_matrix)
phenoData(myRPKM) <- phenoData(data_NOISEQ)
myPCA = dat(myRPKM, type = "PCA")
explo.plot(myPCA, factor = "Group")
#UQUA
myUQUA_matrix <- rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA <- ExpressionSet(assayData = myUQUA_matrix)
phenoData(myUQUA) <- phenoData(data_NOISEQ)
myPCA = dat(myUQUA, type = "PCA")
explo.plot(myPCA, factor = "Group")
#TMM
myTMM_matrix <- rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myTMM <- ExpressionSet(assayData = myTMM_matrix)
phenoData(myTMM) <- phenoData(data_NOISEQ)
myPCA = dat(myTMM, type = "PCA")
explo.plot(myPCA, factor = "Group")
#Recheck for data bias after normalization...
#add the length, gc feactures..
myTMM <- addData(myTMM, length = lengthuse)
featureData(myTMM)$GC <- gc

mylengthbias = dat(myTMM, factor = "Group", type = "lengthbias")
#
myGCbias = dat(myTMM, factor = "Group", type = "GCbias")
#explo.plot(myGCbias, samples = NULL, toplot = "global")
myRNAc = dat(myTMM, type = "cd", norm = FALSE, refColumn = 1)
#####################################################
#3 DATA NORMALIZATION-NOISEQ-YOU MAY REMOVE THIS PART
#####################################################

#############################################
#3 DATA NORMALIZATION & DIFF. ANALYSIS-DESEQ2
#############################################
library(Biobase)
library(DESeq2)
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="Covid19: Acute infection"] <- "Covid19AI"
pDataUSE[pDataUSE=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
pDataUSE[pDataUSE=="Covid19: Recovery 6Mo"] <- "Covid196Mo"
pDataUSE[,1] <- as.factor(pDataUSE[,1])
#read the data into an object that fits the DESeq2 requirments!
GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ -1 + Group)#be mindful of the used model 
#resultsNames(GSE198256_DESeq2)
GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ Group)
#resultsNames(GSE198256_DESeq2)
#removing low or zero counts genes in at least six or more samples...
smallestGroupSize <- 6
keep <- rowSums(counts(GSE198256_DESeq2) >= 10) >= smallestGroupSize
GSE198256_DESeq2_F <- GSE198256_DESeq2[keep,]
# DESeq2 nomrmalize + diff. expression...
GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
GSE198256_res <- results(GSE198256_DESeq2_F)
GSE198256_res
resultsNames(GSE198256_DESeq2_F)

plotMA(GSE198256_res, ylim=c(-2,2))
#do the shrink technique or not based on the plot (look for genes woth low counts or high dispersion )?
#res_shrink <- lfcShrink(GSE198256_DESeq2_F,coef=c("Group_Covid196Mo_vs_Covid193Mo"))
#plotMA(res_shrink, ylim=c(-2,2))

###########################################
#4 PAIRED ANALYSIS & VISUALIZATION-PHEATMAP
###########################################
library("pheatmap")
#paired-analysis (usinf the default mmodel Wald test)...
GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
res <- results(GSE198256_DESeq2_F, contrast=c("Group","Healthy","Covid19AI"))
vsd <- vst(GSE198256_DESeq2_F, blind=FALSE)

valid_indices <- complete.cases(res$log2FoldChange, res$padj)#check for Nan..
filtered_results <- res[valid_indices & abs(res$log2FoldChange) >= 2 & res$padj <= 0.05, ]
top_genes <- filtered_results[order(filtered_results$padj), ][1:30, ]

#Can't sort the heatmap by groups --> changing samples names to group names....
top_genes_rownames <- rownames(top_genes)
selected_indices <- match(top_genes_rownames, rownames(GSE198256_DESeq2_F))
groups <- colData(GSE198256_DESeq2_F)$Group
order_samples <- order(groups)
expression_matrix_ordered <- assay(vsd)[top_genes_rownames, order_samples]
modified_sample_names <- paste(groups[order_samples], 1:sum(groups == "Covid19AI"), sep = "")
colnames(expression_matrix_ordered) <- modified_sample_names
#generate the heatmap...
pheatmap(expression_matrix_ordered,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = FALSE)
#############################
#5 BIOLOGICAL INTERPRETATION
#############################
#There is no gene enrichment analysis here!
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
#5.1 add NCBI gene name and id as column...
ensembl_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl_mart)
attributes <- subset(attributes, name %in% c("entrezgene_id", "external_gene_name"))
gene_mapping <- getBM(attributes = c("entrezgene_id", "external_gene_name"), mart = ensembl_mart)
filtered_results_df <- as.data.frame(filtered_results)
filtered_results_df$NCBI_gene_id <- rownames(filtered_results_df)
filtered_results_with_names <- merge(filtered_results_df, gene_mapping, by.x = "NCBI_gene_id", by.y = "entrezgene_id", all.x = TRUE)
#5.2 prepare the input for clusterProfiler functions...
gene_symbols <- filtered_results_with_names$external_gene_name 
orgdb <- org.Hs.eg.db 
entrez_ids <- select(orgdb, keys = gene_symbols, keytype = "SYMBOL", columns = "SYMBOL")           
colnames(entrez_ids)

ego <- enrichGO(gene = entrez_ids$SYMBOL, 
                OrgDb = org.Hs.eg.db, 
                keyType = "SYMBOL", 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05)
dotplot(ego, showCategory = 20) # show only the top 20!