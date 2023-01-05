### 11-Jul-2019

# Comparative RNAseq analysis of DADA2 and PAN patients using DEseq2

library(dplyr)
library(DESeq2)
library(pheatmap)
library(vsn)
library(hexbin)
library(RColorBrewer)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(tibble)

# Set wd and load all count files 
setwd("/Volumes/users/RNAseq/DADA2_PAN/PedVas/HTSeq_files")

# count files for patients WITH DADA2
CA5002V4 <- read.table("CA5002V4.count", col.names = c("gene_ID", "CA5002V4"))
CA5002V5 <- read.table("CA5002V5.count", col.names = c("gene_ID", "CA5002V5"))
EE1001V1 <- read.table("EE1001V1.count", col.names = c("gene_ID", "EE1001V1"))
EE1001V2 <- read.table("EE1001V2.count", col.names = c("gene_ID", "EE1001V2"))
EE1001V3 <- read.table("EE1001V3.count", col.names = c("gene_ID", "EE1001V3"))
SE1003V5 <- read.table("36-15.count", col.names = c("gene_ID", "SE1003V5"))
SK1003V1 <- read.table("SK1003V1.count", col.names = c("gene_ID", "SK1003V1"))
SK1003V2 <- read.table("SK1003V2.count", col.names = c("gene_ID", "SK1003V2"))
TO5009V4 <- read.table("TO5009V4.count", col.names = c("gene_ID", "TO5009V4"))
TO5009V5 <- read.table("TO5009V5.count", col.names = c("gene_ID", "TO5009V5"))
VA5004V4 <- read.table("VA5004V4.count", col.names = c("gene_ID", "VA5004V4"))
VA5004V5 <- read.table("VA5004V5.count", col.names = c("gene_ID", "VA5004V5"))

# count files for patients with PAN or cPAN
AK1003V3 <- read.table("AK1003V3.count", col.names = c("gene_ID", "AK1003V3"))
BA1001V3 <- read.table("BA1001V3.count", col.names = c("gene_ID", "BA1001V3"))
BR1004V5 <- read.table("BR1004V5.count", col.names = c("gene_ID", "BR1004V5"))
CA1001V2 <- read.table("CA1001V2.count", col.names = c("gene_ID", "CA1001V2"))
CA1010V1 <- read.table("CA1010V1.count", col.names = c("gene_ID", "CA1010V1"))
GA1002V1 <- read.table("GA1002V1.count", col.names = c("gene_ID", "GA1002V1"))
GA1002V3 <- read.table("GA1002V3.count", col.names = c("gene_ID", "GA1002V3"))
VA1006V1 <- read.table("VA1006V1.count", col.names = c("gene_ID", "VA1006V1"))
VA1006V5 <- read.table("VA1006V5.count", col.names = c("gene_ID", "VA1006V5"))

setwd("/Volumes/users/RNAseq/ReSeq_201903/Counts") 
CA1010V2 <- read.table("Pedvas_308.count", col.names = c("gene_ID", "CA1010V2"))

# Make count table by munging individual count files together
counts <- left_join(EE1001V1, EE1001V2, by = "gene_ID")
counts <- left_join(counts, SE1003V5, by = "gene_ID")
counts <- left_join(counts, SK1003V1, by = "gene_ID")
counts <- left_join(counts, SK1003V2, by = "gene_ID")
counts <- left_join(counts, TO5009V4, by = "gene_ID")
counts <- left_join(counts, TO5009V5, by = "gene_ID")
counts <- left_join(counts, VA5004V4, by = "gene_ID")
counts <- left_join(counts, VA5004V5, by = "gene_ID")
counts <- left_join(counts, AK1003V3, by = "gene_ID")
counts <- left_join(counts, BA1001V3, by = "gene_ID")
counts <- left_join(counts, BR1004V5, by = "gene_ID")
counts <- left_join(counts, CA1001V2, by = "gene_ID")
counts <- left_join(counts, CA1010V1, by = "gene_ID")
counts <- left_join(counts, CA1010V2, by = "gene_ID")
counts <- left_join(counts, GA1002V1, by = "gene_ID")
counts <- left_join(counts, GA1002V3, by = "gene_ID")
counts <- left_join(counts, VA1006V1, by = "gene_ID")
counts <- left_join(counts, VA1006V5, by = "gene_ID")
counts <- left_join(counts, CA5002V4, by = "gene_ID")
counts <- left_join(counts, CA5002V5, by = "gene_ID")

counts <- column_to_rownames(counts, var = "gene_ID") # make gene ID column row names

# read in meta-data
setwd("/Volumes/users/RNAseq/DADA2_PAN")
meta <- read.csv("PAN_DADA2_metadata.csv")
meta <- tibble::column_to_rownames(meta, var = "lab_id")

disease_status <- ifelse(meta$PVAS_score>1, "active", "not_active") # make vector to denote active or inactive disease
meta$status <- disease_status

# DADA2 vs PAN paired analysis --------------------------------------------

# keep only DADA2 and PAN samples with paired data (1 active + 1 inactive count file)

counts_paired <- counts[ , c("EE1001V1", "EE1001V2", "SK1003V1","SK1003V2", "TO5009V4", "TO5009V5", "VA5004V4", "VA5004V5",
                             "GA1002V1", "GA1002V3", "VA1006V1", "VA1006V5", "CA1010V1", "CA1010V2")]

# convert NAs to 0s
counts_paired[is.na(counts_paired)] <- 0

meta_paired <- meta[c("EE1001V1", "EE1001V2", "SK1003V1","SK1003V2", "TO5009V4", "TO5009V5", "VA5004V4", "VA5004V5",
                      "GA1002V1", "GA1002V3", "VA1006V1", "VA1006V5", "CA1010V1", "CA1010V2"), ]


Patient <- gl(4, 2, length = 14)
meta_paired$patient <- Patient

# Make DESeqDataSet
colData_paired <- meta_paired[ , c("status", "mut_ADA2", "patient")]
colData_paired$status <- as.factor(colData_paired$status)

colData_paired$status <- relevel(colData_paired$status, ref = "not_active") # relevel factors so not_active is the reference

dds_paired <- DESeqDataSetFromMatrix(countData = counts_paired,
                                     colData = colData_paired,
                                     design = ~ mut_ADA2 + patient + mut_ADA2:status)

# filter weakly expressed and noninformative features
keep <- rowSums(counts(dds_paired)) >= 10
dds_paired <- dds_paired[keep,]

dds_paired <- DESeq(dds_paired)
resultsNames(dds_paired)
res_paired <- results(dds_paired, contrast = list("mut_ADA2negative.statusactive", "mut_ADA2positive.statusactive"))

resOrdered_paired <- res_paired[order(res_paired$pvalue),]

setwd("/Users/kristen/.Volumes/kgibson/RNAseq/DADA2_PAN/DESeq2")
# write.csv(as.data.frame(resOrdered_paired, ), "deg_paired_DESeq2.csv")

plotMA(res_paired, ylim=c(-2,2))   

# log fold change shrinkage for visualization and ranking
resLFC_paired <- lfcShrink(dds_paired, coef="mut_ADA2_positive_vs_negative", type="apeglm")

# MA plots
plotMA(res_paired, ylim = c(-2, 2))
plotMA(resLFC_paired, ylim = c(-2, 2))

### data transformation and visualization
vsd_paired <- vst(dds_paired, blind=FALSE) # variance stabilized normalization

ntd_paired <- normTransform(dds_paired) # normalized counts transformation -- log2(n+1)
meanSdPlot(assay(ntd_paired)) # plot norm-transformed counts
meanSdPlot(assay(vsd_paired)) # plot vsd transformed counts 

# heatmap of count matrix - using vsd normalized count data
select_paired <- order(rowMeans(counts(dds_paired,normalized=TRUE)),
                       decreasing=TRUE)[1:5]
df_paired <- as.data.frame(colData(dds_paired)[,c("status", "mut_ADA2")])
pheatmap(assay(vsd_paired), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_paired)

# heatmap of sample-sample distances
sampleDists_paired <- dist(t(assay(vsd_paired)))
sampleDistMatrix_paired <- as.matrix(sampleDists_paired)
rownames(sampleDistMatrix_paired) <- paste(vsd_paired$mut_ADA2)
colnames(sampleDistMatrix_paired) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_paired,
         clustering_distance_rows = sampleDists_paired,
         clustering_distance_cols = sampleDists_paired,
         col = colors)

# PCA plot of samples
plotPCA(vsd_paired, intgroup=c("mut_ADA2", "status"))
plotPCA(ntd_paired, intgroup = c("mut_ADA2", "status"))

boxplot(log10(assays(dds_paired)[["cooks"]]), range=0, las=2)


# Individual gene plots

plotCounts(dds_paired, "ENSG00000124614", intgroup = c("mut_ADA2", "status")) # RPS10
HSPB1_paired <- plotCounts(dds_paired, "ENSG00000106211", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # HSPB1 *****
C3_paired <- plotCounts(dds_paired, "ENSG00000125730", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # C3 *****
VHL_paired <- plotCounts(dds_paired, "ENSG00000134086", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # VHL *****
plotCounts(dds_paired, "ENSG00000096060", intgroup = c("mut_ADA2", "status")) # FKBP5 
TUBGCP3_paired <- plotCounts(dds_paired, "ENSG00000126216", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # TUBGCP3 *****
AURKB_paired <- plotCounts(dds_paired, "ENSG00000178999", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # AURKB **
plotCounts(dds_paired, "ENSG00000229921", intgroup = c("mut_ADA2", "status")) # KIF25-AF1
plotCounts(dds_paired, "ENSG00000149256", intgroup = c("mut_ADA2", "status")) # TENM4
plotCounts(dds_paired, "ENSG00000152463", intgroup = c("mut_ADA2", "status")) # OLAH
ADAR1_paired <- plotCounts(dds_paired, "ENSG00000160710", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # ADAR1 **
plotCounts(dds_paired, "ENSG00000197381", intgroup = c("mut_ADA2", "status")) # ADAR2 
plotCounts(dds_paired, "ENSG00000185736", intgroup = c("mut_ADA2", "status")) # ADAR3
plotCounts(dds_paired, "ENSG00000055332", intgroup = c("mut_ADA2", "status")) # EIF2AK2
MYO7A_paired <- plotCounts(dds_paired, "ENSG00000137474", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # MYO7A
plotCounts(dds_paired, "ENSG00000102145", intgroup = c("mut_ADA2", "status")) # GATA1
DNMT1_paired <- plotCounts(dds_paired, "ENSG00000130816", intgroup = c("mut_ADA2", "status"), returnData = TRUE) # DNMT1

plotCounts(dds_paired, "ENSG00000112715", intgroup = c("mut_ADA2", "status")) #VEGF
plotCounts(dds_paired, "ENSG00000117394", intgroup = c("mut_ADA2", "status")) #GLUT1
plotCounts(dds_paired, "ENSG00000007171", intgroup = c("mut_ADA2", "status")) #NOS2
plotCounts(dds_paired, "ENSG00000109320", intgroup = c("mut_ADA2", "status")) #NF-KB
plotCounts(dds_paired, "ENSG00000160255", intgroup = c("mut_ADA2", "status")) #ITGB2





library("ggplot2")
library("ggpubr")

groups_paired <- rownames_to_column(meta_paired, "lab_ID")
groups_paired <- mutate(groups_paired, group = ifelse(c(mut_ADA2 == "positive" & status == "active"), "DADA2_active", 
                                                      ifelse(c(mut_ADA2 == "positive" & status == "not_active"), "DADA2_inactive",
                                                             ifelse(c(mut_ADA2 == "negative" & status == "active"), "PAN_active",
                                                                    ifelse(c(mut_ADA2 == "negative" & status == "not_active"), "PAN_inactive", "FAIL")))))
groups_paired <- groups_paired[ , c("lab_ID", "study_id", "group")]

HSPB1_paired <- rownames_to_column(HSPB1_paired, "lab_ID")
HSPB1_paired <- left_join(HSPB1_paired, groups_paired, by = "lab_ID")
ggpaired(HSPB1_paired, x = "group", y = "count", id = "study_id", title = "HSPB1_norm_counts",
         xlab = "", ylab = "normalized counts", label = "study_id")

C3_paired <- rownames_to_column(C3_paired, "lab_ID")
C3_paired <- left_join(C3_paired, groups_paired, by = "lab_ID")
ggpaired(C3_paired, x = "group", y = "count", id = "study_id", title = "C3_norm_counts", 
         xlab = "", ylab = "normalized counts", label = "study_id")

VHL_paired <- rownames_to_column(VHL_paired, "lab_ID")
VHL_paired <- left_join(VHL_paired, groups_paired, by = "lab_ID")
ggpaired(VHL_paired, x = "group", y = "count", id = "study_id", title = "VHL_norm_counts",
         xlab = "", ylab = "normalized counts", label = "study_id")

TUBGCP3_paired <- rownames_to_column(TUBGCP3_paired, "lab_ID")
TUBGCP3_paired <- left_join(TUBGCP3_paired, groups_paired, by = "lab_ID")
ggpaired(TUBGCP3_paired, x = "group", y = "count", id = "study_id", title = "TUBGCP3_norm_counts",
         label = "study_id")

AURKB_paired <- rownames_to_column(AURKB_paired, "lab_ID")
AURKB_paired <- left_join(AURKB_paired, groups_paired, by = "lab_ID")
ggpaired(AURKB_paired, x = "group", y = "count", id = "study_id", label = "study_id")

ADAR1_paired <- rownames_to_column(ADAR1_paired, "lab_ID")
ADAR1_paired <- left_join(ADAR1_paired, groups_paired, by = "lab_ID")
ggpaired(ADAR1_paired, x = "group", y = "count", id = "study_id", title = "ADAR1_norm_counts",
         xlab = "", ylab = "normalized counts", label = "study_id")

MYO7A_paired <- rownames_to_column(MYO7A_paired, "lab_ID")
MYO7A_paired <- left_join(MYO7A_paired, groups_paired, by = "lab_ID")
ggpaired(MYO7A_paired, x = "group", y = "count", id = "study_id")


DNMT1_paired <- rownames_to_column(DNMT1_paired, "lab_ID")
DNMT1_paired <- left_join(DNMT1_paired, groups_paired, by = "lab_ID")
ggpaired(DNMT1_paired, x = "group", y = "count", id = "study_id", title = "DNMT1_norm_counts", 
         xlab = "", ylab = "normalized counts", label = "study_id") ###
# Active analysis ---------------------------------------------------------

# select active count files
counts_active <- counts[ , c("TO5009V4", "VA5004V4", "SK1003V1", "EE1001V1", 
                             "GA1002V1", "VA1006V1", "CA5002V4", "CA1010V1")]        

# select active metadata
meta_active <- meta[c("TO5009V4", "VA5004V4", "SK1003V1", "EE1001V1", 
                      "GA1002V1", "VA1006V1", "CA5002V4", "CA1010V1"),]        

colData_active <- meta_active[ , c("study_id", "mut_ADA2", "status")]

# Make DESeqDataSet
meta_active$mut_ADA2 <- as.factor(meta_active$mut_ADA2)

dds_active <- DESeqDataSetFromMatrix(countData = counts_active,
                                     colData = colData_active,
                                     design = ~ mut_ADA2)

# filter weakly expressed and noninformative features
keep_active <- rowSums(counts(dds_active)) >= 10
dds_active <- dds_active[keep_active,]

keep_active_2 <- rowSums( counts(dds_active, normalized=FALSE) >= 10 ) >= 3
dds_active <- dds_active[keep_active_2,]

# DE analysis
dds_active <- DESeq(dds_active)
res_active <- results(dds_active)

# write csv of results
# setwd("/Users/kristen/.Volumes/kgibson/RNAseq/DADA2_PAN/DESeq2")
# write.csv(as.data.frame(res_active), "deg_active.csv", row.names = TRUE )

# log fold change shrinkage for visualization and ranking
resLFC_active <- lfcShrink(dds_active, coef="mut_ADA2_positive_vs_negative", type="apeglm")

# MA plots
plotMA(res_active, ylim = c(-2, 2))
plotMA(resLFC_active, ylim = c(-2, 2))

### data transformation and visualization
vsd_active <- vst(dds_active, blind=FALSE) # variance stabilized normalization

ntd_active <- normTransform(dds_active) # normalized counts transformation -- log2(n+1)
meanSdPlot(assay(ntd_active)) # plot norm-transformed counts
meanSdPlot(assay(vsd_active)) # plot vsd transformed counts 

# heatmap of count matrix - using vsd normalized count data
select_active <- order(rowMeans(counts(dds_active,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df_active <- as.data.frame(colData(dds_active)[,c("study_id", "mut_ADA2")])
pheatmap(assay(vsd_active), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_active)

# heatmap of sample-sample distances
sampleDists_active <- dist(t(assay(vsd_active)))
sampleDistMatrix_active <- as.matrix(sampleDists_active)
rownames(sampleDistMatrix_active) <- paste(vsd_active$mut_ADA2)
colnames(sampleDistMatrix_active) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_active,
         clustering_distance_rows = sampleDists_active,
         clustering_distance_cols = sampleDists_active,
         col = colors)

# PCA plot of samples
plotPCA(vsd_active, intgroup=c("mut_ADA2"))
plotPCA(ntd_active, intgroup = "mut_ADA2")

boxplot(log10(assays(dds_active)[["cooks"]]), range=0, las=2)

plotCounts(dds_active, "ENSG00000211937", intgroup = c("mut_ADA2")) # IGHV2-5
plotCounts(dds_active, "ENSG00000211893", intgroup = c("mut_ADA2")) # IGHG2
plotCounts(dds_active, "ENSG00000211896", intgroup = c("mut_ADA2")) # 
plotCounts(dds_active, "ENSG00000174944", intgroup = c("mut_ADA2")) # 
plotCounts(dds_active, "ENSG00000079215", intgroup = c("mut_ADA2")) # 
plotCounts(dds_active, "ENSG00000020577", intgroup = c("mut_ADA2")) # 
plotCounts(dds_active, "ENSG00000182534", intgroup = c("mut_ADA2")) # MXRA7 *** NetwrokAnalyst
plotCounts(dds_active, "ENSG00000132465", intgroup = c("mut_ADA2")) # IGJ *** NetworkAnalyst




# Inactive analysis -------------------------------------------------------

# select inactive count files
counts_inactive <- counts[ , c("TO5009V5", "VA5004V5", "EE1001V2","SK1003V2",
                             "BR1004V5", "GA1002V3", "VA1006V5")]        

# select inactive metadata
meta_inactive <- meta[c("TO5009V5", "VA5004V5", "EE1001V2","SK1003V2",
                      "BR1004V5", "GA1002V3", "VA1006V5"),]        

colData_inactive <- meta_inactive[ , c("study_id", "MD_diagnosis", "status")]

# Make DESeqDataSet
meta_inactive$mut_ADA2 <- as.factor(meta_inactive$MD_diagnosis)

dds_inactive <- DESeqDataSetFromMatrix(countData = counts_inactive,
                                     colData = colData_inactive,
                                     design = ~ MD_diagnosis)

# filter weakly expressed and noninformative features
keep_inactive <- rowSums(counts(dds_inactive)) >= 10
dds_inactive <- dds_inactive[keep_inactive,]

# DE analysis
dds_inactive <- DESeq(dds_inactive)
res_inactive <- results(dds_inactive)
res_inactive <- as.data.frame(res_inactive)

# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- res_inactive$ensembleID
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))


res_inactive$ensembleID <- row.names(res_inactive)
res_inactive <- left_join(res_inactive, geneIDs, by = c("ensembleID" = "GENEID"))

res_inactive_sub <- subset(res_inactive, padj < 0.1)

plotCounts(dds_inactive, "ENSG00000126709", intgroup = "MD_diagnosis") # IFI6
plotCounts(dds_inactive, "ENSG00000187608", intgroup = "MD_diagnosis") # ISG15
plotCounts(dds_inactive, "ENSG00000119917", intgroup = "MD_diagnosis") # IFIT3
plotCounts(dds_inactive, "ENSG00000185745", intgroup = "MD_diagnosis") # IFIT1
plotCounts(dds_inactive, "ENSG00000088827", intgroup = "MD_diagnosis") # SIGLEC1
plotCounts(dds_inactive, "ENSG00000125730", intgroup = "mut_ADA2") # C3
plotCounts(dds_inactive, "ENSG00000138646", intgroup = "mut_ADA2") # HERC3
plotCounts(dds_inactive, "ENSG00000126456", intgroup = "mut_ADA2") # IRF3
plotCounts(dds_inactive, "ENSG00000139998", intgroup = "mut_ADA2") # RAB15
plotCounts(dds_inactive, "ENSG00000107201", intgroup = "mut_ADA2") # RIG1***
plotCounts(dds_inactive, "ENSG00000109320", intgroup = "mut_ADA2") # NFKB1
plotCounts(dds_inactive, "ENSG00000123610", intgroup = "mut_ADA2") # TNFAIP6
plotCounts(dds_inactive, "ENSG00000134321", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_inactive, "ENSG00000185507", intgroup = "MD_diagnosis") # IRF7


plotCounts(dds_inactive, "ENSG00000112715", intgroup = "mut_ADA2") # VEGF
plotCounts(dds_inactive, "ENSG00000117394", intgroup = "mut_ADA2") # GLUT1
plotCounts(dds_inactive, "ENSG00000109320", intgroup = "mut_ADA2") # NF-KB
plotCounts(dds_inactive, "ENSG00000160255", intgroup = "mut_ADA2") # ITGB2



plotCounts(dds_active, "ENSG00000174944", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_active, "ENSG00000079215", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_active, "ENSG00000164120", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_active, "ENSG00000103723", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_active, "ENSG00000162772", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_active, "ENSG00000225492", intgroup = "mut_ADA2") # RSAD2
plotCounts(dds_active, "ENSG00000134463", intgroup = "mut_ADA2") # RSAD2


library(EnsDb.Hsapiens.v79)
library(dplyr)

IFN_sig_counts <- as.data.frame(assays(dds_inactive))
IFN_sig_counts$gene <- rownames(IFN_sig_counts)
IFN_sig_counts <- left_join(IFN_sig_counts, geneIDs, by = c("gene" = "GENEID"))

#write table for INF_sig_counts and re-read in - clear envrionemnt and restart as EnsDb package maskes dplyr package funcitons
# setwd("/Volumes/users/RNAseq/DADA2_PAN/DESeq2")
# write.csv(IFN_sig_counts, "IFN_sig_counts.csv", row.names = FALSE)

library(dplyr)
IFN_sig_counts <- read.csv("IFN_sig_counts.csv", header = TRUE)

IFN_sig_counts_sub <- subset(IFN_sig_counts, SYMBOL = c("IFI27", "IFI44L", "IFIT1", "ISG15", "RSAD2" | "SIGLEC1"))
IFN_sig_counts_sub <- subset(IFN_sig_counts, SYMBOL = TSPAN6)


IFI27
IFI44L
IFIT1
ISG15
RSAD2
SIGLEC1

# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- IFN_sig_counts$gene
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))




# write table of hits padj < 0.1 for SB
# setwd("/Volumes/kgibson/RNAseq/DADA2_PAN/DESeq2")
# write.csv(as.data.frame(res_inactive_sub), "deg_inactive.csv", row.names = FALSE )



# Covariate analysis with FactoMineR --------------------------------------

### PAIRED COVARIATE ANALYSIS
# generate covariate tables for each analysis - paired, active, inactive
paired_covar_meta <-  meta[c("EE1001V1", "EE1001V2", "SK1003V1","SK1003V2", "TO5009V4", "TO5009V5", "VA5004V4", "VA5004V5",
                             "GA1002V1", "GA1002V3", "VA1006V1", "VA1006V5", "CA1010V1", "CA1010V2"), ]

paired_covar_meta <- paired_covar_meta[,c("gender", "ANCA_status", "PVAS_score", "mother_ethnicity", "father_ethnicity", "mut_ADA2")]

paired_covar_meta$mut_ADA2 <- if_else(paired_covar_meta$mut_ADA2 == "positive", 1, 2) # change grouping to a qualitative variable
  # 1 = DADA2 ; 2 = PAN

paired_pca <- PCA(paired_covar_meta, quali.sup = 1:5) # perform PCA

paired_scree <- fviz_screeplot(paired_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
paired_eig <- get_eig(paired_pca)


paired_var <- paired_pca$var$contrib


paired_dim <- dimdesc(paired_pca, axes = 1, proba = 0.05)
# paired_dim$Dim.1


### ACTIVE COVARIATE ANALYSIS
active_covar_meta <- meta[c("TO5009V4", "VA5004V4", "SK1003V1", "EE1001V1", 
                      "GA1002V1", "VA1006V1", "CA5002V4", "CA1010V1"),]   

active_covar_meta <- active_covar_meta[,c("gender", "ANCA_status", "PVAS_score", "mother_ethnicity", "father_ethnicity", "mut_ADA2")]

active_covar_meta$mut_ADA2 <- if_else(active_covar_meta$mut_ADA2 == "positive", 1, 2) # change grouping to a qualitative variable
# 1 = DADA2 ; 2 = PAN

active_pca <- PCA(active_covar_meta, quali.sup = 1:5) # perform PCA

active_scree <- fviz_screeplot(active_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
active_eig <- get_eig(active_pca)

active_var <- active_pca$var$contrib

active_dim <- dimdesc(active_pca, axes = 1, proba = 0.05)
# active_dim$Dim.1

### INACTIVE COVARIATE ANALYSIS

inactive_covar_meta <- meta[c("TO5009V5", "VA5004V5", "EE1001V2","SK1003V2",
                        "BR1004V5", "GA1002V3", "VA1006V5"),]        

inactive_covar_meta <- inactive_covar_meta[,c("gender", "ANCA_status", "PVAS_score", "mother_ethnicity", "father_ethnicity", "mut_ADA2")]

inactive_covar_meta$mut_ADA2 <- if_else(inactive_covar_meta$mut_ADA2 == "positive", 1, 2) # change grouping to a qualitative variable
# 1 = DADA2 ; 2 = PAN


inactive_pca <- PCA(inactive_covar_meta, quali.sup = 1:5) # perform PCA

inactive_scree <- fviz_screeplot(inactive_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
inactive_eig <- get_eig(inactive_pca)

inactive_dim <- dimdesc(inactive_pca, axes = 1, proba = 0.05)
# inactive_dim$Dim.1

