### 11-Jun-2019

## Comparative RNA seq analysis between DADA2 and PAN patients

# load packages
library(dplyr)
library(tibble)
library(edgeR)
library(ggplot2)
library(GO.db)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set wd and load all count files 
setwd("/Users/kristen/.Volumes/kgibson/RNAseq/DADA2_PAN/PedVas/HTSeq_files")

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

setwd("/Users/kristen/.Volumes/kgibson/RNAseq/ReSeq_201903/Counts") 
CA1010V2 <- read.table("Pedvas_308.count", col.names = c("gene_ID", "CA1010V2"))

#munge count table together
counts <- left_join(CA5002V4, CA5002V5, by = "gene_ID")
counts <- left_join(counts, EE1001V1, by = "gene_ID")
counts <- left_join(counts, EE1001V2, by = "gene_ID")
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
counts <- left_join(counts, CA1001V2, by = "gene_ID")
counts <- left_join(counts, CA1010V1, by = "gene_ID")
counts <- left_join(counts, CA1010V2, by = "gene_ID")
counts <- left_join(counts, GA1002V1, by = "gene_ID")
counts <- left_join(counts, GA1002V3, by = "gene_ID")
counts <- left_join(counts, VA1006V1, by = "gene_ID")
counts <- left_join(counts, VA1006V5, by = "gene_ID")

counts <- column_to_rownames(counts, var = "gene_ID") # make gene ID column row names

# read in meta-data
setwd("/Users/kristen/.Volumes/kgibson/RNAseq/DADA2_PAN")
meta <- read.csv("PAN_DADA2_metadata.csv")
meta <- tibble::column_to_rownames(meta, var = "lab_id")

disease_status <- ifelse(meta$PVAS_score>1, "active", "not_active") # make vector to denote active or inactive disease
meta$status <- disease_status

# DADA2 vs PAN paired analysis --------------------------------------------

# keep only DADA2 and PAN samples with paired data (1 active + 1 inactive count file)

paired_counts <- counts[ , c("EE1001V1", "EE1001V2", "SK1003V1","SK1003V2", "TO5009V4", "TO5009V5", "VA5004V4", "VA5004V5",
                             "GA1002V1", "GA1002V3", "VA1006V1", "VA1006V5", "CA1010V1", "CA1010V2")]

# convert NAs to 0s
paired_counts[is.na(paired_counts)] <- 0

paired_meta <- meta[c("EE1001V1", "EE1001V2", "SK1003V1","SK1003V2", "TO5009V4", "TO5009V5", "VA5004V4", "VA5004V5",
                      "GA1002V1", "GA1002V3", "VA1006V1", "VA1006V5", "CA1010V1", "CA1010V2"), ]

# filter weakly expressed and noninformative features

#noint <- rownames(paired_counts) %in%
#  c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")

#cpms <- cpm(paired_counts)
#keep <- rowSums(cpms > 1) >= 3 & !noint # remove features without at least 1 read per million in 3 samples
#paired_counts <- paired_counts[keep,]
#paired_counts <- as.matrix(paired_counts)

# create a DGEList object

d <- DGEList(counts = paired_counts, samples = paired_meta)

d$samples$lib.size <- colSums(d$counts) # make column for library size

d$samples$group <- c(1,2,1,2,1,2,1,2,3,4,3,4,3,4) # rename groups for samples 1:DADA2 active, 2:DADA2 inactive, 3:PAN active, 4:PAN inactive
keep <- filterByExpr(d, )
d <- d[keep, ,keep.lib.sizes=FALSE]

d <- calcNormFactors(d)
plotMDS(d, labels = rownames(paired_meta), 
        col = c("red", "blue") [factor(paired_meta$mut_ADA2)])

# make design matrix

Patient <- gl(4, 2, length = 14)
DADA2 <- factor(paired_meta$mut_ADA2, levels = c("negative", "positive"))
Status <- factor(paired_meta$status, levels = c("not_active", "active"))

paired_design <- model.matrix(~DADA2+Patient+DADA2:Status)
#paired_design2 <- paired_design[,-c(7)] # remove the DADA2negative:Patient4 columns, as these contrast cannot be made due to an n=3 for DADA2 negative

# Estimte dispersions
d2 <- estimateDisp(d, paired_design)

fit <- glmQLFit(d2, paired_design)

# find genes that are differentially expresse d during activity in CECR1_var positive vs negative patients
de_paired <- glmQLFTest(fit, contrast = c(0,0,0,0,0,-1,1))
de_toptags <- topTags(de_paired, n = 70000)

summary(decideTests(de_paired))
plotMD(de_paired)

# GO enrichment anlysis
de_paired_sub <- subset(de_paired$table, PValue < 0.05)
gene.ids_paired <- mapIds(org.Hs.eg.db, keys=rownames(de_paired_sub), keytype="ENSEMBL", column = c("ENTREZID"))
de_paired_sub$genes <- gene.ids_paired
go_paired <- goana(de_paired_sub$genes)

kegg_paired <- kegga(de_paired_sub$genes)


# Active DADA2 vs Active PAN (not paired) ----------------------------------------

# select active count files
active_counts <- dplyr::select(counts, ends_with("V1"), ends_with("V4"))

# remove genes with high 0 counts
noint_active <- rownames(active_counts) %in%
  c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")

cpms_active <- cpm(active_counts)
keep_active <- rowSums(cpms_active > 1) >= 2 & !noint_active # remove features without at least 1 read per million in 3 samples
active_counts <- active_counts[keep_active,]

# select active metadata

meta_active <- meta[c("TO5009V4", "VA5004V4", "SK1003V1", "EE1001V1", 
                      "GA1002V1", "VA1006V1", "CA5002V4", "CA1010V1"),]

# make DGEList object
d_active <- DGEList(counts = active_counts, group = meta_active$mut_ADA2)

# estimate normalization factors
d_active <- calcNormFactors(d_active)

# inpect the relationships between samples using multidimensional scaling (MDS) plot
plotMDS(d_active, labels = rownames(meta_active), 
        col = c("red", "blue") [factor(meta_active$mut_ADA2)])


# estimate dispersions
d_active <- estimateCommonDisp(d_active)
d_active <- estimateTagwiseDisp(d_active)

plotBCV(d_active)

# test for differnatially expressed genes

de_active <- exactTest(d_active, pair = c("negative", "positive"))
tophits_active <- topTags(de_active, n = nrow(active_counts), p.value = 0.05)

summary(decideTests(de_active))
plotMD(de_active)

gene.ids <- mapIds(org.Hs.eg.db, keys=rownames(de_active), keytype="ENSEMBL", column = c("ENTREZID"))
de_active$genes <- gene.ids
go <- goana(de_active$genes)

tophits_gene.ds <- mapIds(org.Hs.eg.db, keys=rownames(tophits_active), keytype="ENSEMBL", column = c("ENTREZID"))
tophits_active$table$ID <- tophits_gene.ds
go_active_top <- goana(tophits_active$table$ID)
#tophits_gene.ids <- mapIds(org.Hs.eg.db, keys=rownames(tophits_active), keytype="ENSEMBL", column = c("SYMBOL"))
#tophits_active$table$ID <- tophits_gene.ids

#write.csv(tophits_active$table, "tophits_active.csv", row.names = TRUE)

keg_active <- kegga(de_active$genes)
keg_active_top <- kegga(tophits_active$table$ID)

# make heatmap of DE genes
heatmap()

# Inactive DADA2 vs Inactive PAN (not paired) -----------------------------

# select inactive count files
inactive_counts <- dplyr::select(counts, -ends_with("V1"), -ends_with("V4"))
inactive_counts <- inactive_counts[ , c("TO5009V5", "VA5004V5", "EE1001V2","SK1003V2",
                                        "BR1004V5", "GA1002V3", "VA1006V5")]

# remove genes with high 0 counts
noint_inactive <- rownames(inactive_counts) %in%
  c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")

cpms_inactive <- cpm(inactive_counts)
keep_inactive <- rowSums(cpms_inactive > 1) >= 3 & !noint_inactive # remove features without at least 1 read per million in 3 samples
inactive_counts <- inactive_counts[keep_inactive,]

# select inactive metadata
meta_inactive <- meta[c("TO5009V5", "VA5004V5", "EE1001V2", "SK1003V2",
                          "BR1004V5", "GA1002V3", "VA1006V5"),]

# make DGEList object
d_inactive <- DGEList(counts = inactive_counts, group = meta_inactive$mut_ADA2)

# estimate normalization factors
d_inactive <- calcNormFactors(d_inactive)

# inpect the relationships between samples using multidimensional scaling (MDS) plot
plotMDS(d_inactive, labels = rownames(meta_inactive), 
        col = c("red", "blue") [factor(meta_inactive$mut_ADA2)])

# estimate dispersions
d_inactive <- estimateCommonDisp(d_inactive)
d_inactive <- estimateTagwiseDisp(d_inactive)

plotBCV(d_inactive)

# test for differnatially expressed genes

de_inactive <- exactTest(d_inactive, pair = c("negative", "positive"))
tophits_inactive <- topTags(de_inactive, n = nrow(inactive_counts), p = 0.05)

summary(decideTests(de_inactive))
plotMD(de_inactive)

gene.ids_inactive <- mapIds(org.Hs.eg.db, keys=rownames(de_inactive), keytype="ENSEMBL", column = c("ENTREZID"))
de_inactive$genes <- gene.ids_inactive
go_inactive <- goana(de_inactive$genes)
 
tophits_gene.ids_inactive <- mapIds(org.Hs.eg.db, keys=rownames(tophits_inactive), keytype="ENSEMBL", column = c("ENTREZID"))
tophits_inactive$table$ID <- tophits_gene.ids_inactive
go_inactive_top <- goana(tophits_inactive$table$ID)
#tophits_gene.ids <- mapIds(org.Hs.eg.db, keys=rownames(tophits_inactive), keytype="ENSEMBL", column = c("SYMBOL"))
#tophits_inactive$table$ID <- tophits_gene.ids

# write.csv(tophits_inactive$table, "tophits_inactive.csv", row.names = TRUE)

keg_inactive <- kegga(de_inactive$genes)
kegg_inactive_top <- kegga(tophits_inactive$table$ID)
