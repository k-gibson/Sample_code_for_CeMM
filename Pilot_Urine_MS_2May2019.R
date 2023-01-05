### 02-May-2019

### PedVas Renal Outcomes Pilot - Diesease Activity - Urine MS

# load packages
library(dplyr)
library(plyr)
library(FactoMineR)
library(data.table)
library(DEP)

#set wd and read in master data
setwd("/Volumes/kgibson/Urine_MassSpec/raw_data")
master_samples <- read.csv("UrineProject_Protein_Identifications_02Apr2019.csv")
colnames(master_samples) <- c("Gene", "Protein_Accession", "BE1001V2", "BE1001V3", "HA1002V1", "HA1002V3", "TO5002V4", "TO5002V5") #edit column headers

master_controls <- read.csv("UrineProject_controls_Dec2018.csv")
colnames(master_controls) <- c("Gene", "Protein_Accession", "Ctrl_1", "Ctrl_2", "Ctrl_3", "Ctrl_4", "Ctrl_5", "Ctrl_6", "Ctrl_7", "Ctrl_8", "Ctrl_9", "Ctrl_10")

master <- left_join(master_controls, master_samples, by = "Gene")
master <- subset(master, select = -c(Protein_Accession.y))

master[master == "#N/A"] <- NA # change all #N/A to NA in data.frame (control data)

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan)) # write function to change NaN to NA in data.frame (sample data)
master[is.nan(master)] <- NA # change all NaN to NA in data.frame

# Make unique names using the annotation in the "Gene" column as a primary name and the annotation in "Protein_Accession" as name for those that do not have a gene name
data_unique <- make_unique(master, "Gene", "Protein_Accession.x", delim = ";")

# make all protein quantitaion columns numeric
data_unique$Ctrl_1 <- as.numeric(as.character(data_unique$Ctrl_1))
data_unique$Ctrl_2 <- as.numeric(as.character(data_unique$Ctrl_2))
data_unique$Ctrl_3 <- as.numeric(as.character(data_unique$Ctrl_3))
data_unique$Ctrl_4 <- as.numeric(as.character(data_unique$Ctrl_4))
data_unique$Ctrl_5 <- as.numeric(as.character(data_unique$Ctrl_5))
data_unique$Ctrl_6 <- as.numeric(as.character(data_unique$Ctrl_6))
data_unique$Ctrl_7 <- as.numeric(as.character(data_unique$Ctrl_7))
data_unique$Ctrl_8 <- as.numeric(as.character(data_unique$Ctrl_8))
data_unique$Ctrl_9 <- as.numeric(as.character(data_unique$Ctrl_9))
data_unique$Ctrl_10 <- as.numeric(as.character(data_unique$Ctrl_10))
data_unique$BE1001V2 <- as.numeric(as.character(data_unique$BE1001V2))
data_unique$BE1001V3 <- as.numeric(as.character(data_unique$BE1001V3))
data_unique$HA1002V1 <- as.numeric(as.character(data_unique$HA1002V1))
data_unique$HA1002V3 <- as.numeric(as.character(data_unique$HA1002V3))
data_unique$TO5002V4 <- as.numeric(as.character(data_unique$TO5002V4))
data_unique$TO5002V5 <- as.numeric(as.character(data_unique$TO5002V5))


# Generate a SummarizedExperiment Object
LFQ_columns <- 3:18 # numeric values corresponding to the column in data_unique 

setwd("/Volumes/kgibson/Urine_MassSpec")
exp_design <- read.csv("Urine_MS_experimental_design.csv") # read in experimental design
exp_design$label <- as.character(exp_design$label) # change varables to character vectors - required for make_se function
exp_design$condition <- as.character(exp_design$condition)
exp_design$replicate <- as.character(exp_design$replicate)

data_se <- make_se(data_unique, LFQ_columns, exp_design) # make summarized experiment

# Filter on missing values
protein_frequency_plot <- plot_frequency(data_se) # plot a barplot of the protein identification overlap between samples

data_filt <- filter_missval(data_se, thr = 1) # filter for proteins that are identified in 2 our of 3 replicates of at least one condition

protein_per_sample_plot <- plot_numbers(data_filt) # plot a barplot of the number of identified proteins per sample

plot_coverage(data_filt) # plot a barplot of the protein identification overlp between samples

# Normalization
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm) # visualize normalization by boxplots for all samples before and after normalization

# Impute data for missing values
plot_missval(data_filt) # plot a heatmap of proteins with missing values (only proteins with at least one missing value are visualized)

plot_detect(data_filt) #plot intensity distributions and cumulative fraction of prtoeins with and without missing values (to check whether missing values are biased to lower intence proteins)
# as the data shows proteins with missing values have on average low intensities, the data is "missing not at random" (MNAR) and should be imputed by a left-censored imputaiont method or random draws from a lef-shifted distribution ("MinProb and "man)

data_imp <- impute(data_norm, fun = "MinProb", q = 0.05) # impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) # impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)

plot_imputation(data_norm, data_imp) # plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp_man)

# Differential enrichment analysis
# Differential enrichment analysis based on linear models and empherical Bayes statistics
data_diff <- test_diff(data_imp, type = "control", control = "control") # test every sample vs control (tested contrasts: inactivs_vs_control, active_vs_control)
data_diff_all_contrasts <- test_diff(data_imp, type = "all") # test all possible comparisons of samples (tested contrasts: inactive_vs_active, inactive_vs_control, active_vs_control)

dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5)) # denote significant proteins based on user defined cutoffs
dep_all_contrasts <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5))

# Visualization of results

# PCA PLOT
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4) # plot the first (x) and second (y) principle components

# Correlation matrix

plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds") # plot the Pearson correlation matrix
plot_cor(dep_all_contrasts, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Heatmap of significant proteins

# plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate") )
plot_heatmap(dep_all_contrasts, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate") )

# plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, k = 6, col_limit = 10, show_row_names = FALSE)
plot_heatmap(dep_all_contrasts, type = "contrast", kmeans = TRUE, k = 6, col_limit = 10, show_row_names = FALSE)

# Volcano plots of specific contrasts
plot_volcano(dep_all_contrasts, contrast = "inactive_vs_active", label_size = 4, add_names = TRUE)
plot_volcano(dep_all_contrasts, contrast = "inactive_vs_control", label_size = 4, add_names = TRUE)
plot_volcano(dep_all_contrasts, contrast = "active_vs_control", label_size = 4, add_names = TRUE)

# Barplots of protein of inerest
plot_single(dep_all_contrasts, proteins = c("SIAE", "LDHA", "LRRC15")) # plot barplots of SIAE, LDHA, and LRRC15

# plot barplots for proteins with the data centred
plot_single(dep_all_contrasts, proteins = "SIAE", type = "centered")
plot_single(dep_all_contrasts, proteins = "LDHA", type = "centered")
plot_single(dep_all_contrasts, proteins = "LRRC15", type = "centered")
plot_single(dep_all_contrasts, proteins = "GAPDH", type = "centered")

# Frequency plot of significant proteins and overlap of conditions
plot_cond(dep_all_contrasts)

# Results table
data_results <- get_results(dep) # generaate a results table
data_results_all_contrasts <- get_results(dep_all_contrasts)

# Generate data.frame from the resulting summarized experiment object
df_wide <- get_df_wide(dep_all_contrasts) # genearte a wide data.frame
df_long <- get_df_long(dep_all_contrasts) # generate a long data.frame

# All samples - no imputation -----------------------------------------------------

# generate new se, filtering for proteins that have no missing values
complete_cases <- filter_proteins(data_se, "complete")

plot_frequency(complete_cases)
plot_numbers(complete_cases)
plot_coverage(complete_cases) # plot a barplot of the protein identification overlp between samples

# Normalization
complete_cases_norm <- normalize_vsn(complete_cases)
plot_normalization(complete_cases, complete_cases_norm) # visualize normalization by boxplots for all samples before and after normalization

## skip imputation - we have no missing values

# Differential enrichment analysis
# Differential enrichment analysis based on linear models and empherical Bayes statistics
complete_cases_diff <- test_diff(complete_cases_norm, type = "control", control = "control") # test every sample vs control (tested contrasts: inactivs_vs_control, active_vs_control)
complete_cases_diff_all_contrasts <- test_diff(complete_cases_norm, type = "all") # test all possible comparisons of samples (tested contrasts: inactive_vs_active, inactive_vs_control, active_vs_control)

complete_cases_dep <- add_rejections(complete_cases_diff, alpha = 1, lfc = log2(0)) # denote significant proteins based on user defined cutoffs
complete_cases_dep_all_contrasts <- add_rejections(complete_cases_diff_all_contrasts, alpha = 1, lfc = log2(0))

# PCA PLOT
plot_pca(complete_cases_dep, x = 1, y = 2, n = 119, point_size = 4) # plot the first (x) and second (y) principle components

# Correlation matrix
plot_cor(complete_cases_dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds") # plot the Pearson correlation matrix
plot_cor(complete_cases_dep_all_contrasts, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Heatmap of significant proteins

# plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(complete_cases_dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate") )
plot_heatmap(complete_cases_dep_all_contrasts, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate") )

# plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, k = 6, col_limit = 10, show_row_names = FALSE)
plot_heatmap(dep_all_contrasts, type = "contrast", kmeans = TRUE, k = 6, col_limit = 10, show_row_names = FALSE)

# Volcano plots of specific contrasts
plot_volcano(complete_cases_dep_all_contrasts, contrast = "inactive_vs_active", label_size = 3, add_names = TRUE)
plot_volcano(complete_cases_dep_all_contrasts, contrast = "inactive_vs_control", label_size = 3, add_names = TRUE)
plot_volcano(complete_cases_dep_all_contrasts, contrast = "active_vs_control", label_size = 3, add_names = TRUE)

# Barplots of protein of inerest
plot_single(complete_cases_dep_all_contrasts, proteins = c("L1TD1", "IGKV2-24", "ALB")) 

# plot barplots for proteins with the data centred
plot_single(complete_cases_dep_all_contrasts, proteins = "L1TD1", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "IGKV2-24", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "ALB", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "CD59", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "KNG1", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "APOH", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "COL12A1", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "ORM2", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "IGKV1D-33", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "GC", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "LAMP2", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "HP", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "APOD", type = "centered") # INTERESTING
plot_single(complete_cases_dep_all_contrasts, proteins = "AMY1A", type = "centered") # INTERESITING
plot_single(complete_cases_dep_all_contrasts, proteins = "LEAP2", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "UBB", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "SLURP1", type = "centered")
plot_single(complete_cases_dep_all_contrasts, proteins = "LRP2", type = "centered")

# Frequency plot of significant proteins and overlap of conditions
plot_cond(complete_cases_dep_all_contrasts)

# Results table
complete_cases_results <- get_results(complete_cases_dep) # generaate a results table
complete_cases_results_all_contrasts <- get_results(complete_cases_dep_all_contrasts)

# Generate data.frame from the resulting summarized experiment object
df_wide_complete_cases <- get_df_wide(complete_cases_dep_all_contrasts) # genearte a wide data.frame
df_long_complete_cases <- get_df_long(complete_cases_dep_all_contrasts) # generate a long data.frame



# write file for complete_cases_results_all_contrasts

# write.csv(complete_cases_results_all_contrasts, "complete_cases_results_all_contrasts.csv", row.names = FALSE)


# PedVas samples only -----------------------------------------------------

master_samples[is.nan(master_samples)] <- NA # change all NaN to NA in samples data.frame
  
# Make unique names using the annotation in the "Gene" column as a primary name and the annotation in "Protein_Accession" as name for those that do not have a gene name
samples_unique <- make_unique(master_samples, "Gene", "Protein_Accession", delim = ";")

# make all protein columns numeric
samples_unique$BE1001V2 <- as.numeric(as.character(samples_unique$BE1001V2))
samples_unique$BE1001V3 <- as.numeric(as.character(samples_unique$BE1001V3))
samples_unique$HA1002V1 <- as.numeric(as.character(samples_unique$HA1002V1))
samples_unique$HA1002V3 <- as.numeric(as.character(samples_unique$HA1002V3))
samples_unique$TO5002V4 <- as.numeric(as.character(samples_unique$TO5002V4))
samples_unique$TO5002V5 <- as.numeric(as.character(samples_unique$TO5002V5))

samples_unique <- na.omit(samples_unique) # remove all query's with NAs (removed all due to low sample numbers... dont want to impute on 3 samples)

# Generate a SummarizedExperiment Object
LFQ_columns_samples <- 3:8 # numeric values corresponding to the column in data_unique 

exp_design_samples <- exp_design
exp_design_samples <- subset(exp_design_samples, condition != "control") # remove controls from experimental design

data_se_samples <- make_se(samples_unique, LFQ_columns_samples, exp_design_samples) # make summarized experiment


# Filter on missing values - dont need to filter as removed all protein querys with NAs already
plot_frequency(data_se_samples) # plot a barplot of the protein identification overlap between samples
plot_numbers(data_se_samples) # plot a barplot of the number of identified proteins per sample
plot_coverage(data_se_samples) # plot a barplot of the protein identification overlp between samples

# Normalization
data_norm_samples <- normalize_vsn(data_se_samples)
plot_normalization(data_se_samples, data_norm_samples) # visualize normalization by boxplots for all samples before and after normalization

# Differential enrichment analysis

# Differential enrichment analysis based on linear models and empherical Bayes statistics
data_diff_samples <- test_diff(data_norm_samples, type = "control", control = "inactive") # test every sample vs control (tested contrasts: active_vs_inactive)

dep_samples <- add_rejections(data_diff_samples, alpha = 1, lfc = log2(0)) # denote significant proteins based on user defined cutoffs


# Visualization of results

# PCA PLOT
plot_pca(dep_samples, x = 1, y = 2, n = 129, point_size = 4) # plot the first (x) and second (y) principle components

# Correlation matrix
plot_cor(dep_samples, significant = TRUE, lower = -1, upper = 1, pal = "Reds")

# Heatmap of significant proteins

# plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep_samples, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate") )

# plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep_samples, type = "contrast", kmeans = TRUE, k = 6, col_limit = 10, show_row_names = FALSE)

# Volcano plots of specific contrasts
plot_volcano(dep_samples, contrast = "active_vs_inactive", label_size = 3, add_names = TRUE)

# Barplots of protein of inerest
plot_single(dep_samples, proteins = c("L1TD1", "ALB", "HP", "GC")) # plot barplots of SIAE, LDHA, and LRRC15

# plot barplots for proteins with the data centred
plot_single(dep_samples, proteins = "L1TD1", type = "centered")
plot_single(dep_samples, proteins = "ALB", type = "centered")
plot_single(dep_samples, proteins = "HP", type = "centered")
plot_single(dep_samples, proteins = "LAMP2", type = "centered") # LOOKS INTERESTING - could be significant if done in a paired analysis
plot_single(dep_samples, proteins = "LEAP2", type = "centered")
plot_single(dep_samples, proteins = "IGKV1D-33", type = "centered")
plot_single(dep_samples, proteins = "AMY1A", type = "centered")
plot_single(dep_samples, proteins = "ATP6V0C", type = "centered")
plot_single(dep_samples, proteins = "APOD", type = "centered") # LOOKS INTERESTING - checked paired anlaysis
plot_single(dep_samples, proteins = "KNG1", type = "centered") # LOOKS INTERESTING - checked paired anlaysis
plot_single(dep_samples, proteins = "S100A7", type = "centered")
plot_single(dep_samples, proteins = "TFF2", type = "centered")
plot_single(dep_samples, proteins = "GC", type = "centered")

# Results table
data_results_samples <- get_results(dep_samples) # generaate a results table

# Generate data.frame from the resulting summarized experiment object
df_wide <- get_df_wide(dep_all_contrasts) # genearte a wide data.frame
df_long <- get_df_long(dep_all_contrasts) # generate a long data.frame



