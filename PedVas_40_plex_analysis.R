### 2018/01/04

# 40plex analysis of PedVas samples

library(dplyr)
library(FactoMineR)
library(missMDA)
library(ggplot2)
library(factoextra)

# set wd and read in file for 40plex data
setwd("/Users/kristen/.Volumes/kgibson/40_plex")
master <- read.csv("40plex_data.csv", header = TRUE, row.names = 1) # use pedvas ID as row names

# select for supplementary variables and 40plex variables
meta <- select(master, c(Visit_No, MD_diagnosis, EMA_diagnosis, ANCA_status, Hemoglobin, PVAS_score))

data <- select(master, c(IFNg_conc : VCAM.1_conc))

# remove IL.8 from chemokine panel and VEGF from angiogenesis panel (they are repeated)
data <- select(data, -c(IL.8_HA_conc, VEGF_conc_2)) 

# convert NAs to 0 (as they represent concentrations that were too low to be detected)
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)} # formulate the na.zero function (whcih will change any NAs in the data to 0)

data <- na.zero(data) # apply the na.zero function to the 40plex data table

# I changed the NAs to zero instead of imputing them, as they represent low concentrations
#   undetectable by the analysis. I kept the following three lines of script in case
#   I wanted to change it back the imputation method to deal with the NAs
# impute for incomplete variables using the missMDA package
# nb <- estim_ncpPCA(total_serum_data, ncp.min = 0, ncp.max = 5)
# res.imputed <- imputePCA(total_serum_data, nb$ncp) # ncp of nb is 2

# munge PVAS score and data together in 1 dataframe
data <- merge(select(meta, PVAS_score), data, by = "row.names")
rownames(data) <- data[,1] #add rownames
data <- data[,-1] # remove rownames column

# impute for PVAS score
ncp <- estim_ncpPCA(data, ncp.min = 0, ncp.max = 5) # calculate number of components for imputation
res.impute <- imputePCA(data, ncp=2) # impute for numerical variables

# munge imputed data with metadata
meta <- meta[,-6]
sub_master <- cbind(meta, res.impute$completeObs)

# preform PCA 
total_pca <- PCA(sub_master, quali.sup = 1:5, quanti.sup = 6)
plot(total_pca, shadowtext = TRUE, cex = 0.8, habillage = 3)

total_scree <- fviz_screeplot(total_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
total_eig <- get_eig(total_pca) # amount of variance attributed to each dimension of the PCA
# Plot the correlations of the variables with the components 
# Visualize variables that group together
fviz_pca_var(total_pca, col.var = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 0.5) + theme_minimal()

total_var <- total_pca$var$contrib # describe percentatge variables contribuute to the components
fviz_contrib(total_pca, choice = "var", axes = 1) # visualize contribution of variables to PC1
fviz_contrib(total_pca, choice = "var", axes = 2) # visualize contribution of variables to PC2
fviz_contrib(total_pca, choice = "var", axes = 1:2) # visualize contribution of variables to PC1 and 2


fviz_pca_var(total_pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 50) + theme_minimal()


total_dim <- dimdesc(total_pca, axes = 1:3, proba = 0.05)
# total_dim$Dim.1
# total_dim$Dim.2
# total_dim$Dim.3


# PCA of each panel -------------------------------------------------------

# PRO-INFLAMATORY PANEL
# select for pro-inflammatory panel
pro_inflam <- select(sub_master, c(Visit_No : PVAS_score, IFNg_conc : TNFa_conc))
pro_inflam_pca <- PCA(pro_inflam, quali.sup = 1:5, quanti.sup = 6)
plot(pro_inflam_pca, shadowtext = TRUE, cex = 0.7)

pro_inflam_scree <- fviz_screeplot(pro_inflam_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
pro_inflam_eig <- get_eig(pro_inflam_pca) # amount of variance attributed to each dimension of the PCA
# Plot the correlations of the variables with the components 
# Visualize variables that group together
fviz_pca_var(pro_inflam_pca, col.var = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 0.5) + theme_minimal()

pro_inflam_var <- pro_inflam_pca$var$contrib # describe percentatge variables contribuute to the components
fviz_contrib(pro_inflam_pca, choice = "var", axes = 1) # visualize contribution of variables to PC1
fviz_contrib(pro_inflam_pca, choice = "var", axes = 2) # visualize contribution of variables to PC2
fviz_contrib(pro_inflam_pca, choice = "var", axes = 3) # visualize contribution of variables to PC3
fviz_contrib(pro_inflam_pca, choice = "var", axes = 1:3) # visualize contribution of variables to PC1, 2, and 3
fviz_contrib(pro_inflam_pca, choice = "var", axes = 1:2)

fviz_pca_var(pro_inflam_pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 50) + theme_minimal()

pro_inflam_dim <- dimdesc(pro_inflam_pca)
# pro_inflam_dim$Dim.1
# pro_inflam_dim$Dim.2
# pro_inflam_dim$Dim.3

# CYTOKINE PANEL
#select for cytokine panel
cyto <- select(sub_master, c(Visit_No : PVAS_score, GM.CSF_conc : VEGF_conc_1))
cyto_pca <- PCA(cyto, quali.sup = 1:5, quanti.sup = 6)
plot(cyto_pca, shadowtext = TRUE, cex = 0.7)

# remove outliers and redo PCA
cyto <- cyto[-c(24),]
cyto_pca <- PCA(cyto, quali.sup = 1:5, quanti.sup = 6)
plot(cyto_pca, shadowtext = TRUE, cex = 0.7)

cyto_scree <- fviz_screeplot(cyto_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
cyto_eig <- get_eig(cyto_pca) # amount of variance attributed to each dimension of the PCA
# Plot the correlations of the variables with the components 
# Visualize variables that group together
fviz_pca_var(cyto_pca, col.var = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 0.5) + theme_minimal()

cyto_var <- cyto_pca$var$contrib # describe percentatge variables contribuute to the components
fviz_contrib(cyto_pca, choice = "var", axes = 1) # visualize contribution of variables to PC1
fviz_contrib(cyto_pca, choice = "var", axes = 2) # visualize contribution of variables to PC2
fviz_contrib(cyto_pca, choice = "var", axes = 3) # visualize contribution of variables to PC3
fviz_contrib(cyto_pca, choice = "var", axes = 1:3) # visualize contribution of variables to PC1, 2, and 3
fviz_contrib(cyto_pca, choice = "var", axes = 1:2)

fviz_pca_var(cyto_pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 50) + theme_minimal()

cyto_dim <- dimdesc(cyto_pca)
# cyto_dim$Dim.1
# cyto_dim$Dim.2
# cyto_dim$Dim.3

# CHEMOKINE PANEL
# select for chemokine panel
chemo <- select(sub_master, c(Visit_No : PVAS_score, Eotaxin_conc : TARC_conc))
chemo_pca <- PCA(chemo, quali.sup = 1:5, quanti.sup = 6)
plot(chemo_pca, shadowtext = TRUE, cex = 0.6)

# remove outliers and redo PCA
chemo <- chemo[-c(19),]
chemo_pca <- PCA(chemo, quali.sup = 1:5, quanti.sup = 6)
plot(chemo_pca, shadowtext = TRUE, cex = 0.6)

chemo_scree <- fviz_screeplot(chemo_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
chemo_eig <- get_eig(chemo_pca) # amount of variance attributed to each dimension of the PCA
# Plot the correlations of the variables with the components 
# Visualize variables that group together
fviz_pca_var(chemo_pca, col.var = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 0.5) + theme_minimal()

chemo_var <- chemo_pca$var$contrib # describe percentatge variables contribuute to the components
fviz_contrib(chemo_pca, choice = "var", axes = 1) # visualize contribution of variables to PC1
fviz_contrib(chemo_pca, choice = "var", axes = 2) # visualize contribution of variables to PC2
fviz_contrib(chemo_pca, choice = "var", axes = 3) # visualize contribution of variables to PC3
fviz_contrib(chemo_pca, choice = "var", axes = 1:3) # visualize contribution of variables to PC1, 2, and 3
fviz_contrib(chemo_pca, choice = "var", axes = 1:2)

fviz_pca_var(chemo_pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 50) + theme_minimal()

chemo_dim <- dimdesc(chemo_pca)
# chemo_dim$Dim.1
# chemo_dim$Dim.2
# chemo_dim$Dim.3

# ANGIOGENESIS PANEL
# select for angiogenesis panel
angio <- select(sub_master, c(Visit_No : PVAS_score, bFGF_conc : VEGF.D_conc))
angio_pca <- PCA(angio, quali.sup = 1:5, quanti.sup = 6)
plot(angio_pca, shadowtext = TRUE, cex = 0.6)

# remove outlier and redo PCA
angio <- angio[-c(11),]
angio_pca <- PCA(angio, quali.sup = 1:5, quanti.sup = 6)
plot(angio_pca, shadowtext = TRUE, cex = 0.6)

angio_scree <- fviz_screeplot(angio_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
angio_eig <- get_eig(angio_pca) # amount of variance attributed to each dimension of the PCA
# Plot the correlations of the variables with the components 
# Visualize variables that group together
fviz_pca_var(angio_pca, col.var = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 0.5) + theme_minimal()

angio_var <- angio_pca$var$contrib # describe percentatge variables contribuute to the components
fviz_contrib(angio_pca, choice = "var", axes = 1) # visualize contribution of variables to PC1
fviz_contrib(angio_pca, choice = "var", axes = 2) # visualize contribution of variables to PC2
fviz_contrib(angio_pca, choice = "var", axes = 3) # visualize contribution of variables to PC3
fviz_contrib(angio_pca, choice = "var", axes = 1:3) # visualize contribution of variables to PC1, 2, and 3
fviz_contrib(angio_pca, choice = "var", axes = 1:2)

fviz_pca_var(angio_pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 50) + theme_minimal()

angio_dim <- dimdesc(angio_pca)
# angio_dim$Dim.1
# angio_dim$Dim.2
# angio_dim$Dim.3

# VASCULAR INJURY PANEL
# select for vascular injury panel
vasc_inj <- select(sub_master, c(Visit_No : PVAS_score, CRP_conc : VCAM.1_conc))
vasc_inj <- select(vasc_inj, - SAA_conc) # take out SAA measurement
vasc_inj_pca <- PCA(vasc_inj, quali.sup = 1:5, quanti.sup = 6)
plot(vasc_inj_pca, shadowtext = TRUE, cex = 0.6, habillage = 3)

vasc_inj_scree <- fviz_screeplot(vasc_inj_pca, ncp = 10) # visualize the eigenvalues/variance of the dimensions from the results of the PCA
vasc_inj_eig <- get_eig(vasc_inj_pca) # amount of variance attributed to each dimension of the PCA
# Plot the correlations of the variables with the components 
# Visualize variables that group together
fviz_pca_var(vasc_inj_pca, col.var = "cos2") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 0.5) + theme_minimal()

vasc_inj_var <- vasc_inj_pca$var$contrib # describe percentatge variables contribuute to the components
fviz_contrib(vasc_inj_pca, choice = "var", axes = 1) # visualize contribution of variables to PC1
fviz_contrib(vasc_inj_pca, choice = "var", axes = 2) # visualize contribution of variables to PC2
fviz_contrib(vasc_inj_pca, choice = "var", axes = 3) # visualize contribution of variables to PC3
fviz_contrib(vasc_inj_pca, choice = "var", axes = 1:3) # visualize contribution of variables to PC1, 2, and 3
fviz_contrib(vasc_inj_pca, choice = "var", axes = 1:2)

fviz_pca_var(vasc_inj_pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", mid = "blue", high = "red",
                        midpoint = 50) + theme_minimal()

vasc_inj_dim <- dimdesc(vasc_inj_pca)
# vasc_inj_dim$Dim.1
# vasc_inj_dim$Dim.2
# vasc_inj_dim$Dim.3


# Hierarchical clustering on principle components -------------------------

# total serum data
total_hcpc <- HCPC(total_pca)

# pro-inflammatory panel clustering
pro_inflam_hcpc <- HCPC(pro_inflam_pca)

# cytokine panel clustering
cyto_hcpc <- HCPC(cyto_pca)

# chemokine panel clustering
chemo_hcpc <- HCPC(chemo_pca)

# angiogenesis panel clustering
angio_hcpc <- HCPC(angio_pca)

# vascular injury panel clustering
vasc_inj_hcpc <- HCPC(vasc_inj_pca)


# Select for variables of interest

IL2 <- subset(pro_inflam, select = c(2:6, IL.2_conc))
IL13 <- subset(pro_inflam, select = c(2:4, IL.13_conc))
IL4 <- subset(pro_inflam, select = c(2:4, IL.4_conc)) 
IL12_IL23p40 <- subset(cyto, select = c(2:4, IL.12_IL.23p40_conc))
TARC <- subset(chemo, select = c(2:4, TARC_conc))  
MCP1 <- subset(chemo, select = c(2:4, MCP.1_conc))  
IP10 <- subset(chemo, select = c(2:4, IP.10_conc))
ICAM <- subset(vasc_inj, select = c(1:5, ICAM.1_conc))
MCP4 <- subset(chemo, select = c(1:6, MCP.4_conc))  
Eotaxin3 <- subset(chemo, select = c(1:6, Eotaxin.3_conc))
CRP <- subset(sub_master, select = c(2:6, CRP_conc))  
IL7 <- subset(cyto, select = c(1:6, IL.7_conc))  
VEGFD <- subset(angio, select = c(1:6, VEGF.D_conc))
Heme <- subset(sub_master, select = 1:6)
  