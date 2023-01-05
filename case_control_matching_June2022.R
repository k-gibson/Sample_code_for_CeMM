## PCA matching of cases and crontrols for PedVas GWAS
# 22-Jun-2022

library(PCAmatchR)
library(optmatch)



# Case-Control Matching for PedVas/WTCCC Data -----------------------------

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")
library(dplyr)

# Combined pediatric and adult cases --------------------------------------

#Create PC data frame
PCs_combined <- read.table("WTCCC_Vasculitis_Eu_qc_no_artefacts_PCA.eigenvec")
PCs_combined <- select(PCs_combined, -c(V2))
colnames(PCs_combined) <- c("sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PV10")

#Create eigen value vector
eigenvalues_combined <- read.table("WTCCC_Vasculitis_Eu_qc_no_artefacts_PCA.eigenval")
colnames(eigenvalues_combined) <- "eigen_values"

eigen_vals_combined <- eigenvalues_combined$eigen_values

all_eigenvalues_combined <- read.table("WTCCC_Vasculitis_Eu_qc_no_artefacts_PCA_all.eigenval")
colnames(all_eigenvalues_combined) <- "eigen_values"

all_eigen_vals_combined <- all_eigenvalues_combined$eigen_values


#Asign case/control populations

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019")

meta_combined <- read.table("WTCCC_Vasculitis_Eu_qc_no_artefacts.fam")
colnames(meta_combined) <- c("sample", "fammily", "F", "M", "sex", "case")

meta_combined$case <- ifelse(meta_combined$case == 2, 1, 0)

#Case-Control Matching

match_maker_output_combined<- match_maker(PC = PCs_combined,
                                          eigen_value = eigen_vals_combined,
                                          data = meta_combined,
                                          ids = c("sample"),
                                          case_control="case",
                                          num_controls = 5,
                                          eigen_sum = sum(all_eigen_vals_combined))

PCA_matches_combined <- match_maker_output_combined$matches

### Case-Control Matching Visualization

plot_maker(data=match_maker_output_combined,
           x_var="PC1",
           y_var="PC2",
           case_control="case",
           line=F)


plot_maker(data=match_maker_output_combined,
           x_var="PC1",
           y_var="PC2",
           case_control="case",
           line=T)

# make a plot to show selection of matched cases and control 

library(ggplot2)
library(tidyr)

total_plot <- select(meta_combined, c("sample", "case"))
total_plot <- left_join(total_plot, PCs_combined, by = "sample")

match_ids <- select(PCA_matches_combined, c("sample", "case", "match_final"))
colnames(match_ids) <- c("sample", "matched_controls", "match_key")

total_plot <- left_join(total_plot, match_ids, by = "sample")

total_plot$matched_controls <- replace_na(total_plot$matched_controls, 2)

total_plot$matched_controls <- ifelse(total_plot$matched_controls == 0, "matched_control",
                                      ifelse(total_plot$matched_controls == 1, "case", 
                                             ifelse(total_plot$matched_controls == 2, "control", NA)))

total_plot$matched_controls <- as.factor(total_plot$matched_controls) # change variable to factor and relevel
levels(total_plot$matched_controls) <- c("case", "matched_control", "control")

total_plot$matched_controls <- factor(total_plot$matched_controls, levels = c("case", "matched_control", "control"))


total_plot_sp <- ggplot(total_plot, aes(x=PC1, y=PC2, color=matched_controls)) + geom_point(alpha = 1) +
  scale_color_manual(values=c('red', 'blue', 'gray'))


# save match output file

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

# write.table(PCA_matches_combined, "PCA_matches_combined_1Jul2022.txt", row.names = FALSE, quote = FALSE)

# Pediatric Cases ---------------------------------------------------------

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

PCs_ped <- read.table("PedVas_WTCCC_Eu_qc_no_adult_PCA.eigenvec")
PCs_ped <- select(PCs_ped, -c(V2))
colnames(PCs_ped) <- c("sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PV10")

#Create eigen value vector
eigenvalues_ped <- read.table("PedVas_WTCCC_Eu_qc_no_adult_PCA.eigenval")
colnames(eigenvalues_ped) <- "eigen_values"

eigen_vals_ped <- eigenvalues_ped$eigen_values

all_eigenvalues_ped <- read.table("PedVas_WTCCC_Eu_qc_no_adult_PCA_all.eigenval")
colnames(all_eigenvalues_ped) <- "eigen_values"

all_eigen_vals_ped <- all_eigenvalues_ped$eigen_values


#Asign case/control populations

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019")

meta_ped <- read.table("PedVas_WTCCC_Eu_qc_no_adult.fam")
colnames(meta_ped) <- c("sample", "fammily", "F", "M", "sex", "case")

meta_ped$case <- ifelse(meta_ped$case == 2, 1, 0)

#Case-Control Matching

match_maker_output_ped<- match_maker(PC = PCs_ped,
                                          eigen_value = eigen_vals_ped,
                                          data = meta_ped,
                                          ids = c("sample"),
                                          case_control="case",
                                          num_controls = 5,
                                          eigen_sum = sum(all_eigen_vals_ped))

PCA_matches_ped <- match_maker_output_ped$matches

### Case-Control Matching Visualization

plot_maker(data=match_maker_output_ped,
           x_var="PC1",
           y_var="PC2",
           case_control="case",
           line=F)


plot_maker(data=match_maker_output_ped,
           x_var="PC1",
           y_var="PC2",
           case_control="case",
           line=T)

# save match output file

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

# write.table(PCA_matches_ped, "PCA_matches_pediatric_1Jul2022.txt", row.names = FALSE, quote = FALSE)

# Adult cases -------------------------------------------------------------

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

PCs_adult <- read.table("PedVas_WTCCC_Eu_qc_adult_PCA.eigenvec")
PCs_adult <- select(PCs_adult, -c(V2))
colnames(PCs_adult) <- c("sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PV10")

#Create eigen value vector
eigenvalues_adult <- read.table("PedVas_WTCCC_Eu_qc_adult_PCA.eigenval")
colnames(eigenvalues_adult) <- "eigen_values"

eigen_vals_adult <- eigenvalues_adult$eigen_values

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019")
all_eigenvalues_adult <- read.table("PedVas_WTCCC_Eu_qc_adult_PCA_all.eigenval")
colnames(all_eigenvalues_adult) <- "eigen_values"

all_eigen_vals_adult <- all_eigenvalues_adult$eigen_values


#Asign case/control populations

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019")

meta_adult <- read.table("PedVas_WTCCC_Eu_qc_adult.fam")
colnames(meta_adult) <- c("sample", "fammily", "F", "M", "sex", "case")

meta_adult$case <- ifelse(meta_adult$case == 2, 1, 0)

#Case-Control Matching

match_maker_output_adult <- match_maker(PC = PCs_adult,
                                     eigen_value = eigen_vals_adult,
                                     data = meta_adult,
                                     ids = c("sample"),
                                     case_control="case",
                                     num_controls = 5,
                                     eigen_sum = sum(all_eigen_vals_adult))

PCA_matches_adult <- match_maker_output_adult$matches

### Case-Control Matching Visualization

plot_maker(data=match_maker_output_adult,
           x_var="PC1",
           y_var="PC2",
           case_control="case",
           line=F)


plot_maker(data=match_maker_output_adult,
           x_var="PC1",
           y_var="PC2",
           case_control="case",
           line=T)


# save match output file

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

# write.table(PCA_matches_adult, "PCA_matches_adult_1Jul2022.txt", row.names = FALSE, quote = FALSE)





# Create sub groups based on ANCA  ----------------------------------------
setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

PCA_matches_ped <- read.table("PCA_matches_pediatric_1Jul2022.txt", header = TRUE) # read in PCA match file for peds
PCA_matches_adult <- read.table("PCA_matches_adult_1Jul2022.txt", header = TRUE)  # read in PCA mtch file for adults

setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019")
meta_ANCA <- read.csv("genotyped_samples_2019_master_V3.csv", header = TRUE) # read in meta file of genotyped samples (ped and adult)

ANCA_ped <- select(meta_ANCA, c("Lab_ID", "anca_elisa")) # select ANCA status from metadata

PCA_matches_ped_ANCA <- left_join(PCA_matches_ped, ANCA_ped, by = c("sample" = "Lab_ID")) # munch PCA matches and ANCA status
PCA_matches_ped_ANCA$anca_elisa[is.na(PCA_matches_ped_ANCA$anca_elisa)] = "control" # sub elisa NAs for "control" (all WTCCC controls)

ped_PR3 <- subset(PCA_matches_ped_ANCA, anca_elisa == "PR3") # subset for PR3
ped_MPO_PR3 <- subset(PCA_matches_ped_ANCA, anca_elisa == "MPO_PR3")
ped_PR3 <- rbind(ped_PR3, ped_MPO_PR3)

PCA_matches_ped_PR3 <- left_join(ped_PR3, PCA_matches_ped_ANCA, by = "match_final") 
PCA_matches_ped_PR3 <- select(PCA_matches_ped_PR3, c("sample.y", "fammily.y", "match_final", "anca_elisa.y")) # clean up dataframe
# write.table(PCA_matches_ped_PR3, "PCA_matches_pediatric_PR3_22Sep2022.txt", row.names = FALSE, quote = FALSE)


ped_MPO <- subset(PCA_matches_ped_ANCA, anca_elisa == "MPO")
ped_MPO <- rbind(ped_MPO, ped_MPO_PR3)

PCA_matches_ped_MPO <- left_join(ped_MPO, PCA_matches_ped_ANCA, by = "match_final")
PCA_matches_ped_MPO <- select(PCA_matches_ped_MPO, c("sample.y", "fammily.y", "match_final", "anca_elisa.y"))
# write.table(PCA_matches_ped_MPO, "PCA_matches_pediatric_MPO_22Sep2022.txt", row.names = FALSE, quote = FALSE)





ped_cov <- read.table("PedVas_WTCCC_Eu_no_adult_covar_trimmed.cov", header = TRUE)
ped_HLA_DPB1 <- select(ped_cov, c("FID", "HLA_DPB1_O401"))

ped_PR3 <- left_join(ped_PR3, ped_HLA_DPB1, by = c("sample" = "FID"))
ped_MPO <- left_join(ped_MPO, ped_HLA_DPB1, by = c("sample" = "FID"))

