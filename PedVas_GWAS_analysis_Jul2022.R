# GWAS analysis and plots for AnR revision
# 2-Jul-2022

library(dplyr)
library(qqman)
library(plyr)
library(stringr)
library(ggplot2)

# read in all assoc files
setwd("W:/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")
setwd("/Volumes/31 Colin Ross/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

assoc_log_ped <- read.table("PedVas_WTCCC_Eu_no_adult_trimmed_assoc_corrected.assoc.logistic", header = TRUE)
assoc_model_ped <- read.table("PedVas_WTCCC_Eu_no_adult_trimmed_assoc_corrected.model", header = TRUE)
assoc_HLA_log_ped <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_assoc_corrected.assoc.logistic", header = TRUE)
assoc_HLA_model_ped <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_assoc_corrected.model", header = TRUE)

assoc_log_adult <- read.table("PedVas_WTCCC_Eu_adult_trimmed_assoc_corrected.assoc.logistic", header = TRUE)
assoc_model_adult <- read.table("PedVas_WTCCC_Eu_adult_trimmed_assoc_corrected.model", header = TRUE)
assoc_HLA_log_adult <- read.table("PedVas_WTCCC_Eu_HLA_adult_trimmed_assoc_corrected.assoc.logistic", header = TRUE)
assoc_HLA_model_adult <- read.table("PedVas_WTCCC_Eu_HLA_adult_trimmed_assoc_corrected.model", header = TRUE)

assoc_log_comb <- read.table("WTCCC_Vasculitis_Eu_trimmed_assoc_corrected.assoc.logistic", header = TRUE)
assoc_model_comb <- read.table("WTCCC_Vasculitis_Eu_trimmed_assoc_corrected.model", header = TRUE)
assoc_HLA_log_comb <- read.table("WTCCC_Vasculitis_Eu_HLA_trimmed_assoc_corrected.assoc.logistic", header = TRUE)
assoc_HLA_model_comb <- read.table("WTCCC_Vasculitis_Eu_HLA_trimmed_assoc_corrected.model", header = TRUE)

# for .model files, select for dominant test (TEST == DOM)

assoc_dom_ped <- subset(assoc_model_ped, TEST == "DOM")
assoc_HLA_dom_ped <- subset(assoc_HLA_model_ped, TEST == "DOM")

assoc_dom_adult <- subset(assoc_model_adult, TEST == "DOM")
assoc_HLA_dom_adult <- subset(assoc_HLA_model_adult, TEST == "DOM")

assoc_dom_comb <- subset(assoc_model_comb, TEST == "DOM")
assoc_HLA_dom_comb <- subset(assoc_HLA_model_comb, TEST == "DOM")


# Association analysis of PEDIATRIC cases ---------------------------------

# make assoc table for manhattan and qq plots
res_assoc_log_ped <- subset(assoc_log_ped, TEST == "ADD")
res_assoc_log_ped <- select(res_assoc_log_ped, c("SNP", "CHR", "BP", "P")) # select variables for Manhattan plot
res_assoc_log_ped <- na.omit(res_assoc_log_ped) # remove incomplete cases
res_assoc_log_ped$SNP <- as.character(res_assoc_log_ped$SNP) # change SNP variable to character

man_plot_assoc_log_ped <- manhattan(res_assoc_log_ped, cex = 0.6, cex.axis = 0.8, ylim = c(0, 8), col = c("turquoise3", "black")) # generate manhatten plot

qq_log_ped <- qq(res_assoc_log_ped$P, main = "pediatric Q-Q plot - addative")



res_assoc_dom_ped <- subset(assoc_model_ped, TEST == "DOM")
res_assoc_dom_ped <- select(res_assoc_dom_ped, c("SNP", "CHR", "BP", "P"))
res_assoc_dom_ped <- na.omit(res_assoc_dom_ped)
res_assoc_dom_ped$SNP <- as.character(res_assoc_dom_ped$SNP)

man_plot_assoc_dom_ped 

qq_dom_ped <- qq(res_assoc_dom_ped$P, main = "pediatric Q-Q plot - dominant")



# Association analysis of ADULT cases ---------------------------------


res_assoc_log_adult <- subset(assoc_log_adult, TEST == "ADD")
res_assoc_log_adult <- select(assoc_log_adult, c("SNP", "CHR", "BP", "P")) # select variables for Manhattan plot
res_assoc_log_adult <- na.omit(res_assoc_log_adult) # remove incomplete cases
res_assoc_log_adult$SNP <- as.character(res_assoc_log_adult$SNP) # change SNP variable to character

man_plot_assoc_log_adult <- manhattan(res_assoc_log_adult, cex = 0.6, cex.axis = 0.8, ylim = c(0, 8), col = c("turquoise3", "black")) # generate manhatten plot

qq_log_adult <- qq(res_assoc_log_adult$P, main = "adult Q-Q plot - addative")

BP_adult <- select(res_assoc_log_adult, c("SNP", "BP"))

res_assoc_dom_adult <- subset(assoc_model_adult, TEST == "DOM")
res_assoc_dom_adult <- left_join(res_assoc_dom_adult, BP_adult, by = "SNP")
res_assoc_dom_adult <- select(res_assoc_dom_adult, c("SNP", "CHR", "BP", "P"))
res_assoc_dom_adult <- na.omit(res_assoc_dom_adult)
res_assoc_dom_adult$SNP <- as.character(res_assoc_dom_adult$SNP)

man_plot_assoc_dom_adult <- manhattan(res_assoc_dom_adult, cex = 0.6, cex.axis = 0.8, ylim = c(0, 8), col = c("turquoise3", "black")) # generate manhatten plot

qq_dom_adult <- qq(res_assoc_dom_adult$P, main = "adult Q-Q plot - dominant")


# Association analysis of COMBINED cases ---------------------------------

# make assoc table for manhattan and qq plots
res_assoc_log_comb <- subset(assoc_log_comb, TEST == "ADD")rm()
res_assoc_log_comb <- select(res_assoc_log_comb, c("SNP", "CHR", "BP", "P")) # select variables for Manhattan plot
res_assoc_log_comb <- na.omit(res_assoc_log_comb) # remove incomplete cases
res_assoc_log_comb$SNP <- as.character(res_assoc_log_comb$SNP) # change SNP variable to character

man_plot_assoc_log_comb <- manhattan(res_assoc_log_comb, cex = 0.6, cex.axis = 0.8, ylim = c(0, 10), col = c("turquoise3", "black")) # generate manhatten plot

qq_log_comb <- qq(res_assoc_log_comb$P, main = "combined Q-Q plot - addative")



res_assoc_dom_comb <- subset(assoc_model_comb, TEST == "DOM")
BP_comb <- select(res_assoc_log_comb, c("SNP", "BP"))
res_assoc_dom_comb <- left_join(res_assoc_dom_comb, BP_comb, by = "SNP")
res_assoc_dom_comb <- select(res_assoc_dom_comb, c("SNP", "CHR", "BP", "P"))
res_assoc_dom_comb <- na.omit(res_assoc_dom_comb)
res_assoc_dom_comb$SNP <- as.character(res_assoc_dom_comb$SNP)

man_plot_assoc_dom_comb <- manhattan(res_assoc_dom_comb, cex = 0.6, cex.axis = 0.8, ylim = c(0, 14), col = c("turquoise3", "black")) # generate manhatten plot

qq_dom_comb <- qq(res_assoc_dom_comb$P, main = "combined Q-Q plot - dominant")


# Association analysis of imputed HLA region ------------------------------

# make HLA scatterplot for PEDIATRIC ADDATIVE model
HLA_plot_log_ped <- assoc_HLA_log_ped
HLA_plot_log_ped <- subset(HLA_plot_log_ped, TEST == "ADD")
HLA_plot_log_ped$Type <- HLA_plot_log_ped$SNP # make new column for type of marker (ie. SNP, allele, or AA)
test <- HLA_plot$SNP

SNP <- c("[S][N][P].*", "[r][s].*", "[1][k][g].*")

pattern_SNP <- paste(SNP, collapse = "|") # create a regex pattern
pattern_AA <- paste("[A][A].*")
pattern_HLA <- paste("[H][L][A].*")   

HLA_plot_log_ped$Type <- str_replace(HLA_plot_log_ped$Type, pattern_SNP, "SNP")
HLA_plot_log_ped$Type <- str_replace(HLA_plot_log_ped$Type, pattern_AA, "AA")
HLA_plot_log_ped$Type <- str_replace(HLA_plot_log_ped$Type, pattern_HLA, "HLA")

HLA_plot_log_ped <- subset(HLA_plot_log_ped, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_plot_log_ped$Type <- factor(HLA_plot_log_ped$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_log_ped <- ggplot(HLA_plot_log_ped %>%
                   arrange(Type), 
                 aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))



# make HLA scatterplot for PEDIATRIC DOMINANT model
HLA_plot_dom_ped <- assoc_HLA_model_ped
HLA_plot_dom_ped <- subset(HLA_plot_dom_ped, TEST == "DOM")
HLA_plot_dom_ped$Type <- HLA_plot_dom_ped$SNP # make new column for type of marker (ie. SNP, allele, or AA)

HLA_plot_dom_ped$Type <- str_replace(HLA_plot_dom_ped$Type, pattern_SNP, "SNP")
HLA_plot_dom_ped$Type <- str_replace(HLA_plot_dom_ped$Type, pattern_AA, "AA")
HLA_plot_dom_ped$Type <- str_replace(HLA_plot_dom_ped$Type, pattern_HLA, "HLA")

HLA_plot_dom_ped <- subset(HLA_plot_dom_ped, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_BP <- select(HLA_plot_log_ped, c("SNP", "BP"))
HLA_plot_dom_ped <- left_join(HLA_plot_dom_ped, HLA_BP, by = "SNP")

HLA_plot_dom_ped$Type <- factor(HLA_plot_dom_ped$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_dom_ped <- ggplot(HLA_plot_dom_ped %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))



# make HLA scatterplot for ADULT ADDATIVE model
HLA_plot_log_adult <- assoc_HLA_log_adult
HLA_plot_log_adult <- subset(HLA_plot_log_adult, TEST == "ADD")
HLA_plot_log_adult$Type <- HLA_plot_log_adult$SNP # make new column for type of marker (ie. SNP, allele, or AA)

HLA_plot_log_adult$Type <- str_replace(HLA_plot_log_adult$Type, pattern_SNP, "SNP")
HLA_plot_log_adult$Type <- str_replace(HLA_plot_log_adult$Type, pattern_AA, "AA")
HLA_plot_log_adult$Type <- str_replace(HLA_plot_log_adult$Type, pattern_HLA, "HLA")

HLA_plot_log_adult <- subset(HLA_plot_log_adult, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_plot_log_adult$Type <- factor(HLA_plot_log_adult$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_log_adult <- ggplot(HLA_plot_log_adult %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))



# make HLA scatterplot for ADULT DOMINANT model
HLA_plot_dom_adult <- assoc_HLA_model_adult
HLA_plot_dom_adult <- subset(HLA_plot_dom_adult, TEST == "DOM")
HLA_plot_dom_adult$Type <- HLA_plot_dom_adult$SNP # make new column for type of marker (ie. SNP, allele, or AA)

HLA_plot_dom_adult$Type <- str_replace(HLA_plot_dom_adult$Type, pattern_SNP, "SNP")
HLA_plot_dom_adult$Type <- str_replace(HLA_plot_dom_adult$Type, pattern_AA, "AA")
HLA_plot_dom_adult$Type <- str_replace(HLA_plot_dom_adult$Type, pattern_HLA, "HLA")

HLA_plot_dom_adult <- subset(HLA_plot_dom_adult, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_BP <- select(HLA_plot_log_adult, c("SNP", "BP"))
HLA_plot_dom_adult <- left_join(HLA_plot_dom_adult, HLA_BP, by = "SNP")

HLA_plot_dom_adult$Type <- factor(HLA_plot_dom_adult$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_dom_adult <- ggplot(HLA_plot_dom_adult %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))


# make HLA scatterplot for COMBINED ADDATIVE model
HLA_plot_log_comb <- assoc_HLA_log_comb
HLA_plot_log_comb <- subset(HLA_plot_log_comb, TEST == "ADD")
HLA_plot_log_comb$Type <- HLA_plot_log_comb$SNP # make new column for type of marker (ie. SNP, allele, or AA)

HLA_plot_log_comb$Type <- str_replace(HLA_plot_log_comb$Type, pattern_SNP, "SNP")
HLA_plot_log_comb$Type <- str_replace(HLA_plot_log_comb$Type, pattern_AA, "AA")
HLA_plot_log_comb$Type <- str_replace(HLA_plot_log_comb$Type, pattern_HLA, "HLA")

HLA_plot_log_comb <- subset(HLA_plot_log_comb, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_plot_log_comb$Type <- factor(HLA_plot_log_comb$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_log_comb <- ggplot(HLA_plot_log_comb %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))



# make HLA scatterplot for COMBINED DOMINANT model
HLA_plot_dom_comb <- assoc_HLA_model_comb
HLA_plot_dom_comb <- subset(HLA_plot_dom_comb, TEST == "DOM")
HLA_plot_dom_comb$Type <- HLA_plot_dom_comb$SNP # make new column for type of marker (ie. SNP, allele, or AA)

HLA_plot_dom_comb$Type <- str_replace(HLA_plot_dom_comb$Type, pattern_SNP, "SNP")
HLA_plot_dom_comb$Type <- str_replace(HLA_plot_dom_comb$Type, pattern_AA, "AA")
HLA_plot_dom_comb$Type <- str_replace(HLA_plot_dom_comb$Type, pattern_HLA, "HLA")

HLA_plot_dom_comb <- subset(HLA_plot_dom_comb, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_BP <- select(HLA_plot_log_comb, c("SNP", "BP"))
HLA_plot_dom_comb <- left_join(HLA_plot_dom_comb, HLA_BP, by = "SNP")

HLA_plot_dom_comb$Type <- factor(HLA_plot_dom_comb$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_dom_comb <- ggplot(HLA_plot_dom_comb %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))


# Association analyses correcteding for HLA-DPB1 allelic dosage -----------

setwd("W:Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

assoc_log_ped_HLADPB1 <- read.table("PedVas_WTCCC_Eu_no_adult_trimmed_assoc_HLADPB1_corrected.assoc.logistic", header = TRUE) 
res_assoc_log_ped_HLADPB1 <- subset(assoc_log_ped_HLADPB1, TEST == "ADD")
res_assoc_log_ped_HLADPB1 <- select(res_assoc_log_ped_HLADPB1, c("SNP", "CHR", "BP", "P")) # select variables for Manhattan plot
res_assoc_log_ped_HLADPB1 <- na.omit(res_assoc_log_ped_HLADPB1) # remove incomplete cases
res_assoc_log_ped_HLADPB1$SNP <- as.character(res_assoc_log_ped_HLADPB1$SNP) # change SNP variable to character

man_plot_assoc_log_ped_HLADPB1 <- manhattan(res_assoc_log_ped_HLADPB1, cex = 0.6, cex.axis = 0.8, ylim = c(0, 8), col = c("turquoise3", "black")) # generate manhatten plot

qq_log_ped_HLADPB1 <- qq(res_assoc_log_ped_HLADPB1$P, main = "pediatric Q-Q plot - addative corrected for HLA_DPB1")



assoc_log_adult_HLADPB1 <- read.table("PedVas_WTCCC_Eu_adult_trimmed_assoc_HLADPB1_corrected.assoc.logistic", header = TRUE)
res_assoc_log_adult_HLADPB1 <- subset(assoc_log_adult_HLADPB1, TEST == "ADD")
res_assoc_log_adult_HLADPB1 <- select(res_assoc_log_adult_HLADPB1, c("SNP", "CHR", "BP", "P")) # select variables for Manhattan plot
res_assoc_log_adult_HLADPB1 <- na.omit(res_assoc_log_adult_HLADPB1) # remove incomplete cases
res_assoc_log_adult_HLADPB1$SNP <- as.character(res_assoc_log_adult_HLADPB1$SNP) # change SNP variable to character

man_plot_assoc_log_adult_HLADPB1 <- manhattan(res_assoc_log_adult_HLADPB1, cex = 0.6, cex.axis = 0.8, ylim = c(0, 8), col = c("turquoise3", "black")) # generate manhatten plot

qq_log_adult_HLADPB1 <- qq(res_assoc_log_adult_HLADPB1$P, main = "adult Q-Q plot - addative corrected for HLA_DPB1")


assoc_log_comb_HLADPB1 <- read.table("WTCCC_Vasculitis_Eu_trimmed_assoc_HLA_DPB1_corrected.assoc.logistic", header = TRUE)
res_assoc_log_comb_HLADPB1 <- subset(assoc_log_comb_HLADPB1, TEST == "ADD")
res_assoc_log_comb_HLADPB1 <- select(res_assoc_log_comb_HLADPB1, c("SNP", "CHR", "BP", "P")) # select variables for Manhattan plot
res_assoc_log_comb_HLADPB1 <- na.omit(res_assoc_log_comb_HLADPB1) # remove incomplete cases
res_assoc_log_comb_HLADPB1$SNP <- as.character(res_assoc_log_comb_HLADPB1$SNP) # change SNP variable to character

man_plot_assoc_log_comb_HLADPB1 <- manhattan(res_assoc_log_comb_HLADPB1, cex = 0.6, cex.axis = 0.8, ylim = c(0, 8), col = c("turquoise3", "black")) # generate manhatten plot

qq_log_comb_HLADPB1 <- qq(res_assoc_log_comb_HLADPB1$P, main = "combined Q-Q plot - addative corrected for HLA_DPB1")


# anslysis of HLA region after correcting for HLA_DPB1 allelic dosage 

assoc_log_ped_HLA_HLADPB1 <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_assoc_HLADPB1_corrected.assoc.logistic", header = TRUE)
assoc_log_adult_HLA_HLADPB1 <- read.table("PedVas_WTCCC_Eu_HLA_adult_trimmed_assoc_HLADPB1_corrected.assoc.logistic", header = TRUE)
assoc_log_comb_HLA_HLADPB1 <- read.table("WTCCC_Vasculitis_Eu_HLA_trimmed_assoc_HLADPB1_corrected.assoc.logistic", header = TRUE)


# make HLA scatterplot for PEDIATRIC ADDATIVE model
HLA_plot_log_ped_HLADPB1 <- assoc_log_ped_HLA_HLADPB1
HLA_plot_log_ped_HLADPB1 <- subset(HLA_plot_log_ped_HLADPB1, TEST == "ADD")
HLA_plot_log_ped_HLADPB1$Type <- HLA_plot_log_ped_HLADPB1$SNP # make new column for type of marker (ie. SNP, allele, or AA)

SNP <- c("[S][N][P].*", "[r][s].*", "[1][k][g].*")

pattern_SNP <- paste(SNP, collapse = "|") # create a regex pattern
pattern_AA <- paste("[A][A].*")
pattern_HLA <- paste("[H][L][A].*")   

HLA_plot_log_ped_HLADPB1$Type <- str_replace(HLA_plot_log_ped_HLADPB1$Type, pattern_SNP, "SNP")
HLA_plot_log_ped_HLADPB1$Type <- str_replace(HLA_plot_log_ped_HLADPB1$Type, pattern_AA, "AA")
HLA_plot_log_ped_HLADPB1$Type <- str_replace(HLA_plot_log_ped_HLADPB1$Type, pattern_HLA, "HLA")

HLA_plot_log_ped_HLADPB1 <- subset(HLA_plot_log_ped_HLADPB1, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_plot_log_ped_HLADPB1$Type <- factor(HLA_plot_log_ped_HLADPB1$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_log_ped_HLADPB1 <- ggplot(HLA_plot_log_ped_HLADPB1 %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))







### MAKE TOP HITS LIST FOR DOWNSTREM ANALYSES

setwd("W:/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

tophits_log_ped <- subset(res_assoc_log_ped, P < 0.001)
tophits_log_ped <- select(tophits_log_ped, -c("A1"))
# write.table(tophits_log_ped, "tophits_log_ped.txt", row.names = FALSE, quote = FALSE)

tophits_log_adult <- subset(res_assoc_log_adult, P < 0.001)
tophits_log_adult <- select(tophits_log_adult, -c("A1"))
# write.table(tophits_log_adult, "tophits_log_adult.txt", row.names = FALSE, quote = FALSE)

tophits_log_comb <- subset(res_assoc_log_comb, P < 0.001)
tophits_log_comb <- select(tophits_log_comb, -c("A1"))
# write.table(tophits_log_comb, "tophits_log_comb.txt", row.names = FALSE, quote = FALSE)






### Read in frequency files
setwd( "W:/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022/Freq_files")

freq_ped <- read.table("PedVas_WTCCC_Eu_qc_no_adult_trimmed_freq.frq.cc", header = TRUE)
freq_ped_HLA <- read.table("PedVas_WTCCC_Eu_qc_HLA_no_adult_trimmed_freq.frq.cc", header = TRUE)

freq_adult <- read.table("PedVas_WTCCC_Eu_adult_trimmed_freq.frq.cc", header = TRUE)
freq_adult_HLA <- read.table("PedVas_WTCCC_Eu_qc_HLA_adult_trimmed_freq.frq.cc", header = TRUE)

freq_comb <- read.table("WTCCC_Vasculitis_Eu_trimmed_freq.frq.cc", header = TRUE)
freq_comb_HLA <- read.table("WTCCC_Vasculitis_Eu_HLA_trimmed_freq.frq.cc", header = TRUE)


# make CHR tables for locuszoom
setwd("W:/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022/LocusZoom")

ped_log_CHR1 <- subset(res_assoc_log_ped, CHR == 1)
# write.table(ped_log_CHR1, "ped_log_CHR1.txt", row.names = FALSE, sep = "\t")
ped_log_CHR5 <- subset(res_assoc_log_ped, CHR == 5)
# write.table(ped_log_CHR5, "ped_log_CHR5.txt", row.names = FALSE, sep = "\t")
ped_log_CHR14 <- subset(res_assoc_log_ped, CHR == 14)
# write.table(ped_log_CHR14, "ped_log_CHR14.txt", row.names = FALSE, sep = "\t")
ped_log_CHR19 <- subset(res_assoc_log_ped, CHR == 19)
# write.table(ped_log_CHR19, "ped_log_CHR19.txt", row.names = FALSE, sep = "\t")







CHR_19_SNPs_vec <- c("rs62132279", "rs72982187", "rs149482902", "rs147113697", "rs76537326", "rs76111962", "rs74913538",
                 "rs77733763", "rs62132294", "rs3826947", "rs62132295", "rs62132296", "rs62132297", "rs2301879")

CHR_19_SNPs <- as.data.frame(CHR_19_SNPs_vec)
colnames(CHR_19_SNPs) <- "SNP"

CHR_19_qry <- left_join(CHR_19_SNPs, ped_log_CHR19, by = "SNP")


CHR_14_SNPs_vec <- c("rs5007068", "rs5007067", "rs5007065", "rs5007064", "rs5007063", "rs5007062")
CHR_14_SNPs <- as.data.frame(CHR_14_SNPs_vec)
colnames(CHR_14_SNPs) <- "SNP"

CHR_14_qry <- left_join(CHR_14_SNPs, ped_log_CHR14, by = "SNP")






# ANCA subset ANALYSIS -----------------------------------------------------------

setwd("W:/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022")

assoc_HLA_log_ped_PR3 <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_PR3_assoc.assoc.logistic", header = TRUE)

HLA_plot_log_ped_PR3 <- assoc_HLA_log_ped_PR3
HLA_plot_log_ped_PR3 <- subset(HLA_plot_log_ped_PR3, TEST == "ADD")
HLA_plot_log_ped_PR3$Type <- HLA_plot_log_ped_PR3$SNP # make new column for type of marker (ie. SNP, allele, or AA)

SNP <- c("[S][N][P].*", "[r][s].*", "[1][k][g].*")

pattern_SNP <- paste(SNP, collapse = "|") # create a regex pattern
pattern_AA <- paste("[A][A].*")
pattern_HLA <- paste("[H][L][A].*")   

HLA_plot_log_ped_PR3$Type <- str_replace(HLA_plot_log_ped_PR3$Type, pattern_SNP, "SNP")
HLA_plot_log_ped_PR3$Type <- str_replace(HLA_plot_log_ped_PR3$Type, pattern_AA, "AA")
HLA_plot_log_ped_PR3$Type <- str_replace(HLA_plot_log_ped_PR3$Type, pattern_HLA, "HLA")

HLA_plot_log_ped_PR3 <- subset(HLA_plot_log_ped_PR3, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_plot_log_ped_PR3$Type <- factor(HLA_plot_log_ped_PR3$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_log_ped_PR3 <- ggplot(HLA_plot_log_ped_PR3 %>%
                           arrange(Type), 
                         aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))




assoc_HLA_log_ped_MPO <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_MPO_assoc.assoc.logistic", header = TRUE)

HLA_plot_log_ped_MPO <- assoc_HLA_log_ped_MPO
HLA_plot_log_ped_MPO <- subset(HLA_plot_log_ped_MPO, TEST == "ADD")
HLA_plot_log_ped_MPO$Type <- HLA_plot_log_ped_MPO$SNP # make new column for type of marker (ie. SNP, allele, or AA)


HLA_plot_log_ped_MPO$Type <- str_replace(HLA_plot_log_ped_MPO$Type, pattern_SNP, "SNP")
HLA_plot_log_ped_MPO$Type <- str_replace(HLA_plot_log_ped_MPO$Type, pattern_AA, "AA")
HLA_plot_log_ped_MPO$Type <- str_replace(HLA_plot_log_ped_MPO$Type, pattern_HLA, "HLA")

HLA_plot_log_ped_MPO <- subset(HLA_plot_log_ped_MPO, Type == "AA" | Type == "HLA" | Type == "SNP")

HLA_plot_log_ped_MPO$Type <- factor(HLA_plot_log_ped_MPO$Type, levels = c("SNP", "AA", "HLA"))

HLA_sp_log_ped_MPO <- ggplot(HLA_plot_log_ped_MPO %>%
                               arrange(Type), 
                             aes(BP/1000000, -log(P, base = 10), colour = Type)) + 
  geom_point() + 
  geom_hline(yintercept = -log(5*10^(-8), base = 10), linetype = "dashed", color = "black") +
  scale_color_manual(values = c( "#A6ACAF", "#2E86C1", "#C0392B")) +
  theme_classic() +
  labs(x = "Position on Chr 6 (Mb)", y = "-log[10]p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10))

### Read in frequency files of ANCA subgroups
setwd( "W:/Colin Ross Lab/2. Individual Data/Kristen Gibson/PedVas_GWAS_2019/AnR_revision_data_2022/Freq_files")

freq_ped_HLA_PR3 <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_PR3_freq.frq.cc", header = TRUE)
freq_ped_HLA_MPO <- read.table("PedVas_WTCCC_Eu_HLA_no_adult_trimmed_MPO_freq.frq.cc", header = TRUE)



