# Sample_code_for_CeMM


### DADA2_PAN_RNAseq_DESq2.R
Whole blood transcriptome analysis comparing patients with DADA2 and PAN
    
Quality control of Fastq files was done using FastQC and MultiQC then aligned to the human genome (GRCh38.93) using STAR. Read count tables were      generated with HTSeq-count. Read counts of globin genes were removed bioinformatically. Raw RNA-Seq counts were normalized for library size and heteroskedasticity by variance stabilizing transformation with the vst function available in the DESeq2 package. Differential expression analysis was performed using DESeq2
    

### PedVas_40plex_analysis
Multiple factor analysis of 40 cytokines and clinical metadata of pediatric patients with ANCA-associated vasculitis 


### PedVas_GWAS_analysis.R
Genome-wide association study of pediatric ANCA-associated vasculitis

Association analyses were done using .plink, the attached code is primarily for the visualization of the data
    
    
### PedVas_pilot_urine_massspec.R
Proteomic anaysis of urine from pediatric vasculits patients with varying stages of renal disease
    
