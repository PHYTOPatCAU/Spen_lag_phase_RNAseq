# Perform GSEA analysis for multiple genotypes
# repeat for 48 hpi etc. 
# Severin Einspanier

rm(list=ls())
library(tidyverse)
library(clusterProfiler)

# Load GO term databases (MF only as in original code)
TERM2GENE <- read.csv("term2gene_BP.csv")
TERM2NAME <- read.csv("term2name_BP.csv")

# Genotypes to analyze
genotypes <- c("1282", "1809", "1941")

# Loop through each genotype
for (geno in genotypes) {
  # Read and preprocess data
  input_file <- sprintf("expr_dat_%s_inf_mock.csv", geno)
  gene_data <- read.csv(input_file) %>% 
    filter(timepoint == "24hpi") %>% 
    dplyr::select(GeneID, log2FoldChange) %>% 
    arrange(-log2FoldChange)
  
  # Create ranked gene list
  gene_list <- gene_data$log2FoldChange
  names(gene_list) <- as.character(gene_data$GeneID)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA
  gsea_result <- GSEA(
    geneList = gene_list,
    pvalueCutoff = 0.05,
    eps = 1e-6,
    gson=NULL,
    minGSSize=10,
    exponent=.5,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    pAdjustMethod = "BH",
    scoreType = "pos",
    by = "fgsea"
  )
  
  # Create and store results dataframe
  results_df <- as.data.frame(gsea_result) %>%
    arrange(-NES) %>% 
    dplyr::select(ID, Description, NES, p.adjust, setSize)
  
  assign(paste0("gsea_results", geno), results_df)
}

# Results are now in:
# gsea_results1282, gsea_results1809, and gsea_results1941

# print for table

write.csv(gsea_results1282, "tables/gsea_1282_24hpi_BP.csv")
write.csv(gsea_results1809, "tables/gsea_1809_24hpi_BP.csv")
write.csv(gsea_results1941, "tables/gsea_1941_24hpi_BP.csv")
