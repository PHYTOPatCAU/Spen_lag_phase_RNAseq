# Get DEGs between genotypes at each timepoint 
# in MOCK conditions
# Severin Einspanier
# 2025_02_11#
rm(list=ls())
library(DESeq2)
library(PCAtools)
library(tidyverse)
options(scipen = 999, big.mark = ",")

cts <- read.delim("salmon.merged.gene_counts.tsv", sep = "\t")

cts_matrix <- cts %>% 
  select(!gene_name) %>% 
  column_to_rownames("gene_id")

colnames(cts_matrix) <- gsub("^X", "", colnames(cts_matrix))
cts_matrix <- as.matrix(cts_matrix)

coldata <- read.csv("sample_sheet.csv")
rownames(coldata)=coldata$sample

if (!all(rownames(coldata) == colnames(cts_matrix))) {
  # Reorder the columns of cts_matrix to match the order of rownames in coldata
  cts_matrix <- cts_matrix[, rownames(coldata)]
}

all(rownames(coldata) == colnames(cts_matrix))

coldata$inoculum <- as.factor(coldata$inoculum)
coldata$accession <- as.factor(coldata$accession)
coldata$timepoint <- as.factor(coldata$timepoint)


# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(cts_matrix),
                              colData = coldata,
                              design = ~ accession + inoculum + timepoint)




dds$contrasts <- factor(paste0(dds$accession, dds$inoculum, dds$timepoint))
design(dds) <- ~0+contrasts
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 4

dds <- DESeq2::DESeq(dds)

# Prefilter the dataset
# smallestGroupSize <- 4
# keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
# dds <- dds[keep, ]
DESeq2::resultsNames(dds)

# ORDER: contrasts1mock contrasts1sclero contrasts2mock contrasts2sclero
TreatmentMatrix <- rbind(
  # 24 hpi
  "LA1282_vs_LA1809_24" = c(  1, 0, 0, 0, 0, 0, 0, 0,  # LA1282ss24 - LA1809ss24
                             -1, 0, 0, 0, 0, 0, 0, 0, 
                              0, 0, 0, 0, 0, 0, 0, 0),
  "LA1282_vs_LA1941_24" = c( 1, 0, 0, 0, 0, 0, 0, 0,  # LA1282ss24 - LA1941ss24
                             0, 0, 0, 0, 0, 0, 0, 0, 
                             -1, 0, 0, 0,0, 0, 0, 0),
  "LA1941_vs_LA1809_24" = c( 0, 0, 0, 0, 0, 0, 0, 0,  # LA1809ss24 - LA1941ss24
                             -1, 0, 0, 0,0, 0, 0, 0, 
                             1, 0, 0, 0, 0, 0, 0, 0),
  # 48 hpi
  "LA1282_vs_LA1809_48" = c(0, 1, 0, 0, 0, 0, 0, 0,  
                            0,-1, 0, 0, 0,0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 0),
  "LA1282_vs_LA1941_48" = c(0, 1, 0, 0, 0, 0, 0, 0,  
                            0, 0, 0, 0, 0, 0, 0, 0, 
                            0,-1, 0, 0, 0, 0, 0, 0),
  "LA1941_vs_LA1809_48" = c(0, 0, 0, 0, 0, 0, 0, 0,  
                            0,-1, 0, 0, 0, 0, 0, 0, 
                            0, 1, 0, 0, 0, 0, 0, 0),
  # 72 hpi
  "LA1282_vs_LA1809_72" = c(0, 0, 1, 0, 0, 0, 0, 0,  
                            0, 0,-1, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 0),
  "LA1282_vs_LA1941_72" = c(0, 0, 1, 0, 0, 0, 0, 0,  
                            0, 0, 0, 0, 0, 0, 0, 0, 
                            0, 0,-1, 0, 0, 0, 0, 0),
  "LA1941_vs_LA1809_72" = c(0, 0, 0, 0, 0, 0, 0, 0,  
                            0, 0,-1, 0, 0, 0, 0, 0, 
                            0, 0, 1, 0, 0, 0, 0, 0),
  # 96 hpi
  "LA1282_vs_LA1809_96" = c(0, 0, 0, 1, 0, 0, 0, 0,  
                            0, 0, 0,-1, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 0),
  "LA1282_vs_LA1941_96" = c(0, 0, 0, 1, 0, 0, 0, 0,  
                            0, 0, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0,-1, 0, 0, 0, 0),
  "LA1941_vs_LA1809_96" = c(0, 0, 0, 0, 0, 0, 0, 0,  
                            0, 0, 0,-1, 0, 0, 0, 0, 
                            0, 0, 0, 1, 0, 0, 0, 0)
)

DEG_counts <- list()

# Iterate over each contrast in the TreatmentMatrix
for (contrast_name in rownames(TreatmentMatrix)) {
  selected_contrast <- TreatmentMatrix[contrast_name, ]
  DEG_counts[[contrast_name]] <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type = "ashr")
}

DEG_long_df <- map_dfr(names(DEG_counts), function(contrast_name) {
  deg_result <- DEG_counts[[contrast_name]]  # Extract DESeq2 results
  
  # Convert to dataframe and add contrast name
  df <- as.data.frame(deg_result) %>%
    rownames_to_column(var = "GeneID") %>%
    mutate(Contrast = contrast_name)
  
  return(df)
})


write.csv(DEG_long_df, "expr_contrasts_genotypes_mock.csv", 
          row.names = FALSE)

p_tot <- ggplot(DEG_long_df, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Contrast) +
  theme_minimal() +
  labs(title = "Volcano Plots for All Contrasts")


ggsave(p_tot, filename="2025_02_11_DEGs_between_genotypes_mock.png",
       bg="white", dpi=300)
