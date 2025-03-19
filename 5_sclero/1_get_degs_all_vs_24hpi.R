# Overview on fungal reads per timepoint 
# Lets make DEG contrast to 24 hpi 
# per genotype 

rm(list=ls())
library(tidyverse)
library(DESeq2)
options(scipen = 999, big.mark = ",")

cts <- read.delim("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/nf_core/sclero/star_salmon/salmon.merged.gene_counts.tsv", 
                  sep="\t")

cts_matrix <- cts %>% 
  select(!gene_name) %>% 
  column_to_rownames("gene_id")

colnames(cts_matrix) <- gsub("^X", "", colnames(cts_matrix))
cts_matrix <- as.matrix(cts_matrix)

coldata <- read.csv("C:/Users/suaph281/Desktop/GitLab/spen_lag_phase_rnaseq/data/sample_sheet.csv")
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
DESeq2::resultsNames(dds)

# ORDER: contrasts1mock contrasts1sclero contrasts2mock contrasts2sclero

TreatmentMatrix <- rbind(
  # 1282 time vs 24 hpi INFECTED only
  "LA1282_48_vs_24" = c(0, 0, 0, 0, -1, 1, 0, 0,  
                        0, 0, 0, 0, 0,  0, 0, 0, 
                        0, 0, 0, 0, 0,  0, 0, 0),
  "LA1282_72_vs_24"  = c(0, 0, 0, 0, -1, 0, 1, 0,  
                         0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0),
  "LA1282_96_vs_24"  = c(0, 0, 0, 0, -1, 0, 0, 1,  
                         0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0, 0, 0, 0, 0),
  # 1809 time vs 24 hpi INFECTED only
  "LA1809_48_vs_24" = c(0, 0, 0, 0, 0, 0, 0, 0,  
                        0, 0, 0, 0, -1,1, 0, 0, 
                        0, 0, 0, 0,  0,0, 0, 0),
  "LA1809_72_vs_24"  = c(0, 0, 0, 0, 0,0, 0, 0,  
                         0, 0, 0, 0,-1,0, 1, 0, 
                         0, 0, 0, 0, 0,0, 0, 0),
  "LA1809_96_vs_24"  = c(0, 0, 0, 0, 0,0, 0, 0,  
                         0, 0, 0, 0,-1,0, 0, 1, 
                         0, 0, 0, 0, 0,0, 0, 0),
  # 72 hpi
  "LA1941_48_vs_24" = c(0, 0, 0, 0, 0,  0, 0, 0,  
                        0, 0, 0, 0, 0,  0, 0, 0, 
                        0, 0, 0, 0, -1, 1, 0, 0),
  "LA1941_72_vs_24"  = c(0, 0, 0, 0, 0, 0, 0, 0,  
                         0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0,-1, 0, 1, 0),
  "LA1941_96_vs_24"  = c(0, 0, 0, 0, 0, 0, 0, 0,  
                         0, 0, 0, 0, 0, 0, 0, 0, 
                         0, 0, 0, 0,-1, 0, 0, 1)

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

write.csv(DEG_long_df, "C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/sclero/expr_all_vs_24hpi_sclero.csv", 
          row.names = FALSE)

(p_tot <- ggplot(DEG_long_df, aes(x = log2FoldChange, y = -log10(pvalue), color = padj < 0.05)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ Contrast) +
  theme_minimal() +
  labs(title = "Volcano Plots for All Contrasts")
)

ggsave(p_tot, filename="C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing_2/figures/2025_02_12_DEGs_sclero_against24hpi.png",
       bg="white", dpi=300)
