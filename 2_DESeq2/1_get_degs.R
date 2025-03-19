# Determine DEGs for each time point
# The genotypes are tested seperately
# Contrast: Inf-Mock.
# Severin Einspanier 2025_01_23
# 0. Load libraries
rm(list=ls())
pacman::p_load(DESeq2, tidyverse)
options(scipen = 999, big.mark = ",")

# 1. get filtered reads
# input: Proteome-filtered counts (see 0_filter_read_table.r)
counts_data <- read.csv("LA1809_filtered_counts.csv",
                        row.names = 1)
colnames(counts_data) <- gsub("^X", "", colnames(counts_data))
cts_matrix <- as.matrix(counts_data)

coldata <- read.csv("sample_sheet.csv") %>% 
  filter(accession %in% c("LA1809")) %>% 
  select(sample, accession, timepoint, inoculum)
rownames(coldata)=coldata$sample

if (!all(rownames(coldata) == colnames(counts_data))) {
  # Reorder the columns of cts_matrix to match the order of rownames in coldata
  counts_data <- counts_data[, rownames(coldata)]
}


coldata$inoculum <- as.factor(coldata$inoculum)
coldata$accession <- as.factor(coldata$accession)
coldata$timepoint <- as.factor(coldata$timepoint)

dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = coldata,
                              design = ~  inoculum + timepoint)
dds$contrasts <- factor(paste0(dds$accession, dds$inoculum, dds$timepoint))
design(dds) <- ~0+contrasts
dds <- DESeq2::DESeq(dds)

# Prefilter the dataset
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep, ]
DESeq2::resultsNames(dds)

# "contrastsLA1809pdb_mock24" "contrastsLA1809pdb_mock48""contrastsLA1809pdb_mock72""contrastsLA1809pdb_mock96"
#"contrastsLA1809ss24""contrastsLA1809ss48" "contrastsLA1809ss72""contrastsLA1809ss96" 


TreatmentMatrix <- rbind(
  "inf_mock_24" = c(-1, 0, 0, 0, 1, 0, 0,0),
  "inf_mock_48" = c(0, -1, 0, 0, 0, 1, 0,0),
  "inf_mock_72" = c(0, 0, -1, 0, 0, 0, 1,0),
  "inf_mock_96" = c(0, 0, 0, -1, 0, 0, 0,1)
)




# Now, get the list of TFs.
TFs <- read.csv("spen_TFs.txt", header=T)
colnames(TFs)=c("Gene", "TF")

# Return the selected contrast
selected_contrast <- TreatmentMatrix["inf_mock_96", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_96 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "96hpi")%>% 
  rownames_to_column("GeneID")

# repeat for 24, 48, 72 hpi

selected_contrast <- TreatmentMatrix["inf_mock_72", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_72 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "72hpi")%>% 
  rownames_to_column("GeneID")

selected_contrast <- TreatmentMatrix["inf_mock_48", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_48 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "48hpi")%>% 
  rownames_to_column("GeneID")

selected_contrast <- TreatmentMatrix["inf_mock_24", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_24 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "24hpi") %>% 
  rownames_to_column("GeneID")

merged_data_TFs1809 <- rbind(df_24, df_48, df_72, df_96) %>% 
  drop_na(padj)





# repeat for 1282

# 1. get filtered reads

counts_data <- read.csv("LA1282_filtered_counts.csv")
colnames(counts_data) <- gsub("^X", "", colnames(counts_data))
rownames(counts_data) <- counts_data[,1]
cts_matrix <- as.matrix(counts_data)

coldata <- read.csv("sample_sheet.csv") %>% 
  filter(accession %in% c("LA1282")) %>% 
  select(sample, accession, timepoint, inoculum)
rownames(coldata)=coldata$sample

if (!all(rownames(coldata) == colnames(counts_data))) {
  # Reorder the columns of cts_matrix to match the order of rownames in coldata
  counts_data <- counts_data[, rownames(coldata)]
}


coldata$inoculum <- as.factor(coldata$inoculum)
coldata$accession <- as.factor(coldata$accession)
coldata$timepoint <- as.factor(coldata$timepoint)


dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = coldata,
                              design = ~  inoculum + timepoint)
dds$contrasts <- factor(paste0(dds$accession, dds$inoculum, dds$timepoint))
design(dds) <- ~0+contrasts
dds <- DESeq2::DESeq(dds)

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep, ]
DESeq2::resultsNames(dds)


# Return the selected contrast
selected_contrast <- TreatmentMatrix["inf_mock_96", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")


# Now, get the list of TFs.
TFs <- read.csv("spen_TFs.txt", header=T)
colnames(TFs)=c("Gene", "TF")

df_96 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "96hpi")%>% 
  rownames_to_column("GeneID")

# repeat for 24, 48, 72 hpi

selected_contrast <- TreatmentMatrix["inf_mock_72", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_72 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "72hpi")%>% 
  rownames_to_column("GeneID")

selected_contrast <- TreatmentMatrix["inf_mock_48", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_48 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "48hpi")%>% 
  rownames_to_column("GeneID")

selected_contrast <- TreatmentMatrix["inf_mock_24", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_24 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "24hpi") %>% 
  rownames_to_column("GeneID")

# merge 
merged_data_TFs1282 <- df_24 %>% 
  bind_rows(df_48, df_72, df_96) %>% 
  drop_na(padj)


# 1941

# 1. get filtered reads

counts_data <- read.csv("LA1941_filtered_counts.csv",
                        row.names = 1)
colnames(counts_data) <- gsub("^X", "", colnames(counts_data))
cts_matrix <- as.matrix(counts_data)

coldata <- read.csv("sample_sheet.csv") %>% 
  filter(accession %in% c("LA1941")) %>% 
  select(sample, accession, timepoint, inoculum)
rownames(coldata)=coldata$sample

if (!all(rownames(coldata) == colnames(counts_data))) {
  # Reorder the columns of cts_matrix to match the order of rownames in coldata
  counts_data <- counts_data[, rownames(coldata)]
}


coldata$inoculum <- as.factor(coldata$inoculum)
coldata$accession <- as.factor(coldata$accession)
coldata$timepoint <- as.factor(coldata$timepoint)

dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = coldata,
                              design = ~  inoculum + timepoint)
dds$contrasts <- factor(paste0(dds$accession, dds$inoculum, dds$timepoint))
design(dds) <- ~0+contrasts
dds <- DESeq2::DESeq(dds)

# Prefilter the dataset
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep, ]
DESeq2::resultsNames(dds)


# Return the selected contrast
selected_contrast <- TreatmentMatrix["inf_mock_96", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")


# Now, get the list of TFs.
TFs <- read.csv("data/spen_TFs.txt", header=T)
colnames(TFs)=c("Gene", "TF")

df_96 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "96hpi")%>% 
  rownames_to_column("GeneID")

# repeat for 24, 48, 72 hpi

selected_contrast <- TreatmentMatrix["inf_mock_72", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_72 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "72hpi")%>% 
  rownames_to_column("GeneID")

selected_contrast <- TreatmentMatrix["inf_mock_48", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_48 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "48hpi")%>% 
  rownames_to_column("GeneID")

selected_contrast <- TreatmentMatrix["inf_mock_24", ]
return_contrast <- DESeq2::lfcShrink(dds, contrast = selected_contrast, type="ashr")
df_24 <- as.data.frame(return_contrast) %>% 
  #filter(row.names(return_contrast) %in% TFs$Gene) %>%
  mutate(timepoint = "24hpi") %>% 
  rownames_to_column("GeneID")

# merge 
merged_data_TFs1941 <- rbind(df_24, df_48, df_72, df_96) %>% 
  drop_na(padj)


# write 

write.csv(merged_data_TFs1809, "expr_dat_1809_inf_mock.csv")
write.csv(merged_data_TFs1941, "data/expr_dat_1941_inf_mock.csv")
write.csv(merged_data_TFs1282, "expr_dat_1282_inf_mock.csv")
