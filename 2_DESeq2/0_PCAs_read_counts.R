#BiocManager::install("DESeq2")
rm(list=ls())
pacman::p_load(PCAtools, tidyverse, DESeq2)

cts <- read.delim("salmon.merged.gene_counts.tsv", sep = "\t")

cts_matrix <- cts %>% 
  select(!gene_name) %>% 
  column_to_rownames("gene_id")
  
colnames(cts_matrix) <- gsub("^X", "", colnames(cts_matrix))
cts_matrix <- as.matrix(cts_matrix)

coldata <- read.csv("data/sample_sheet.csv")
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

rld <- vst(dds)

plotPCA(rld, intgroup="timepoint", ntop=20000)

# deseqdataset to matrix

#pca <- prcomp(t(assay(rld)))
coldata2 <- coldata %>% 
  select(sample, timepoint, accession, inoculum) %>% 
  mutate(inoculum=ifelse(inoculum=="pdb_mock", "Mock", "SS")) %>% 
  rename("inoculum"="Inoculum", "accession"="Accession")

p <- pca(assay(rld), metadata = coldata2, removeVar = 0.1)

screeplot(p, axisLabSize = 18, titleLabSize = 22)

p1 <- biplot(p,
       x="PC1", 
       y="PC2",
       showLoadings = F, 
       colby="Accession", 
       shape="Inoculum",
       labSize = NA, pointSize = 5, sizeLoadingsNames = 5,
       # ellipse config
       lab=NA,
       ellipse = F,
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)+
  theme(legend.position = "top", 
        axis.text=element_text(size=14, color="black"))

p2 <- biplot(p,
       x="PC2", 
       y="PC3",
       showLoadings = F, 
       colby="Accession", 
       shape="Inoculum",
       labSize = NA, pointSize = 5, sizeLoadingsNames = 5,
       # ellipse config
       lab=NA,
       ellipse = F,
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)+
  theme(legend.position = "top", 
        axis.text=element_text(size=14, color="black"))
svg("figures/fig_2/PCA.svg",
    width = 11, height = 7)
(p_tot <- ggpubr::ggarrange(p1, p2 , nrow=1, common.legend = T))
dev.off()

ggsave(p_tot, filename="ig_2/PCA.png",
       width = 10, height = 11, units = "in", dpi = 900, bg="white")
