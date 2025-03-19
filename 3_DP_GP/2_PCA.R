# Make PCA for the focal DEGs, cluster tempraol trajectories
# Severin Einspanier

pacman::p_load(PCAtools, tidyverse, DESeq2)
# for this, each gene is a sample and the time point measurements are the variables



counts <- read.table( "in_data/CONSENSUS_48hpi_DEGs.txt",
                      header=T) %>% 
  column_to_rownames("GeneID") %>% 
  arrange(rownames(.)) %>% 
  t()

# deseqdataset to matrix

# col= Genotype 
coldata <- read.table( "in_data/CONSENSUS_48hpi_DEGs.txt",
                      header=T) %>%
  mutate(genotype = sub("_.*", "", GeneID)) %>% 
  select(GeneID, genotype) %>%
  column_to_rownames("GeneID") %>% 
  arrange(rownames(.))

# Ensure coldata is a data frame
#coldata <- as.data.frame(coldata)


# Sort coldata by row names
#coldata <- coldata[order(rownames(coldata)), ]

# Check column names of counts matrix
colnames(counts)

# Check row names of coldata data frame
rownames(coldata)

# Now you can run the PCA
pca_result <- pca(counts, coldata, removeVar = 0.1)


screeplot(p, axisLabSize = 18, titleLabSize = 22)

svg(paste0("figures/fig_4/", Sys.Date(), "_PCA.svg"), 
    width = 8, height = 7)
(p1 <- biplot(p,
             x="PC1", 
             y="PC2",
             showLoadings = F, 
             colby="genotype", 
             #shape="Inoculum",
             labSize = NA, pointSize = 5, sizeLoadingsNames = 5,
             # ellipse config
             lab=NA,
             ellipse = T,
             legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)+
  theme(legend.position = "top", 
        axis.text=element_text(size=14, color="black"))
)

dev.off()
