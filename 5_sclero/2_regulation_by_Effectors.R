# Test whether my genes are regulated by effectors

focal_genes <- c("Sopen01g026160", "Sopen03g033250",  "Sopen03g041120",  "Sopen12g023430")
data_edges <- read.table("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/GRN/sclero/LA1282_infected_with_effectors_grn_network_edges.txt",
                         header=T, sep="\t") %>% 
  filter(targetGene %in%focal_genes)



data_edges_target_filtered_more <- data_edges
# summary 
data_edges_target_filtered_more %>% 
  group_by(targetGene) %>% 
  summarise(n=n(), sum=sum(weight), mean=mean(weight))

# get how many regulatory genes 
data_edges_target_filtered_more %>% 
  group_by(regulatoryGene) %>% 
  summarise(n=n()) %>% 
  group_by(n) %>% 
  summarise(n=n())

# #1    88
# 2    30
# 3    20
# 4     5

put_regulators <- data_edges_target_filtered_more %>% 
  group_by(regulatoryGene) %>% 
  summarise(n=n()) %>% 
  group_by(n) %>% 
  filter(n>3)



# only plant is regulating those target genes 

# are there any interaction between effectors and host genes?
effectors <- read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/sclero/effectors.csv") %>% 
  mutate(gene = tolower(Version.Two.ID)) %>% 
  dplyr::select(gene) 

data_edgeseffects <- read.table("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/GRN/sclero/LA1282_infected_with_effectors_grn_network_edges.txt",
                         header=T, sep="\t") %>% 
  filter(regulatoryGene %in%effectors$gene) 

# Are the tomato genes in the LRTs?
lrt1282 <- read.delim("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/DP_GP/in_data/LA1282_LFC_LRT.txt") %>% 
  filter(gene %in%data_edgeseffects$targetGene) 
# nor really interesting 


lrt1282_regulators <- read.delim("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/DP_GP/in_data/LA1282_LFC_LRT.txt") %>% 
  filter(gene %in%put_regulators$regulatoryGene)

genotype_mock <- read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/DeSeq/expr_contrasts_genotypes_mock.csv") %>% 
  filter(GeneID %in%put_regulators$regulatoryGene)

genotype_infected <- read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/DeSeq/expr_contrasts_genotypes.csv") %>% 
  filter(GeneID %in%put_regulators$regulatoryGene)
  

# 2) Extract time from the 'Contrast' column (e.g. "24", "48", "72", "96")
df <- genotype_mock %>%
  mutate(Time = stringr::str_extract(Contrast, "\\d+$")) %>%  # gets the digits at the end
  mutate(Comparison = stringr::str_replace(Contrast, "_\\d+$", ""))  # e.g. "LA1282_vs_LA1809"

# Now you have columns: GeneID, log2FoldChange, Time, Comparison, etc.


df_wide <- df %>%
  pivot_wider(
    id_cols = GeneID,
    names_from = c(Comparison, Time),  # Make a combined column name
    values_from = log2FoldChange,
    names_sep = "_"
  )
# Turn the wide dataframe into a matrix
mat <- df_wide %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

# Basic heatmap
pheatmap(
  mat,
  scale = "row",            # Normalizes each gene's expression across columns
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = T,    # Turn off rownames for large gene lists
  show_colnames = TRUE,
  main = "Heatmap of log2FoldChanges"
)
