# Extract major regulators for focal genes
# Severin Einspanier
rm(list=ls())
library(tidyverse)


genes <- read.table("2025_02_15_48erDEGs_optimal_clustering.txt",
                    header=T) %>% 
  #filter(cluster=="2") %>% 
  pull(gene)

hub_nodes <- read.table("2025_02_17_infected_with_effectors_grn_network_nodes_07.txt",
                        header=T, sep="\t") %>% 
  filter(hub_grn =="hub") 

all_edges <- read.table("2025_02_17_infected_with_effectors_grn_network_edges.txt",
                        header=T, sep="\t") %>% 
  filter(targetGene %in% genes & regulatoryGene %in% hub_nodes$gene | targetGene %in% genes & regulatoryGene %in% genes) 

hist(all_edges$weight, breaks=100)

# summary 
all_edges %>% 
  group_by(targetGene) %>% 
  summarise(n=n(), sum=sum(weight), mean=mean(weight))

# get how many regulatory genes 
all_edges_data <-all_edges %>% 
  group_by(regulatoryGene) %>% 
  summarise(n=n()) 

hist(all_edges$weight, breaks=100)

all_edges %>%  filter(regulatoryGene == "Sopen02g011350")
## This is it. 
all_edges %>%  filter(regulatoryGene == "Sopen02g011350") %>% 
  write.csv2("2025_02_17_shared_TFGRN.csv")


# get within network regulation 

all_edges_within <- read.table("2025_02_17_infected_with_effectors_grn_network_edges.txt",
                        header=T, sep="\t") %>% 
  filter(targetGene %in% genes & regulatoryGene %in% genes & regulatoryGene %in% hub_nodes$genes) 

edges_out <- all_edges %>%  filter(regulatoryGene == "Sopen02g011350")

merge_out <- edges_out %>% 
  bind_rows(all_edges_within) 

write.table(merge_out, "2025_02_18_grn_focal_genes.tsv",
            sep="\t", row.names = F)
