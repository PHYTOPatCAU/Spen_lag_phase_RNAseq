# Check basal RLP and RLK expression
# basal, mock dEGs genotypes
# Severin Einspanier 

library(tidyverse)
# CRLPs
crps <- read.table("rlps/Solpeni_CRPs.txt") %>% 
  mutate(GeneID = sub("\\.1$", "", V1)) %>% 
  mutate(GeneID = sub("\\.2$", "", GeneID)) %>% 
  select(GeneID)

clks <- read.table("rlps/solpeni_kinase.txt") %>% 
  mutate(GeneID = sub("\\.1$", "", V1)) %>% 
  mutate(GeneID = sub("\\.2$", "", GeneID)) %>% 
  select(GeneID) %>% 
  unique()


# How different are they regulated between the genotypes in mock?

mock_1282_1809_24hpi <- read.csv("DeSeq/expr_contrasts_genotypes_mock.csv") %>% 
  filter(Contrast=="LA1282_vs_LA1809_24")


pti_genes_degs <- crps %>% 
  left_join(mock_1282_1809_24hpi, by=c("GeneID"="GeneID")) %>% 
  mutate(DEG = ifelse(padj<0.05 & abs(log2FoldChange) > 1, "DEG", "Not DEG")) 
write.csv(pti_genes_degs, "crlps_genes_1282_1809.csv", 
          row.names=FALSE)

kinases_degs    <- clks %>% 
  left_join(mock_1282_1809_24hpi, by=c("GeneID"="GeneID")) %>% 
  mutate(DEG = ifelse(padj<0.05 & abs(log2FoldChange) > 1, "DEG", "Not DEG"))
write.csv(kinases_degs, "kinases_genes_1282_1809.csv", 
          row.names=FALSE)
