# Get basal expression DEGs
# PTI genes
# Severin Einspanier

rm(list=ls())
library(tidyverse)
# CERK1 homologs: 
cerk1 <- c("Sopen02g025710", "Sopen02g025720", "Sopen07g025160")

pepr1 <- "Sopen03g041600" # ath has 2 but tomato only 1?

bak <- c("Sopen01g047490", "Sopen04g028770", "Sopen10g017310")

lyk3 <- "Sopen03g039880"
lyk4 <- "Sopen02g034600"
lym1 <- "Sopen11g007610"
sobir1 <- c("Sopen03g030890", "Sopen06g028100")


pti_genes <- data.frame(
  gene = c("cerk1.1", "cerk1.2", "cerk1.3", "pepr1", "bak.1","bak.2", "bak.3",
           "lyk3", "lyk4", "lym1", "sobir1.1", "sobir1.2"),
  geneID=c("Sopen02g025710", "Sopen02g025720", "Sopen07g025160", "Sopen03g041600",
           "Sopen01g047490", "Sopen04g028770", "Sopen10g017310", "Sopen03g039880",
           "Sopen02g034600", "Sopen11g007610", "Sopen03g030890", "Sopen06g028100")
)


# How different are they regulated between the genotypes in mock?

mock_1282_1809_24hpi <- read.csv("expr_contrasts_genotypes_mock.csv") %>% 
  filter(Contrast=="LA1282_vs_LA1809_24")


pti_genes_degs <- pti_genes %>% 
  left_join(mock_1282_1809_24hpi, by=c("geneID"="GeneID")) %>% 
  mutate(DEG = ifelse(padj<0.05 & abs(log2FoldChange) > 1, "DEG", "Not DEG"))

write.csv(pti_genes_degs, "pti_genes_1282_1809.csv", row.names=FALSE)
