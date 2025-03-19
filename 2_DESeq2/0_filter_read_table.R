# get filtered reads for DeSeq:

rm(list=ls())
pacman::p_load(tidyverse)
options(scipen = 999, big.mark = ",")

cts <- read.delim("salmon.merged.gene_counts.tsv", sep = "\t")


# Filter genes first
proteome_ids <- read.table("spen_curated_proteome_OG_pannzer_dedub_ids.txt")%>%
  mutate(gene=gsub(">", "", V1)) %>%
  mutate(gene=gsub("GeneExt~", "", gene))%>% 
  mutate(gene=gsub("mRNA:", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
  select(gene) %>% 
  unique()

# 1282 
matrix_1282 <- cts %>% 
  filter(gene_id %in%proteome_ids$gene) %>% 
  select(!gene_name) %>% 
  column_to_rownames("gene_id")%>% 
  as.data.frame() %>% 
  select(contains("1282"))

colnames(matrix_1282) <- gsub("^X", "", colnames(matrix_1282))

write.csv(matrix_1282, "LA1282_filtered_counts.csv")

# 1809
matrix_1809 <- cts %>% 
  filter(gene_id %in%proteome_ids$gene) %>% 
  select(!gene_name) %>% 
  column_to_rownames("gene_id")%>% 
  as.data.frame() %>% 
  select(contains("1809"))

colnames(matrix_1809) <- gsub("^X", "", colnames(matrix_1809))

write.csv(matrix_1809, "LA1809_filtered_counts.csv")


#1941
matrix_1941 <- cts %>% 
  filter(gene_id %in%proteome_ids$gene) %>% 
  select(!gene_name) %>% 
  column_to_rownames("gene_id")%>% 
  as.data.frame() %>% 
  select(contains("1941"))

colnames(matrix_1941) <- gsub("^X", "", colnames(matrix_1941))

write.csv(matrix_1941, "LA1941_filtered_counts.csv")


# Filter TF - file

tfs <- read.table("data/spen_TFs.txt", header=F) %>% 
  mutate(gene=gsub(">", "", V1)) %>%
  mutate(gene=gsub("GeneExt~", "", gene))%>% 
  mutate(gene=gsub("mRNA:", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
  select(gene, V2) %>% 
  unique() %>% 
  write.csv(., "data/spen_TFs.txt", row.names = F)
