rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Get IDs of the focal DEGs
focal_degs <-  read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/DeSeq/1282_48_ss_DEGs_ID.csv",
                        row.names = 1)
# Load Expression Data of each Mock (VST and pre-filtered)

data_1282 <- read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/vst_counts_filtered_1282.csv") %>% 
  #filter(grepl("_pdb_mock_24hpi", X)) %>% 
  column_to_rownames(var = "X") %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% focal_degs$GeneID) %>% 
  rownames_to_column(var = "GeneID")

data_1809 <- read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/vst_counts_filtered_1809.csv") %>% 
  #filter(grepl("_pdb_mock_24hpi", X)) %>% 
  column_to_rownames(var = "X") %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% focal_degs$GeneID)%>% 
  rownames_to_column(var = "GeneID")

data_1941 <- read.csv("C:/Users/suaph281/Desktop/nesh_local/spen_lag_phase_rnaseq/vst_counts_filtered_1941.csv") %>% 
  #filter(grepl("_pdb_mock_24hpi", X)) %>% 
  column_to_rownames(var = "X") %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% focal_degs$GeneID)%>% 
  rownames_to_column(var = "GeneID")


data_1282_mako <- data_1282 %>% 
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>% 
  mutate(
    genotype = "LA1282",
    # Extract the last two digits before "hpi" for timepoint
    timepoint = sub(".*_(\\d{2})hpi$", "\\1", Sample),
    # Identify and label the inoculum
    inoculum = case_when(
      grepl("pdb_mock", Sample) ~ "pdb_mock",
      grepl("ss", Sample)       ~ "ss",
      TRUE                      ~ NA_character_
    )
  ) %>% 
  group_by(GeneID,genotype,timepoint, inoculum) %>% 
  summarise(expression=mean(Expression))


data_1809_mako <- data_1809 %>% 
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>% 
  mutate(
    genotype = "LA1809",
    # Extract the last two digits before "hpi" for timepoint
    timepoint = sub(".*_(\\d{2})hpi$", "\\1", Sample),
    # Identify and label the inoculum
    inoculum = case_when(
      grepl("pdb_mock", Sample) ~ "pdb_mock",
      grepl("ss", Sample)       ~ "ss",
      TRUE                      ~ NA_character_
    )
  ) %>% 
  group_by(GeneID,genotype,timepoint, inoculum) %>% 
  summarise(expression=mean(Expression))

data_1941_mako <- data_1941 %>% 
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>% 
  mutate(
    genotype = "LA1941",
    # Extract the last two digits before "hpi" for timepoint
    timepoint = sub(".*_(\\d{2})hpi$", "\\1", Sample),
    # Identify and label the inoculum
    inoculum = case_when(
      grepl("pdb_mock", Sample) ~ "pdb_mock",
      grepl("ss", Sample)       ~ "ss",
      TRUE                      ~ NA_character_
    )
  ) %>% 
  group_by(GeneID, genotype, timepoint, inoculum) %>% 
  summarise(expression=mean(Expression))


merged <- data_1282_mako %>% 
  bind_rows(data_1809_mako) %>%
  bind_rows(data_1941_mako) 

# Plot 
(merged %>% 
    ggplot(aes(x=as.numeric(timepoint),
                       y=expression,
           group = interaction(GeneID, inoculum, genotype),
           col=inoculum)
           )+
    geom_point()+
    geom_line()+
    facet_grid(.~genotype)

)
