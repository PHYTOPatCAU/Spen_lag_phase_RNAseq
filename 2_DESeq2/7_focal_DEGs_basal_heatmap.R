# Visualize the 71 DEGs of interest in mock conditions between the genotypes
# counts 
# three genotypes 
# Severin Einspanier
# 71 DEGs
rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Get IDs of the focal DEGs
focal_degs <-  read.csv("1282_48_ss_DEGs_ID.csv",
                                    row.names = 1)
# Load Expression Data of each Mock (VST and pre-filtered)

data_1282 <- read.csv("vst_counts_filtered_1282.csv") %>% 
  filter(grepl("_pdb_mock_24hpi", X)) %>% 
  column_to_rownames(var = "X") %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% focal_degs$GeneID) %>% 
  rownames_to_column(var = "GeneID")

data_1809 <- read.csv("vst_counts_filtered_1809.csv") %>% 
  filter(grepl("_pdb_mock_24hpi", X)) %>% 
  column_to_rownames(var = "X") %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% focal_degs$GeneID)%>% 
  rownames_to_column(var = "GeneID")

data_1941 <- read.csv("vst_counts_filtered_1941.csv") %>% 
  filter(grepl("_pdb_mock_24hpi", X)) %>% 
  column_to_rownames(var = "X") %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rownames(.) %in% focal_degs$GeneID)%>% 
  rownames_to_column(var = "GeneID")


# ---- 1. Calculate Mean Expression Per Genotype ----
mean_1282 <- data_1282%>% filter(GeneID != "Sopen11g009900") %>% column_to_rownames("GeneID") %>% rowMeans() 
mean_1809 <- data_1809%>% filter(GeneID != "Sopen11g009900") %>% column_to_rownames("GeneID") %>% rowMeans()
mean_1941 <- data_1941 %>% column_to_rownames("GeneID") %>% rowMeans()

# Combine into a single data frame
heatmap_matrix <- data.frame(
  GeneID = names(mean_1282),
  `LA1282` = mean_1282,
  `LA1809` = mean_1809,
  `LA1941` = mean_1941
  ) %>%
  select(-GeneID) %>%
  as.matrix()


###
# ---- Normalize the Whole Dataframe (Global Min-Max Normalization) ----
normalize_minmax_global <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  (x - min_val) / (max_val - min_val)
}

normalized_matrix_global <- normalize_minmax_global(heatmap_matrix)

data1282 <- read.csv("expr_dat_1282_inf_mock.csv",
                     row.name=1) %>% 
  filter(timepoint=="48hpi" & abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(GeneID) %>% 
  left_join(read.csv("genid2goid_spen.csv",
                     row.name=1), by=c("GeneID"="GeneID") ) 

Gene_Info <- read.delim("ITAG_LA1282_DEGs.txt", header = FALSE) %>% 
  dplyr::rename("OG_ID" = V1, "firstline" = V2, "ITAG" = V3, "itag_function" = V4, "ath_function" = V5) %>% 
  left_join(data1282, by = c("OG_ID" = "OG")) %>% 
  mutate(itag_function = sub("^[^ ]+ ", "", itag_function)) %>% 
  mutate(itag_function = sub("\\(AH.*", "", itag_function)) %>% 
  mutate(ath_function = gsub(".*Symbols:([^|]+)\\|.*", "\\1", ath_function)) %>%
  unique() %>% 
  mutate(merged_function = ifelse(itag_function == "", ath_function, itag_function)) %>% 
  select(GeneID, OG_ID, firstline, ITAG, merged_function) %>%
  mutate(funct_unique = ifelse(merged_function == "", GeneID, 
                               ifelse(duplicated(merged_function) | duplicated(merged_function, fromLast = TRUE),
                                      paste0(merged_function, " (", GeneID, ")"), merged_function))) %>% 
  mutate(funct_unique = sub(" ", "", funct_unique)) %>% 
  select(-merged_function, -ITAG, -firstline, -OG_ID)


normalized_matrix_global_df <- as.data.frame(heatmap_matrix) %>% 
  rownames_to_column(var = "GeneID") %>% 
  left_join(Gene_Info, by="GeneID") %>% 
  mutate(funct_unique = ifelse(is.na(funct_unique), GeneID, funct_unique)) %>%
  select(-GeneID) %>% 
  column_to_rownames(var = "funct_unique") %>% 
  as.matrix()





heatmap <- Heatmap(normalized_matrix_global_df, 
                   name = "Variance-Stabilized Counts",
                   cluster_rows = TRUE, 
                   cluster_columns = TRUE, 
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   row_dend_width = unit(.5, "cm"),
                   column_dend_height = unit(.5, "cm"),
                   row_names_gp = grid::gpar(fontsize = 8),
                   width             = ncol(normalized_matrix_global_df) * unit(5, "mm"),
                   height            = nrow(normalized_matrix_global_df) * unit(3, "mm"),
                    col = colorRamp2(
                      breaks = c(1, 2,3,4, 
                                 5,6,7, 8, 
                                 9,10,11),
                      colors = c("#000000", "#29241a", "#3b331e", "#504321",
                                 "#655324","#7c6327", "#94742a", "#ad842c",
                                 "#c7962e", "#e3a731", "#ffb833")
                    ),
                   #row_dend_width    = unit(1, "cm"),
                   rect_gp           = gpar(col = "black", lwd = 2),
                   heatmap_legend_param = list(
                     # Choose where you want the break lines, e.g. at 1, 4, 8, 12
                     at     = c(1, 4, 8, 12),  
                     # Corresponding labels
                     labels = c("1", "4", "8", "12")
                   ))

svg(paste0( Sys.Date(), 
           "_basal_regulation_focal_genes.svg"),
    width=6, height=10)
# Draw the heatmap
draw(heatmap)
dev.off()


png(paste0(Sys.Date(), 
           "_basal_regulation_focal_genes.png"),
    width=6, height=10, units="in", res=900)
# Draw the heatmap
draw(heatmap)
dev.off()
