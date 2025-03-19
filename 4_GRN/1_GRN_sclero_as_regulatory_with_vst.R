# Make GRN with effetors as Regulatory Genes 
library(DESeq2)
library(GENIE3)
library(segmented)
library(tidyverse)

# 1. VST discretly 

## 1.1 plants
setwd("C:/Users/suaph281/Desktop/GitLab/spen_lag_phase_rnaseq/")
setwd("/gxfs_home/cau/suaph281/spen_lag_phase_rnaseq")

cts <- read.delim("/gxfs_work/cau/suaph281/RNAseq/spen_lag_timeseries/data/2024_11_12_nfcore_out/star_salmon/salmon.merged.gene_counts.tsv", sep = "\t")
# Filter genes first
proteome_ids <- read.table("data/spen_curated_proteome_OG_pannzer_dedub_ids.txt")%>%
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
  select(contains("1282")) %>% 
  select(contains("ss"))

colnames(matrix_1282) <- gsub("^X", "", colnames(matrix_1282))

coldata <- read.csv("data/sample_sheet.csv") %>% 
  filter(accession %in% c("LA1282")&inoculum=="ss") %>% 
  select(sample, accession, timepoint, inoculum)
rownames(coldata)=coldata$sample

if (!all(rownames(coldata) == colnames(matrix_1282))) {
  # Reorder the columns of cts_matrix to match the order of rownames in coldata
  matrix_1282 <- matrix_1282[, rownames(coldata)]
}

all(rownames(coldata) == colnames(matrix_1282))

coldata$inoculum <- as.factor(coldata$inoculum)
coldata$accession <- as.factor(coldata$accession)
coldata$timepoint <- as.factor(coldata$timepoint)



# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(matrix_1282),
                              colData = coldata,
                              design = ~ 1)
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize 
dds <- dds[keep, ]
dds_norm <- vst(dds)

normalized_counts_1282 <- assay(dds_norm) %>% 
  as.data.frame() 
## 1.2 VST fungus 

cts_sc<- read.delim("/gxfs_work/cau/suaph281/RNAseq/spen_lag_timeseries/data/2025_02_10_scleromap/star_salmon/salmon.merged.gene_counts.tsv", 
                  sep="\t")

matrix_sc <- cts_sc %>% 
  select(!gene_name) %>% 
  mutate(gene_id=gsub("gene-", "", gene_id)) %>%
  mutate(gene_id = gsub("_", "", gene_id)) %>%
  column_to_rownames("gene_id")%>% 
  as.data.frame() %>% 
  select(contains("1282")) %>% 
  select(contains("ss"))

colnames(matrix_sc) <- gsub("^X", "", colnames(matrix_sc))


if (!all(rownames(coldata) == colnames(matrix_sc))) {
  # Reorder the columns of cts_matrix to match the order of rownames in coldata
  matrix_sc <- matrix_sc[, rownames(coldata)]
}

all(rownames(coldata) == colnames(matrix_sc))

# Create DESeqDataSet object
dds_S <- DESeqDataSetFromMatrix(countData = round(matrix_sc),
                              colData = coldata,
                              design = ~ 1)
# keep_s <- rowSums(counts(dds_S) >= 5) >= smallestGroupSize 
# dds_S <- dds_S[keep_s, ]
# I skip filtering as there are so few left...

dds_S_norm <- varianceStabilizingTransformation(dds_S)

normalized_counts_sclero <- assay(dds_S_norm) %>% 
  as.matrix()

dim(normalized_counts_1282)
dim(normalized_counts_sclero)

merged_counts <- as.data.frame(normalized_counts_1282) %>% 
  bind_rows(as.data.frame(normalized_counts_sclero)) %>% 
  as.matrix()


# filter the gene expression:
#colnames(merged_counts)[1] <- "sample"

datExpr <-merged_counts

# MERGE TFs and Effectors

effectors <- read.csv("data/effectors.csv") %>% 
  mutate(gene = tolower(Version.Two.ID)) %>% 
  dplyr::select(gene ) 
  
TF <- read.delim("data/spen_TFs.txt", 
                 sep=",", header=T) %>%
  mutate(gene=gsub(">", "", gene)) %>%
  mutate(gene=gsub("GeneExt~", "", gene))%>% 
  mutate(gene=gsub("mRNA_", "", gene))%>%  
  mutate(gene=gsub("t\\.peak", "g.peak", gene )) %>%
  mutate(gene=gsub("t\\.minus", "g.minus", gene)) %>% 
  mutate(gene=gsub("t\\.plus", "g.plus", gene)) %>% 
  mutate(gene=gsub("\\.[1-9].*|\\.p[1-9].*", "",gene))%>%
  dplyr::select(gene) %>%
  unique()

TF_ids <- TF %>%
  bind_rows(effectors) %>% 
  dplyr::select(gene) %>%
  unique()

head(TF_ids)
dim(TF_ids)
tail(TF_ids)

tail(datExpr)

expressed_TFs <- TF_ids %>% 
  filter(gene %in% rownames(datExpr))

dim(datExpr)

summary(datExpr)

weightMat <- GENIE3(datExpr, regulators=expressed_TFs$gene, nCores=30)

# Created "linked list" from a weighted adjacency_matrix

wam_linked_list <- getLinkList(weightMatrix = weightMat) %>%
  mutate_if(is.factor, as.character)


# calculate edge weight threshold

brkpnt_fun <- function(x, y, output_file = "breakpoint_plot.png") {
  if (length(x) > 10000) {
    set.seed(54321)
    x <- sample(x, size = 10000)
  }
  df <- data.frame(  # generate data frame
    vec = sort(x), num = seq_along(sort(x))
  )
  glm_ <- glm(vec ~ num, data = df)  # create linear model
  seg <- segmented(glm_, seg.Z =  ~ num, npsi = y)  # calculate breakpoint(s)
  brkpnt <- seg[["psi"]][, 2] |> round() # return breakpoint(s)
  cat(  # print breakpoint(s)
    "Breakpoints are \nx:\t",
    paste0(brkpnt[1:y], collapse = ", "),
    "\ny:\t",
    paste0(df$vec[brkpnt[1:y]], collapse = ", "),
    "\n"
  )
  
  # Generate plot
  p <- ggplot() +  
    geom_point(data = df, aes(x = num, y = vec)) +
    geom_vline(xintercept = as.numeric(brkpnt), color='red') +
    geom_hline(yintercept = as.numeric(df$vec[brkpnt]), color='forestgreen') +
    theme_classic()
  
  # Save plot to file
  ggsave(output_file, plot = p)
}

breakpont <- brkpnt_fun(wam_linked_list$weight,9,paste0("documentation/pics/", Sys.Date(), "_1282_SCLERO_grn_breakpoints.png"))

# picked 0.0043565137634722
# filter dataframe with eigencentrality threshold
grn_edges_all <- wam_linked_list[wam_linked_list$weight > 0.00435705522345596 ,]

genes <- append(wam_linked_list$regulatoryGene, wam_linked_list$targetGene) %>% as.character() %>% unique() %>% sort()

 genes_names <- as.data.frame(genes) #%>% 
#   left_join(GO_TERM, by=c("genes"="gene"))%>%
#   mutate(protein=ifelse(is.na(desc), "unknown", desc)) %>%
#   dplyr::select(!desc)
dim(genes_names)

#### define hubs 

grn_nodes_all_eigen <- igraph::graph_from_data_frame(grn_edges_all)  # generate igraph onject
grn_all_eigen <- igraph::eigen_centrality(grn_nodes_all_eigen)  # calculate eigenvector centrality

grn_all_eigen <- data.frame(egnvctr = unname(grn_all_eigen[["vector"]]),  # generate dataframe of genes
                            gene = names(grn_all_eigen[["vector"]]))     # and eigenvector centrality


breakpont <- brkpnt_fun((grn_all_eigen$egnvctr),5, paste0("documentation/pics/", Sys.Date(), "_1282_SCLERO_breakpoints_hubs.png"))
hubs <- grn_all_eigen[grn_all_eigen$egnvctr > 0.132547712716916 ,]

genes_names_wth_hubs <- genes_names %>% 
  left_join(hubs, by=c("genes"="gene")) %>%
  mutate(hub_grn=ifelse(is.na(egnvctr), "non-hub", "hub")) %>%
  dplyr::select(!egnvctr)  
#left_join(modules, by=c("genes"="GeneID")) %>%
#left_join(wgcna_hub, by=c("genes"="gene")) 
# also append WGRC hubs 


dim(genes_names_wth_hubs)
head(genes_names_wth_hubs)

write_delim(genes_names_wth_hubs, file = "/gxfs_work/cau/suaph281/RNAseq/spen_lag_timeseries/data/6_GRN/sclero/LA1282_infected_with_effectors_grn_network_nodes.txt", delim = "\t")

# output GRN network
write_delim(grn_edges_all, file = "/gxfs_work/cau/suaph281/RNAseq/spen_lag_timeseries/data/6_GRN/sclero/LA1282_infected_with_effectors_grn_network_edges.txt", delim = "\t")
