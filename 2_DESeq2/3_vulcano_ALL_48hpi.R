# Draw annotated vulcano plot for 48 hpi 
# LA1809, LA1282, LA1941
# Severin Einspanier

rm(list=ls())

library(tidyverse)
library(ggrepel)


# get information about those genes:

data1282 <- read.csv("expr_dat_1282_inf_mock.csv",
                     row.name=1) %>% 
  filter(timepoint=="48hpi" & abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(GeneID) %>% 
  left_join(read.csv("genid2goid_spen.csv",
                     row.name=1), by=c("GeneID"="GeneID") ) 

write.csv(data1282, "1282_48_ss_DEGs_ID.csv")

Gene_Info <- read.delim("ITAG_LA1282_DEGs.txt",
                        header=F) %>% 
  dplyr::rename("OG_ID"=V1, "firstline"=V2, "ITAG"=V3, "itag_function"=V4, "ath_function"=V5) %>% 
  left_join(data1282, by=c("OG_ID"="OG")) %>% 
  mutate(itag_function=sub("^[^ ]+ ", "", itag_function)) %>% 
  mutate(itag_function=sub("\\(AH.*", "",itag_function)) %>% 
  mutate(ath_function=gsub(".*Symbols:([^|]+)\\|.*", "\\1", ath_function)) %>%
  unique() %>% 
  mutate(merged_function=ifelse(itag_function == "", ath_function, itag_function)) %>% 
  select(GeneID, OG_ID, firstline, ITAG, merged_function)



# select interesting genes #Sopen03g029450, Sopen11g001970,Sopen04g026840

interesting <- c("Sopen01g040940", "Sopen01g040950","Sopen05g011470",
                 "Sopen09g028060","Sopen09g036010",  "Sopen07g001220",
                 "Sopen08g028950", "Sopen09g034580", "Sopen11g027650",
                 "Sopen01g048970", "Sopen01g048960", "Sopen01g026170", "Sopen01g026150",
                 "Sopen08g009740", "Sopen02g031830", "Sopen07g024910", "Sopen07g024880", 
                 "Sopen04g028240")

annotation <- Gene_Info %>% 
  filter(GeneID %in% interesting) %>% 
  select(GeneID, merged_function) #%>% 

# Vulcano

data1282_2 <- read.csv("expr_dat_1282_inf_mock.csv",
                       row.name=1) %>% 
  filter(timepoint=="48hpi") %>% 
  left_join(annotation, by=c("GeneID"="GeneID")) 

(volcano_plot1282 <- data1282_2 %>%
    # Remove NA p-values if any
    filter(!is.na(pvalue)) %>%
    mutate(
      Expression = case_when(
        log2FoldChange >= 1 & pvalue <= 0.05 ~ "Up-regulated",
        log2FoldChange <= -1 & pvalue <= 0.05 ~ "Down-regulated",
        TRUE ~ "Unchanged"
      )
    ) %>%
    mutate(genotype="LA1282") %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
    scale_fill_manual(values = c("Up-regulated" = "#2085A2", 
                                 "Unchanged" = "gray90", 
                                 "Down-regulated" = "#F06060")) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    ylim(0,10)+
    xlim(-5,5)+
    geom_hline(yintercept = -log10(0.05), size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = c(-1, 1), size = 1, linetype = "dashed", color = "grey70") +
    geom_point(aes(fill = Expression), size = 4, col="black", 
               shape=21, stroke=1) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12, color = "Black"),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, color = "black"),
      axis.title = element_text(size = 12, color = "Black"),
      legend.title = element_text(size = 12, margin = margin(t = 0, r = 10, b = 5, l = 0), face = "bold"),
      legend.text = element_text(size = 12), 
      legend.spacing.x = unit(0.5, 'cm'), 
      strip.text = element_text(size=13, face="bold", color="black"),
      strip.background = element_rect(fill="#bae0e7ff"),
      panel.grid.minor.y = element_blank()
    ) +
    geom_label_repel(aes(label=merged_function),max.overlaps = Inf,
                     seed = 51983, force = 2, force_pull = -0.1,
                     size=4, segment.size=.51, max.time = 2,
                     box.padding = 1,
                     label.padding = .2)+
    guides(fill="none")+
    facet_grid(.~genotype)
)



# LA1809

data1809 <- read.csv("expr_dat_1809_inf_mock.csv",
                     row.name=1) %>% 
  filter(timepoint=="48hpi" & abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(GeneID) %>% 
  left_join(read.csv("genid2goid_spen.csv",
                     row.name=1), by=c("GeneID"="GeneID") ) 


data1809_2 <- read.csv("expr_dat_1809_inf_mock.csv",
                       row.name=1) %>% 
  filter(timepoint=="48hpi") %>% 
  left_join(annotation, by=c("GeneID"="GeneID")) 

(volcano_plot1809 <- data1809_2 %>%
    # Remove NA p-values if any
    filter(!is.na(pvalue)) %>%
    mutate(
      Expression = case_when(
        log2FoldChange >= 1 & pvalue <= 0.05 ~ "Up-regulated",
        log2FoldChange <= -1 & pvalue <= 0.05 ~ "Down-regulated",
        TRUE ~ "Unchanged"
      )
    ) %>%
    mutate(genotype="LA1809") %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
    scale_fill_manual(values = c("Up-regulated" = "#2085A2", 
                                 "Unchanged" = "gray90", 
                                 "Down-regulated" = "#F06060")) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    ylim(0,10)+
    xlim(-5,5)+
    geom_hline(yintercept = -log10(0.05), size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = c(-1, 1), size = 1, linetype = "dashed", color = "grey70") +
    geom_point(aes(fill = Expression), size = 4, col="black", 
               shape=21, stroke=1) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12, color = "Black"),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, color = "black"),
      axis.title = element_text(size = 12, color = "Black"),
      legend.title = element_text(size = 12, margin = margin(t = 0, r = 10, b = 5, l = 0), face = "bold"),
      legend.text = element_text(size = 12), 
      strip.text = element_text(size=13, face="bold", color="black"),
      strip.background = element_rect(fill="#efe9d1ff"), 
      legend.spacing.x = unit(0.5, 'cm'),
      panel.grid.minor.y = element_blank()
    )+ 
    guides(fill="none")+
    facet_grid(.~genotype)
)

# LA1941

data1941 <- read.csv("expr_dat_1941_inf_mock.csv",
                     row.name=1) %>% 
  filter(timepoint=="48hpi" & abs(log2FoldChange) > 1 & padj < 0.05) %>% 
  select(GeneID) %>% 
  left_join(read.csv("genid2goid_spen.csv",
                     row.name=1), by=c("GeneID"="GeneID") )


Gene_Info <- read.delim("ITAG_LA1941_DEGs.txt",
                        header=F) %>% 
  dplyr::rename("OG_ID"=V1, "firstline"=V2, "ITAG"=V3, "itag_function"=V4, "ath_function"=V5) %>% 
  left_join(data1941, by=c("OG_ID"="OG")) %>% 
  mutate(itag_function=sub("^[^ ]+ ", "", itag_function)) %>% 
  mutate(itag_function=sub("\\(AH.*", "",itag_function)) %>% 
  mutate(ath_function=gsub(".*Symbols:([^|]+)\\|.*", "\\1", ath_function)) %>%
  unique() %>% 
  mutate(merged_function=ifelse(itag_function == "", ath_function, itag_function)) %>% 
  select(GeneID, OG_ID, firstline, ITAG, merged_function)



# select interesting genes 

interesting <- c("Sopen09g032860", "Sopen04g035160","Sopen04g001120","Sopen10g001910",
                 "Sopen06g024710","Sopen03g040350", "Sopen04g024830")

annotation <- Gene_Info %>% 
  filter(GeneID %in% interesting) %>% 
  select(GeneID, merged_function) #%>% 

# Vulcano


data1941_2 <- read.csv("expr_dat_1941_inf_mock.csv",
                       row.name=1) %>% 
  filter(timepoint=="48hpi") %>% 
  left_join(annotation, by=c("GeneID"="GeneID")) 

(volcano_plot1941 <- data1941_2 %>%
    # Remove NA p-values if any
    filter(!is.na(pvalue)) %>%
    mutate(
      Expression = case_when(
        log2FoldChange >= 1 & pvalue <= 0.05 ~ "Up-regulated",
        log2FoldChange <= -1 & pvalue <= 0.05 ~ "Down-regulated",
        TRUE ~ "Unchanged"
      )
    ) %>%
    mutate(genotype="LA1941") %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
    scale_fill_manual(values = c("Up-regulated" = "#2085A2", 
                                 "Unchanged" = "gray90", 
                                 "Down-regulated" = "#F06060")) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    ylim(0,10)+
    xlim(-5,5)+
    geom_hline(yintercept = -log10(0.05), size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = c(-1, 1), size = 1, linetype = "dashed", color = "grey70") +
    geom_point(aes(fill = Expression), size = 4, col="black", 
               shape=21, stroke=1) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12, color = "Black"),
      axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 1, color = "black"),
      axis.title = element_text(size = 12, color = "Black"),
      legend.title = element_text(size = 12, margin = margin(t = 0, r = 10, b = 5, l = 0), face = "bold"),
      legend.text = element_text(size = 12), 
      strip.text = element_text(size=13, face="bold", color="black"),
      strip.background = element_rect(fill="#d9a09fff"),
      legend.spacing.x = unit(0.5, 'cm'),
      panel.grid.minor.y = element_blank()
    ) +
    geom_label_repel(aes(label=merged_function),max.overlaps = Inf,
                     seed = 3464, force = 1, force_pull = -.1,
                     size=5, segment.size=.51,
                     box.padding = 2)+
    facet_grid(.~genotype)
)

png(paste0(Sys.Date(), "_combined_vulcanos.png"),
    width=15, height=27, unit="cm", res=900)
ggpubr::ggarrange(volcano_plot1282, volcano_plot1809, volcano_plot1941,  nrow = 3,
                  common.legend = T, 
                  labels = c("A", "B", "C"),
                  heights = c(3,2,3),
                  font.label = list(size = 20))
dev.off()

