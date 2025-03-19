# Lineplot of DEGs per Genotype and Time point
# Severin Einspanier
rm(list=ls())
options(stringsAsFactors = FALSE,
        OutDec = ",")


pacman::p_load(tidyverse, DESeq2)



data1282 <- read.csv("expr_dat_1282_inf_mock.csv",
                     row.names = 1) %>% 
  group_by(timepoint) %>% 
  summarise(
    UP = sum(log2FoldChange > 1 & padj < 0.05, na.rm = TRUE),  # Count UP-regulated genes
    DOWN = sum(log2FoldChange < -1 & padj < 0.05, na.rm = TRUE) # Count DOWN-regulated genes
  ) %>% 
  mutate(genotype="LA1282")

data1809 <- read.csv("expr_dat_1809_inf_mock.csv",
                     row.names = 1) %>% 
  group_by(timepoint) %>% 
  summarise(
    UP = sum(log2FoldChange > 1 & padj < 0.05, na.rm = TRUE),  # Count UP-regulated genes
    DOWN = sum(log2FoldChange < -1 & padj < 0.05, na.rm = TRUE) # Count DOWN-regulated genes
  )%>% 
  mutate(genotype="LA1809")
data1941 <- read.csv("expr_dat_1941_inf_mock.csv",
                     row.names = 1) %>% 
  group_by(timepoint) %>% 
  summarise(
    UP = sum(log2FoldChange > 1 & padj < 0.05, na.rm = TRUE),  # Count UP-regulated genes
    DOWN = sum(log2FoldChange < -1 & padj < 0.05, na.rm = TRUE) # Count DOWN-regulated genes
  )%>% 
  mutate(genotype="LA1941")

combined <- data1282 %>% 
  bind_rows(data1809) %>% 
  bind_rows(data1941) %>% 
  mutate(DOWN = -DOWN) %>% 
  pivot_longer(cols = c(UP, DOWN), names_to = "type", values_to = "count")

svg(paste0("fig_2/", Sys.Date(), "DEGs_time.svg"),
    width = 8, height = 5)
(# Plot the data
  p <- combined %>%
    mutate(type = as.factor(type)) %>%
    ggplot(aes(x=timepoint,y=count, alpha=type, color=genotype, fill=genotype))+
    # geom_histogram(stat="count",
    #                position="dodge")
    geom_line(aes(group = interaction(type, genotype)),
              size=2) +
    geom_point(size=5)+
    theme_bw()+
    scale_fill_manual(values = c("LA1282" = "#007F94", "LA1809" = "#EED78D", "LA1941" = "#C22B26")) +
    scale_color_manual(values = c("LA1282" = "#007F94", "LA1809" = "#EED78D", "LA1941" = "#C22B26")) +
    scale_alpha_manual(values = c("UP" = 1, "DOWN" = 0.75))+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(vjust=1, 
                                     size=11, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank()
    )+
    geom_hline(yintercept = 0, linetype="dashed")+
    scale_y_continuous(labels = function(x) abs(x)) +
    ggrepel::geom_label_repel(    data = combined %>% filter(count != 0),  # Exclude zero values
                                  aes(label = abs(count)),
                                  fill = "white",
                                  position = position_nudge(y = 0.5),
                                  force = 2,
                                  max.overlaps = Inf,
                                  size = 4
    ) )
dev.off()
  