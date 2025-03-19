# Visualization of lag phase duration
# Based on this data:DOI: 10.34133/plantphenomics.0214
# Severin Einspanier

rm(list=ls())
library(tidyverse)

# Add the phenotype 
phenotypes <- read.csv("AllSlopes_Ann_wthMock_noContaminants_threeReps_filtered.csv") %>% 
  filter(genotype %in% c("LA1282", "LA1809", "LA1941")) 
# Calculate the mean lag for each genotype

mean_lag <- aggregate(lag ~ genotype, data = phenotypes, FUN = mean)

# Order the genotypes based on the mean lag
phenotypes$genotype <- factor(phenotypes$genotype, levels = mean_lag$genotype[order(mean_lag$lag)])
svg("fig_1/lag_boxplots.svg", 
    width = 6, height = 4)
# Create the boxplot with jitter
(p1 <- phenotypes %>% 
    filter(comm=="ok") %>% 
    mutate(lag_days = lag/(60*24)) %>%
    ggplot(aes(x = genotype, y = lag_days, fill=genotype)) +
      geom_boxplot() +
    #geom_jitter(width = 0.2, alpha=.2) +
    labs(x = "Genotype", y = "Lag phase [days]", 
         fill="Genotype") +
    theme_bw()+
    scale_fill_manual(values = c("LA1282" = "#007F94", "LA1809" = "#EED78D", "LA1941" = "#C22B26")) +
  scale_y_continuous(limits=c(0,8.4), breaks = c(2,4,6))+
  theme(axis.text = element_text(size=11, color="Black"),
        axis.text.x = element_text(vjust=1, 
                                    size=11, color="black"),
        axis.title = element_text(size=11, color="Black"),
        legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=11), 
        legend.spacing.x = unit(.5, 'cm'),
        panel.grid.major.x =  element_blank()
  )+
    ggpubr::stat_compare_means(comparisons = list(c("LA1809", "LA1282"),
                                                  c("LA1809", "LA1941"),
                                                  c("LA1282", "LA1941")),
                               method = "wilcox.test", label = "p.signif")
)
