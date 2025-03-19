# Mapping statistics
# Severin Einspanier

pacman::p_load(tidyverse)

map_dat <- read.table("mapping_stats/star_alignment_plot.tsv",
                      sep="\t", header=TRUE) %>% 
  left_join(read.csv("sample_sheet.csv"), 
            by=c("Sample"="sample")) %>% 
  select(-fastq_1, -fastq_2, -strandedness, -ID) %>% 
  mutate(inoculum=ifelse(inoculum=="pdb_mock", "mock", "ss"))

# Define a position_dodge object with a specified width
pd <- position_dodge(width = 0.4)

svg("fig_2/mapping.svg", 
    width=7, height=5)
# Create the modified plot with position dodging
(p1 <- map_dat %>% 
    ggplot(aes(x = as.factor(timepoint), 
               y = Uniquely.mapped + Mapped.to.multiple.loci, 
               color = inoculum,
               fill = inoculum)) +

    # Add error bars for ±1 SD with dodge
    stat_summary(fun.data = mean_sdl, 
                 fun.args = list(mult = 1), 
                 linewidth=1,
                 col="black",
                 geom = "errorbar", 
                 width = 0.3,         # Width of error bars
                 position = pd) +
    # Add points for the mean with dodge
    stat_summary(fun = mean, 
                 geom = "point", 
                 size = 5,
                 stroke=2,
                 col="black",
                 shape = 23,         # Diamond shape
                 position = pd) +
    # Add individual data points with jitter and dodge
    geom_point(position=position_jitterdodge(jitter.width=.2, 
                                             jitter.height=0, 
                                             dodge.width = 0.4),
               alpha = 1,        # Transparency
               size = 4,
               stroke=1, # Point size
               shape = 21,         # Shape with fill
               color = "black") +
    # Labels and theme adjustments
    labs(x = "Time [h]", 
         y = "Mapping Rate [%]", 
         color = "Inoculum", 
         fill = "Inoculum") +
    
    theme_bw() +
    
    # Manual color and fill scale for 'inoculum'
    #scale_color_manual(values = c("mock" = "#007F94", "ss" = "#C22B26")) +
    scale_fill_manual(values = c("mock" = "#007F94", "ss" = "#C22B26")) +
    
    # Theme customizations for text, legend, and facets
    theme(
      axis.text = element_text(size = 11, color = "Black"),
      axis.text.x = element_text(vjust = 1, size = 11, color = "black"),
      axis.title = element_text(size = 13, color = "Black"),
      legend.title = element_text(size = 13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      legend.text = element_text(size = 11), 
      legend.spacing.x = unit(0.5, 'cm'),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = 13, face = "bold"),
      strip.background = element_rect(fill = "white")
    ) +
    ylim(50,100)+
    # Facet by 'accession'
    facet_wrap(. ~ accession)
)
dev.off()
ggsave()

# Sclero mapping 

sclero_data <- read.table("star_salmon/multiqc_report_data/star_alignment_plot.txt",
                          sep="\t", header=TRUE) %>% 
  left_join(read.csv("sample_sheet.csv"), 
            by=c("Sample"="sample")) %>% 
  select(-fastq_1, -fastq_2, -strandedness, -ID) %>% 
  mutate(inoculum=ifelse(inoculum=="pdb_mock", "mock", "ss")) %>% 
  mutate(all=Uniquely.mapped + Mapped.to.multiple.loci + Unmapped..too.many.mismatches+
                 Unmapped..too.short+ Unmapped..other,
         mappingrate = (Uniquely.mapped + Mapped.to.multiple.loci)/all*100)


# Define a position_dodge object with a specified width
pd <- position_dodge(width = 0.4)

svg("mapping_sclero.svg", 
    width=7, height=5)
# Create the modified plot with position dodging
(p1 <- sclero_data %>% 
    ggplot(aes(x = as.factor(timepoint), 
               y = mappingrate, 
               color = inoculum,
               fill = inoculum)) +
    
    # Add error bars for ±1 SD with dodge
    stat_summary(fun.data = mean_sdl, 
                 fun.args = list(mult = 1), 
                 linewidth=1,
                 col="black",
                 geom = "errorbar", 
                 width = 0.3,         # Width of error bars
                 position = pd) +
    # Add points for the mean with dodge
    stat_summary(fun = mean, 
                 geom = "point", 
                 size = 5,
                 stroke=2,
                 col="black",
                 shape = 23,         # Diamond shape
                 position = pd) +
    # Add individual data points with jitter and dodge
    geom_point(position=position_jitterdodge(jitter.width=.2, 
                                             jitter.height=0, 
                                             dodge.width = 0.4),
               alpha = 1,        # Transparency
               size = 4,
               stroke=1, # Point size
               shape = 21,         # Shape with fill
               color = "black") +
    # Labels and theme adjustments
    labs(x = "Time [h]", 
         y = "Mapping Rate [%]", 
         color = "Inoculum", 
         fill = "Inoculum") +
    
    theme_bw() +
    
    # Manual color and fill scale for 'inoculum'
    #scale_color_manual(values = c("mock" = "#007F94", "ss" = "#C22B26")) +
    scale_fill_manual(values = c("mock" = "#007F94", "ss" = "#C22B26")) +
    
    # Theme customizations for text, legend, and facets
    theme(
      axis.text = element_text(size = 11, color = "Black"),
      axis.text.x = element_text(vjust = 1, size = 11, color = "black"),
      axis.title = element_text(size = 13, color = "Black"),
      legend.title = element_text(size = 13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      legend.text = element_text(size = 11), 
      legend.spacing.x = unit(0.5, 'cm'),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = 13, face = "bold"),
      strip.background = element_rect(fill = "white")
    ) +
    #ylim(50,100)+
    # Facet by 'accession'
    facet_wrap(. ~ accession)
)
dev.off()


