# Plot Cluster - Membership


cluster_asignment <- read.table("CONSENSUS_48degs/2025_02_15_48erDEGs_optimal_clustering.txt",
                                header=T)%>%
  mutate(genotype = sub("_.*", "", gene)
  ) 


# idea: number of genes per genotype and cluster
# cluster as facet, share (of total number in this cluster)

cluster_sum<- cluster_asignment %>%
  mutate(cluster = factor(cluster, levels = 1:5),
         genotype = as.factor(genotype)) %>%
  group_by(cluster, genotype) %>%
  summarize(share = n() / 71 * 100, .groups = 'drop') %>%
  complete(cluster, genotype, fill = list(share = 0))

(p1 <- cluster_sum %>% 
    ggplot(aes(x=genotype, y=share, fill=genotype)) +
    geom_col()+theme_bw()+
    scale_fill_manual(values = c("LA1282" = "#007F94", "LA1809" = "#EED78D", "LA1941" = "#C22B26")) +
    scale_y_continuous(limits=c(0,100), breaks = c(25,50,75))+
    theme(axis.text = element_text(size=11, color="Black"),
          axis.text.x = element_text(vjust=1, 
                                     size=11, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          strip.background = element_rect(fill="white"),
          panel.grid.major.x =  element_blank(),
          strip.
    )+
    facet_grid(.~cluster))
