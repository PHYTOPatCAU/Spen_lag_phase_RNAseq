# Evaluate PE-data 
# Print figures
# Severin Einspanier
# 2025_01_25

rm(list=ls())
pacman::p_load(tidyverse)

data <- read.csv2("Data_Analysis_Report_20250125_1746.csv", 
                 header=T)

# summary of object size

hist(data$Obj.Size, main="Object size", 
     xlab="Object size [px]", col="lightblue")
# first, we're only interested in the fv/fm categories

data_fvfm <- data %>% 
  select(File, Date, Time, Obj.No, Obj.Size, Fv.fm.V....,
         Fv.fm.I...., Fv.fm.II...., Fv.fm.III...., Fv.fm.IV....,
         Mask.Border) %>% 
  mutate(
    timepoint=as.numeric(sub(".*P(\\d+)N.*", "\\1", .$File)),  # Extracts value after P,
    plate = as.numeric(sub(".*N(\\d+).*", "\\1", .$File)),
  # forgot to reset plate. 
  # for timepoint=6, do plate-5
  plate = ifelse(timepoint == 6, plate-5, plate),
  #timepoint = ifelse(timepoint == 6 & plate > 5, 7, timepoint),
  rep = ifelse(plate < 3, 1, 2),
  treatment = ifelse(plate == 1 | plate == 3, "control", "infected")
  ) 
colnames(data_fvfm) <- c("File", "Date", "Time", "Obj.No", "Obj.Size", "Fv.fm.V",
                         "Fv.fm.I", "Fv.fm.II", "Fv.fm.III", "Fv.fm.IV", "Mask.Border",
                         "timepoint", "plate", "rep", "treatment")
  
summary(as.factor(data_fvfm$timepoint))

# plot fvfm test wise 

(data_fvfm %>%
    ggplot(aes(x=treatment, fill=treatment, y=Fv.fm.IV)) +
    geom_boxplot()+
    facet_grid(~timepoint) 
)


hist(data_fvfm$Obj.Size, main="Object size", 
     xlab="Object size [px]", col="lightblue")

# plot number of objects per timepoint & plate
(data_fvfm %>%
    ggplot(aes(x=plate, fill=treatment)) +
    geom_bar()+
    facet_grid(~timepoint) 
)

# only one outlier in plate 2 of date 1
hist(data_fvfm$Mask.Border)


data_fvfm_genotypes <- data_fvfm %>% 
  filter(Mask.Border> 15 )

(data_fvfm_genotypes %>%
    ggplot(aes(x=plate, fill=treatment)) +
    geom_bar()+
    facet_grid(~timepoint) 
)

# classify the genotypes with object-number 
# plate 1: 1-18 1282, 19-36 1809, 37-54 1941
# plate 2: 1-18 1282, 19- 38, 1941: 39-58  waste 40 
# plate 3: 1-18 1282, 19-36 1809, 37-54 1941
# plate 4: 1-18 1282, 19-36 1809, 37-54 1941

data_fvfm_genotypes <- data_fvfm %>% 
  filter(Mask.Border> 15 ) %>% 
  mutate(genotype = ifelse(Obj.No <= 18, "LA1282", 
                           ifelse(Obj.No > 18 & Obj.No <= 36, "LA1809", 
                                  "LA1941")
  ))

(data_fvfm_genotypes %>%
    ggplot(aes(x=plate, fill=genotype)) +
    geom_bar()+
    facet_grid(~timepoint) 
)                                                                                               
# remove technical artifacts
# plate 4: time point 9 remove obj. 15
# plate 4: time point 8: remove 33 & 11
# plate 4: time point 7: remove 11

data_fvfm_genotypes <- data_fvfm_genotypes %>% 
  filter(!(plate == 4 & timepoint == 9 & Obj.No == 15)) %>%
  filter(!(plate == 4 & timepoint == 8 & (Obj.No == 33 | Obj.No == 11))) %>%
  filter(!(plate == 4 & timepoint == 7 & Obj.No == 11))

# vis stacked bar plot

summary_data <- data_fvfm_genotypes %>%
  filter(treatment=="infected") %>% 
  pivot_longer(cols = starts_with("Fv.fm"), names_to = "category", values_to = "FvFm") %>%
  group_by(genotype, timepoint, category, treatment) %>%
  summarize(Count = mean(FvFm, na.rm = TRUE)) %>% 
  mutate(timepoint=(timepoint-1)*12)

# Stacked bar plot


(plot_cat <- ggplot(summary_data, aes(x = as.factor(timepoint), 
                         y = Count, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~genotype) +
  labs(title = "Proportions of Fv/Fm Categories Over Time",
       y = "Proportion of the leaf [%]",
       x = "Timepoint [hpi]") +
  scale_fill_manual(
    values = c("#000000", "#3b1c24", "#74303c", "#b14651", "#f06060"),
    labels = c("I   [0-0.160 fv/fm]", "II  [0.160-0.320 fv/fm]", 
               "III [0.320-0.480 fv/fm]", "IV [0.480-0.640 fv/fm]", 
               "V  [0.640-1 fv/fm]"),   
    name = "Category"  # Optional: Rename legend title
  ) +
  theme_bw()+
  theme(axis.text = element_text(size=11, color="Black"),
        axis.text.x = element_text(vjust=1, 
                                   size=11, color="black"),
        axis.title = element_text(size=13, color="Black"),
        legend.title = element_text(size=13, 
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size=11), 
        legend.spacing.x = unit(.5, 'cm'),
        panel.grid.major.x =  element_blank(),
        strip.text = element_text(size=13, face="bold"),
        strip.background = element_rect(fill="white")
  )
)

ggsave(plot_cat, 
       filename = "fig_1/fvfm.png", 
       width = 8, height =3, dpi = 900)
