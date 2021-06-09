## R script for visualising PW distances between chosen tree structures
## Tamika Lunn, Griffith University
## Version 1 - December 2020

library(plyr) #for ddply 
library(dplyr) #for renaming collumns
library(ggplot2)

setwd("D:/PhD/Project/Manuscripts/First author/Planned_Roost structure #3 (density and disease)/Data/Processed/NN distances")
PWdist <- read.csv("PW distances.csv")

site.labs <- c("Sparse", "Intermediate", "Dense")
names(site.labs) <- c("01-DTOW", "02-DCLU", "03-DLIS")

PWdist %>%
  #filter(site.code == "DCLU") %>%
  ggplot(aes(x=PWdistance)) +
  geom_histogram(alpha=0.5, position="identity") +
  #scale_fill_manual(values = c("#8FD744FF","#29AF7FFF", "#39568CFF","#482677FF",   "#000000","#6A6C71FF","#969797","#D2D3D2")) +
  guides(fill = guide_legend(reverse = TRUE))+
  facet_grid(site.code~.,labeller = labeller(site.code = site.labs), scales = "free") +
  theme(axis.text.x = element_text(size=20), #angle = 70, hjust = 1
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(x="Distance between tree pairs", y = "Count of tree pairs")


ddply(PWdist, c("site.code", "subplot"), summarise,
      median = median(PWdistance), 
      mean = mean(PWdistance))

