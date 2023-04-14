
## check for normal distribution
require(ggplot2)
library(rcompanion)


plotNormalHistogram(list_meta$meta_all$NO3_umol.l,  main = "NO3")
plotNormalHistogram(list_meta$meta_all$PO4_umol.l,  main = "PO4")
plotNormalHistogram(list_meta$meta_all$Si_umol.l,  main = "Si")
plotNormalHistogram(list_meta$meta_all$temperature...C.,  main = "temperature")
plotNormalHistogram(list_meta$meta_all$salinity..psu.,  main = "salinity")
plotNormalHistogram(list_meta$meta_all$O2umol.l,  main = "oxygen")
plotNormalHistogram(list_meta$meta_all$Fluorometer,  main = "fluorescence")


## boxplots
box_plot_H <- ggplot(data = list_meta$meta_all, aes(x= Region, y=NO3_umol.l)) + 
  geom_point() + 
  geom_boxplot()

print(box_plot_H)


box_plot_G <- ggplot(data = list_meta$meta_all, 
                     aes(x= Region, y=PO4_umol.l)) +
  geom_point() +
  geom_boxplot()

print(box_plot_G)



box_plot_J <- ggplot(data = list_meta$meta_all, 
                     aes(x= Region, y=Si_umol.l)) +
  geom_point()+
  geom_boxplot()

print(box_plot_J)



#meta_all$Glacial.influence <- as.factor(meta_all$Glacial.influence)
print(kruskal.test(temperature...C. + salinity..psu. + O2umol.l + Fluorometer + PO4_umol.l + NO3_umol.l + Si_umol.l ~ Bioclimatic_subzone, list_meta$meta_all))


## T-S plot

h <- ggplot(list_meta$meta_all, aes(salinity..psu., temperature...C., color=Fjord))+
         geom_point()+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
print(h)

h <- ggplot(list_meta$meta_all, aes(salinity..psu., temperature...C., color=Region))+
  geom_point()+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(h)


 #clean up!

rm(box_plot_J, box_plot_H, box_plot_G, box_plot_D, box_plot_B, box_plot_AB, h, box_plot_A)
rm(meta_all_glacier, meta_all_Nglacier)
       