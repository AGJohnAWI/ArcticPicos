
#subset metadata with and without marine-terminating glaciers

meta_all_glacier <- list_meta$meta_all%>%dplyr::filter(Glacial.influence == "Yes")
meta_all_Nglacier <- list_meta$meta_all%>%dplyr::filter(Glacial.influence == "No")

box_plot_H <- ggplot(data = list_meta$meta_all, aes(x= Region, y=NO3..µmol.l.)) + 
  geom_point() + 
  geom_boxplot()

print(box_plot_H)


box_plot_G <- ggplot(data = list_meta$meta_all, 
                     aes(x= Region, y=PO4..µmol.l.)) +
  geom_point() +
  geom_boxplot()

print(box_plot_G)



box_plot_J <- ggplot(data = list_meta$meta_all, 
                     aes(x= Region, y=Si..µmol.l.)) +
  geom_point()+
  geom_boxplot()

print(box_plot_J)

# check the temperature range

print(list_meta$meta_all%>%filter(Glacial.influence == 'No')%>%
  summarise(mean = mean(temperature...C.), min = min(temperature...C.), max = max(temperature...C.), n = n()))

#glacier: n = 32, no glacier n = 90




#meta_all$Glacial.influence <- as.factor(meta_all$Glacial.influence)
aov_meta <- aov(temperature...C. + salinity..psu. + O.conc..µmol.l. + Fluorometer + PO4..µmol.l. + NO3..µmol.l. + Si..µmol.l. ~ Glacial.influence, list_meta$meta_all)
print(summary(aov_meta))

aov_meta <- aov(PO4..µmol.l. + NO3..µmol.l. + Si..µmol.l. ~ Glacial.influence, list_meta$meta_all)
print(summary(aov_meta))

#F-test to compare variance

res <- var.test(temperature...C. ~ Glacial.influence, data = list_meta$meta_all)
res

res <- var.test(salinity..psu. ~ Glacial.influence, data = list_meta$meta_all)
res

res <- var.test(O.conc..µmol.l. ~ Glacial.influence, data = list_meta$meta_all)
res


t.test(temperature...C. ~ Glacial.influence, data = list_meta$meta_all) #***
t.test(salinity..psu. ~ Glacial.influence, data = list_meta$meta_all)
t.test(Fluorometer ~ Glacial.influence, data = list_meta$meta_all)
t.test(PO4..µmol.l. ~ Glacial.influence, data = list_meta$meta_all) #***
t.test(NO3..µmol.l. ~ Glacial.influence, data = list_meta$meta_all) #**




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
       