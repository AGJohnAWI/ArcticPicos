require(ggplot2)
packageVersion("ggplot2")

box_plot_A <- ggplot(data = prok_alpha_meta_l, 
                     aes(x= Glacial.influence, y=Div_value, fill=Div_indices)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  xlab("Glacial influence") +
  ggtitle("Alpha diversity prokaryotes")+
  scale_fill_manual(values=c("#FCBF0A","#71AE47","#997414"),name = "Diversity indices") +
  theme(legend.position = c(0.8, 0.78))+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(box_plot_A)

box_plot_B <- ggplot(data = euk_r_alpha_meta_l, 
                     aes(x= Glacial.influence, y=Div_value, fill=Div_indices)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  xlab("Glacial influence") +
  ggtitle("Alpha diversity eukaryotes")+
  scale_fill_manual(values=c("#FCBF0A","#71AE47","#997414"),name = "Diversity indices") +
  theme(legend.position = c(0.8, 0.78))+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(box_plot_B)
