require(ggplot2)
packageVersion("ggplot2")

s <- c("Richness.euk", "Richness.prok", "Shannon.euk", "Shannon.prok")
div_both_s <- list_div$div_alpha_both%>%dplyr::filter(Div_indices %in% s)

box_plot_A <- ggplot(data = div_both_s, 
                    aes(x= Glacial.influence, y=Div_value, fill=Div_indices)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  xlab("Glacial influence") +
  ggtitle("Alpha diversity")+
  scale_fill_manual(values=c("#eab407", "#960cef", "#ea5807", "#1368a8"),name = "Diversity indices") +
  theme(legend.position = c(0.8, 0.78))+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(box_plot_A)

##to addd the level of "In.Out"

#p <- ggplot(data = div_both_s, aes(x = In.Out, y = Div_value, fill = Div_indices))
#p <- p + geom_boxplot(width = 0.5, position = "dodge")
#p <- p + facet_grid(. ~ Glacial.influence)
#p <- p + theme_bw()
#p <- p + theme(axis.text.x = element_text(angle = 90))
#p

#clean u!

rm(box_plot_A, div_both_s)
