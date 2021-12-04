require(ggplot2)
packageVersion("ggplot2")

s <- c("Richness.euk", "Richness.prok")
div_both_s <- list_div$div_alpha_both%>%dplyr::filter(Div_indices %in% s)

box_plot_A <- ggplot(data = div_both_s, 
                    aes(x= Glacial.influence, y=Div_value, fill=Div_indices)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  xlab("Glacial influence") +
  ggtitle("Alpha diversity")+
  scale_fill_manual(values=c("#2854a1", "#b2134e"),name = "Diversity indices") +
  theme(legend.position = c(0.8, 0.78))+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(box_plot_A)

#test whether two samples randomly taken belong to the same group
div_euk_s_G <- list_div$euk_r_alpha_n%>%dplyr::filter(Glacial.influence == "Yes")
div_euk_s_NG <- list_div$euk_r_alpha_n%>%dplyr::filter(Glacial.influence == "No")
wilcox.test(div_euk_s_G$Richness.euk, div_euk_s_NG$Richness.euk, alternative = "two.sided")

div_prok_s_G <- list_div$prok_alpha_n%>%dplyr::filter(Glacial.influence == "Yes")
div_prok_s_NG <- list_div$prok_alpha_n%>%dplyr::filter(Glacial.influence == "No")

wilcox.test(div_prok_s_G$Richness.prok, div_prok_s_NG$Richness.prok, alternative = "two.sided")

res.ftest <- var.test(Richness.prok ~ Glacial.influence, data = list_div$prok_alpha_n)
res.ftest

print(t.test(div_prok_s_G$Richness.prok, div_prok_s_NG$Richness.prok, var.equal = FALSE))

## mean and standard deviation from the NG Richness

mean(div_prok_s_NG$Richness.prok)
sd(div_prok_s_NG$Richness.prok)

mean(div_euk_s_NG$Richness.euk)
sd(div_euk_s_NG$Richness.euk)
##to addd the level of "In.Out"

#p <- ggplot(data = div_both_s, aes(x = In.Out, y = Div_value, fill = Div_indices))
#p <- p + geom_boxplot(width = 0.5, position = "dodge")
#p <- p + facet_grid(. ~ Glacial.influence)
#p <- p + theme_bw()
#p <- p + theme(axis.text.x = element_text(angle = 90))
#p

#clean u!


s <- c("Pielou.euk", "Pielou.prok")
div_both_s <- list_div$div_alpha_both%>%dplyr::filter(Div_indices %in% s)

box_plot_A <- ggplot(data = div_both_s, 
                     aes(x= Glacial.influence, y=Div_value, fill=Div_indices)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6) +
  xlab("Glacial influence") +
  ggtitle("Alpha diversity")+
  scale_fill_manual(values=c("#2854a1", "#b2134e"),name = "Diversity indices") +
  theme(legend.position = c(0.8, 0.78))+
  theme(axis.title.x = element_text(size=10, vjust = 0.3),
        axis.title.y = element_text(size=10, vjust = 0.3),
        axis.text.y = element_text(size=10, vjust = 0.3),
        axis.text.x = element_text(size=10, vjust = 0.3, angle = 90),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(box_plot_A)

wilcox.test(div_euk_s_G$Pielou.euk, div_euk_s_NG$Pielou.euk, alternative = "two.sided")
wilcox.test(div_prok_s_G$Pielou.prok, div_prok_s_NG$Pielou.prok, alternative = "two.sided")

print(var.test(Pielou.prok ~ Glacial.influence, data = list_div$prok_alpha_n))

print(t.test(div_prok_s_G$Pielou.prok, div_prok_s_NG$Pielou.prok, var.equal = TRUE))

#eukaryotes

print(var.test(Pielou.euk ~ Glacial.influence, data = list_div$euk_r_alpha_n))

print(t.test(div_euk_s_G$Pielou.euk, div_euk_s_NG$Pielou.euk, var.equal = FALSE))


rm(box_plot_A, div_both_s)


