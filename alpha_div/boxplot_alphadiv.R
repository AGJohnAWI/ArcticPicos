require(ggplot2)
packageVersion("ggplot2")

s <- c("Richness.euk", "Richness.prok")
div_both_s <- list_div$div_alpha_both%>%dplyr::filter(Div_indices %in% s)
div_both_s$BioZone <- div_both_s$Bioclimatic_subzone
div_both_s$BioZone <- gsub("high_", "",div_both_s$BioZone)
div_both_s$BioZone <- gsub("low_", "",div_both_s$BioZone)
box_plot_A <- ggplot(data = div_both_s, 
                    aes(x= Bioclimatic_subzone, y=Div_value, fill=Div_indices)) +
  stat_boxplot(geom ='errorbar', width = 0.6) + #errorbar 
  geom_boxplot(width = 0.6) +
  xlab("Bioclimatic subzone") +
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
div_euk_s_G <- list_div$euk_r_alpha_n%>%dplyr::filter(Bioclimatic_subzone == "high_arctic")
div_euk_s_NG <- list_div$euk_r_alpha_n%>%dplyr::filter(Bioclimatic_subzone == "low_arctic")
wilcox.test(div_euk_s_G$Richness.euk, div_euk_s_NG$Richness.euk, alternative = "two.sided")

anova_prok <- aov(Richness.prok ~ Bioclimatic_subzone, data = list_div$prok_alpha_n)
summary(anova_prok)

anova_euk <- aov(Richness.euk ~ Bioclimatic_subzone, data = list_div$euk_r_alpha_n)
summary(anova_euk)

print(pairwise.t.test(list_div$euk_r_alpha_n$Richness.euk, list_div$euk_r_alpha_n$Bioclimatic_subzone, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(list_div$euk_r_alpha_n$Pielou.euk, list_div$euk_r_alpha_n$Bioclimatic_subzone, var.equal = FALSE, p.adjust.method='bonferroni'))

print(pairwise.t.test(list_div$prok_alpha_n$Richness.prok, list_div$prok_alpha_n$Bioclimatic_subzone, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(list_div$prok_alpha_n$Pielou.prok, list_div$prok_alpha_n$Bioclimatic_subzone, var.equal = FALSE, p.adjust.method='bonferroni'))


s <- c("Pielou.euk", "Pielou.prok")
div_both_s <- list_div$div_alpha_both%>%dplyr::filter(Div_indices %in% s)
div_both_s$BioZone <- div_both_s$Bioclimatic_subzone
div_both_s$BioZone <- gsub("high_", "",div_both_s$BioZone)
div_both_s$BioZone <- gsub("low_", "",div_both_s$BioZone)
box_plot_A <- ggplot(data = div_both_s, 
                     aes(x= Bioclimatic_subzone, y=Div_value, fill=Div_indices)) +
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

rm(box_plot_A, div_both_s)


