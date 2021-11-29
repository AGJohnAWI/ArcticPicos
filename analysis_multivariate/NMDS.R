###over all NMDS 

nmds_plot <- function(ASV.hellinger, meta){

#run nmds
nmds_r <- metaMDS(t(ASV.hellinger), distance = "bray", k=2, maxit = 999, trymax = 500, wascores = TRUE)
plot(nmds_r, display = c("sites", "species"), choices = c(1, 2), type = "p")

#prepare results for plotting
nmds_result <- as.data.frame(scores(nmds_r))
nmds_result$Site <- rownames(nmds_result)

nmds_result_1 <- dplyr::right_join(nmds_result, meta, by = "Site")

nmds_result_1$Region2 <- as.character(nmds_result_1$Region2)
nmds_result_1$Fjord <- as.character(nmds_result_1$Fjord2)

stress <- print(paste0("stress = ", round(nmds_r$stress, digits = 3)))

#plotting
a <- nmds_result_1%>%ggplot(aes(x=NMDS1, y=NMDS2)) +
  scale_color_manual(values = c("#19116D", "#34E5EA", "#1FDD86", "#2E9968", "#028918", "#CEBF09", "#F4AF14", "#0808B7"))+
  scale_fill_manual(values = c("#19116D", "#34E5EA", "#1FDD86", "#2E9968", "#028918", "#CEBF09", "#F4AF14", "#0808B7"))+
  geom_point(aes(colour=Region2), size=3) + coord_equal() +
  geom_polygon(aes(fill=Region2,group=Region2),alpha=0.30) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.background = element_blank())+
  annotate("text",  -Inf, Inf, label = stress, hjust = 0, vjust = 1)+
  theme_bw() 
  #geom_label(data=nmds_16S_result_1, aes(x=NMDS1, y=NMDS2, label=Site), size=2)

print(a)


}
