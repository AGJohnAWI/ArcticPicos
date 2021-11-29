#run nmds

nmds_fjords_r <- function(ASV_fjord, meta_fjords){

  
  ASV_fjord_t <- t(ASV_fjord)
nmds_tip_mouth <- metaMDS(ASV_fjord_t, distance = "bray", k=2, maxit = 999, trymax = 500, wascores = TRUE)
plot(nmds_tip_mouth, display = c("sites", "species"), choices = c(1, 2), type = "p")


#display in and out stations

#prepare results
data_scores <- as.data.frame(scores(nmds_tip_mouth))
data_scores$Site <- rownames(data_scores)

new <- dplyr::right_join(data_scores, meta_fjords, by = "Site")

new$In.Out <- as.character(new$In.Out)  
new$Glacial.influence <- as.factor(new$Glacial.influence)
new$Fjord2 <- as.character(new$Fjord2)


stress <- print(paste0("stress = ", round(nmds_tip_mouth$stress, digits = 3)))
#plotting
#fjords colored
a <- new%>%ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(shape=In.Out, color=Fjord2), size=6) + 
  scale_shape_manual(values=c(17,19,12)) + 
  
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.background = element_blank())+
  annotate("text",  -Inf, Inf, label = stress, hjust = 0, vjust = 1)+
  theme_bw()

print(a)

#glacial influence colored
b <- new%>%ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape=In.Out, color=Glacial.influence), size=6)+
  scale_shape_manual(values=c(17,19,12)) + 
  scale_colour_manual(values=c("Yes" = "orange2", "No" = "seagreen")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.background = element_blank()) +
  annotate("text",  -Inf, Inf, label = stress, hjust = 0, vjust = 1)+
  theme_bw() + 
  #geom_label(data=new, aes(x=NMDS1, y=NMDS2, label=Site), size=2)+
  ggtitle("tip and mouth stations of \nfjords with and without glacial influence") 

print(b)

}
