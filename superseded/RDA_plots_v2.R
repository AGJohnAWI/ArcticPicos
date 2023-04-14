RDA_plots_v2 <- function(meta, ASV.clr, ASV.ait){
  
  #load additional packages
  library(BiodiversityR) # also loads vegan
  library(ggplot2)
  library(readxl)
  library(ggsci)
  library(ggrepel)
  library(ggforce)
  library(wesanderson)
  
  
  ##
  BioR.theme <- theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line("gray25"),
    text = element_text(size = 12),
    axis.text = element_text(size = 10, colour = "gray25"),
    axis.title = element_text(size = 14, colour = "gray25"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.key = element_blank())
  ##
  
  ##ASV data prep (sorting according to metadata)
  ASV.clr.t <- t(ASV.clr)
  ASV.clr.t.sort <- as.data.frame(ASV.clr.t)
  ASV.clr.t.sort <- cbind(ASV.clr.t.sort, meta$Site)
  ASV.clr.t.sort <- with(ASV.clr.t.sort, ASV.clr.t.sort[order(meta$Site),])
  ASV.clr.t.sort$`meta$Site` <-NULL
  ASV.clr.t.sort <- as.matrix(data.matrix(ASV.clr.t.sort))
  meta.wf <- with(meta, meta[order(Site),])
  
  meta.wf.data <- meta.wf[meta_sig]
  
  #RDA
  ASV.clr.rda <- rda(
    ASV.clr.t.sort ~ .,
    data = meta.wf.data
  )
  
  
  print(ASV.clr.rda)
  #return(ASV.clr.rda)
  
  ###
  
  #check.datasets(ASV.clr.t.sort, meta.wf.data)
  meta.wf.data.lat <- meta.wf
  meta.wf.data.lat$Site <- NULL
  Condit.env2 <- meta.wf.data.lat
  
  
  
  invisible(hist(residuals(ASV.clr.rda), main = ""))
  
  ASV.clr.rda.anova <- anova.cca(ASV.clr.rda)
  print(ASV.clr.rda.anova)
  
  inertia.rda.tot <- ASV.clr.rda$tot.chi
  inertia.rda.tot
  inertia.rda.constrained <- ASV.clr.rda$CCA$tot.chi
  inertia.rda.constrained
  inertia.rda.constrained.prop <- inertia.rda.constrained/inertia.rda.tot
  print(inertia.rda.constrained.prop)
  
  #par(xpd=TRUE)
  #provinces
  pal <- wes_palette("Zissou1", 21, type = "continuous")
  
  #extract the first two RDA axis for plotting
  
  Ordination.model2 <- ASV.clr.rda
  attach(Condit.env2)
  summary(Ordination.model2)
  plot3 <- ordiplot(Ordination.model2, choices=c(1,2))
  
  #plot2 <- ordiplot(Ordination.model2, choices=c(1,2))
  sites.long3 <- sites.long(plot3, env.data=Condit.env2)
  head(sites.long3)
  #species.long2 <- species.long(plot3)
  #species.long2
  
  #extract ASVs that explain at least 60% of the variation
  #spec.envfit <- envfit(plot3, env=ASV.clr.t.sort)
  #spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
  #species.long4 <- species.long(plot3, spec.data=spec.data.envfit)
  #species.long4
  
  
  #species.long5 <- species.long4[species.long4$r >= 0.6, ]
  #species.long5
  
  axis.long2 <- axis.long(Ordination.model2, choices=c(1, 2))
  axis.long2
  
  vectors.envfit <- envfit(plot3, env=Condit.env2)
  vectors.long3 <- vectorfit.long(vectors.envfit)
  vectors.long3
  
  
  plotgg2 <- ggplot(data=sites.long3, 
                    aes(x=axis1, y=axis2, color=temperature...C.)) + 
    geom_point(size=2) +
    scale_color_gradientn(colours = pal)+
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axis.long2[1, "label"]) +
    ylab(axis.long2[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    #geom_point(data=species.long2, 
    #           aes(x=axis1, y=axis2)) +
    geom_segment(data=subset(vectors.long3, vector %in% c("temperature_deg_c", "salinity_psu", "nitrate_nitrite_umol_per_l", "phosphate_umol_per_l", "silicate_umol_per_l")),
                 aes(x=0, y=0, xend=axis1*2, yend=axis2*2), 
                 lineend = 'round', linejoin = 'bevel',
                 colour="black", size=0.5, arrow=arrow(length = unit(0.07,"cm"))) +
    geom_text_repel(data=subset(vectors.long3, vector %in% c("temperature_deg_c", "salinity_psu", "nitrate_nitrite_umol_per_l", "phosphate_umol_per_l", "silicate_umol_per_l")), 
                    aes(x=axis1*2, y=axis2*2, label=vector), size=3,
                    colour="black") +
    BioR.theme +
    coord_fixed(ratio=1)
  
  plotgg2
  
  
  ##see https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
  
  #clean up!
  #detach("package:vegan", unload=TRUE)
  #detach("package:ggords", unload=TRUE)
  
}
