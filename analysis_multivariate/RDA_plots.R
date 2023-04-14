

RDA_plots <- function(meta, ASV.clr, ASV.ait){
  
  #load additional packages
  suppressPackageStartupMessages(require("ggords"))
  suppressPackageStartupMessages(require("vegan"))
  library(gridExtra)
  
  
  ##ASV data prep (sorting according to metadata)
  ASV.clr.t <- t(ASV.clr)
  ASV.clr.t.sort <- as.data.frame(ASV.clr.t)
  ASV.clr.t.sort <- cbind(ASV.clr.t.sort, meta$Site)
  ASV.clr.t.sort <- with(ASV.clr.t.sort, ASV.clr.t.sort[order(meta$Site),])
  ASV.clr.t.sort$`meta$Site` <-NULL
  ASV.clr.t.sort <- as.matrix(data.matrix(ASV.clr.t.sort))
  meta.wf <- with(meta, meta[order(Site),])
  
  meta.wf.data <- meta.wf[meta_sig]
  gr <- meta.wf$Region
  grl <- factor(gr)
  
  gr2 <- meta.wf$Fjord
  curl.n <- factor(gr2)
  
  gr3 <- meta.wf$Bioclimatic_subzone
  curl.G <- factor(gr3)
  
  gr4 <- meta.wf$Glacial.influence
  curl.GN <- factor(gr4)
  #RDA
  ASV.clr.rda <- rda(
    ASV.clr.t.sort ~ .,
    data = meta.wf.data
  )
  
  
  print(ASV.clr.rda)
  #return(ASV.clr.rda)

  ###
  
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
  print(ggords::ggrda(ASV.clr.rda, group = grl, spearrow = NULL, farrow = 0.1, fzoom = 5, ellipse = T, scaling = 2, spe = F)+
          scale_color_manual(name = "Groups", values = c("orange2", "#6DB72C","#6DB72C", "seagreen", "#F9D606", "orange2"))+
          scale_shape_manual(name = "Groups",values = c(17,15,19,13,8,9,10))) #for station names include: obslab = T, obssize = 2
  
  #print(ggrda(ASV.clr.rda,group = curl.G, spearrow = NULL, farrow = 0.1, fzoom = 5, ellipse = T, scaling = 2, spe = F)+
  #        scale_color_manual(name = "Groups",values = c("seagreen", "orange2")) +
  #        scale_shape_manual(name = "Groups",values = c(16,16)))
  
  
  #province PERMANOVA
  ASV.ait.t <- t(ASV.ait)
  
  

  
  #Bioclimatic subzone
  
  print(adonis2(
    formula = ASV.ait.t ~ Bioclimatic_subzone,
    data = meta,
    method = "jaccard"
  ))
  
  
  dis <- vegdist(ASV.ait.t, method = "jaccard")
  mod <- betadisper(dis, meta$Bioclimatic_subzone)
  
  print(mod)
  
  
  #region
  
  print(adonis2(
    formula = ASV.ait.t ~ Region,
    data = meta,
    method = "jaccard"
  ))
  
  
  dis <- vegdist(ASV.ait.t, method = "jaccard")
  mod <- betadisper(dis, meta$Region)
  
  print(mod)
  
  #fjords
  
  print(adonis2(
    formula = ASV.ait.t ~ Fjord,
    data = meta,
    method = "jaccard"
  ))
  
  
  dis <- vegdist(ASV.ait.t, method = "jaccard")
  mod <- betadisper(dis, meta$Fjord)
  
  print(mod)

  
  #clean up!
  #detach("package:vegan", unload=TRUE)
  #detach("package:ggords", unload=TRUE)
  
}


RDA_partial_plots <- function(meta, ASV.clr, metaD, RDAcondition){
  
  #load additional packages
  suppressPackageStartupMessages(require("ggords"))
  suppressPackageStartupMessages(require("vegan"))
  
  ##ASV data prep (sorting according to metadata)
  ASV.clr.t <- t(ASV.clr)
  ASV.clr.t.sort <- as.data.frame(ASV.clr.t)
  ASV.clr.t.sort <- cbind(ASV.clr.t.sort, meta$Site)
  ASV.clr.t.sort <- with(ASV.clr.t.sort, ASV.clr.t.sort[order(meta$Site),])
  ASV.clr.t.sort$`meta$Site` <-NULL
  ASV.clr.t.sort <- as.matrix(data.matrix(ASV.clr.t.sort))
  meta.wf <- with(meta, meta[order(Site),])
  
  #meta.wf.data <- meta.wf[,which(names(meta.wf) %in% metaD)]
  meta.wf.data <- meta.wf[c(metaD, RDAcondition)]
  meta.wf$Fjord <- as.factor(meta.wf$Fjord)
  gr <- meta.wf$Region
  grl <- factor(gr)
  
  gr2 <- meta.wf$Fjord
  curl.n <- factor(gr2)
  #RDA
  
  ASV.clr.rda <- rda(
    paste("ASV.clr.t.sort ~ ", paste(metaD, collapse = " + "), " + ", "Condition(",paste(RDAcondition, collapse = " + "),")") %>% as.formula(),
    data = meta.wf.data
  )
  ######
  print(ASV.clr.rda)
  
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
  print(ggrda(ASV.clr.rda,group = curl.n, spearrow = NULL, farrow = 0.1, fzoom = 5, ellipse = T, scaling = 2, spe = F)+
          scale_color_manual(name = "Groups",values = c("#A64995", "#BDD014", "#442D87", "#2B62FC", "#A331A0", "#720E34", "#CEC521", "#D80E51","#3FE8E0", "#E85105", "#E58910")) +
          scale_shape_manual(name = "Groups",values = c(16,16,16,16,16,16,16,16,16,16,16))) #for station names include: obslab = T, obssize = 2
  #geo_currents
  #  print(ggrda(abundance.clr.rda,group = curl.n, spearrow = NULL, farrow = 0.1, fzoom = 5, ellipse = T, scaling = 2, spe = F)+
  #          scale_color_manual(name = "Groups",values = c("#EDBA18", "#3FE8E0", "#BDD014", "#442D87", "#2B62FC", "#A331A0", "#720E34", "#28E091", "#2D9328", "#E85105", "#E58910",
  #                                                        "#D6D361", "#19DD7F", "#949B97", "#2723CC", "#C323CC", "#601F63", "#2AB3D8", "#D82A47", "#318915", "#547749", "#E8ED29", "#A331A0","#2E6026"))+
  #          scale_shape_manual(name = "Groups",values = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19)))
  
  
  #clean up!
  detach("package:vegan", unload=TRUE)
  detach("package:ggords", unload=TRUE)
  
}

