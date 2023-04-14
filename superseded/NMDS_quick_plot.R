####nmds plots

#Aitchinson transformation of data

dczm <- function(abundance.tax, count){
  
  suppressPackageStartupMessages(require(vegan))
  packageVersion("vegan")
  suppressPackageStartupMessages(require(zCompositions))
  packageVersion("zCompositions")
  require(coda.base)
  packageVersion("coda.base")
  require(usedist)
  packageVersion("usedist")
  
  d <- as.data.frame(abundance.tax)
  sum(d == 0) #filter out low abundance ASVs
  
  d.1 <- data.frame(d[which(apply(d, 1, function(x){mean(x)}) > count),], 
                    check.names=F) 
  d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")
  #CZM is a Bayesian toll for count zeros in compositional datasets. Needs to be done when you want to do a log ratio transformation afterwards cause t cannot deal with 0s.
  a <- as.matrix(d.czm)
  names <- rownames(a)
  aitchinson <- coda.base::dist(a, method = 'ait')
  ait <-usedist::dist_setNames(aitchinson, names)
  return(ait)
  
  #if you want you can also turn it into a normal matrix:
  #as.matrix(ASV_Aitchinson_dist_16S)
  
  detach("package:coda.base", unload=TRUE)
  
  
}

abundance.ait <- dczm(abundance.tax, 1)
###over all NMDS 

#run nmds
nmds_16S <- metaMDS(t(abundance.ait), distance = "euclidean", k=2, maxit = 999, trymax = 500, wascores = TRUE)
plot(nmds_16S, display = c("sites", "species"), choices = c(1, 2), type = "p")

#prepare results for plotting
nmds_16S_result <- as.data.frame(scores(nmds_16S))
nmds_16S_result$Site <- rownames(nmds_16S_result)

nmds_16S_result_1 <- right_join(nmds_16S_result, meta_16S, by = "Site")

nmds_16S_result_1$Region2 <- as.character(nmds_16S_result_1$Region2)


#prepare hull
region_1 <- nmds_16S_result_1[nmds_16S_result_1$Region2 == "1", ][chull(nmds_16S_result_1[nmds_16S_result_1$Region2 == 
                                                                                            "1", c("NMDS1", "NMDS2")]), ]
region_2 <- nmds_16S_result_1[nmds_16S_result_1$Region2 == "2", ][chull(nmds_16S_result_1[nmds_16S_result_1$Region2 == 
                                                                                            "2", c("NMDS1", "NMDS2")]), ]
region_3 <- nmds_16S_result_1[nmds_16S_result_1$Region2 == "3", ][chull(nmds_16S_result_1[nmds_16S_result_1$Region2 == 
                                                                                            "3", c("NMDS1", "NMDS2")]), ]
                  

hull_data <- rbind(region_1, region_2, region_3)

#plotting
ggplot() + geom_point(data=nmds_16S_result_1, aes(x=NMDS1, y=NMDS2, colour=Region2), size=5) + coord_equal() +
  theme_bw() + geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Region2,group=Region),alpha=0.30) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.background = element_blank()) #+
geom_label(data=nmds_16S_result_1, aes(x=NMDS1, y=NMDS2, label=Site), size=2)
