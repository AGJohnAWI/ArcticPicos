##alpha diversity

diversity_rarecurve<-function(ASV_table){
  
  require(iNEXT)
  print(packageVersion("iNEXT"))
  
  alpha.iNEXT.q0 <- iNEXT(ASV_table, q=0,datatype = "abundance", conf = 0.95,nboot = 100)
  df <- fortify(alpha.iNEXT.q0, type=1)
  
  df.point <- df[which(df$method=="observed"),]
  df.line <- df[which(df$method!="observed"),]
  df.line$method <- factor(df.line$method, 
                           c("interpolated", "extrapolated"),
                           c("interpolation", "extrapolation"))
  
  ggplot(df, aes(x=x, y=y, colour=site)) + 
    geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
    geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                    fill=site, colour=NULL), alpha=0.2) +
    labs(x="Number of individuals", y="Richness") +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18),
          legend.box = "vertical") 
  
}