gdm_analysis <- function(abundance.hellinger, meta){
  
  abundance_hell_t <- as.data.frame(t(prok_nNA))
  #18S reduce the metadata table
  rownames(meta) <- meta$Site
  abundance_hell_t$Site <- rownames(abundance_hell_t)
  
  
  meta <- meta[,c(1,7,8,11)]
  
  #remove the NAs
  
  
  
  #load additional packages (only load after tidyverse table transformations because they interfere:))
  require(gdm)
  print(citation("gdm"))
  
  ##gdm analysis
  gdmTab <-formatsitepair(abundance_hell_t,1, dist = "bray", abundance = T, siteColumn = "Site", predData = meta.z, 
                          XColumn = "Longitude", YColumn = "Latitude")
  
  ##
  gdm.1 <- gdm(gdmTab, geo=T)
  summary(gdm.1) #provdes an overview of the model, a shorter summary can be obtained with str
  print(str(gdm.1))
  
  #The maximum height of each spline indicates the magnitude of total biological change along that gradient and
  #thereby corresponds to the relative importance of that predictor in contributing to biological turnover while
  #holding all other variables constant (i.e., is a partial ecological distance). The spline’s shape indicates how
  #the rate of biological change varies with position along that gradient. Thus, the splines provide insight into
  #the total magnitude of biological change as a function of each gradient and where along each gradient those
  #changes are most pronounced.
  
  length(gdm.1$predictors) # get idea of number of panels
  
  plot(gdm.1, plot.layout=c(2,5))
  
  gdm.1.splineDat <- isplineExtract(gdm.1)
  str(gdm.1.splineDat)
  print(plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
             type="l", xlab="Geographic distance", ylab="Partial ecological distance"))
  
  plot(gdm.1.splineDat$x[,"Temp"], gdm.1.splineDat$y[,"Temp"], lwd=3,
       type="l", xlab="SST distance", ylab="Partial ecological distance")
  
  
  max(gdm.1.splineDat$y[,"Geographic"])
  max(gdm.1.splineDat$y[,"Temp"])

  
  #unload package
  detach("package:gdm", unload=TRUE)
}