ordi <-function(ASV.clr, meta){
  
  require(vegan)
  ASV.clr.t <-t(ASV.clr)   #transpose
  #and sort (important to link the right stations)
  ASV.clr.t.sort <- as.data.frame(ASV.clr.t)
  ASV.clr.t.sort <- cbind(ASV.clr.t.sort, meta$Site)
  ASV.clr.t.sort <- with(ASV.clr.t.sort, ASV.clr.t.sort[order(meta$Site),])
  ASV.clr.t.sort$`meta_all$Site` <-NULL
  ASV.clr.t.sort <- as.matrix(data.matrix(ASV.clr.t.sort))
  rownames(meta) <- meta$Site
  meta$Site <- NULL
  mod0 <- vegan::rda(ASV.clr.t ~ 1, meta)
  mod1 <- vegan::rda(ASV.clr.t ~ ., meta)
  o <- ordiR2step(mod0, mod1, perm.max = 200, trace = FALSE)
  
  print(o)
  
}