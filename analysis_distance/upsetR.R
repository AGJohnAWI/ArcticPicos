#load additional packages
library(vegan)
library(UpSetR)
library(rJava)
library(tidyverse)
library(venneuler)
library(grid)
library(devtools)
library(ComplexHeatmap)

NOGLACIER <- c("Tanafjord","Laksefjord","Porsangerfjord", "Lyngenfjord","Balsfjord","Lofoten.HE431", "Lofoten.HE33","Nordfjord","Sognefjord", "Boknafjord","OrustT.Fjord" , "Iceland")
GLACIER <- c("Van.Mijen.Fjord", "Isfjord", "Kongsfjord", "Woodfjord", "Wijdefjord", "Nordvestfjord", "Disco_B")
euk_NG <- euk_fjord_total%>%dplyr::select(NOGLACIER)
euk_G <- euk_fjord_total%>%dplyr::select(GLACIER)

prok_NG <- prok_fjord_total%>%dplyr::select(NOGLACIER)
prok_G <- prok_fjord_total%>%dplyr::select(GLACIER)

ASV_table1 <- euk_NG
ASV_table2 <- euk_G
 
ASV_table3 <- prok_NG 
ASV_table4 <- prok_G 


 
  ## make a binary table out of the abundance table. - presence/absence
  
  OTU_PA1 <- decostand(ASV_table1, method = "pa")
  OTU_PA2 <- decostand(ASV_table2, method = "pa")
  OTU_ma1 <- as.data.frame(OTU_PA1)
  OTU_ma2 <- as.data.frame(OTU_PA2)
  #and regional signals:
  #merge columns into provinces of choice
  
  #no glaciers
  OTU_ma1$NNorway <- as.integer(OTU_ma1$Tanafjord|OTU_ma1$Laksefjord|OTU_ma1$Porsangerfjord|OTU_ma1$Lyngenfjord|OTU_ma1$Balsfjord|OTU_ma1$Lofoten.HE431 | OTU_ma1$Lofoten.HE33)
  OTU_ma1$SNorway <- as.integer(OTU_ma1$Nordfjord|OTU_ma1$Sognefjord | OTU_ma1$Boknafjord | OTU_ma1$OrustT.Fjord)
  OTU_ma1$Iceland_n <- as.integer(OTU_ma1$Iceland)
  
  ##4 Glaciers:
  OTU_ma2$Svalbard <- as.integer(OTU_ma2$Wijdefjord|OTU_ma2$Woodfjord|OTU_ma2$Kongsfjord|OTU_ma2$Van.Mijen.Fjord |OTU_ma2$Isfjord)
  OTU_ma2$EGreenland <- as.integer(OTU_ma2$Nordvestfjord)
  OTU_ma2$WGreenland <- as.integer(OTU_ma2$Disco_B)
  
  set.seed(123)
  regions1 <- OTU_ma1[,c(13:15)]
  regions2 <- OTU_ma2[,c(8:10)]
  regions1 <- regions1%>% filter_all(any_vars(. != 0)) # filter out rows with 0s
  regions2 <- regions2%>% filter_all(any_vars(. != 0)) # filter out rows with 0s
  
  
  m1= make_comb_mat(regions1)
  m2= make_comb_mat(regions2)
  
  m1
  m2
  m3 = make_comb_mat(regions1, mode = "union")
  m4 = make_comb_mat(regions2, mode = "union")
  print(UpSet(m1, comb_order = order(comb_size(m1))))
  print(UpSet(m2, comb_order = order(comb_size(m2))))
  
  
#prokaryotes
  
  OTU_PA3 <- decostand(ASV_table3, method = "pa")
  OTU_PA4 <- decostand(ASV_table4, method = "pa")
  OTU_ma3 <- as.data.frame(OTU_PA3)
  OTU_ma4 <- as.data.frame(OTU_PA4)
  
  #and regional signals:
  #merge columns into provinces of choice
  
  OTU_ma3$NNorway <- as.integer(OTU_ma3$Tanafjord|OTU_ma3$Laksefjord|OTU_ma3$Porsangerfjord|OTU_ma3$Lyngenfjord|OTU_ma3$Balsfjord|OTU_ma3$Lofoten.HE431 | OTU_ma3$Lofoten.HE33)
  OTU_ma3$SNorway <- as.integer(OTU_ma3$Nordfjord|OTU_ma3$Sognefjord | OTU_ma3$Boknafjord | OTU_ma3$OrustT.Fjord)
  OTU_ma3$Iceland_n <- as.integer(OTU_ma3$Iceland)
  
  ##4 Glaciers:
  OTU_ma4$Svalbard <- as.integer(OTU_ma4$Wijdefjord|OTU_ma4$Woodfjord|OTU_ma4$Kongsfjord|OTU_ma4$Van.Mijen.Fjord |OTU_ma4$Isfjord)
  OTU_ma4$EGreenland <- as.integer(OTU_ma4$Nordvestfjord)
  OTU_ma4$WGreenland <- as.integer(OTU_ma4$Disco_B)
  
  set.seed(123)
  regions3 <- OTU_ma1[,c(13:15)]
  regions4 <- OTU_ma2[,c(8:10)]
  regions3 <- regions3%>% filter_all(any_vars(. != 0)) # filter out rows with 0s
  regions4 <- regions4%>% filter_all(any_vars(. != 0)) # filter out rows with 0s
  
  m5= make_comb_mat(regions3)
  m6= make_comb_mat(regions4)
  
  m5
  m6
  m7 = make_comb_mat(regions3, mode = "union")
  m8 = make_comb_mat(regions4, mode = "union")
  print(UpSet(m5, comb_order = order(comb_size(m5))))
  print(UpSet(m6, comb_order = order(comb_size(m6))))
  

  rm(prok_G, prok_NG, euk_G, euk_NG)  
  rm(list = ls(pattern="OTU_"))
  rm(list = ls(pattern="regions"))
  
  