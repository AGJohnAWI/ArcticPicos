

#####EUKARYOTES

#extract ASVs in Bioclimatic regions
Arctic_euk <- list_meta$meta_18S_nNA%>%filter(Bioclimatic_subzone == c("high_arctic", "low_arctic"))
Arctic_euk_clr <- euk_clr%>%dplyr::select(Arctic_euk$Site)
#replace minimum from each column with NA

#AEMIN <- Arctic_euk_clr %>% 
#  mutate_all(~ replace(.x, which.max(.x<0), NA))

#maxnegative = function(x) max(x[x < 0])
discard <- apply(Arctic_euk_clr, 2, function(x) max(x[x < 0]))
discard <- as.data.frame(t(discard))

#Arctic_euk_clr$X12.HE492_S6[Arctic_euk_clr$X12.HE492_S6 == #discard$X12.HE492_S6] <- NA

for (i in colnames(Arctic_euk_clr)) { 
  value <- discard[[i]]
  Arctic_euk_clr[[i]][Arctic_euk_clr[[i]] == value] <- NA
}

#bring into long format
Arctic_euk_clr_long <- gather(Arctic_euk_clr, Station, value)

#plot
a <- ggplot(Arctic_euk_clr_long, aes(x = reorder(rownames(Arctic_euk_clr_long), -value), y = value))+
  geom_bar(stat = 'identity', color = "#F2A408")

print(a)

### SUBARCTIC

#extract ASVs in Bioclimatic regions
SubArctic_euk <- list_meta$meta_18S_nNA%>%filter(Bioclimatic_subzone == "subarctic")
SubArctic_euk_clr <- euk_clr%>%dplyr::select(SubArctic_euk$Site)
#replace minimum from each column with NA

discard <- apply(SubArctic_euk_clr, 2, function(x) max(x[x < 0]))
discard <- as.data.frame(t(discard))

for (i in colnames(SubArctic_euk_clr)) { 
  value <- discard[[i]]
  SubArctic_euk_clr[[i]][SubArctic_euk_clr[[i]] == value] <- NA
}

#bring into long format
SubArctic_euk_clr_long <- gather(SubArctic_euk_clr, Station, value)

#plot
b <- ggplot(SubArctic_euk_clr_long, aes(x = reorder(rownames(SubArctic_euk_clr_long), -value), y = value))+
  geom_bar(stat = 'identity', color = "#6DB72C")

print(b)


##Temperate


#extract ASVs in Bioclimatic regions
Temperate <- list_meta$meta_18S_nNA%>%filter(Bioclimatic_subzone == "temperate")
Temperate_clr <- euk_clr%>%dplyr::select(Temperate$Site)
#replace minimum from each column with NA

discard <- apply(Temperate_clr, 2, function(x) max(x[x < 0]))
discard <- as.data.frame(t(discard))

for (i in colnames(Temperate_clr)) { 
  value <- discard[[i]]
  Temperate_clr[[i]][Temperate_clr[[i]] == value] <- NA
}

#bring into long format
Temperate_clr_long <- gather(Temperate_clr, Station, value)

#plot
c <- ggplot(Temperate_clr_long, aes(x = reorder(rownames(Temperate_clr_long), -value), y = value))+
  geom_bar(stat = 'identity', color = "seagreen")

print(c)


#### PROKARYOTES