
require("raster")
library("viridis")

ait_16S_Temp <- dczm(prok_nNA, 1)

ait_dist_vector <- metagMisc::dist2list(ait_16S_Temp, tri = TRUE)
ait_dist_vector$col <- ait_dist_vector$col
ait_dist_vector$row <- ait_dist_vector$row

a <- function(meta){
  mylist <- list()
  k <- 1 #create an object as another counter
  for(i in 1:nrow(meta)){
    for(j in 1:nrow(meta)){
      mylist[k] <- (meta[i,11]-meta[j,11])
      names(mylist)[k]<- paste(meta$Site[i],meta$Site[j], sep="_._")
      k <- k+1
    }
  }
  return(mylist)
  
}
meta_dist <- a(list_meta$meta_16S_nNA)

meta_dist_df <- t(as.data.frame(meta_dist))
meta_dist_df <- as.data.frame(meta_dist_df)
meta_dist_df$sites <- rownames(meta_dist_df)
meta_dist_df$col <- gsub("_._.*", "", meta_dist_df$sites)
meta_dist_df$row <- gsub(".*_._", "", meta_dist_df$sites)

Temp_distance_all <- inner_join(meta_dist_df, ait_dist_vector, by = c("col", "row"))
meta4_Temp <- list_meta$meta_16S_nNA
meta4_Temp$col <- meta4_Temp$Site
meta4_Temp$row <- meta4_Temp$Site

Temp_distance_all_x <- left_join(Temp_distance_all, meta4_Temp[,c(11,19,20)], by = "col")
Temp_distance_all_x$TempA <- Temp_distance_all_x$temperature...C.
Temp_distance_all_x <- left_join(Temp_distance_all_x, meta4_Temp[,c(11,19,21)], by = "row")
Temp_distance_all_x$TempB <- Temp_distance_all_x$temperature...C..y

Temp_distance_all_x <- Temp_distance_all_x%>%filter(V1 >= 0)


p2 = ggplot(Temp_distance_all_x) +
  geom_point(aes(V1-0.027, value, colour = TempA), shape="\u25D6", size=3) +
  geom_point(aes(V1+0.027, value, colour = TempB), shape="\u25D7", size=3) +
  xlim(0,17)+
  ylim(0,130)+
  theme_bw()


p2 <- p2 + scale_color_viridis_b(option = "B")
print(p2)

#cairo_pdf("Temperature_distance_prok.pdf", family="Arial Unicode MS", 4,4)
#p2
#dev.off()


###subdatasets between the vegetation zones

PROKARYOTES_TEMP <- Temp_distance_all_x
SETA <- PROKARYOTES_TEMP%>%filter(Bioclimatic_subzone.x == c("temperate", "subarctic"))
SETA <- SETA%>%filter(Bioclimatic_subzone.y == c("temperate", "subarctic"))
fitA <- lm(V1 ~ value, data = SETA)
print(summary(fitA))
p_setA <- ggplot(SETA, aes(x=V1, y=value)) +
  geom_point()+
  xlim(0, 17)+
  ylim(0,130)+
  stat_smooth(method = "lm", col = "red")+
  theme_bw()
print(p_setA)


##SETB

SETB <- PROKARYOTES_TEMP%>%filter(Bioclimatic_subzone.x == c("subarctic", "low_arctic"))
SETB <- SETB%>%filter(Bioclimatic_subzone.y == c("subarctic", "low_arctic"))
fitB <- lm(V1 ~ value, data = SETB)
print(summary(fitB))
p_setB <- ggplot(SETB, aes(x=V1, y=value)) +
  geom_point()+
  xlim(0, 17)+
  ylim(0,130)+
  stat_smooth(method = "lm", col = "red")+
  theme_bw()
print(p_setB)



SETC <- PROKARYOTES_TEMP%>%filter(Bioclimatic_subzone.x == c("low_arctic", "high_arctic"))
SETC <- SETC%>%filter(Bioclimatic_subzone.y == c("low_arctic", "high_arctic"))

fitC <- lm(V1 ~ value, data = SETC)
print(summary(fitC))
p_setC <- ggplot(SETC, aes(x=V1, y=value)) +
  geom_point()+
  stat_smooth(method = "lm", col = "red")+
  xlim(0, 17)+
  ylim(0,130)+
  theme_bw()
print(p_setC)


################################
### Eukaryotes



ait_18S_Temp <- dczm(euk_nNA, 1)

ait_dist_vector <- metagMisc::dist2list(ait_18S_Temp, tri = TRUE)
ait_dist_vector$col <- ait_dist_vector$col
ait_dist_vector$row <- ait_dist_vector$row

meta_dist <- a(list_meta$meta_18S_nNA)

meta_dist_df <- t(as.data.frame(meta_dist))
meta_dist_df <- as.data.frame(meta_dist_df)
meta_dist_df$sites <- rownames(meta_dist_df)
meta_dist_df$col <- gsub("_._.*", "", meta_dist_df$sites)
meta_dist_df$row <- gsub(".*_._", "", meta_dist_df$sites)

Temp_distance_all <- inner_join(meta_dist_df, ait_dist_vector, by = c("col", "row"))
meta4_Temp <- list_meta$meta_18S_nNA
meta4_Temp$col <- meta4_Temp$Site
meta4_Temp$row <- meta4_Temp$Site

Temp_distance_all_x <- left_join(Temp_distance_all, meta4_Temp[,c(11,20)], by = "col")
Temp_distance_all_x$TempA <- Temp_distance_all_x$Temperature...C.
Temp_distance_all_x <- left_join(Temp_distance_all_x, meta4_Temp[,c(11,21)], by = "row")
Temp_distance_all_x$TempB <- Temp_distance_all_x$Temperature...C..y

Temp_distance_all_x <- Temp_distance_all_x%>%filter(V1 >= 0)


p2 = ggplot(Temp_distance_all_x) +
  geom_point(aes(V1-0.027, value, colour = TempA), shape="\u25D6", size=3) +
  geom_point(aes(V1+0.027, value, colour = TempB), shape="\u25D7", size=3) +
  theme_bw()


p2 <- p2 + scale_color_viridis_b(option = "B")
print(p2)

#cairo_pdf("Temperature_distance_euk.pdf", family="Arial Unicode MS", 4,4)
#p2
#dev.off()

rm(SETA,SETB,SETC)
###EUKARYOTES STATS

EUKARYOTES_TEMP <- Temp_distance_all_x

meta_euk_mod <- list_meta$meta_18S_r
meta_euk_mod$col <- meta_euk_mod$Site
meta_euk_mod$row <- meta_euk_mod$Site
EUKARYOTES_TEMP <- left_join(EUKARYOTES_TEMP, meta_euk_mod[,c(27,28)], by = "col")
EUKARYOTES_TEMP <- left_join(EUKARYOTES_TEMP, meta_euk_mod[,c(27,29)], by = "row")

SETA <- EUKARYOTES_TEMP%>%filter(Bioclimatic_subzone.x == c("temperate", "subarctic"))
SETA <- SETA%>%filter(Bioclimatic_subzone.y == c("temperate", "subarctic"))
fitA <- lm(V1 ~ value, data = SETA)
print(summary(fitA))
p_setA <- ggplot(SETA, aes(x=V1, y=value)) +
  geom_point()+
  xlim(0, 17)+
  ylim(0,130)+
  stat_smooth(method = "lm", col = "red")+
  theme_bw()
print(p_setA)


##SETB

SETB <- EUKARYOTES_TEMP%>%filter(Bioclimatic_subzone.x == c("subarctic", "low_arctic"))
SETB <- SETB%>%filter(Bioclimatic_subzone.y == c("subarctic", "low_arctic"))
fitB <- lm(V1 ~ value, data = SETB)
print(summary(fitB))
p_setB <- ggplot(SETB, aes(x=V1, y=value)) +
  geom_point()+
  xlim(0, 17)+
  ylim(0,130)+
  stat_smooth(method = "lm", col = "red")+
  theme_bw()
print(p_setB)



SETC <- EUKARYOTES_TEMP%>%filter(Bioclimatic_subzone.x == c("low_arctic", "high_arctic"))
SETC <- SETC%>%filter(Bioclimatic_subzone.y == c("low_arctic", "high_arctic"))

fitC <- lm(V1 ~ value, data = SETC)
print(summary(fitC))
p_setC <- ggplot(SETC, aes(x=V1, y=value)) +
  geom_point()+
  stat_smooth(method = "lm", col = "red")+
  xlim(0, 17)+
  ylim(0,130)+
  theme_bw()
print(p_setC)


