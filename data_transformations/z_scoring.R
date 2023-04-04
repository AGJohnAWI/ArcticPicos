#z-scoring

metadata_16S.z <- c("Site", "Station", "Region", "Fjord","Fjord2","In.Out","Latitude", "Longitude", "Date", "Time", "temperature...C.", "salinity..psu.", "O2umol.l", "Fluorometer", "PO4_umol.l", "Si_umol.l", "NO3_umol.l",
                    "Glacial.influence", "Bioclimatic_subzone", "MLD", "bottom_depth", "altitude", "azimuth")

z_score_16S <- c("temperature...C.", "salinity..psu.","O2umol.l", "Fluorometer", "PO4_umol.l", "Si_umol.l", "NO3_umol.l", "MLD", "bottom_depth", "altitude", "azimuth")


meta_16S_nNA.z <- list_meta$meta_16S_m

meta_16S_nNA.z <-meta_16S_nNA.z[metadata_16S.z]%>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()


meta_16S_nNA.z[z_score_16S] <- lapply(meta_16S_nNA.z[z_score_16S], function(x) {y <-scale(x, center = TRUE, scale = TRUE)})

#z-scoring

z_score_18S <- c("Temperature...C.", "Salinity..psu.","O2umol.l", "Fluorometer", "PO4_umol.l", "Si_umol.l", "NO3_umol.l", "MLD", "bottom_depth", "altitude", "azimuth")

metadata_18S.z <- c("Site","Station", "Region", "Fjord","Fjord2","In.Out", "Latitude", "Longitude", "Date", "Time", "Temperature...C.", "Salinity..psu.", "O2umol.l", "Fluorometer", "PO4_umol.l", "Si_umol.l", "NO3_umol.l",
                    "Glacial.influence","Bioclimatic_subzone", "MLD", "bottom_depth", "altitude", "azimuth")

meta_18S_nNA.z <- list_meta$meta_18S_m

meta_18S_nNA.z <-meta_18S_nNA.z[metadata_18S.z]%>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()



meta_18S_nNA.z[z_score_18S] <- lapply(meta_18S_nNA.z[z_score_18S], function(x) {y <-scale(x, center = TRUE, scale = TRUE)})


meta_all.z <- list_meta$meta_all
meta_all.z$Fjord

meta_all_fjord <- meta_all.z%>%group_by(Fjord)%>% 
  #tally()%>%
  summarise(mean_temp = mean(temperature...C., na.rm = TRUE), 
            mean_sal = mean(salinity..psu., na.rm = TRUE),
            mean_Fluo = mean(Fluorometer, na.rm = TRUE),
            mean_PO4 = mean(PO4_umol.l, na.rm = TRUE),
            mean_NO3 = mean(NO3_umol.l, na.rm = TRUE),
            mean_Si = mean(Si_umol.l, na.rm = TRUE),
            mean_MLD = mean(MLD, na.rm = TRUE),
            mean_bottomD = mean(bottom_depth, na.rm = TRUE)
  )
meta_all_fjord_n <- inner_join(meta_all_fjord, meta_all.z[,c(3,4)], by = "Fjord")

meta_all_fjord <- meta_all_fjord_n[!duplicated(meta_all_fjord_n),]

meta_all.z[z_score_16S] <- lapply(meta_all.z[z_score_16S], function(x) {y <-scale(x, center = TRUE, scale = TRUE)})
