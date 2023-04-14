##prokaryotes

#replace below detection with actual number
meta_16S$PO4_umol.l <- gsub("<LOQ","0.001", meta_16S$PO4_umol.l)
meta_16S$Si_umol.l <- gsub("<LOQ","0.001", meta_16S$Si_umol.l)

meta_16S$PO4_umol.l <- as.numeric(meta_16S$PO4_umol.l)
meta_16S$Si_umol.l <- as.numeric(meta_16S$Si_umol.l)


#remove NAs and reduce to relevant metadata

metadata_16S <- c("Site", "Station", "Region", "Fjord","Fjord2","In.Out", "Latitude", "Longitude", "Date", "Time",
                  "temperature...C.", "salinity..psu.","O2umol.l", "Fluorometer",
                  "PO4_umol.l", "Si_umol.l", "NO3_umol.l",
                  "Glacial.influence", "Bioclimatic_subzone")

meta_16S_nNA <-meta_16S[metadata_16S]%>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()


##eukaryotes

#replace below detection with actual number
meta_18S_r$PO4_umol.l <- gsub("<LOQ","0.001", meta_18S_r$PO4_umol.l)
meta_18S_r$Si_umol.l <- gsub("<LOQ","0.001", meta_18S_r$Si_umol.l)

meta_18S_r$PO4_umol.l <- as.numeric(meta_18S_r$PO4_umol.l)
meta_18S_r$Si_umol.l <- as.numeric(meta_18S_r$Si_umol.l)


#remove NAs and reduce to relevant metadata

metadata_18S <- c("Site","Station", "Region", "Fjord","Fjord2", "In.Out", "Latitude", "Longitude", "Date", "Time",
                  "Temperature...C.", "Salinity..psu.","O2umol.l", "Fluorometer",
                  "PO4_umol.l", "Si_umol.l", "NO3_umol.l",
                  "Glacial.influence", "Bioclimatic_subzone")

meta_18S_nNA <-meta_18S_r[metadata_18S]%>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()

meta_map_1 <- dplyr::full_join(meta_16S_nNA[,c("Station", "Glacial.influence", "Latitude", "Longitude")], meta_18S_nNA[,c("Station", "Glacial.influence", "Latitude", "Longitude")])
meta_map <- meta_map_1

