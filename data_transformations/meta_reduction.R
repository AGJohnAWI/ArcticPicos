##prokaryotes

#replace below detection with actual number
meta_16S$PO4..Âµmol.l. <- gsub("<LOQ","0.001", meta_16S$PO4..Âµmol.l.)
meta_16S$Si..Âµmol.l. <- gsub("<LOQ","0.001", meta_16S$Si..Âµmol.l.)

meta_16S$PO4..Âµmol.l. <- as.numeric(meta_16S$PO4..Âµmol.l.)
meta_16S$Si..Âµmol.l. <- as.numeric(meta_16S$Si..Âµmol.l.)


#remove NAs and reduce to relevant metadata

metadata_16S <- c("Site", "Station", "Region", "Fjord","In.Out", "Sill", "Latitude", "Longitude", "Date", "Time",
                  "temperature...C.", "salinity..psu.","O.conc..Âµmol.l.", "Fluorometer",
                  "PO4..Âµmol.l.", "Si..Âµmol.l.", "NO3..Âµmol.l.",
                  "Glacial.influence", "Current_flow", "Current_score")

meta_16S_nNA <-meta_16S[metadata_16S]%>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()


##eukaryotes

#replace below detection with actual number
meta_18S_r$PO4..Âµmol.l. <- gsub("<LOQ","0.001", meta_18S_r$PO4..Âµmol.l.)
meta_18S_r$Si..Âµmol.l. <- gsub("<LOQ","0.001", meta_18S_r$Si..Âµmol.l.)

meta_18S_r$PO4..Âµmol.l. <- as.numeric(meta_18S_r$PO4..Âµmol.l.)
meta_18S_r$Si..Âµmol.l. <- as.numeric(meta_18S_r$Si..Âµmol.l.)


#remove NAs and reduce to relevant metadata

metadata_18S <- c("Site","Station", "Region", "Fjord","In.Out", "Sill", "Latitude", "Longitude", "Date", "Time",
                  "Temperature...C.", "Salinity..psu.","O.conc..Âµmol.l.", "Fluorometer",
                  "PO4..Âµmol.l.", "Si..Âµmol.l.", "NO3..Âµmol.l.",
                  "Glacial.influence", "Current_flow", "Current_score")

meta_18S_nNA <-meta_18S_r[metadata_18S]%>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()

meta_map_1 <- dplyr::full_join(meta_16S[,c("Station", "Glacial.influence", "Latitude", "Longitude")], meta_18S[,c("Station", "Glacial.influence", "Latitude", "Longitude")])
meta_map <- meta_map_1