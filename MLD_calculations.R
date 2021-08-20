#MLD calculations
cTD_HE533 <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/CTD/HE533_CTD.txt", sep = "\t",  encoding = "UTF-8", stringsAsFactors=FALSE)
CTD_HE492 <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/CTD/HE492_CTD.txt", sep = "\t",  encoding = "UTF-8", stringsAsFactors=FALSE)
CTD_HE431 <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/CTD/HE431_CTD.txt", sep = "\t",  encoding = "UTF-8", stringsAsFactors=FALSE)
CTD_MSM56 <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/CTD/MSM56_CTD.txt", sep = "\t",  encoding = "UTF-8", stringsAsFactors=FALSE)
CTD_MSM21 <- read.csv("C:/Users/choerstm/Documents/Studenten/arctic_picos/CTD/MSM21_CTD.txt", sep = "\t",  encoding = "UTF-8", stringsAsFactors=FALSE)
#HE533_1 <- cTD_HE533%>%group_split(Event)

HE533 <- dplyr::inner_join(meta_map_1, cTD_HE533, by = c("Latitude", "Longitude"))
HE533_MLD <- HE533%>%group_by(Event)%>%
  summarise(MLD = MLD_DK(depth = Depth.water..m., variable = Temp....C., threshold = 0.5, 
                   print.info = F, plot = F, depth.max = 1000), bottom_depth = max(Elevation..m.), 
            Latitude = mean(Latitude), Longitude = mean(Longitude))


HE492 <- dplyr::inner_join(meta_map_1, CTD_HE492, by = c("Latitude", "Longitude"))
HE492_MLD <- HE492%>%group_by(Event)%>%
  summarise(MLD = MLD_DK(depth = Depth.water..m., variable = Temp....C., threshold = 0.5, 
                         print.info = F, plot = F, depth.max = 1000), bottom_depth = max(Elevation..m.),
            Latitude = mean(Latitude), Longitude = mean(Longitude))

HE431 <- dplyr::inner_join(meta_map_1, CTD_HE431, by = c("Latitude", "Longitude"))
HE431_MLD <- HE431%>%group_by(Event)%>%
  summarise(MLD = MLD_DK(depth = Depth.water..m., variable = Temp....C., threshold = 0.5, 
                         print.info = F, plot = F, depth.max = 1000), bottom_depth = max(Elevation..m.),
            Latitude = mean(Latitude), Longitude = mean(Longitude))

#remove doubles

HE431_MLD_s <- HE431_MLD[!(HE431_MLD$Event == "HE431/28-1"|HE431_MLD$Event == "HE431/29-1" | HE431_MLD$Event == "HE431/31-1"),] 

CTD_MSM21$Station <- gsub("/3_", ".F02.", CTD_MSM21$Event)
CTD_MSM21$Station <- gsub("-.*", "", CTD_MSM21$Station)
MSM21 <- dplyr::inner_join(meta_map_1, CTD_MSM21, by = "Station")
MSM21$Temp....C. <- MSM21$Temp...C...ITS.90..Sensor.1..CTD..SEA.BI....
MSM21_MLD <- MSM21%>%group_by(Event)%>%
  summarise(MLD = MLD_DK(depth = Depth.water..m., variable = Temp....C., threshold = 0.5, 
                         print.info = F, plot = F, depth.max = 1000), bottom_depth = max(Elevation..m.),
            Latitude = mean(Latitude.x), Longitude = mean(Longitude.x))

MSM56 <- dplyr::inner_join(meta_map_1, CTD_MSM56, by = c("Latitude", "Longitude"))
MSM56$Temp....C. <- MSM56$Temp...C...Sensor.1..ITS.90..CTD..SEA.BI....
MSM56_MLD <- MSM56%>%group_by(Event)%>%
  summarise(MLD = MLD_DK(depth = Depth.water..m., variable = Temp....C., threshold = 0.5, 
                         print.info = F, plot = F, depth.max = 1000), bottom_depth = max(Elevation..m.),
            Latitude = mean(Latitude), Longitude = mean(Longitude))

all_MLD <- rbind(HE533_MLD, HE492_MLD, MSM21_MLD, MSM56_MLD, HE431_MLD_s)

meta_map_MLD <- dplyr::left_join(meta_map_1, all_MLD, by = c("Latitude", "Longitude"))
#now the calculate MLD needs to be added to metadata tables

#remove the samples with NA Event

meta_map_MLD <- meta_map_MLD[!is.na(meta_map_MLD$Event), ]


meta_16S_m <- dplyr::left_join(meta_16S_nNA, meta_map_MLD[, c("Station", "MLD", "bottom_depth")], by = "Station")
meta_18S_m <- dplyr::left_join(meta_18S_nNA, meta_map_MLD[, c("Station", "MLD", "bottom_depth")], by = "Station")


#those samples with well mixed need to be turned into numbers

meta_16S_m$MLD <- ifelse(is.na(meta_16S_m$MLD), -(meta_16S_m$bottom_depth), meta_16S_m$MLD)
meta_18S_m$MLD <- ifelse(is.na(meta_18S_m$MLD), -(meta_18S_m$bottom_depth), meta_18S_m$MLD)


rm(all_MLD, MSM21, MSM21_MLD, MSM56, MSM56_MLD, HE431, HE431_MLD, HE431_MLD_s, HE533, HE533_MLD, HE492, HE492_MLD)
rm(CTD_HE431, CTD_HE492, CTD_MSM21, CTD_MSM56, cTD_HE533)
