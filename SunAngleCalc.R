
Sunanglecalc <- function(meta){
  require(suncalc)
  meta$date <- as.POSIXct(paste(meta$Date, meta$Time),format="%d.%m.%Y %H:%M", tz = "GMT")
  meta$lat <- meta$Latitude
  meta$lon <- meta$Longitude
  meta_sun_subset <- c("Site","date", "lat", "lon")
  
  meta_sun <- meta[,meta_sun_subset]
  
  meta_sun2 <-getSunlightPosition(data = meta_sun, keep = c("altitude", "azimuth"))
  
  meta_sun2 <- dplyr::inner_join(meta_sun, meta_sun2, by = c("date", "lat", "lon"))
  
  #one could also calculate sunrise/sunset but most stations don't have that this time of thehe year thus the angle makes more sense. 
  return(meta_sun2)
}