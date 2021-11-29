
#source: https://hansenjohnson.org/post/bathymetric-maps-in-r/

#library(marmap)

library(oce)
packageVersion("oce")
library(ocedata)
packageVersion("ocedata")
#library(sp)
library(rgdal)
packageVersion("rgdal")
data("coastlineWorldFine")


# transform datapoints in polar projection
# transform the points data.frame into a spatial object

meta_map_1 <- dplyr::full_join(meta_16S[,c("Station", "Glacial.influence", "Latitude", "Longitude")], meta_18S[,c("Station", "Glacial.influence", "Latitude", "Longitude")])
meta_map <- meta_map_1

coordinates(meta_map) <- c("Longitude", "Latitude")
proj4string(meta_map) <- CRS("+proj=longlat")

# project the points with Mercator projection (maybe not the wright one?)

meta_map <- spTransform(meta_map, CRS("+proj=stere +lat_0=90"))

# plot coastline (no projection)
mapPlot(coastlineWorldFine, projection="+proj=stere +lat_0=90",
     longitudelim = c(-30, 5),
     latitudelim = c(55, 85), col='grey')
points(meta_map, col=ifelse(meta_map$Glacial.influence=="Yes", "orange2","seagreen"), pch = 18, cex = .8)


detach("package:marmap", unload = T)
detach("package:oce", unload = T)
detach("package:ocedata", unload = T)
detach("package:rgdal", unload = T)
#detach("package:sp", unload = T)

rm(coastlineWorldFine, meta_map)
