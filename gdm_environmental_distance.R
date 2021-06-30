
#gdm analysis

gdm_analysis <- function(abundance.hellinger, meta.z){
  
  require(vegan)
  
abundance_hell_t <- as.data.frame(t(abundance.hellinger))
#18S reduce the metadata table
rownames(meta.z) <- meta.z$Site
meta.z <- meta.z[c(12:15,21:25)]
df <-meta.z 
#remove NAs

meta.v <-df %>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()


#replace below detection with actual number
meta.v$PO4..µmol.l. <- gsub("<LOQ","0.001", meta.v$PO4..µmol.l.)
meta.v$Si..µmol.l. <- gsub("<LOQ","0.001", meta.v$Si..µmol.l.)
meta.v$Event <- rownames(meta.v)

meta.v$PO4..µmol.l. <- as.numeric(meta.v$PO4..µmol.l.)
meta.v$Si..µmol.l. <- as.numeric(meta.v$Si..µmol.l.)

abundance_hell_r <- abundance.hellinger%>%dplyr::select(rownames(meta.v))

abundance_hell_t_r <- abundance_hell_t%>%dplyr::filter(rownames(abundance_hell_t) %in% rownames(meta.v))
abundance_hell_t_r <- abundance_hell_t_r[,colSums(abundance_hell_t_r)>0]
abundance_hell_t_r <- as.data.frame(abundance_hell_t_r)

abundance_hell_t_r$Event <- rownames(abundance_hell_t_r)

lat_lon <- meta.v[, c(1,2)]

#all samples
spatial_dist <- geodist(lat_lon, paired = F, sequential = F, pad = FALSE, measure = "geodesic")
spatial_dist <- as.data.frame(spatial_dist)
rownames(spatial_dist) <- rownames(lat_lon)
colnames(spatial_dist) <- rownames(lat_lon)
spatial_dist_matrix <- as.matrix(spatial_dist)
spatial_dist <- melt(spatial_dist_matrix)[melt(upper.tri(spatial_dist_matrix))$value,]
names(spatial_dist) <- c("c1", "c2", "distance")

#environmental distance (euclidean) (without PP and N2 fix)
phy_dist <- vegdist(meta.v[3:4], method = "euclidean", na.rm = FALSE)
env_dist <- vegdist(meta.v[4:9], method = "euclidean", na.rm = FALSE)


#PCoA of environmental metadata - look at the eigenvalues, if they are negative then they don't make sense and shouldn't be included in further analysis
##the problem with PCoA is that it sorts the eigenvalues and you don't know which is what
### perform the PCoA
myPcoa <- cmdscale(env_dist,k=2,eig = TRUE)

x<-myPcoa$points[,1]
y<-myPcoa$points[,2]

plot(x, y, pch = 19, xlim=range(x) + c(-0.1,0.1), ylim=range(y)+c(-0.1,0.1))
text(x, y, pos = 4, labels = rownames(meta.z), cex = 0.8)

### see variation captured in the PCoA axes
barplot(myPcoa$eig)
myPcoa$eig

##rate distance
#PP_dis <- vegdist(meta.z.PP, method = "euclidean", na.rm = FALSE)

#Bray-Curtis distance of abundance table


bray_dist <- vegdist(t(abundance_hell_t_r), method = "bray", na.rm = FALSE)

#Mantel test
#mantel r gives you something like a correlation coefficient and sgígnificance by permutation

mantel(spatial_dist_matrix, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)
mantel(env_dist, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)
mantel(phy_dist, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)

#mantel(PP_dis, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)


mantel(env_dist, spatial_dist_matrix, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)


##plot distances
phy_dist_vector <- as.vector(phy_dist)
env_dist_vector <- as.vector(env_dist)
bray_dist_vector <- as.vector(bray_dist)

#PP_dis_vector <- as.vector(PP_dis)

all <- cbind(spatial_dist, phy_dist_vector, env_dist_vector, bray_dist_vector)


a <- ggplot(all, aes(y=bray_dist_vector, x=phy_dist_vector))+
  geom_point(shape=1, size=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Bray-Curtis-Distance") + 
  xlab("T/S Dissimilarity)") + 
  ggtitle("Effect of temperature/salinity dissimilarity on ASV turnover")
print(a)

b <- ggplot(all, aes(y=bray_dist_vector, x=env_dist_vector))+
  geom_point(shape=1, size=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Bray-Curtis-Distance") + 
  xlab("Environmental Dissimilarity)") + 
  ggtitle("Effect of environmental dissimilarity on ASV turnover")
print(b)

c <-ggplot(all, aes(x=distance, y=bray_dist_vector)) + 
  geom_point(shape=1, size=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Geographic distance") + 
  ylab("Bray-Curtis-Distance") + 
  ggtitle("Effect of geographic distance on ASV turnover")


print(c)

ggplot(all, aes(y=env_dist_vector, x=distance)) + geom_point(shape=1, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Environmental Dissimilarity") + 
  xlab("Geographic Distance") + 
  ggtitle("Effect of geographic distance on environmental dissimilarity")


##no glacier


meta_euk_glacier <- meta_18S_r%>%dplyr::filter(Glacial.influence == "No")
rownames(meta_euk_glacier) <- meta_euk_glacier$Site
meta_euk_glacier <- meta_euk_glacier[c(12:15,21:25)]
df <-meta_euk_glacier
#remove NAs

meta.v <-df %>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit()


#replace below detection with actual number
meta.v$PO4..µmol.l. <- gsub("<LOQ","0.001", meta.v$PO4..µmol.l.)
meta.v$Si..µmol.l. <- gsub("<LOQ","0.001", meta.v$Si..µmol.l.)

meta.v$PO4..µmol.l. <- as.numeric(meta.v$PO4..µmol.l.)
meta.v$Si..µmol.l. <- as.numeric(meta.v$Si..µmol.l.)

meta_euk_glacier_c <- meta.v

hell_euk_glacier <- euk_hellinger_df%>%dplyr::select(rownames(meta_euk_glacier_c))
hell_euk_glacier <- hell_euk_glacier[rowSums(hell_euk_glacier)>0,]

lat_lon <- meta_euk_glacier_c[, c(1,2)]

#all samples
spatial_dist <- geodist(lat_lon, paired = F, sequential = F, pad = FALSE, measure = "geodesic")
spatial_dist <- as.data.frame(spatial_dist)
rownames(spatial_dist) <- rownames(lat_lon)
colnames(spatial_dist) <- rownames(lat_lon)
spatial_dist_matrix <- as.matrix(spatial_dist)
spatial_dist <- melt(spatial_dist_matrix)[melt(upper.tri(spatial_dist_matrix))$value,]
names(spatial_dist) <- c("c1", "c2", "distance")

#environmental distance (euclidean) (without PP and N2 fix)
phy_dist <- vegdist(meta_euk_glacier_c[3:4], method = "euclidean", na.rm = FALSE)
env_dist <- vegdist(meta_euk_glacier_c[4:9], method = "euclidean", na.rm = FALSE)


#PCoA of environmental metadata - look at the eigenvalues, if they are negative then they don't make sense and shouldn't be included in further analysis
##the problem with PCoA is that it sorts the eigenvalues and you don't know which is what
### perform the PCoA
myPcoa <- cmdscale(env_dist,k=2,eig = TRUE)

x<-myPcoa$points[,1]
y<-myPcoa$points[,2]

plot(x, y, pch = 19, xlim=range(x) + c(-0.1,0.1), ylim=range(y)+c(-0.1,0.1))
text(x, y, pos = 4, labels = rownames(meta.z), cex = 0.8)

### see variation captured in the PCoA axes
barplot(myPcoa$eig)
myPcoa$eig

##rate distance
#PP_dis <- vegdist(meta.z.PP, method = "euclidean", na.rm = FALSE)

#Bray-Curtis distance of abundance table

bray_dist <- vegdist(t(hell_euk_glacier), method = "bray", na.rm = FALSE)

#Mantel test
#mantel r gives you something like a correlation coefficient and sgígnificance by permutation

mantel(spatial_dist_matrix, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)
mantel(env_dist, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)
mantel(phy_dist, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)

#mantel(PP_dis, bray_dist, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)


mantel(env_dist, spatial_dist_matrix, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1)


##plot distances
phy_dist_vector <- as.vector(phy_dist)
env_dist_vector <- as.vector(env_dist)
bray_dist_vector <- as.vector(bray_dist)

#PP_dis_vector <- as.vector(PP_dis)

all <- cbind(spatial_dist, phy_dist_vector, env_dist_vector, bray_dist_vector)


a <- ggplot(all, aes(y=bray_dist_vector, x=phy_dist_vector))+
  geom_point(shape=1, size=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Bray-Curtis-Distance") + 
  xlab("T/S Dissimilarity)") + 
  ggtitle("Effect of temperature/salinity dissimilarity on ASV turnover")
print(a)

b <- ggplot(all, aes(y=bray_dist_vector, x=env_dist_vector))+
  geom_point(shape=1, size=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Bray-Curtis-Distance") + 
  xlab("Environmental Dissimilarity)") + 
  ggtitle("Effect of environmental dissimilarity on ASV turnover")
print(b)

c <-ggplot(all, aes(x=distance, y=bray_dist_vector)) + 
  geom_point(shape=1, size=1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Geographic distance") + 
  ylab("Bray-Curtis-Distance") + 
  ggtitle("Effect of geographic distance on ASV turnover")


print(c)

ggplot(all, aes(y=env_dist_vector, x=distance)) + geom_point(shape=1, size=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Environmental Dissimilarity") + 
  xlab("Geographic Distance") + 
  ggtitle("Effect of geographic distance on environmental dissimilarity")











##################







#load additionla packages (only load after tidyverse table transformations because they interfere:))
suppressPackageStartupMessages(require(gdm))
print(citation("gdm"))

##gdm analysis
gdmTab <-formatsitepair(abundance_hell_t_r,1, dist = "bray", abundance = T, siteColumn = "Event", predData = meta.v, 
               XColumn = "Longitude", YColumn = "Latitude")

##
gdm.1 <- gdm(gdmTab, geo=T)
summary(gdm.1) #provdes an overview of the model, a shorter summary can be obtained with str
print(str(gdm.1))

#The maximum height of each spline indicates the magnitude of total biological change along that gradient and
#thereby corresponds to the relative importance of that predictor in contributing to biological turnover while
#holding all other variables constant (i.e., is a partial ecological distance). The spline’s shape indicates how
#the rate of biological change varies with position along that gradient. Thus, the splines provide insight into
#the total magnitude of biological change as a function of each gradient and where along each gradient those
#changes are most pronounced.

length(gdm.1$predictors) # get idea of number of panels

plot(gdm.1, plot.layout=c(2,2))

gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
print(plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance"))

plot(gdm.1.splineDat$x[,"Temp"], gdm.1.splineDat$y[,"Temp"], lwd=3,
     type="l", xlab="SST distance", ylab="Partial ecological distance")


print(plot(gdm.1.splineDat$x[,"mean_c.fix"], gdm.1.splineDat$y[,"mean_c.fix"], lwd=3,
     type="l", xlab="PP distance", ylab="Partial ecological distance"))

print(plot(gdm.1.splineDat$x[,"mean_N2fix"], gdm.1.splineDat$y[,"mean_N2fix"], lwd=3,
     type="l", xlab="N2 fix distance", ylab="Partial ecological distance"))

max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"Sal"])
max(gdm.1.splineDat$y[,"Temp"])
max(gdm.1.splineDat$y[,"oxygen"])
max(gdm.1.splineDat$y[,"MLD..m."])
max(gdm.1.splineDat$y[,"P_umol.l"])
max(gdm.1.splineDat$y[,"Si_umol.l"])
max(gdm.1.splineDat$y[,"chl.a_ug.l"])
max(gdm.1.splineDat$y[,"POCPN_ratio"])
max(gdm.1.splineDat$y[,"mean_c.fix"])
max(gdm.1.splineDat$y[,"mean_PB"])
max(gdm.1.splineDat$y[,"mean_N2fix"])


#unload package
detach("package:gdm", unload=TRUE)
}
