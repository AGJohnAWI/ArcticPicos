require(metagMisc)
require(kmed)
require(energy)
require(vegan)
require(rcompanion)

#subset those where we have drifter model
distance_calc_prok <- function(prok_r, info.sort, drifter_matrix, meta){
  
prok_drifter <- prok_r%>%dplyr::select(info.sort$Site_prok)
ait_dist_16S_drifter <- dczm(prok_drifter, 1)

ait_dist_vector <- dist2list(ait_dist_16S_drifter, tri = TRUE)
ait_dist_vector$col <- ait_dist_vector$col
ait_dist_vector$row <- ait_dist_vector$row


#normalize the drifter matrix

drifter_matrix_p <- 1/drifter_matrix[]*(-1)

#take information from both directions of the connectivity matrix

drifter_matrix1 <- drifter_matrix_p
drifter_matrix2 <- drifter_matrix_p
drifter_matrix1$col <- rownames(drifter_matrix1)
drift_long1 <- drifter_matrix1 %>% gather(row, value, -col)
drifter_matrix2$row <- rownames(drifter_matrix2)
drift_long2 <- drifter_matrix2 %>% gather(col, value, -row)

all_16S <- inner_join(ait_dist_vector, drift_long1, by = c('col', 'row'))
all_16S <- inner_join(all_16S, drift_long2, by = c('col', 'row'))
all_16S$value.y <- as.numeric(all_16S$value.y)
all_16S$total_particle <- all_16S$value.y + all_16S$value

#normalize the particle concentration

all_16S$particle_normalized <- normalize(all_16S$total_particle, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")

#all_16S$pmax <- pmax(all_16S$value.y, all_16S$value)
#all_16S$pmin <- pmin(all_16S$value.y, all_16S$value)

#all_16S_oneDirectional <- all_16S%>%dplyr::filter(!pmin > 0.1)
all_16S_c <- all_16S %>%dplyr::filter(!particle_normalized == 0)

meta$col <- meta$Site
meta$row <- meta$Site
all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,31)], by = "col")
all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,32)], by = "row")

all_16S_c <- unite(all_16S_c, var, c(Region2.x, Region2.y), sep = "_")

#scatterplot

f <- ggplot(all_16S_c, aes(x = particle_normalized, y = value.x, color =var))+
  geom_point()+
  scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
  scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))
  #geom_point() 

print(f)

g <- ggplot(all_16S_c, aes(x = particle_normalized, y = value.x))+
  geom_smooth(aes(y=value.x, x=particle_normalized), method='lm') +
  geom_point()+
  scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
  scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))

print(g)

##with 0

h <- ggplot(all_16S, aes(x = particle_normalized, y = value.x))+
  geom_smooth(aes(y=value.x, x=particle_normalized), method='lm') +
  geom_point()+
  scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
  scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))

print(h)

print(cor.test(all_16S_c$particle_normalized, all_16S_c$value.x,  method = "pearson"))
a <-lm(particle_normalized ~ value.x, data = all_16S_c)
plotNormalHistogram(a$residual)


}


# eukaryotes


distance_calc_euk <- function(euk_r, info.sort, drifter_matrix, meta){
  
  euk_drifter <- euk_r%>%dplyr::select(info.sort$Site_euk)
  ait_dist_18S_drifter <- dczm(euk_drifter, 1)
  
  ait_dist_vector <- dist2list(ait_dist_18S_drifter, tri = TRUE)
  ait_dist_vector$col <- ait_dist_vector$col
  ait_dist_vector$row <- ait_dist_vector$row
  
  #normalize the drifter matrix
  
  drifter_matrix_p <- 1/drifter_matrix[]*(-1)
  
  #take information from both directions of the connectivity matrix
  
  drifter_matrix1 <- drifter_matrix_p
  drifter_matrix2 <- drifter_matrix_p
  drifter_matrix1$col <- rownames(drifter_matrix1)
  drift_long1 <- drifter_matrix1 %>% gather(row, value, -col)
  drifter_matrix2$row <- rownames(drifter_matrix2)
  drift_long2 <- drifter_matrix2 %>% gather(col, value, -row)
  
  all_16S <- inner_join(ait_dist_vector, drift_long1, by = c('col', 'row'))
  all_16S <- inner_join(all_16S, drift_long2, by = c('col', 'row'))
  all_16S$value.y <- as.numeric(all_16S$value.y)
  all_16S$total_particle <- all_16S$value.y + all_16S$value
  
  #normalize the particle concentration
  
  all_16S$particle_normalized <- normalize(all_16S$total_particle, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
  
  #all_16S$pmax <- pmax(all_16S$value.y, all_16S$value)
  #all_16S$pmin <- pmin(all_16S$value.y, all_16S$value)
  
  #all_16S_oneDirectional <- all_16S%>%dplyr::filter(!pmin > 0.1)
  all_16S_c <- all_16S %>%dplyr::filter(!particle_normalized == 0)
  
  meta$col <- meta$Site
  meta$row <- meta$Site
  all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,30)], by = "col")
  all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,31)], by = "row")
  
  all_16S_c <- unite(all_16S_c, var, c(Region2.x, Region2.y), sep = "_")
  
  #scatterplot
  
  f <- ggplot(all_16S_c, aes(x = particle_normalized, y = value.x, color =var))+
    geom_point()+
    scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
    scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))
  #geom_point() 
  
  print(f)
  
  g <- ggplot(all_16S_c, aes(x = particle_normalized, y = value.x))+
    geom_smooth(aes(y=value.x, x=particle_normalized), method='lm') +
    geom_point()+
    scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
    scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))
  
  print(g)
  
  ##with 0
  
  h <- ggplot(all_16S, aes(x = particle_normalized, y = value.x))+
    geom_smooth(aes(y=value.x, x=particle_normalized), method='lm') +
    geom_point()+
    scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
    scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))
  
  print(h)
  
  print(cor.test(all_16S_c$particle_normalized, all_16S_c$value.x,  method = "pearson"))
  a <-lm(particle_normalized ~ value.x, data = all_16S_c)
  plotNormalHistogram(a$residual)
  
  
}
