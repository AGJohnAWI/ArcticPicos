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

#drifter_matrix_P <- drifter_matrix[] *(-1)
drifter_matrix_P[drifter_matrix_P < 0.0011] <- 0
env_dist_current <- vegan::vegdist(drifter_matrix_P, method = "jaccard")
env_dist_current_df <- dist2list(env_dist_current, tri = TRUE)#jaccard because it's rank-based
drifter_matrix1 <- drifter_matrix
drifter_matrix2 <- drifter_matrix
drifter_matrix1$col <- rownames(drifter_matrix1)
drift_long1 <- drifter_matrix1 %>% gather(row, value, -col)
drifter_matrix2$row <- rownames(drifter_matrix2)
drift_long2 <- drifter_matrix2 %>% gather(col, value, -row)

env_dist_current_df$col <- as.character(env_dist_current_df$col)
env_dist_current_df$row <- as.character(env_dist_current_df$row)

#all_16S <- inner_join(ait_dist_vector, env_dist_current_df, by = c('col', 'row'))
all_16S <- inner_join(ait_dist_vector, drift_long1, by = c('col', 'row'))
all_16S <- inner_join(all_16S, drift_long2, by = c('col', 'row'))
all_16S$value.y <- as.numeric(all_16S$value.y)
all_16S$total_particle <- all_16S$value.y + all_16S$value

all_16S$pmax <- pmax(all_16S$value.y, all_16S$value)
all_16S$pmin <- pmin(all_16S$value.y, all_16S$value)

#all_16S_oneDirectional <- all_16S%>%dplyr::filter(!pmin > 0.1)
all_16S_c <- all_16S %>%dplyr::filter(total_particle < -0.003)

meta$col <- meta$Site
meta$row <- meta$Site
all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,31)], by = "col")
all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,32)], by = "row")

all_16S_c <- unite(all_16S_c, var, c(Region2.x, Region2.y), sep = "_")
#scatterplot

#f <- ggplot(all_16S, aes(x = pmax, y = pmin, color = value.x))+
#  geom_point()

#print(f)

all_16S_c$total_particle_p <- all_16S_c$total_particle * (-1)

f <- ggplot(all_16S_c, aes(x = total_particle, y = value.x, color =var))+
  geom_point()
  #geom_point() 

print(f)

g <- ggplot(all_16S_c, aes(x = total_particle_p, y = value.x))+
  geom_smooth(aes(y=value.x, x=total_particle_p), method='lm') +
  geom_point()+
  scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
  scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))

print(g)

print(cor.test(all_16S_c$total_particle_p, all_16S_c$value.x,  method = "pearson"))
a <-lm(total_particle_p ~ value.x, data = all_16S_c)
plotNormalHistogram(a$residual)


#dcor(all_16S_c$total_particle, all_16S_c$value.x, index = 1.0)

## now with the similarity 

all_16S_si <- inner_join(ait_dist_vector, env_dist_current_df, by = c('col', 'row'))


a <- ggplot(all_16S_si, aes(y=value.x, x=value.y), color = "black")+
  geom_point(shape=1, size=1) + 
  geom_smooth(aes(y=value.x, x=value.y), method='lm') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")+
  ylab("Aitchinson distance 16S") + 
  xlab("oceanographic conncetivity") + 
  ggtitle("current influence on ASV distance")

print(a)


all_16S_t_connceted <- all_16S_si%>%dplyr::filter(!value.y == 1)
all_16S_t_connceted$value.yz <- log(1/all_16S_t_connceted$value.y)

b <- ggplot(all_16S_t_connceted, aes(y=value.x, x=value.y), color = "black")+
  geom_point(shape=1, size=1) + 
  geom_smooth(aes(y=value.x, x=value.y), method='lm') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none")+
  scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))#+ 
#scale_x_continuous(name="oceanographic conncetivity", limits=c(-1,0))

print(b)
##mantel test


print(vegan::mantel(env_dist_current, ait_dist_16S_drifter, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1))


}
#eukaryotes


distance_calc_euk <- function(euk_r, info.sort, drifter_matrix, meta){
  
  euk_drifter <- euk_r%>%dplyr::select(info.sort$Site_euk)
  ait_dist_18S_drifter <- dczm(euk_drifter, 1)
  
  ait_dist_vector <- dist2list(ait_dist_18S_drifter, tri = TRUE)
  ait_dist_vector$col <- ait_dist_vector$col
  ait_dist_vector$row <- ait_dist_vector$row
  
  #drifter_matrix_P <- drifter_matrix[] *(-1)
  drifter_matrix_P[drifter_matrix_P < 0.0011] <- 0
  env_dist_current <- vegan::vegdist(drifter_matrix_P, method = "jaccard")
  env_dist_current_df <- dist2list(env_dist_current, tri = TRUE)#jaccard because it's rank-based
  drifter_matrix1 <- drifter_matrix
  drifter_matrix2 <- drifter_matrix
  drifter_matrix1$col <- rownames(drifter_matrix1)
  drift_long1 <- drifter_matrix1 %>% gather(row, value, -col)
  drifter_matrix2$row <- rownames(drifter_matrix2)
  drift_long2 <- drifter_matrix2 %>% gather(col, value, -row)
  
  env_dist_current_df$col <- as.character(env_dist_current_df$col)
  env_dist_current_df$row <- as.character(env_dist_current_df$row)
  
  #all_16S <- inner_join(ait_dist_vector, env_dist_current_df, by = c('col', 'row'))
  all_16S <- inner_join(ait_dist_vector, drift_long1, by = c('col', 'row'))
  all_16S <- inner_join(all_16S, drift_long2, by = c('col', 'row'))
  all_16S$value.y <- as.numeric(all_16S$value.y)
  all_16S$total_particle <- all_16S$value.y + all_16S$value
  
  all_16S$pmax <- pmax(all_16S$value.y, all_16S$value)
  all_16S$pmin <- pmin(all_16S$value.y, all_16S$value)
  
  #all_16S_oneDirectional <- all_16S%>%dplyr::filter(!pmin > 0.1)
  all_16S_c <- all_16S%>%dplyr::filter(total_particle < -0.003)
  
  meta$col <- meta$Site
  meta$row <- meta$Site
  all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,30)], by = "col")
  all_16S_c <- dplyr::left_join(all_16S_c, meta[,c(4,31)], by = "row")
  
  all_16S_c <- unite(all_16S_c, var, c(Region2.x, Region2.y), sep = "_")
  #scatterplot
  
  #f <- ggplot(all_16S, aes(x = pmax, y = pmin, color = value.x))+
  #  geom_point()
  
  #print(f)
  
  all_16S_c$total_particle_p <- all_16S_c$total_particle * (-1)
  
  f <- ggplot(all_16S_c, aes(x = total_particle_p, y = value.x, color = var))+
    geom_point()+
    scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
    scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))
  
  print(f)
  
  g <- ggplot(all_16S_c, aes(x = total_particle_p, y = value.x))+
    geom_smooth(aes(y=value.x, x=total_particle_p), method='lm') +
    geom_point()+
    scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))+ 
    scale_x_continuous(name="oceanographic conncetivity", limits=c(0, 1))
  
  
  print(g)
  
  #only line
  
  
  print(cor.test(all_16S_c$total_particle_p, all_16S_c$value.x,  method = "pearson"))
  a <-lm(total_particle_p ~ value.x, data = all_16S_c)
  plotNormalHistogram(a$residual)
  
  
  #dcor(all_16S_c$total_particle_p, all_16S_c$value.x, index = 1.0)
  
  ## now with the similarity 
  
  all_16S_si <- inner_join(ait_dist_vector, env_dist_current_df, by = c('col', 'row'))
  
  
  a <- ggplot(all_16S_si, aes(y=value.x, x=value.y), color = "black")+
    geom_point(shape=1, size=1) + 
    geom_smooth(aes(y=value.x, x=value.y), method='lm') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none")+
    ylab("Aitchinson distance 18S") + 
    xlab("oceanographic conncetivity") + 
    ggtitle("current influence on ASV distance")
  
  print(a)
  
  
  all_16S_t_connceted <- all_16S_si%>%dplyr::filter(!value.y == 1)
  all_16S_t_connceted$value.yz <- log(1/all_16S_t_connceted$value.y)
  
  b <- ggplot(all_16S_t_connceted, aes(y=value.x, x=value.y), color = "black")+
    geom_point(shape=1, size=1) + 
    geom_smooth(aes(y=value.x, x=value.y), method='lm') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none")+
    scale_y_continuous(name="Aitchinson distance", limits=c(0, 150))#+ 
    #scale_x_continuous(name="oceanographic conncetivity", limits=c(-1,0))
  
  print(b)
  
  ##mantel test
  
  
  print(vegan::mantel(env_dist_current, ait_dist_18S_drifter, method="pearson", permutations=999, strata = NULL, na.rm = FALSE, parallel = 1))
  
  
  
}
