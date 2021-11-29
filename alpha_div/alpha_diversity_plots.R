#boxplot


Richness.BP <- euk_r_alpha_meta %>% ggplot(aes(x=Glacial.influence, y=Richness)) + 
  geom_boxplot() +
  geom_point()+
  labs(subtitle="Boxplot Richness")+
  theme_bw(base_size=16)

Shannon.BP <- euk_r_alpha_meta %>% ggplot(aes(x=Glacial.influence, y=Shannon)) + 
  geom_boxplot() +
  geom_point()+
  labs(subtitle="Boxplot Shannon")+
  theme_bw(base_size=16)

Simpson.BP <- euk_r_alpha_meta %>% ggplot(aes(x=Glacial.influence, y=Simpson)) + 
  geom_boxplot() +
  geom_point()+
  labs(subtitle="Boxplot Simpson")+
  theme_bw(base_size=16)

library(cowplot)

print(plot_grid(Richness.BP, Shannon.BP, Simpson.BP))