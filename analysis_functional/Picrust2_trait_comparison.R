library(ggpmisc)
KO_Qiime <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Datasets/16S_functional/KO_pathabun_exported/feature-table.biom.tsv", skip = 1, sep = "\t")

#energy pathway KO

#see https://www.genome.jp/kegg/pathway.html
#from raes et al. 2021

KO_number_pathway <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Datasets/16S_functional/KO_pathabun_exported/KO_numbers_with_pathways_names.csv")
KO_energy <- KO_number_pathway%>%filter(TAX2 == " Energy metabolism")
KO_energy_Cfix <- KO_energy%>%filter(str_detect(TAX3a, "Carbon fixation pathways in prokaryotes"))
#00195 Photosyhtesis

#Energy_pathway_C <- c("K00190","K00195", "K00196", "K00710", "K00720")

KO_C_energy <- KO_Qiime%>%filter(X.OTU.ID %in% KO_energy_Cfix$KO)
rownames(KO_C_energy) <- KO_C_energy$X.OTU.ID
KO_C_energy$X.OTU.ID <- NULL
KO_C_energy_t <- as.data.frame(t(KO_C_energy))
KO_C_energy_t$sum_C_KOs <- rowSums(KO_C_energy_t)
KO_C_energy_t$Site <- rownames(KO_C_energy_t)


Tax_traits <- as.data.frame(Trophic_annotation)
Tax_traits_Bac <- Tax_traits%>%filter(Tax_1 == " Bacteria")
Tax_traits_Bac_auto <- Tax_traits_Bac%>%filter(trophy_simple == "auto")
Abundance_traits_all <- as.data.frame(merged_clr_all_norm_P_a)

#the clr-transfomred abundance need to be shifted to positive values to calculate sums
prok_clr_norm2 <- apply(prok_clr, 2, function(x) x - min(x[x < 0]))
prok_clr_norm2 <- as.data.frame(prok_clr_norm2)
Bac_auto_abun_norm <- prok_clr_norm2%>%filter(rownames(prok_clr_norm2) %in% rownames(Tax_traits_Bac_auto))

Bac_auto_abun_norm_t <- as.data.frame(t(Bac_auto_abun_norm))
Bac_auto_abun_norm_t$SUM_abundance <- rowSums(Bac_auto_abun_norm_t)
Bac_auto_abun_norm_t$Site <- rownames(Bac_auto_abun_norm_t)
#KO_C_energy_t_R <- KO_C_energy_t%>%filter(Site %in% rownames(Bac_auto_abun_norm_t))

KO_Trait_merged <- left_join(KO_C_energy_t, Bac_auto_abun_norm_t, by = "Site")

ggplot(KO_Trait_merged,aes(x=sum_C_KOs, y = SUM_abundance))+
  geom_point()+
  stat_poly_line()+
ggpmisc::stat_poly_eq(use_label(c("eq", "R2", "p", "n")))

