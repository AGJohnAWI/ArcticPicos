
meta_all_A <- meta_all%>%filter(Bioclimatic_subzone == c("high_arctic", "low_arctic"))%>%
  cbind(Bioclimatic_subzone_B = paste0("Arctic"))
meta_all_B <- meta_all%>%filter(Bioclimatic_subzone == "temperate")%>%
  cbind(Bioclimatic_subzone_B = paste0("Temperate"))
meta_all_C <- meta_all%>%filter(Bioclimatic_subzone == "subarctic")%>%
  cbind(Bioclimatic_subzone_B = paste0("Subarctic"))
meta_all <- rbind(meta_all_A, meta_all_B, meta_all_C)

BZA <- meta_all_A$Site.x
merged_CLR_A <- merged_clr_all%>%dplyr::select(meta_all_A$Station)
merged_CLR_A <- merged_clr_all%>%dplyr::select(meta_all_B$Station)

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

merged_CLR_A = amgut1.filt.cs
d <- ncol(merged_CLR_A)
n <- nrow(merged_CLR_A)
e <- d

set.seed(10010)
graph <- SpiecEasi::make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(merged_CLR_A, mar=2, distr='zinegbin', Sigma=Cor, n=n)
