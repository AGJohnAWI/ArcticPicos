#bring the ASV table into long format

c1 <- c(">ASV_6", ">ASV_16", ">ASV_7")
c2 <- c(">ASV_205", ">ASV_12", ">ASV_21",">ASV_17", ">ASV_24", ">ASV_30", ">ASV_50", ">ASV_77", ">ASV_66", ">ASV_139", ">ASV_147")

Cluster1 <- euk_clr%>%dplyr::filter(rownames(euk_clr) %in% c1)
Cluster1$ASV <- rownames(Cluster1)

Cluster1_l <- reshape2::melt(Cluster1, id.vars="ASV", measure.vars=1:156)



#create phyloseq
taxonomy_18S_r <- taxonomy_18S%>%dplyr::filter(ASV %in% rownames(euk_clr))
rownames(taxonomy_18S_r) <- taxonomy_18S_r$ASV

OTU = otu_table(euk_clr, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy_18S_r)
TAX = tax_table(tax.matrix)
rownames(meta_18S) <- meta_18S$Site
map = sample_data(meta_18S)

phyloseq_merged = phyloseq(OTU, TAX)
all = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all

####


rank_names(all)

erie_all <- all %>%
  #tax_glom(taxrank = "genus") %>%   # agglomerate at phylum level
  
  psmelt() %>% # Melt to long format
  #filter(Abundance > 0) %>%    #be careful- here we removed everything what is relatively depleted in the sample!
  arrange(phylum) # Sort data frame alphabetically by phylum
# Set colors for plotting


cluster1 <- erie_all%>%filter(OTU %in% c1)


p <- ggplot(data=cluster1, aes(x=Fjord2, y=Abundance, fill=OTU)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
print(p)


cluster2 <- erie_all%>%filter(OTU %in% c2)


p <- ggplot(data=cluster2, aes(x=Fjord2, y=Abundance, fill=OTU)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
print(p)



