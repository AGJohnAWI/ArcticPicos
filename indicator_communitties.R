#simplyfy the communtities to indicator communitties

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")


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
##

rank_names(all)

# Agglomerate to genus level to extract Micromonas

erie_genus <- all %>%
  #tax_glom(taxrank = "genus") %>%   # agglomerate at phylum level
  
  psmelt() %>% # Melt to long format
  #filter(Abundance > 0) %>%    #be careful- here we removed everything what is relatively depleted in the sample!
  arrange(phylum) # Sort data frame alphabetically by phylum
# Set colors for plotting

Top10 <- erie_genus%>%dplyr::filter(OTU %in% indicator_10_euks_NOglacier$ASV) #change here

p <- ggplot(data=Top10, aes(x=Sample, y=Abundance, fill=Glacial.influence)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
print(p)

print(indicator_10_euks_NOglacier$ASV) #change here, too

Top10$Abundance_indi <- Top10$Abundance

#Top10_1 <- Top10%>%dplyr::filter(ASV == ">ASV_382")
Top10_2 <- Top10%>%dplyr::filter(ASV == ">ASV_4")

##atach the outcome to the 16S metadata

#meta_top10 <- inner_join(meta_16S, Top10_1[c(5,41)], by = "Station")
meta_top10 <- inner_join(meta_16S, Top10_2[c(5,41)], by = "Station")


##which of the prokaryotes correlate with Micromonas?

proks <- prokaryotes%>%dplyr::select(meta_top10$Site)
ASV_prok <- t(proks)

wt <- sapply(seq.int(dim(meta_top10)[1]), function(i) cor.test( meta_top10[,28], ASV_prok[,i])$p.value)
ww <- which(wt <0.001, arr.ind = TRUE)
z <- wt[wt <0.001] #extract the values of w for which cor>.6
cb=cbind(ww,z)
df=data.frame(cbind(unlist(cb), #this creates the desired data frame
                    colnames(ASV_prok)[ww]
))

Top10_tax <- taxonomy_16S%>%dplyr::filter(rownames(taxonomy_16S) %in% df$V3)

#NON-GLACIER

#1 >ASV_101 = 59 ASVs
#2 >ASV_69 = 54 ASVs
#3 >ASV_375 = 42 ASVs
#4 >ASV_130 = 43 ASVs
#5 >ASV_1620 = 27 ASVs
#6 >ASV_35 = 56 ASVs
#7 >ASV_135 = 37 ASVs
#8 >ASV_216 = 28ASVs
#9 >ASV_1770 = 36ASVs
#10 >>ASV_4 = 1ASV

#GLACIER
#1 >ASV_382 = 44 ASVs
#2 >ASV_819 = 44
#3 >ASV_1277 = 31ASVs
#4 >ASV_2 = 49 ASVs
#5 >ASV_229 = 22 ASVs
#6 >ASV_536 = 41
#7 >ASV_330 = 61
#8 >ASV_248 = 20
#9 >ASV_760 = 27
#10>ASV_226 = 39 
##########

#PHAEOCYSTIS (genus)


sub <- erie_genus%>%dplyr::filter(genus == "Phaeocystis")

p <- ggplot(data=sub, aes(x=Sample, y=Abundance, fill=Glacial.influence)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
print(p)



sub$Phaeocystis_abundance <- sub$Abundance

##atach the outcome to the 16S metadata

meta_16S_euks <- left_join(meta_16S, sub[c(5,37)], by = "Station")


##which of the prokaryotes correlate with Phaeocystis

ASV_prok <- t(prokaryotes)
wt <- sapply(seq.int(dim(meta_16S_Micromonas)[1]), function(i) cor.test( meta_16S_euks[,28], ASV_prok[,i])$p.value)
ww <- which(wt <0.001, arr.ind = TRUE)
z <- wt[wt <0.001] #extract the values of w for which cor>.6
cb=cbind(ww,z)
df=data.frame(cbind(unlist(cb), #this creates the desired data frame
                    colnames(ASV_prok)[ww]
))

Phaeocystis_prok <- taxonomy_16S%>%filter(rownames(taxonomy_16S) %in% df$V3)




sub <- erie_all%>%dplyr::filter(species == "Katablepharis_japonica")


p <- ggplot(data=sub, aes(x=Sample, y=Abundance, fill=Glacial.influence)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()
print(p)
