#simplyfy the communtities to indicator communitties

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")


rownames(taxonomy_18S) <- taxonomy_18S$ASV
rownames(eukaryotes) <- rownames(eukaryotes)

###create phyloseq

OTU = otu_table(eukaryotes, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy_18S)
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
  tax_glom(taxrank = "genus") %>%   # agglomerate at phylum level
  
  psmelt() %>% # Melt to long format
  #filter(Abundance > 0) %>%    #be careful- here we removed everything what is relatively depleted in the sample!
  arrange(phylum) # Sort data frame alphabetically by phylum
# Set colors for plotting

Micromonas <- erie_genus%>%dplyr::filter(genus == "Micromonas")

Micromonas$Micromonas_abundance <- Micromonas$Abundance

##atach the outcome to the 16S metadata

meta_16S_Micromonas <- left_join(meta_16S, Micromonas[c(5,37)], by = "Station")


##which of the prokaryotes correlate with Micromonas?


ASV_prok <- t(prokaryotes)
wt <- sapply(seq.int(dim(meta_16S_Micromonas)[1]), function(i) cor.test( meta_16S_Micromonas[,28], ASV_prok[,i])$p.value)
ww <- which(wt <0.001, arr.ind = TRUE)
z <- wt[wt <0.001] #extract the values of w for which cor>.6
cb=cbind(ww,z)
df=data.frame(cbind(unlist(cb), #this creates the desired data frame
                    colnames(ASV_prok)[ww]
))

Micromonas_prok <- taxonomy_16S%>%filter(rownames(taxonomy_16S) %in% df$V3)
##########

#PHAEOCYSTIS (genus)


sub <- erie_genus%>%dplyr::filter(genus == "Phaeocystis")

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
