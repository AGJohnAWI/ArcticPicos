
#soure: https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap

require(phyloseq)

#attach the 18 Data to 16S metadata
indicator_glacier <- indicator_euks_GLACIER$ASV
ASVTable_euk_glacier_indicator <- Euk_r_clr_rename[indicator_glacier,]
taxonomy_indicator_glacier <- indicator_euks_GLACIER[,c(6:16)]
rownames(taxonomy_indicator_glacier) <- taxonomy_indicator_glacier$ASV
taxonomy_indicator_glacier$ASV <- NULL

###

#create phyloseq to melt on species level or similar


OTU = otu_table(ASVTable_euk_glacier_indicator, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy_indicator_glacier)
TAX = tax_table(tax.matrix)
rownames(meta_all) <- meta_all$Station
map = sample_data(meta_all)

phyloseq_merged = phyloseq(OTU, TAX)
all = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all

colnames(tax_table(all))

taxtab <- all@tax_table@.Data
# Find undefined taxa (in this data set, unknowns occur only up to Rank5)
miss_f <- which(taxtab[, "genus"] == "")
miss_g <- which(taxtab[, "family"] == "")

# indicate unspecified genera
taxtab[miss_f, "genus"] <- paste0("f__")
taxtab[miss_g, "family"] <- paste0("g__")

# Find duplicate genera
#dupl_g <- which(duplicated(taxtab[, "genus"]) |
 #                 duplicated(taxtab[, "genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "order"], ")")
  } else if(i %in% miss_f){
    taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "family"], ")")
  }
  
  # Family names are added to duplicate genera
  #if(i %in% dupl_g){
  #  taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "genus"], ")")
  }


all@tax_table@.Data <- taxtab

erie_family <- all %>%
  tax_glom(taxrank = "genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t", "f__()")) %>%   # agglomerate at phylum level
  
  psmelt() %>% # Melt to long format
  #filter(Abundance > 0) %>%    #be careful- here we removed everything what is relatively depleted in the sample!
  arrange(phylum)

erie_genus <- erie_family[, c(3,5,38)]

erie_spread <- spread(erie_genus, genus, Abundance)

rownames(erie_spread) <- erie_spread$Station
erie_spread$Station <- NULL

#################
#################

#attach the 18 Data to 16S metadata
indicator_NOglacier <- indicator_euks_NO_GLACIER$ASV
ASVTable_euk_NOglacier_indicator <- Euk_r_clr_rename[indicator_NOglacier,]
taxonomy_indicator_NOglacier <- indicator_euks_NO_GLACIER[,c(6:16)]
rownames(taxonomy_indicator_NOglacier) <- taxonomy_indicator_NOglacier$ASV
taxonomy_indicator_NOglacier$ASV <- NULL

###

#create phyloseq to melt on species level or similar


OTU = otu_table(ASVTable_euk_NOglacier_indicator, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy_indicator_NOglacier)
TAX = tax_table(tax.matrix)
rownames(meta_all) <- meta_all$Station
map = sample_data(meta_all)

phyloseq_merged = phyloseq(OTU, TAX)
all_NG = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all_NG

colnames(tax_table(all_NG))

taxtab <- all_NG@tax_table@.Data
# Find undefined taxa (in this data set, unknowns occur only up to Rank5)
miss_f <- which(taxtab[, "genus"] == "")
miss_g <- which(taxtab[, "family"] == "")

# indicate unspecified genera
taxtab[miss_f, "genus"] <- paste0("f__")
taxtab[miss_g, "family"] <- paste0("g__")

# Find duplicate genera
#dupl_g <- which(duplicated(taxtab[, "genus"]) |
#                 duplicated(taxtab[, "genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "order"], ")")
  } else if(i %in% miss_f){
    taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "family"], ")")
  }
  
  # Family names are added to duplicate genera
  #if(i %in% dupl_g){
  #  taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "genus"], ")")
}


all_NG@tax_table@.Data <- taxtab

erie_NG <- all_NG %>%
  tax_glom(taxrank = "genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t", "f__()")) %>%   # agglomerate at phylum level
  
  psmelt() %>% # Melt to long format
  #filter(Abundance > 0) %>%    #be careful- here we removed everything what is relatively depleted in the sample!
  arrange(phylum)

erie_NG_genus <- erie_NG[, c(3,5,38)]

erie_NG_spread <- spread(erie_NG_genus, genus, Abundance)

rownames(erie_NG_spread) <- erie_NG_spread$Station
erie_NG_spread$Station <- NULL


#####################
###prokaryotes phyloseq

OTU = otu_table(Prok_r_clr_rename, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy_16S_r)
TAX = tax_table(tax.matrix)
rownames(meta_all) <- meta_all$Station
map = sample_data(meta_all)

phyloseq_merged = phyloseq(OTU, TAX)
all_p = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all_p

colnames(tax_table(all_p))

taxtab <- all_p@tax_table@.Data
# Find undefined taxa (in this data set, unknowns occur only up to Rank5)

missing <- c("", " ", "uncultured", "\t")

miss_f <- which(taxtab[, "species"] %in% missing)
miss_g <- which(taxtab[, "family"] %in% missing)

# indicate unspecified genera
taxtab[miss_f, "species"] <- paste0("f__")
taxtab[miss_g, "family"] <- paste0("g__")

# Find duplicate genera
#dupl_g <- which(duplicated(taxtab[, "species"]) |
#                 duplicated(taxtab[, "species"], fromLast = TRUE))

for(i in seq_along(taxtab)){
  # The next higher non-missing rank is assigned to unspecified genera
  if(i %in% miss_f && i %in% miss_g){
    taxtab[i, "species"] <- paste0(taxtab[i, "species"], "(", taxtab[i, "genus"], ")")
  } else if(i %in% miss_f){
    taxtab[i, "species"] <- paste0(taxtab[i, "species"], "(", taxtab[i, "family"], ")")
  }
  
  # Family names are added to duplicate genera
  #if(i %in% dupl_g){
  #  taxtab[i, "species"] <- paste0(taxtab[i, "species"], "(", taxtab[i, "species"], ")")
}


all_p@tax_table@.Data <- taxtab

erie_family_p <- all_p %>%
  tax_glom(taxrank = "species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t","f__()", "f__( )")) %>%   # agglomerate at phylum level
  
  psmelt() %>% # Melt to long format
  #filter(Abundance > 0) %>%    #be careful- here we removed everything what is relatively depleted in the sample!
  arrange(family)

erie_species <- erie_family_p[, c(3,5,37)]

erie_spread_prokaryotes <- spread(erie_species, species, Abundance)

rownames(erie_spread_prokaryotes) <- erie_spread_prokaryotes$Station
erie_spread_prokaryotes$Station <- NULL


##correlation with prokaryotes

cor_glacier_genus <- cor(x = erie_spread_prokaryotes, y = erie_spread)

cor_NOglacier_genus <- cor(x = erie_spread_prokaryotes, y = erie_NG_spread)

dendro_heat <- function(correlation_matrix){


# Obtain the dendrogram

library(plyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)

sample_names <- colnames(correlation_matrix)

dend <- as.dendrogram(hclust(dist(correlation_matrix)))
dend_data <- dendro_data(dend)


# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names) %>%
  dplyr::mutate(x_center = (1:n()), 
         width = 1)

# Neglecting the gap parameters
heatmap_data <- correlation_matrix %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

# Heatmap plot
plt_hmap <- ggplot(heatmap_data, 
                   aes(x = x_center, y = y_center, fill = expr, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_fill_gradient2("expr", high = "darkblue", low = "darkred") +
  scale_x_continuous(breaks = sample_pos_table$x_center, 
                     labels = sample_pos_table$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Sample", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

library(cowplot)
plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(1, 1))



}

