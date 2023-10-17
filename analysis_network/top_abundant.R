## top abundant picos


#                                   BiocManager::repositories()))

##Import the Trophic annotation
Trophic_annotation <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Submission_CommBio/Arctic_fjords_taxonomy_annotated.csv", sep = ",")
Trophic_annotation_auto <- Trophic_annotation%>%filter(str_detect(trophy, "auto"))%>%
  cbind(trophy_simple = paste0("auto"))
Trophic_annotation_het <- Trophic_annotation%>%filter(str_detect(trophy, "het"))%>%
  cbind(trophy_simple = paste0("het"))
Trophic_annotation_mixo <- Trophic_annotation%>%filter(str_detect(trophy, "mixo"))%>%
  cbind(trophy_simple = paste0("mixo"))
Trophic_annotation_unknown <- Trophic_annotation%>%filter(str_detect(trophy, "unknown"))%>%
  cbind(trophy_simple = paste0("unknown"))
Trophic_annotation <- rbind(Trophic_annotation_auto, Trophic_annotation_het, Trophic_annotation_mixo,
                            Trophic_annotation_unknown)

rm(Trophic_annotation_auto, Trophic_annotation_het, Trophic_annotation_mixo, Trophic_annotation_unknown)
##

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")
require(metagMisc)
packageVersion("metagMisc")
#require("NetCoMi")
#packageVersion("NetCoMi")
#require("igraph")
#packageVersion("igraph")
#require("qgraph")
#packageVersion("qgraph")
#require(limma)
#install.packages("igraph")
#install.packages("qgraph")
###create phyloseq


#devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic"); packageVersion("fantaxtic")

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")

##

OTU = otu_table(merged_clr_all, taxa_are_rows = TRUE)

#merge high and low arctic into one bioclimatic subzone

meta_all_A <- meta_all%>%filter(Bioclimatic_subzone == c("high_arctic", "low_arctic"))%>%
  cbind(Bioclimatic_subzone_B = paste0("Arctic"))
meta_all_B <- meta_all%>%filter(Bioclimatic_subzone == "temperate")%>%
  cbind(Bioclimatic_subzone_B = paste0("Temperate"))
meta_all_C <- meta_all%>%filter(Bioclimatic_subzone == "subarctic")%>%
  cbind(Bioclimatic_subzone_B = paste0("Subarctic"))
meta_all <- rbind(meta_all_A, meta_all_B, meta_all_C)

tax.matrix<- as.matrix(taxonomy_merged_clr)
TAX = tax_table(tax.matrix)
rownames(meta_all) <- meta_all$Station
map = sample_data(meta_all)

phyloseq_merged = phyloseq(OTU, TAX)
all = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all  


###

top <- fantaxtic::top_taxa(all, 
                tax_level = "Tax_5", 
                n_taxa = 10,
                grouping = "Bioclimatic_subzone_B")

top_taxa <- top$top_taxa

top_taxa <- dplyr::left_join(top_taxa, Trophic_annotation, by = "taxid")

print(top_taxa)


## Figure with % clr-transformed data

Prok_clr_norm <- apply(Prok_r_clr_rename, 2, function(x) x - min(x[x < 0]))
Euk_clr_norm <- apply(Euk_r_clr_rename, 2, function(x) x - min(x[x < 0]))

merged_clr_all_norm <- rbind(Prok_clr_norm, Euk_clr_norm)

merged_clr_all_norm_P <- apply(merged_clr_all_norm, 2, function(x) {x/sum(x)})

merged_clr_all_norm_P_a <- as.matrix(merged_clr_all_norm_P)

## Turn into other phyloseq object

OTU = otu_table(merged_clr_all_norm_P_a, taxa_are_rows = TRUE)
rownames(Trophic_annotation) <- Trophic_annotation$taxid
tax.matrix<- as.matrix(Trophic_annotation)
TAX = tax_table(tax.matrix)
map = sample_data(meta_all)

phyloseq_merged = phyloseq::phyloseq(OTU, TAX)
all_2 = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all_2 

##rearrange the stack of the barplot


all_3 <- psmelt(all_2)
all_3$Tax_1 <- gsub(" ", "", all_3$Tax_1)
all_3$function_Tax <- paste(all_3$Tax_1, all_3$trophy_simple)

all_3$function_Tax2 <- factor(all_3$function_Tax, levels = c("Archaea auto","Bacteria auto", "Eukaryota auto","Bacteria mixo", "Eukaryota mixo","Bacteria het","Archaea het","Eukaryota het", "Archaea unknown", "Bacteria unknown","Eukaryota unknown"))
all3_nNA <- all_3%>%filter(!Abundance == 0)
##

d <- ggplot(all3_nNA, aes(x = Sample, y = Abundance, color = function_Tax2))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values= c("#4ca583","#60E27B","#99BC08","#db6adb", "#9324EA", "#A3610A", "#824B04","#eaeaea","#818281", "#4e4f4f"))+
  scale_color_manual(values= c("#4ca583","#60E27B", "#99BC08","#db6adb","#9324EA", "#A3610A", "#824B04","#eaeaea","#818281", "#4e4f4f"))+
  theme(axis.text.x = element_text(size=4, vjust = 0.3, angle = 90))+
  facet_wrap(~Bioclimatic_subzone_B,
             scales = "free_x")
print(d)

#STATS
#bring into wide format according to the function_tax

all_4 <- spread(all_3, function_Tax, Abundance)

pairwise.wilcox.test(all_4$`Bacteria auto`, all_4$Bioclimatic_subzone_B,
                     p.adjust.method = "BH")

sTATS <- all_3%>%
        group_by(Sample,function_Tax)%>%
  summarise(sum=sum(Abundance))
sTATS$Station <- sTATS$Sample
sTATS <- left_join(sTATS, meta_all[,c("Station", "Bioclimatic_subzone_B")], by="Station")

sTATS2 <- spread(sTATS, function_Tax, sum)

print(ggplot(data = sTATS, mapping = aes(x = Bioclimatic_subzone_B, y = sum, color = function_Tax)) +
  geom_boxplot())


#STATS
print(summary(aov(`Bacteria auto` ~ Bioclimatic_subzone_B, data = sTATS2)))
print(summary(aov(`Eukaryota auto` ~ Bioclimatic_subzone_B, data = sTATS2)))
print(summary(aov(`Eukaryota mixo` ~ Bioclimatic_subzone_B, data = sTATS2)))
print(summary(aov(`Bacteria het` ~ Bioclimatic_subzone_B, data = sTATS2)))
print(summary(aov(`Archaea het` ~ Bioclimatic_subzone_B, data = sTATS2)))
print(summary(aov(`Bacteria auto` ~ Bioclimatic_subzone_B, data = sTATS2)))

#No of annotated ASVs
trophy_count <- as.tibble(Trophic_annotation)

#trophy_count <- as.factor(trophy_count$trophy_simple)
trophy_count%>%group_by(Tax_1)%>%dplyr::count(trophy_simple)
#sum(Trophic_annotation$trophy_simple == "auto") 
##t-test

print(pairwise.t.test(sTATS2$`Bacteria auto`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(sTATS2$`Archaea auto`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(sTATS2$`Bacteria mixo`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(sTATS2$`Eukaryota auto`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(sTATS2$`Eukaryota mixo`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(sTATS2$`Bacteria het`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
#print(pairwise.t.test(sTATS2$`Archaea het`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))
print(pairwise.t.test(sTATS2$`Eukaryota het`, sTATS2$Bioclimatic_subzone_B, var.equal = FALSE, p.adjust.method='bonferroni'))


f <- ggplot(all_3, aes(x = Sample, y = Abundance, color = Tax_1))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values= c("#012347", "#003F7D", "#FF8F01"))+
  scale_color_manual(values= c("#012347", "#003F7D", "#FF8F01"))+
  theme(axis.text.x = element_text(size=4, vjust = 0.3, angle = 90))+
  facet_wrap(~Bioclimatic_subzone_B,
             scales = "free_x")
print(f)


rm(top, all, top_taxa, tax.matrix, meta_all_A, meta_all_B, meta_all_C)
