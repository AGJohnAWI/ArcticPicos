
setwd("/Users/corahoerstmann/Documents/AWI_ArcticFjords/")

#eukaryotes

eukaryotes <- read.csv("./Datasets/arctic_picos_18S_final_merge.csv", sep = ";", row.names = 1, header = T, stringsAsFactors = FALSE)
meta_18S <- read.csv("./Datasets/meta_18s_arcticpicos.csv", sep = ",", header = T)
taxonomy_18S <- read.csv("./Datasets/arcticPicos18S_taxonomy.csv", sep = ";", row.names = 1)
translation <- read.csv("./Datasets/arcticPicos18S_translation.txt", sep = "\t")
#asta file 4 taxonomy

taxonomy_18S$ASV <- row.names(taxonomy_18S)

taxonomy_18S <- dplyr::left_join(taxonomy_18S, translation, by = "ASV")

eukaryotes$sequence <- rownames(eukaryotes)
eukaryotes <- dplyr::left_join(eukaryotes, translation, by = "sequence")
rownames(eukaryotes) <- eukaryotes$ASV
eukaryotes$ASV <- NULL
eukaryotes$sequence <- NULL

rm(translation)
#euk <- t(eukaryotes)

#asv_seqs <- colnames(euk)
#asv_headers <- vector(dim(euk)[2], mode="character")
#for (i in 1:dim(euk)[2]) {
#  asv_headers[i] <- paste(">ASV", i, sep="_")
#}
# making and writing out a fasta of our final ASV seqs:
#asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, "arctic_picos_18S_ASVs.fasta")

#taxonomy filtering

taxonomy_18S <- taxonomy_18S%>%dplyr::filter(phylum!= "Metazoa") #11.926 ASVs anotated as metazoans

#prepare ASV table

eukaryotes <- eukaryotes%>%dplyr::filter(rownames(eukaryotes) %in% taxonomy_18S$ASV)
eukaryotes <- eukaryotes[rowSums(eukaryotes)>0,]

taxonomy_18S <- taxonomy_18S%>%dplyr::filter(ASV %in% rownames(eukaryotes))

rownames(taxonomy_18S) <- taxonomy_18S$ASV

#rownames(taxonomy_18S) <- taxonomy_18S$X

#eukaryotes <- eukaryotes%>%dplyr::filter(rownames(eukaryotes) %in% taxonomy_18S$X)

#prokaryotes
#prokaryotes <- read.csv("./Datasets/16S_functional/assigned_all_16S_silvanr99_v138.1_speciesR.csv", row.names = 1)

prokaryotes <- read.csv("./Datasets/numeric_cleaned_prokaryotes_ASV1-X.csv", sep = ";", row.names = 1, stringsAsFactors = FALSE)
meta_16S <- read.csv("./Datasets/meta_16s_arcticpicos.csv", sep = ",", header = T)
rep_16S <- c("HE533.Prok.F02.25B_S90",  "HE533.Prok.F02.25C_S91", "HE533.Prok.F02.26B_S94", "HE533.Prok.F02.26C_S95", "HE533.Prok.F02.27B_S98", "HE533.Prok.F02.27C_S99", "HE533.Prok.F02.28B_S102", "HE533.Prok.F02.28C_S103", "HE533.Prok.F02.21B_S78",  "HE533.Prok.F02.21C_S79", "HE533.Prok.F02.22B_S82",  "HE533.Prok.F02.22C_S83", "HE533.Prok.F02.23B_S86", "HE533.Prok.F02.23C_S87", "HE533.Prok.F02.17B_S62", "HE533.Prok.F02.17C_S63", "HE533.Prok.F02.18B_S66",  "HE533.Prok.F02.18C_S67", "HE533.Prok.F02.19B_S70",  "HE533.Prok.F02.19C_S71", "HE533.Prok.F02.20B_S74", "HE533.Prok.F02.20C_S75", "HE533.Prok.F02.11B_S38",  "HE533.Prok.F02.11C_S39", "HE533.Prok.F02.12B_S42",  "HE533.Prok.F02.12C_S43", "HE533.Prok.F02.13B_S46", "HE533.Prok.F02.13C_S47", "HE533.Prok.F02.14B_S50", "HE533.Prok.F02.14C_S51", "HE533.Prok.F02.15B_S54",  "HE533.Prok.F02.15C_S55", "HE533.Prok.F02.16B_S58",  "HE533.Prok.F02.16C_S59", "HE533.Prok.F02.10B_S22",  "HE533.Prok.F02.7B_S34",  "HE533.Prok.F02.7C_S23", "HE533.Prok.F02.8B_S26",   "HE533.Prok.F02.8C_S27", "HE533.Prok.F02.2B_S2", "HE533.Prok.F02.2C_S3", "HE533.Prok.F02.3B_S6", "HE533.Prok.F02.3C_S7", "HE533.Prok.F02.4B_S10", "HE533.Prok.F02.4C_S11", "HE533.Prok.F02.5B_S14", "HE533.Prok.F02.5C_S15", "HE533.Prok.F02.6B_S18", "HE533.Prok.F02.6C_S19", "HE533.Prok.F02.9B_S30", "HE533.Prok.F02.9C_S31")

#meta_16S <- meta_16S%>%dplyr::filter(!Site %in% rep_16S)
taxonomy_16S <- read.csv("./Datasets/silvaNGS_taxonomy.csv", sep = ";", row.names = 1)
#taxonomy_16S <- prokaryotes%>%dplyr::select(!meta_16S$Site)
#prokaryotes <- prokaryotes%>%dplyr::select(meta_16S$Site)
taxa_16S_fasta <- read.csv("./Datasets/arctic_picos---ssu---otus.txt", sep = "\t")
taxa_16S_fasta$rownames <- gsub("_", "",taxa_16S_fasta$cluster.acc)
#taxonomy_16S$true_name <- rownames(taxonomy_16S)
taxonomy_16S$rownames <-  gsub("\\..*","", rownames(taxonomy_16S))
taxonomy_16S$rownames <- gsub(">", "", taxonomy_16S$rownames)
taxonomy_16S$rownames <- gsub("_", "", taxonomy_16S$rownames)

taxonomy_16S <- left_join(taxonomy_16S[,c(2:7,11)], taxa_16S_fasta[,c(6,7,8,11)], by = "rownames")
rownames(taxonomy_16S) <- taxonomy_16S$rownames
rm(taxa_16S_fasta)
#ASVTranslation16S <- read.csv("./Submission/taxonomy/ASV_translation_16S.txt", sep = "\t")
#ASVTranslation16S$ASV <- gsub("_", "", ASVTranslation16S$ASV)

#filter out mitochondria and chloroplasts
taxonomy_16S <-taxonomy_16S%>%filter(family!= "Mitochondria") #300 ASVs removed #2186
taxonomy_16S <-taxonomy_16S%>%filter(genus!="Chloroplast") #534ASVs removed #1372
# 2092 ASVs were removed

#reduce abundance table to what we have left from the taxonomy table

prokaryotes <- prokaryotes%>%dplyr::filter(rownames(prokaryotes) %in% rownames(taxonomy_16S)) 


#meta_16S in/out: 1=in, 0=out

prokaryotes$HE533.Prok.F02.10C_S35 <- NULL #richness too high
prokaryotes$MSM21.Prok.F02.540_S38 <- NULL #richness too low
prokaryotes$HE431.Prok.F02.02_S2 <- NULL #richness too low
prokaryotes$HE431.Prok.F02.18_S15 <- NULL #richness too low
prokaryotes$HE492_88_S31 <- NULL #no metadata

prokaryotes <- prokaryotes[rowSums(prokaryotes)>0,]

taxonomy_16S <- taxonomy_16S%>%dplyr::filter(rownames(taxonomy_16S) %in% rownames(prokaryotes)) 

