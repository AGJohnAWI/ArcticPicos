
setwd("C:/Users/choerstm/Documents/Studenten/arctic_picos/")

#eukaryotes

eukaryotes <- read.csv("./arctic_picos_18S_final_merge.csv", sep = ";", row.names = 1, header = T, stringsAsFactors = FALSE)
meta_18S <- read.csv("./Submission/meta_18s_arcticpicos.csv", sep = ";", header = T)
taxonomy_18S <- read.csv("./Submission/arcticPicos18S_taxonomy.csv", sep = ";", row.names = 1)
translation <- read.csv("./arcticPicos18S_translation.txt", sep = "\t")
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

#rownames(taxonomy_18S) <- taxonomy_18S$X

#eukaryotes <- eukaryotes%>%dplyr::filter(rownames(eukaryotes) %in% taxonomy_18S$X)

#prokaryotes
prokaryotes <- read.csv("./numeric_cleaned_prokaryotes_ASV1-X.csv", sep = ";", row.names = 1, stringsAsFactors = FALSE)
meta_16S <- read.csv("./Submission/meta_16s_arcticpicos.csv", sep = ";", header = T)
taxonomy_16S <- read.csv("./Submission/silvaNGS_taxonomy.csv", sep = ";", row.names = 1)

prokaryotes <- prokaryotes%>%dplyr::select(meta_16S$Site)

taxonomy_16S$rownames <-  gsub("\\..*","", rownames(taxonomy_16S))
taxonomy_16S$rownames <- gsub(">", "", taxonomy_16S$rownames)
taxonomy_16S$rownames <- gsub("_", "", taxonomy_16S$rownames)

rownames(taxonomy_16S) <- taxonomy_16S$rownames

#ASVTranslation16S <- read.csv("./Submission/taxonomy/ASV_translation_16S.txt", sep = "\t")
#ASVTranslation16S$ASV <- gsub("_", "", ASVTranslation16S$ASV)


#filter out mitochondria and chloroplasts
taxonomy_16S <-taxonomy_16S%>%filter(family!= "Mitochondria") #300 ASVs removed
taxonomy_16S <-taxonomy_16S%>%filter(genus!="Chloroplast") #534ASVs removed
# 2092 ASVs were removed

#reduce abundance table to what we have left from the taxonomy table

prokaryotes <- prokaryotes%>%dplyr::filter(rownames(prokaryotes) %in% taxonomy_16S$rownames) 


#meta_16S in/out: 1=in, 0=out

#prokaryotes$HE533.Prok.F02.10C_S35 <- NULL #richness too high
#prokaryotes$MSM21.Prok.F02.540_S38 <- NULL #richness too low
#prokaryotes$HE431.Prok.F02.02_S2 <- NULL #richness too low
#prokaryotes$HE431.Prok.F02.18_S15 <- NULL #richness too low

prokaryotes <- prokaryotes[rowSums(prokaryotes)>0,]

taxonomy_16S <- taxonomy_16S%>%dplyr::filter(rownames %in% rownames(prokaryotes)) 
