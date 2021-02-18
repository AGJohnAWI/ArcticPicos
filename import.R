
setwd("C:/Users/choerstm/Documents/Studenten/arctic_picos/")

#eukaryotes

eukaryotes <- read.csv("./arctic_picos_18S_final_merge.csv", sep = ";", row.names = 1, header = T, stringsAsFactors = FALSE)
meta_18S <- read.csv("./Submission/meta_18s_arcticpicos.csv", sep = ";", header = T)
rownames(meta_18S) <- meta_18S$Site

#prepare ASV table
eukaryotes <- eukaryotes[rowSums(eukaryotes)>0,]

#prokaryotes
prokaryotes <- read.csv("./numeric_cleaned_prokaryotes_ASV1-X.csv", sep = ";", row.names = 1, stringsAsFactors = FALSE)
meta_16S <- read.csv("./Submission/meta_16s_arcticpicos.csv", sep = ";", header = T)
#meta_16S in/out: 1=in, 0=out

prokaryotes$HE533.Prok.F02.10C_S35 <- NULL #richness too high
prokaryotes$MSM21.Prok.F02.540_S38 <- NULL #richness too low
prokaryotes$HE431.Prok.F02.02_S2 <- NULL #richness too low
prokaryotes$HE431.Prok.F02.18_S15 <- NULL #richness too low

prokaryotes <- prokaryotes[rowSums(prokaryotes)>0,]
