Arctic <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Submission_CommBio/Supplementary/Table_S5_Trophic_annotation.csv", sep = ",")
SubArctic <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Submission_CommBio/Supplementary/Table_S5_Subarctic.csv", sep = ",")
Temperate <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Submission_CommBio/Supplementary/Table_S5_Trophic_annotation_Temperate.csv", sep = ",")

Trophic_annotation <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Submission_CommBio/Supplementary/Arctic_fjords_taxonomy_network.csv", sep = ",")

Arctic$Tax_6 <- gsub("−", "-", Arctic$Tax_6)
SubArctic$Tax_6 <- gsub("−", "-", SubArctic$Tax_6)
Temperate$Tax_6 <- gsub("−", "-", Temperate$Tax_6)

Arctic <- dplyr::left_join(Arctic, Trophic_annotation, by = "Tax_6")
 Arctic_stats <- Arctic%>%group_by(Trophy)%>%count()
 Arctic_stats_eukprok <- Arctic%>%group_by(Tax_1)%>%count()
 
SubArctic <- dplyr::left_join(SubArctic, Trophic_annotation, by = "Tax_6")
SubArctic_stats <- SubArctic%>%group_by(Trophy)%>%count()
SubArctic_stats_eukprok <- SubArctic%>%group_by(Tax_1)%>%count()
 
Temperate <- dplyr::left_join(Temperate, Trophic_annotation, by = "Tax_6")
Temperate_stats <- Temperate%>%group_by(Trophy)%>%count()
Temperate_stats_eukprok <- Temperate%>%group_by(Tax_1)%>%count()

##only look at het components
het_anno <- c("het", "prelim het")
Arctic_het <- Arctic%>%filter(Trophy == het_anno)
Arctic_het_stats <- Arctic_het%>%group_by(Tax_1)%>%count()

Temperate_het <- Temperate%>%dplyr::filter(Trophy == het_anno)
Temperate_het_stats <- Temperate_het%>%group_by(Tax_1)%>%count()


##Temperate

# f__2(SAR86 clade)             3738 het prol
#f__5(Flavobacteriaceae)       2626 prelim het prok
#Tontoniidae_A                 1841 mixo euk
#Dino-Group-II-Clade-2         1804 het euk
#Picozoa_XXX                   1781 prelim unknown
#TAGIRI1-lineage               1613 het euk

#f__2(SAR86 clade)             41.72975 het prok
#Ascidiaceihabitans            37.97399 prelim het prok
#f__5(Flavobacteriaceae)       37.47814 het prok
#f__99(NA)                     35.99627 het < Rhizaria euk
#f__8(AEGEAN-169 marine group) 35.77350 prelim het prok
#Pseudohongiella               34.36061 prelim het prok
#f__6(Arenicellaceae)          33.93672 prelim het prok
#Planktomarina                 33.39229 prelim hhet prok
#Dino-Group-II-Clade-10-and-11 32.40841 het euk
#f__98(MAST-4)                 32.25261 het euk

 ##Subarctic

 #SAR92 prelim het
 #f9 prelim het
 #MAST-3L het
 #f4 prelim het
 #f11 prelim het
 #NS5 prelim het
 #Candidatus Actinomarina prelim het
 
 ##ARCTIC
 #Marinoscillum prelim het
 #f__1(Nitrincolaceae) prelim het
 #Dictyochophyceae_XX auto
 #Candidatus Aquiluna het
 #f7 prelim het
 #Sulfitobacter prelim het
 #f__17(Kordiimonadales) prelim het
 #Parmales_env_1 auto
#Ulvibacter het
#IS-44 prelim het (prok)
#SAR92 het
#Brevundimonas het (prok)
#Loktanella prelim het (prok)


 