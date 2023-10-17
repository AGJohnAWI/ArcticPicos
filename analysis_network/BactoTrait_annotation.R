##TRAIT
#? https://doi.org/10.1371/journal.pone.0265756
#https://www.sciencedirect.com/science/article/pii/S1470160X21007123
Bacteria_traits <- all_3%>%filter(Tax_1 == "Bacteria")
Archaea_traits <- all_3%>%filter(Tax_1 == "Archaea")
prok_traits <- rbind(Bacteria_traits, Archaea_traits)
trait_16S_db <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Datasets/16S_functional/Traits/BactoTraits_databaseV2_Jun2022.csv", sep = ";")
trait_16S_db_red <- trait_16S_db[,c(1:11, 87:105)]
colnames(trait_16S_db_red) <- trait_16S_db_red[2,]
trait_16S_db_red <- trait_16S_db_red[-c(1:2),]

Tax_16S_FUNCT_annotated <- taxonomy_16S_r
Tax_16S_FUNCT_annotated$sequence <- rownames(Tax_16S_FUNCT_annotated)
#Tax_16S_FUNCT_annotated$ASV <- paste0("ASV_", 1:length(Tax_16S_FUNCT_annotated$sequence))
Tax_16S_FUNCT_annotated$Genus = Tax_16S_FUNCT_annotated$species
Tax_16S_FUNCT_annotated <- left_join(Tax_16S_FUNCT_annotated, trait_16S_db_red[,8:19], by = "Genus")

#remove rows with no trophic information


Tax_16S_FUNCT_annotated <- Tax_16S_FUNCT_annotated[complete.cases(Tax_16S_FUNCT_annotated[,12:13]), ]
