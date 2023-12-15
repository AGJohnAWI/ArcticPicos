#combined 16S and 18S

#meta_all, merged_clr_all, taxonomy_merged_clr

#devtools::install_github("stefpeschel/NetCoMi", 
#                         dependencies = c("Depends", "Imports", "LinkingTo"),
#                         repos = c("https://cloud.r-project.org/",
#                                   BiocManager::repositories()))

##Import the Trophic annotation
Trophic_annotation <- read.csv("/Users/corahoerstmann/Documents/AWI_ArcticFjords/Submission_CommBio/Supplementary/Arctic_fjords_taxonomy_network.csv", sep = ",")
##

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")
require(metagMisc)
packageVersion("metagMisc")
require("NetCoMi")
packageVersion("NetCoMi")
require("igraph")
packageVersion("igraph")
require("qgraph")
packageVersion("qgraph")
require(limma)
#install.packages("igraph")
#install.packages("qgraph")
###create phyloseq

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


  # Agglomerate to order level
  
  all_T6 <- phyloseq::tax_glom(all, taxrank = "Tax_6")
  taxtab <- all_T6@tax_table@.Data
  
  miss_f <- which(taxtab[, "Tax_6"] == "f__")
  miss_g <- which(taxtab[, "Tax_5"] == "g__")
  
  # Number unspecified genera
  taxtab[miss_f, "Tax_6"] <- paste0("f__", 1:length(miss_f))
  taxtab[miss_g, "Tax_5"] <- paste0("g__", 1:length(miss_g))
  
  # Find duplicate genera
  dupl_g <- which(duplicated(taxtab[, "Tax_6"]) |
                    duplicated(taxtab[, "Tax_6"], fromLast = TRUE))
  
  for(i in seq_along(taxtab)){
    # The next higher non-missing rank is assigned to unspecified genera
    if(i %in% miss_f && i %in% miss_g){
      taxtab[i, "Tax_6"] <- paste0(taxtab[i, "Tax_6"], "(", taxtab[i, "Tax_4"], ")")
    } else if(i %in% miss_f){
      taxtab[i, "Tax_6"] <- paste0(taxtab[i, "Tax_6"], "(", taxtab[i, "Tax_5"], ")")
    }
    
    # Family names are added to duplicate genera
    if(i %in% dupl_g){
      taxtab[i, "Tax_6"] <- paste0(taxtab[i, "Tax_6"], "(", taxtab[i, "Tax_6"], ")")
    }
  }
  
  all_T6@tax_table@.Data <- taxtab
  
  ###
  
  own_split <- metagMisc::phyloseq_sep_variable(all_T6, 
                                                "Bioclimatic_subzone_B") #Before was Glacial.influence
  
  # Find undefined taxa (in this data set, unknowns occur only up to Rank5)
  
  rownames(own_split$Arctic@otu_table@.Data) <- own_split$Arctic@tax_table@.Data[, "Tax_6"]
  rownames(own_split$Subarctic@otu_table@.Data) <- own_split$Subarctic@tax_table@.Data[, "Tax_6"]
  rownames(own_split$Temperate@otu_table@.Data) <- own_split$Temperate@tax_table@.Data[, "Tax_6"]
  
###TEMPERATE
  
  
  net_season <- netConstruct(data = own_split$Temperate, 
                             measure = "spieceasi", ######
                             normMethod = "none", 
                             zeroMethod = "multRepl",
                             sparsMethod = "threshold", 
                             thresh = 0.3, 
                             dissFunc = "signed",
                             verbose = 3,
                             seed = 123456)
  props_season <- netAnalyze(net_season, 
                             centrLCC = FALSE,
                             avDissIgnoreInf = TRUE,
                             sPathNorm = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = c("degree", "between", "closeness"),
                             hubQuant = 0.9,
                             lnormFit = TRUE,
                             normDeg = FALSE,
                             normBetw = FALSE,
                             normClose = FALSE,
                             normEigen = FALSE)
  #aa <-net_season[["simMat1"]]
  #max(aa[aa != max(aa)])
  summary(props_season)
  
  # Get phyla names from the taxonomic table created before
  phyla <- as.factor(taxtab[, "Tax_1"])
  names(phyla) <- taxtab[, "Tax_6"]
  #table(phyla)
  
  # Define phylum colors
  phylcol <- c( "#ae1816", "#d49dc7","#23378d")
  
  
  plot(props_season, 
       sameLayout = FALSE, 
       nodeColor = "feature", #alternative here: "cluster"
       featVecCol = phyla, 
       colorVec =  phylcol,
       nodeSize = "clr",
       shortenLabels = "none",
       rmSingles = TRUE,
       labelScale = FALSE,
       posCol = "darkgreen", 
       negCol = "darkgrey",
       cexNodes = 0.7, 
       cexLabels = 0.5,
       cexHubLabels = 0.7,
       cexTitle = 1.7,
       #groupNames = c("Arctic", "Temperate"),
       hubBorderCol  = "gray40")
  
  
  
  Temperate_network_species <- as.data.frame(sort(colSums(net_season$normCounts1), decreasing = TRUE))
  Temperate_network_species$Tax_6 <- rownames(Temperate_network_species)
  
  #Stats on Trophic modes
  
  
  Temperate_trophic <- dplyr::left_join(Temperate_network_species, Trophic_annotation, by = "Tax_6")
  Temperate_stats <- Temperate_trophic%>%group_by(Tax_1, Trophy)%>%count()
  Temperate_stats_eukprok <- Temperate_trophic%>%group_by(Tax_1)%>%count()
  
  print(Temperate_stats)
  print(Temperate_stats_eukprok)
  
  rm(Temperate_trophic, Temperate_network_species, Temperate_stats, Temperate_stats_eukprok)
  
  #####
  
  print(summary(props_season, numbNodes = 5L))
  
  

  
###SUBARCTIC
  net_season <- netConstruct(data = own_split$Subarctic, 
                             measure = "spieceasi", ######
                             normMethod = "none", 
                             zeroMethod = "multRepl",
                             sparsMethod = "threshold", 
                             thresh = 0.3, 
                             dissFunc = "signed",
                             verbose = 3,
                             seed = 123456)
  props_season <- netAnalyze(net_season, 
                             centrLCC = FALSE,
                             avDissIgnoreInf = TRUE,
                             sPathNorm = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = c("degree", "between", "closeness"),
                             hubQuant = 0.9,
                             lnormFit = TRUE,
                             normDeg = FALSE,
                             normBetw = FALSE,
                             normClose = FALSE,
                             normEigen = FALSE)
  #aa <-net_season[["simMat1"]]
  #max(aa[aa != max(aa)])
  summary(props_season)

  # Get phyla names from the taxonomic table created before
  phyla <- as.factor(taxtab[, "Tax_1"])
  names(phyla) <- taxtab[, "Tax_6"]
  #table(phyla)
  
  # Define phylum colors
  phylcol <- c( "#ae1816", "#d49dc7","#23378d")
  
  
  plot(props_season, 
       sameLayout = FALSE, 
       nodeColor = "feature", #alternative here: "cluster"
       featVecCol = phyla, 
       colorVec =  phylcol,
       nodeSize = "clr",
       shortenLabels = "none",
       rmSingles = TRUE,
       labelScale = FALSE,
       posCol = "darkgreen", 
       negCol = "darkgrey",
       cexNodes = 0.7, 
       cexLabels = 0.5,
       cexHubLabels = 0.7,
       cexTitle = 1.7,
       #groupNames = c("Arctic", "Temperate"),
       hubBorderCol  = "gray40")


  
  print(sort(colSums(net_season$normCounts1), decreasing = TRUE)[1:10])
  
  Subarctic_network_species <- as.data.frame(sort(colSums(net_season$normCounts1), decreasing = TRUE))
  Subarctic_network_species$Tax_6 <- rownames(Subarctic_network_species)
  
  #Stats on Trophic modes
  
  
  Subarctic_trophic <- dplyr::left_join(Subarctic_network_species, Trophic_annotation, by = "Tax_6")
  Subarctic_stats <- Subarctic_trophic%>%group_by(Tax_1, Trophy)%>%count()
  Subarctic_stats_eukprok <- Subarctic_trophic%>%group_by(Tax_1)%>%count()
  
  print(Subarctic_stats)
  print(Subarctic_stats_eukprok)
  
  rm(Subarctic_trophic, Subarctic_network_species, Subarctic_stats, Subarctic_stats_eukprok)
  
  
  #####
  
  print(summary(props_season, numbNodes = 5L))

  
  
  ###ARCTIC
  
  
  net_season <- netConstruct(data = own_split$Arctic, 
                             measure = "spieceasi", ######
                             normMethod = "none", 
                             zeroMethod = "multRepl",
                             sparsMethod = "threshold", 
                             thresh = 0.3, 
                             dissFunc = "signed",
                             verbose = 3,
                             seed = 123456)
  props_season <- netAnalyze(net_season, 
                             centrLCC = FALSE,
                             avDissIgnoreInf = TRUE,
                             sPathNorm = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = c("degree", "between", "closeness"),
                             hubQuant = 0.9,
                             lnormFit = TRUE,
                             normDeg = FALSE,
                             normBetw = FALSE,
                             normClose = FALSE,
                             normEigen = FALSE)
  #aa <-net_season[["simMat1"]]
  #max(aa[aa != max(aa)])
  summary(props_season)
  
  # Get phyla names from the taxonomic table created before
  phyla <- as.factor(taxtab[, "Tax_1"])
  names(phyla) <- taxtab[, "Tax_6"]
  #table(phyla)
  
  # Define phylum colors
  phylcol <- c( "#ae1816", "#d49dc7","#23378d")
  
  
  plot(props_season, 
       sameLayout = FALSE, 
       nodeColor = "feature", #alternative here: "cluster"
       featVecCol = phyla, 
       colorVec =  phylcol,
       nodeSize = "clr",
       shortenLabels = "none",
       rmSingles = TRUE,
       labelScale = FALSE,
       posCol = "darkgreen", 
       negCol = "darkgrey",
       cexNodes = 0.7, 
       cexLabels = 0.5,
       cexHubLabels = 0.7,
       cexTitle = 1.7,
       #groupNames = c("Arctic", "Temperate"),
       hubBorderCol  = "gray40")
  
  
  
  print(sort(colSums(net_season$normCounts1), decreasing = TRUE)[1:10])
  
  
  
  arctic_network_species <- as.data.frame(sort(colSums(net_season$normCounts1), decreasing = TRUE))
  arctic_network_species$Tax_6 <- rownames(arctic_network_species)
  
  #Stats on Trophic modes
  
  
  arctic_trophic <- dplyr::left_join(arctic_network_species, Trophic_annotation, by = "Tax_6")
  arctic_stats <- arctic_trophic%>%group_by(Tax_1, Trophy)%>%count()
  arctic_stats_eukprok <- arctic_trophic%>%group_by(Tax_1)%>%count()
  
  print(arctic_stats)
  print(arctic_stats_eukprok)
  
  rm(arctic_trophic, arctic_network_species, arctic_stats, arctic_stats_eukprok)
  
  
  
  #####
  
  print(summary(props_season, numbNodes = 5L))