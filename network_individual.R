#only eukaryotes

taxonomy_18S_r_rename$Tax_6 <- gsub("uncultured","",taxonomy_18S_r_rename$Tax_6)
taxonomy_18S_r_rename$Tax_5 <- gsub("uncultured","",taxonomy_18S_r_rename$Tax_5)
taxonomy_18S_r_rename$Tax_4 <- gsub("uncultured","",taxonomy_18S_r_rename$Tax_4)


taxonomy_18S_r_rename[taxonomy_18S_r_rename==" "]<-NA
taxonomy_18S_r_rename[taxonomy_18S_r_rename==""]<-NA

taxonomy_18S_r_rename$Tax_6[is.na(taxonomy_18S_r_rename$Tax_6)] <- "f__"
taxonomy_18S_r_rename$Tax_5[is.na(taxonomy_18S_r_rename$Tax_5)] <- "g__"


OTU = otu_table(Euk_r_clr_rename, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy_18S_r_rename)
TAX = tax_table(tax.matrix)
rownames(meta_all) <- meta_all$Station
map = sample_data(meta_all)

phyloseq_merged = phyloseq(OTU, TAX)
all_e = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all_e


  # Agglomerate to order level
  
  all_T6 <- all_e
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
                                                "Glacial.influence")
  
  # Find undefined taxa (in this data set, unknowns occur only up to Rank5)
  
  rownames(own_split$No@otu_table@.Data) <- own_split$No@tax_table@.Data[, "Tax_6"]
  rownames(own_split$Yes@otu_table@.Data) <- own_split$Yes@tax_table@.Data[, "Tax_6"]
  
  net_season <- netConstruct(data = own_split$No, 
                             data2 = own_split$Yes,  
                             measure = "spieceasi",
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
  
  summary(props_season)
  
  plot(props_season, 
       sameLayout = FALSE, 
       nodeColor = "cluster",
       nodeSize = "mclr",
       shortenLabels = "none",
       labelScale = FALSE,
       cexNodes = 0.7, 
       cexLabels = 0.5,
       cexHubLabels = 0.7,
       cexTitle = 1.7,
       groupNames = c("no glacial influence", "glacial influence"),
       hubBorderCol  = "gray40")
  
  
  
  print(sort(colSums(net_season$normCounts1), decreasing = TRUE)[1:10])
  #####
  
  print(summary(props_season, numbNodes = 5L))
  

##prokaryotes
  
  taxonomy_16S_r_rename$Tax_6 <- gsub("uncultured","",taxonomy_16S_r_rename$Tax_6)
  taxonomy_16S_r_rename$Tax_5 <- gsub("uncultured","",taxonomy_16S_r_rename$Tax_5)
  taxonomy_16S_r_rename$Tax_4 <- gsub("uncultured","",taxonomy_16S_r_rename$Tax_4)
  
  
  taxonomy_16S_r_rename[taxonomy_16S_r_rename==" "]<-NA
  taxonomy_16S_r_rename[taxonomy_16S_r_rename==""]<-NA
  
  taxonomy_16S_r_rename$Tax_6[is.na(taxonomy_16S_r_rename$Tax_6)] <- "f__"
  taxonomy_16S_r_rename$Tax_5[is.na(taxonomy_16S_r_rename$Tax_5)] <- "g__"
  
  OTU = otu_table(Prok_r_clr_rename, taxa_are_rows = TRUE)
  
  tax.matrix<- as.matrix(taxonomy_16S_r_rename)
  TAX = tax_table(tax.matrix)
  rownames(meta_all) <- meta_all$Station
  map = sample_data(meta_all)
  
  phyloseq_merged = phyloseq(OTU, TAX)
  all_p = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
  all_p

all_T6 <- all_p
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
                                              "Glacial.influence")

# Find undefined taxa (in this data set, unknowns occur only up to Rank5)

rownames(own_split$No@otu_table@.Data) <- own_split$No@tax_table@.Data[, "Tax_6"]
rownames(own_split$Yes@otu_table@.Data) <- own_split$Yes@tax_table@.Data[, "Tax_6"]

net_season <- netConstruct(data = own_split$No, 
                           data2 = own_split$Yes,  
                           measure = "spieceasi",
                           normMethod = "none", 
                           zeroMethod = "multRepl",
                           sparsMethod = "threshold", 
                           thresh = 0.3,  #see this discussion: https://github.com/zdk123/SpiecEasi/issues/85 
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

summary(props_season)

plot(props_season, 
     sameLayout = FALSE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     shortenLabels = "none",
     labelScale = FALSE,
     cexNodes = 0.7, 
     cexLabels = 0.5,
     cexHubLabels = 0.7,
     cexTitle = 1.7,
     groupNames = c("no glacial influence", "glacial influence"),
     hubBorderCol  = "gray40")



print(sort(colSums(net_season$normCounts1), decreasing = TRUE)[1:10])
#####

print(summary(props_season, numbNodes = 5L))



