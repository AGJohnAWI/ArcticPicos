#combined 16S and 18S

#meta_all, merged_clr_all, taxonomy_merged_clr

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
#install.packages("igraph")
#install.packages("qgraph")
###create phyloseq

OTU = otu_table(merged_clr_all, taxa_are_rows = TRUE)

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
