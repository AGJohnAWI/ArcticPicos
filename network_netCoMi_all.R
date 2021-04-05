#combined 16S and 18S

#meta_all, merged_clr_all, taxonomy_merged_clr

  
  # Agglomerate to order level
  
a <- function(amgut_genus){taxtab <- amgut_genus@tax_table@.Data
  
  # Find undefined taxa (in this data set, unknowns occur only up to Rank5)
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
    } else if(i %in% miss_g){
      taxtab[i, "Tax_6"] <- paste0(taxtab[i, "Tax_6"], "(", taxtab[i, "Tax_5"], ")")
    }
    
    # Family names are added to duplicate genera
    if(i %in% dupl_g){
      taxtab[i, "Tax_6"] <- paste0(taxtab[i, "Tax_6"], "(", taxtab[i, "Tax_6"], ")")
    }
  }
  
  amgut_genus@tax_table@.Data <- taxtab
  rownames(amgut_genus@otu_table@.Data) <- taxtab[, "Tax_6"]
  
  # Network construction and analysis
  net_single3 <- netConstruct(amgut_genus, 
                              measure = "spieceasi",
                              zeroMethod = "multRepl",
                              normMethod = "none", 
                              sparsMethod = "threshold", 
                              thresh = 0.3, 
                              verbose = 3)
  
  
  props_single3 <- netAnalyze(net_single3, clustMethod = "cluster_fast_greedy")
  
  # Compute layout
  graph3 <- igraph::graph_from_adjacency_matrix(net_single3$adjaMat1, weighted = TRUE)
  lay_fr <- igraph::layout_with_fr(graph3)
  # Note that row names of the layout matrix must match the node names
  rownames(lay_fr) <- rownames(net_single3$adjaMat1)
  
  set.seed(123456)
  graph3 <- igraph::graph_from_adjacency_matrix(net_single3$adjaMat1, weighted = TRUE)
  lay_fr <- igraph::layout_with_fr(graph3)
  rownames(lay_fr) <- rownames(net_single3$adjaMat1)
  
  plot(props_single3,
       layout = "layout_with_fr",
       shortenLabels = "simple",
       labelLength = 10,
       charToRm = "g__",
       labelScale = FALSE,
       rmSingles = TRUE,
       nodeSize = "clr",
       nodeColor = "cluster",
       hubBorderCol = "darkgray",
       cexNodes = 0.7,
       cexLabels = 0.5,
       cexHubLabels = 0.7,
       title1 = "Network on genus level with Pearson correlations", 
       showTitle = TRUE,
       cexTitle = 0.8)
  
  legend(0.01, 0.9, cex = 0.7, title = "estimated correlation:",
         legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
         bty = "n", horiz = TRUE)
  
  print(summary(props_single3, numbNodes = 5L))
  
}
