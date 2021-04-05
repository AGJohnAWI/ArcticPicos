## create a phyloseq
#euk_r, tax_18S_r, meta_18S_r
network_netCoMi <- function(ASVTable, taxonomy, meta){

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")
require("NetCoMi")

#netCoMi seems to only allow 6 characters for names so we have to cut out the ">ASV_"

rownames(taxonomy) <- gsub(">ASV_", "", taxonomy$ASV)
rownames(ASVTable) <- gsub(">ASV_", "", rownames(ASVTable))

###create phyloseq

OTU = otu_table(ASVTable, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(taxonomy)
TAX = tax_table(tax.matrix)
rownames(meta) <- meta$Site
map = sample_data(meta)

phyloseq_merged = phyloseq(OTU, TAX)
all = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all

## single network

ASV_t <- t(ASVTable)


net_single_own <- netConstruct(ASV_t,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 5),
                           measure = "spring",
                           measurePar = list(nlambda=20, 
                                             rep.num=20),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 123456)

props_single_own <- netAnalyze(net_single_own, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
print(summary(props_single_own, numbNodes = 5L))

pi <- plot(props_single_own, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network", 
          showTitle = TRUE,
          cexTitle = 0.7)

legend(0.7, 1.1, cex = 0.7, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)



# group network analysis

require(metagMisc)

own_split <- metagMisc::phyloseq_sep_variable(all, 
                                              "Glacial.influence") 

#######################################################
########with taxonomy

#check the taxonomy ranks

rank_names(all)

# Agglomerate to order level
amgut_genus <- all
taxtab <- amgut_genus@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Rank5)
miss_f <- which(taxtab[, "genus"] == "f__")
miss_g <- which(taxtab[, "class"] == "g__")

# Number unspecified genera
taxtab[miss_f, "genus"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "class"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "genus"]) |
                        duplicated(taxtab[, "genus"], fromLast = TRUE))

for(i in seq_along(taxtab)){
        # The next higher non-missing rank is assigned to unspecified genera
        if(i %in% miss_f && i %in% miss_g){
                taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "phylum"], ")")
        } else if(i %in% miss_g){
                taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "class"], ")")
        }
        
        # Family names are added to duplicate genera
        if(i %in% dupl_g){
                taxtab[i, "genus"] <- paste0(taxtab[i, "genus"], "(", taxtab[i, "genus"], ")")
        }
}

amgut_genus@tax_table@.Data <- taxtab
rownames(amgut_genus@otu_table@.Data) <- taxtab[, "genus"]

# Network construction and analysis
net_single3 <- netConstruct(amgut_genus, 
                            measure = "pearson",
                            zeroMethod = "multRepl",
                            normMethod = "clr", 
                            sparsMethod = "threshold", 
                            thresh = 0.3, 
                            verbose = 3)


props_single3 <- netAnalyze(net_single3, clustMethod = "cluster_fast_greedy")

# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_single3$adjaMat1, weighted = TRUE)
lay_fr <- igraph::layout_with_fr(graph3)
# Note that row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_single3$adjaMat1)

plot(props_single3,
     layout = lay_fr,
     shortenLabels = "simple",
     labelLength = 10,
     nodeSize = "fix",
     nodeColor = "gray",
     cexNodes = 0.8,
     cexHubs = 0.7,
     cexLabels = 0.7,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 0.7)

legend(0.5, 0.5, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)


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

#######################
amgut_genus <- phyloseq::tax_glom(own_split$Yes, taxrank = "order")
taxtab <- amgut_genus@tax_table@.Data

# Find undefined taxa (in this data set, unknowns occur only up to Rank5)
miss_f <- which(taxtab[, "order"] == "f__")
miss_g <- which(taxtab[, "class"] == "g__")

# Number unspecified genera
taxtab[miss_f, "order"] <- paste0("f__", 1:length(miss_f))
taxtab[miss_g, "class"] <- paste0("g__", 1:length(miss_g))

# Find duplicate genera
dupl_g <- which(duplicated(taxtab[, "order"]) |
                        duplicated(taxtab[, "order"], fromLast = TRUE))

for(i in seq_along(taxtab)){
        # The next higher non-missing rank is assigned to unspecified genera
        if(i %in% miss_f && i %in% miss_g){
                taxtab[i, "order"] <- paste0(taxtab[i, "order"], "(", taxtab[i, "phylum"], ")")
        } else if(i %in% miss_g){
                taxtab[i, "order"] <- paste0(taxtab[i, "order"], "(", taxtab[i, "class"], ")")
        }
        
        # Family names are added to duplicate genera
        if(i %in% dupl_g){
                taxtab[i, "order"] <- paste0(taxtab[i, "order"], "(", taxtab[i, "order"], ")")
        }
}

amgut_genus@tax_table@.Data <- taxtab
rownames(amgut_genus@otu_table@.Data) <- taxtab[, "order"]

# Network construction and analysis
net_single3 <- netConstruct(amgut_genus, 
                            measure = "pearson",
                            zeroMethod = "multRepl",
                            normMethod = "clr", 
                            sparsMethod = "threshold", 
                            thresh = 0.3, 
                            verbose = 3)


props_single3 <- netAnalyze(net_single3, clustMethod = "cluster_fast_greedy")

# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_single3$adjaMat1, weighted = TRUE)
lay_fr <- igraph::layout_with_fr(graph3)
# Note that row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_single3$adjaMat1)

plot(props_single3,
     layout = lay_fr,
     shortenLabels = "simple",
     labelLength = 10,
     nodeSize = "fix",
     nodeColor = "gray",
     cexNodes = 0.8,
     cexHubs = 0.7,
     cexLabels = 0.7,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 0.7)

legend(0.5, 0.5, cex = 0.7, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)


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


}
