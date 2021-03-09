## create a phyloseq

network_netCoMi <- function(ASVTable, taxonomy, meta){

suppressPackageStartupMessages(require(phyloseq))
packageVersion("phyloseq")

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
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
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
summary(props_single_own, numbNodes = 5L)

pi <- plot(props_single_own, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network", 
          showTitle = TRUE,
          cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)



# group network analysis

require(metagMisc)

own_split <- metagMisc::phyloseq_sep_variable(all, 
                                              "Glacial.influence")


# Network construction
net_glacial <- netConstruct(data = own_split$No, 
                            data2 = own_split$Yes,  
                            filtTax = "highestVar",
                            filtTaxPar = list(highestVar = 100),
                            measure = "spring",
                            measurePar = list(nlambda=10, 
                                              rep.num=10),
                            normMethod = "none", 
                            zeroMethod = "none",
                            sparsMethod = "none", 
                            dissFunc = "signed",
                            verbose = 3,
                            seed = 123456)


props_glacier <- netAnalyze(net_glacial, 
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

summary(props_glacier)


###visualize


plot(props_glacier, 
     sameLayout = FALSE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     cexNodes = 0.5, 
     cexLabels = 0.5,
     cexHubLabels = 0.7,
     cexTitle = 1.2,
     groupNames = c("No glacial influence", "Glacial influence"),
     hubBorderCol  = "gray40")

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 1, lty = 1, lwd = 1, 
       bty = "n", horiz = TRUE)


plot(props_glacier, 
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "mclr", 
     labelScale = FALSE,
     cexNodes = 0.7, 
     cexLabels = 0.7,
     cexHubLabels = 0.7,
     cexTitle = 0.8,
     groupNames = c("No glacier", "glacier"),
     hubBorderCol  = "gray40")

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)

#######################################################
########with taxonomy

#check the taxonomy ranks

rank_names(all)

# Agglomerate to order level
amgut_genus <- phyloseq::tax_glom(all, taxrank = "order")
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
     cexHubs = 1.1,
     cexLabels = 1.2,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
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

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)

}
