
#create phyloseq

rownames(tax_18S_tipmouth) <- tax_18S_tipmouth$ASV

OTU = otu_table(asv_tip_mouth_euk_a, taxa_are_rows = TRUE)

tax.matrix<- as.matrix(tax_18S_tipmouth)
TAX = tax_table(tax.matrix)
rownames(meta_tip_mouth_18S) <- meta_tip_mouth_18S$Site
map = sample_data(meta_tip_mouth_18S)

phyloseq_merged = phyloseq(OTU, TAX)
all_tm = merge_phyloseq(phyloseq_merged, map) ####merge data into phyloseq
all_tm



amgut_split <- metagMisc::phyloseq_sep_variable(all_tm, 
                                                "Glacial.influence")

# Network construction
net_season <- netConstruct(data = amgut_split$No, 
                           data2 = amgut_split$Yes,  
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 500),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
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
