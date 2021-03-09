#Indicator species extraction

#install.packages("indicspecies")

indicator <- function(ASV_hellinger, grouping_var, taxonomy){
library(indicspecies)
abund = t(ASV_hellinger)
var_d = grouping_var
inv = multipatt(abund, var_d, func = "r.g", control = how(nperm=9999))
print(summary(inv))

Indicator_summary <- inv$sign
Indicator_summary$ASV <- rownames(Indicator_summary)
Indicator_summary_t <<- dplyr::left_join(Indicator_summary, taxonomy, by="ASV") # this needs to be adjusted

print(Indicator_summary_t)

rm(abund, var_d)
}