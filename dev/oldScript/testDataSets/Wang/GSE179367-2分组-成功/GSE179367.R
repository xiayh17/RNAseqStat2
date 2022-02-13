rm(list = ls())
devtools::load_all()
load(file = 'input.Rdata')
counts_input[1:4,1:4]
table(group_list)
library(RNAseqStat)
runAll(count_data = counts_input,
       group_list = group_list,
       case_group = "MMP9", control_group = "Ctrl",
       OrgDb = "org.Hs.eg.db",
       dir ="GSE179367_results/")

