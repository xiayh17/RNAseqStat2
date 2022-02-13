rm(list = ls())
load(file = 'testDataSets/Wang/GSE186581-2分组2重复-ensembl的ID/input.Rdata')
counts_input[1:4,1:4]
table(group_list)



# library(RNAseqStat)
devtools::load_all()
runAll(count_data = counts_input,
       group_list = group_list,
       case_group = "autophagy", control_group = "control",
       OrgDb = "org.Hs.eg.db",
       dir ="GSE186581_results/")

results_dir = "GSE186581_results/"
pre_check(counts_data = counts_input, group_list = group_list, dir = results_dir)
run_deg(counts_input,group_list,case_group = "T", control_group = "C",dir = results_dir)
run_hyperGO(DEG_df, x = "log2FoldChange", y = "pvalue", dir = results_dir)
run_hyperKEGG(DEG_df, x = "log2FoldChange", y = "pvalue", dir = results_dir)
run_gseKEGG(DEG_df, x = "log2FoldChange", dir = results_dir)
