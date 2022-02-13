count_data = exp
group_list = c("C","C","T")
control_group = "C"
case_group = "T"
dir = "results"

library(RNAseqStat)
runAll(count_data = count_data,
       group_list = group_list,
       case_group = case_group, control_group = control_group,
       OrgDb = "org.Mm.eg.db",
       dir = dir)

#message(glue("Step3: EnrichGO analysis"))
load("results/2-DEG_results.Rdata")
run_hyperGO(deg_results@deg_df_DESeq2,
             x = "log2FoldChange", y = "pvalue",
             cut_FC = 3.328,
             cut_FDR = 0.05,
             OrgDb = org.Mm.eg.db, dir = dir,prefix = "3-EnrichGO-DESeq2")
