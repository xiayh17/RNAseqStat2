devtools::load_all()
load(file = 'testDataSets/Wang/GSE184648-2分组3重复-ensembl的ID/input.Rdata')
counts_input[1:4,1:4]
table(group_list)
colnames(counts_input)=LETTERS[1:6]
counts_input <- as.data.frame(counts_input)
# library(RNAseqStat)

counts_input2 <- counts_input[,c(1,2,3,4)]
group_list2 <- group_list[c(1,2,3,4)]
runAll(count_data = counts_input2,
       group_list = group_list2,
       case_group = "p71",
       control_group = "p66",
       OrgDb = "org.Hs.eg.db",
       dir ="GSE184648_results/")



# 修改keytype为ENSEMBL
load("GSE184648_results/2-DEG_results.Rdata")
go_results = hyper_go(deg_data = deg_results@deg_df_limma,
                       x = "logFC" , y = "P.Value",
                       keyType = "ENSEMBL",
                       OrgDb = "org.Hs.eg.db",
                       cut_FC = 2,
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.2)


save(go_results,file = "./GSE184648_results/go_results")

load("./GSE184648_results/go_results")
library(ggplot2)
go_down = go_results[["Down"]]@result
rownames(go_down)<-NULL
Description_order=factor(as.integer(rownames(go_down)),labels=go_down$Description)
pdf(file = "./GSE184648_results/GO.pdf")
ggplot(data = go_down, aes(x=Description_order,y=Count,fill=ONTOLOGY)) +
  geom_bar(stat="identity", width=0.8)+
  theme_bw()+
  xlab("GO term") + ylab("Num of Genes")+
  coord_flip()

dev.off()

