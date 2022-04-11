library(AnnoProbe)
suppressPackageStartupMessages(library(GEOquery))

## 使用AnnoProbe 获取数据
gset=AnnoProbe::geoChina('GSE1009')

# check the ExpressionSet
eSet=gset[[1]]
# extract the expression matrix and phenotype data
probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
probes_expr=limma::normalizeBetweenArrays(probes_expr)
boxplot(probes_expr,las=2)
## pheno info
phenoDat <- pData(eSet)
head(phenoDat[,1:4])

## check GPL and annotate the probes to genes.
gpl=eSet@annotation
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
genes_expr <- filterEM(probes_expr,probe2gene )
head(genes_expr)

# do DEG
## define the group
group_list= c(rep('Control',3),rep('Diabetes',3))
table(group_list)

## Create DEGContainer
data_i <- Create_DEGContainer(
  species = "Human",
  dataType = "Array",
  expMatrix = genes_expr,
  groupInfo = group_list,
  caseGroup = "Diabetes",idType = "SYMBOL",filterMethod = NULL
)

## Run workflow
data_o <- runALL(data_i,dir = "testArray",top = 10,parallel = T,GO = T,KEGG = T)

