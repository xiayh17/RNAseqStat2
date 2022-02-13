library(msigdbr)
library(GSEABase)
library(GSVA)
msigdbr::msigdbr_species()

c("Homo sapiens","Mus musculus")

## 一次下载全部会很大 4263110行
c1_gene_sets = msigdbr(species = "Homo sapiens",category = "C8")

# 查看各个子集的数量
table(all_gene_sets$gs_cat)

# 提取一个子集分析
H_gene_sets = all_gene_sets[all_gene_sets$gs_cat=='H',]

# 只保留基因名和通路并去重
H_gene_sets=H_gene_sets[,3:4]
H_gene_sets=unique(H_gene_sets)

# 做GSEA
geneList <- GSEA_GS(data_g,"limma",type = "SYMBOL")

gsea_res_limma <- GSEA(geneList = geneList,TERM2GENE = H_gene_sets)

is.null(gsea_res_limma)

#

# 做GSVA
# 需要转化功能基因表为list
H_gene_sets_list = split(H_gene_sets$gene_symbol,H_gene_sets$gs_name)

# 构建genecollection
geneset <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, H_gene_sets_list ,  names(H_gene_sets_list)))
geneset
# gsva 是针对矩阵 ，需要的是 geneset这样的 GeneSetCollection对象
# 使用原始count 矩阵需要添加kcdf="Poisson"
# 或者使用logcpm
es.max <- gsva( as.matrix(matrixFiltered(data_g)) , geneset,
                mx.diff=FALSE,
                verbose=FALSE,
                kcdf="Poisson",
                parallel.sz= 2)

# 对GSVA 差异结果
gsva_limma_resolve <- function(gsva_data, group_list, case_group) {
  control_group = setdiff(group_list,case_group)
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(gsva_data)

  # dge <- DGEList(counts=gsva_data)
  # dge <- calcNormFactors(dge)
  # logCPM <- cpm(dge, log=TRUE, prior.count=3)

  # v <- voom(dge,design,plot=TRUE, normalize.method="quantile")
  fit <- lmFit(gsva_data, design)

  con=paste0(case_group,'-',control_group)

  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)

  tempOutput = topTable(fit2, ,coef=con,adjust='BH', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)

  return(DEG_limma_voom)

}

t <- gsva_limma_resolve(es.max,group_list = groupInfo(data_g),case_group = caseGroup(data_g))

pathway <- gsub("HALLMARK_", "",row.names(t))
df <- data.frame(ID = pathway, score = t$t)

