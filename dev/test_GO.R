list_ok <- list()
a=c(1,2,3,4,5)
b=c(6,7,8,9,10)
c=c(11,12,13,14,15)
d=c(16, 17, 18, 19, 20)
e=c(21,22,23,24,25)

list_ok[[1]]=a
list_ok[[2]]=b
list_ok[[3]]=c
list_ok[[4]]=d
list_ok[[5]]=e

new_list <- vector(mode = "list", length=5)

for (i in 1:5) {
  for(k in 1:5) {
    new_list[[i]][k] <- list_ok[[i]][k]
  }
}

new_list

library(purrr)

map(list_ok, function(x)
  map_dbl(seq_along(list_ok), function(y) x[[y]])
  )
DEG_df_g <- cut_much(DEG_df,x = "log2FoldChange",y = "pvalue",cut_FC = 2,cut_FDR = 0.01)
gene_list <- list(
  Up = row.names(DEG_df_g[which(DEG_df_g$group == "Up"),]),
  Down = row.names(DEG_df_g[which(DEG_df_g$group == "Down"),])
)

DEG_df_g <- cut_much(DEG_df,x = "log2FoldChange",y = "pvalue",cut_FC = 2,cut_FDR = 0.01)
ll <- DEG_df_g[which(DEG_df_g$group %in% c("Up","Down")),]

ont = c("BP","CC","MF","ALL")

ego <- map(gene_list,function(x)
  enhance_enrichGO(gene =x, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont="BP", pvalueCutoff=0.05)
)

GO_DATA <- enhance_get_GO_data(OrgDb = 'org.Hs.eg.db', ont = "ALL", keytype = "SYMBOL")

microbenchmark::microbenchmark(
  test1 <- enhance_enricher_internal(gene = gene_list[[1]],
                                     pvalueCutoff=0.05,
                                     pAdjustMethod="BH",
                                     universe = NULL,
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500,
                                     USER_DATA = GO_DATA
  ),

  test2 <- DOSE:::enricher_internal(
    gene = gene_list[[1]],
    pvalueCutoff=0.05,
    pAdjustMethod="BH",
    universe = NULL,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    USER_DATA = GO_DATA
  )
)

library(furrr)
plan(multisession, workers = 1)

tt <- enhance_enrichGO(gene =gene_list[[1]], OrgDb = 'org.Hs.eg.db',
                       keyType = "SYMBOL", ont="ALL", pvalueCutoff=0.05, simplify = TRUE)

tt2 <- enrichGO(gene =gene_list[[1]], OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont="ALL", pvalueCutoff=0.05)

lres <- furrr::future_map(c("BP", "CC", "MF"), function(ont)
  suppressMessages(enhance_enrichGO(gene, OrgDb, keyType, ont,
                                    pvalueCutoff, pAdjustMethod, universe,
                                    qvalueCutoff, minGSSize, maxGSSize
  )),.options = furrr_options(seed = TRUE)
)

lres2 <- lres[!vapply(lres, is.null, logical(1))]
lres3 <- furrr::future_map(lres2, function(x) clusterProfiler::simplify(x))

df <- do.call('rbind', future_map(lres, as.data.frame))

wt <- as.data.frame(lres[[1]])

data(geneList, package = "DOSE")
de <- names(geneList)[1:100]

library(furrr)
plan(multisession, workers = 1)
# create cluster object
cl <- makeCluster(3)

# test each number in sample_numbers for primality
microbenchmark::microbenchmark(
  ego <- map(gene_list,function(x)
    RNAseqStat::enhance_enrichGO(gene =x, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont="CC", pvalueCutoff=0.05)
  ),
  # yy1 <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01),
  ego2 <- furrr::future_map(gene_list,function(x)
    RNAseqStat::enhance_enrichGO(gene =x, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont="CC", pvalueCutoff=0.05)
  ),
  ego3 <- parLapply(cl,gene_list,function(x)
    RNAseqStat::enhance_enrichGO(gene =x, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont="CC", pvalueCutoff=0.05)),
  times = 1L
)

# close cluster object
stopCluster(cl)

