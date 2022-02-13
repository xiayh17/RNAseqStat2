

rm(list = "i")

usethis::ui_info("Enrich KEGG analysis Start. This process will take a few minutes.")

test <- lapply(gene_list, function(x)
  suppressMessages(enrichKEGG(
    gene = x,
    organism = organism,
    keyType = keyType,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    qvalueCutoff = qvalueCutoff,
    use_internal_data = use_internal_data)))

my.cluster <- parallel::makeCluster(
  spec = length(gene_list),
  type = "PSOCK"
)
doParallel::registerDoParallel(cores = length(gene_list))
foreach::getDoParWorkers()
res <- list()
`%dopar%` <- foreach::`%dopar%`
require(foreach)
res[[names(gene_list)[i]]] <- foreach::foreach(i = 1:length(gene_list)) %dopar% {

  suppressMessages(clusterProfiler::enrichKEGG(
    gene = gene_list[[i]],
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 50,
    maxGSSize = 500,
    qvalueCutoff = 1,
    use_internal_data = T))

}
parallel::stopCluster(cl = my.cluster)
