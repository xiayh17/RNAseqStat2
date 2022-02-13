hyperResolve <- function(object) {

  ## 获取GeneList
  test <- deg_here(object)
  ok <- names(test)[which(test == TRUE)] ## 取有效数据

  hyperKEGG_GeneSets = list()
  for (i in ok) {
    hyperKEGG_GeneSets[[i]] <- suppressWarnings(hyper_GS(object,which = i,type = "ENTREZID"))
  }

  ## 富集分析
  keggParams <- hyperKEGGparam(object)

  hyperKEGG_res <- lapply(hyperKEGG_GeneSets, function(x){
    hyper_keggResolve(geneSet_list = x,keggParams = keggParams)
  })

  ## 保存结果
  tmp <- hyperRes(object)
  tmp["hyperKEGG_res"] <- list(hyperKEGG_res)
  hyperRes(object) <- tmp

  return(object)
}

# goResolve <- function() {
#
#   gene_list <- hyperKEGG_GS(object)
#   test <- mclapply(gene_list, function(x)
#     suppressMessages(enhance_enrichGO(gene = x,OrgDb = OrgDb,keyType = keyType,ont = ont,simplify = simplify,
#                                       pvalueCutoff = pvalueCutoff,
#                                       pAdjustMethod = pAdjustMethod,
#                                       qvalueCutoff = qvalueCutoff,
#                                       minGSSize = minGSSize,
#                                       maxGSSize = maxGSSize,
#                                       readable = FALSE,
#                                       pool = FALSE)),mc.cores = getOption("mc.cores", mc.cores))
# }

#' @importFrom clusterProfiler enrichKEGG
hyper_keggResolve <- function(...,geneSet_list,keggParams) {

  usethis::ui_info("Enrich KEGG analysis Start. This process will take a few minutes.")

  keggres <- lapply(geneSet_list, function(x){

    tryCatch(
      expr = {
        hyper_keggCore(gene=x,keggParams=keggParams)
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured. try again.")
        hyper_keggCore(gene=x,keggParams=keggParams)
      },
      finally = {
        usethis::ui_line("Enrich KEGG analysis done")
      }
    )

  })

}

hyper_keggCore <- function(...,keggParams){

  params <- list(...)
  keggParams <- modifyList(params, keggParams)
  kegg_core <- suppressMessages(do.call("enrichKEGG", modifyList(
    list(),
    keggParams)
  ))

  return(kegg_core)

}
