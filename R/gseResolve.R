#' @export
gseResolve <- function(object) {

  ## 获取GeneList
  test <- deg_here(object)
  ok <- names(test)[which(test == TRUE)] ## 取有效数据

  gseKEGG_GeneSets = list()
  for (i in ok) {
    gseKEGG_GeneSets[[i]] <- suppressWarnings(GSEA_GS(object,which = i,type = "ENTREZID"))
  }

  ## 富集分析
  keggParams <- gseKEGGparam(object)

  gseKEGG_res <- gse_keggResolve(geneSet_list = gseKEGG_GeneSets,keggParams = keggParams)


  ## 保存结果
  tmp <- gseRes(object)
  tmp["gseKEGG_res"] <- list(gseKEGG_res)
  gseRes(object) <- tmp

  return(object)
}

#' @importFrom clusterProfiler gseKEGG
#' @export
gse_keggResolve <- function(...,geneSet_list,keggParams) {

  usethis::ui_info("Enrich KEGG analysis Start. This process will take a few minutes.")

  keggres <- lapply(geneSet_list, function(x){

    tryCatch(
      expr = {
        gse_keggCore(gene=x,keggParams=keggParams)
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured. try again.")
        gse_keggCore(gene=x,keggParams=keggParams)
      },
      finally = {
        usethis::ui_line("Enrich KEGG analysis done")
      }
    )

  })

}

gse_keggCore <- function(...,keggParams){

  params <- list(...)
  keggParams <- modifyList(params, keggParams)
  kegg_core <- suppressMessages(do.call("gseKEGG", modifyList(
    list(),
    keggParams)
  ))

  return(kegg_core)

}

# gse_goCore <- function(...,goParams){
#
#   params <- list(...)
#   goParams <- modifyList(params, goParams)
#   go_core <- suppressMessages(do.call("gseGO", modifyList(
#     list(),
#     goParams)
#   ))
#
#   return(go_core)
#
# }
