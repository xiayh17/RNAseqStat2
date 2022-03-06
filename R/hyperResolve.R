#' @export
hyperResolve <- function(object, GO = FALSE, KEGG = TRUE) {

  ## 获取GeneList
  test <- deg_here(object)
  ok <- names(test)[which(test == TRUE)] ## 取有效数据
  OrgDb <- hyperGOparam(object)[["OrgDb"]]

  if (KEGG) {

    usethis::ui_info("Enrich KEGG analysis Start. This process will take a few minutes.")

    ## get parameters
    keggParams <- hyperKEGGparam(object)
    hyperKEGG_GeneSets = list()

    if("TERM2GENE" %in% names(keggParams)){

      for (i in ok) {
        hyperKEGG_GeneSets[[i]] <- suppressWarnings(hyper_GS(object,which = i,type = "SYMBOL"))
      }

    } else {

      for (i in ok) {
        hyperKEGG_GeneSets[[i]] <- suppressWarnings(hyper_GS(object,which = i,type = "ENTREZID",OrgDb = OrgDb))
      }

    }

    ## 富集分析
    hyperKEGG_res <- lapply(seq_along(hyperKEGG_GeneSets), function(x){
      geneSet_list = hyperKEGG_GeneSets[[x]]
      res <- hyper_keggResolve(geneSet_list = geneSet_list,keggParams = keggParams)
      ui_done("Enrich KEGG {names(hyperKEGG_GeneSets)[x]} analysis done")

      return(res)
    })

    names(hyperKEGG_res) <- names(hyperKEGG_GeneSets)

  } else {

    hyperKEGG_res <- NULL

  }

  if (GO) {

    usethis::ui_info("Enrich GO analysis Start. This process will take a few minutes.")

    goParams <- hyperGOparam(object)
    hyperGO_GeneSets = list()

    for (i in ok) {
      hyperGO_GeneSets[[i]] <- suppressWarnings(hyper_GS(object,which = i,type = "SYMBOL"))
    }

    ## 富集分析
    hyperGO_res <- lapply(seq_along(hyperGO_GeneSets), function(x){
      geneSet_list = hyperGO_GeneSets[[x]]
      res <- hyper_goResolve(geneSet_list = x,goParams = goParams)
      ui_done("Enrich KEGG {names(hyperGO_GeneSets)[x]} analysis done")
      return(res)
    })

    names(hyperGO_res) <- names(hyperGO_GeneSets)

  } else {

    hyperGO_res <- NULL

  }

  ## 保存结果
  tmp <- hyperRes(object)
  tmp["hyperKEGG_res"] <- list(hyperKEGG_res)
  tmp["hyperGO_res"] <- list(hyperGO_res)
  hyperRes(object) <- tmp

  return(object)
}

#' @export
hyper_keggResolve <- function(...,geneSet_list,keggParams) {

  keggres <- lapply(seq_along(geneSet_list), function(x){

    gene = geneSet_list[[x]]

    tryCatch(
      expr = {
        hyperCore(gene=gene,fparams = keggParams,f = "enrichKEGG")
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured. try again.")
        hyperCore(gene=gene,fparams = keggParams,f = "enrichKEGG")
      },
      finally = {
        usethis::ui_line("Enrich KEGG {names(geneSet_list)[x]} analysis done")
      }
    )

  })

  names(keggres) <- names(geneSet_list)

  return(keggres)

}

#' @importFrom clusterProfiler enrichGO
#' @export
hyper_goResolve <- function(...,geneSet_list,goParams) {

  gores <- lapply(seq_along(geneSet_list), function(x){

    gene = geneSet_list[[x]]

    tryCatch(
      expr = {
        hyperCore(gene=gene,fparams = goParams,f = "enrichGO")
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured. try again.")
        hyperCore(gene=gene,fparams = goParams,f = "enrichGO")
      },
      finally = {
        usethis::ui_line("Enrich GO {names(geneSet_list)[x]} analysis done")
      }
    )

  })

  names(gores) <- names(geneSet_list)

  return(gores)

}

#' @importFrom clusterProfiler enrichKEGG enricher
hyperCore <- function(..., fparams, f = "enrichKEGG") {

  params <- list(...)
  fparams <- modifyList(params, fparams)

  if("TERM2GENE" %in% names(fparams)){

    core <- suppressMessages(do.call("enricher", modifyList(
      list(),
      fparams)
    ))

  } else {

    core <- suppressMessages(do.call(f, modifyList(
      list(),
      fparams)
    ))

  }

  return(core)

}
