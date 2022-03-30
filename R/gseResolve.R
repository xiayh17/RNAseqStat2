#' GSEA analysis
#'
#' @param object a DEGContainer
#' @param GO run GO database
#' @param KEGG run KEGG database
#'
#' @return a DEGContainer
#' @export
#'
#' @examples
#' gseResolve(DEGContainer)
gseResolve <- function(object, GO = FALSE, KEGG = TRUE) {

  ## 获取GeneList
  test <- deg_here(object)
  ok <- names(test)[which(test == TRUE)] ## 取有效数据

  if (is.null(gseGOparam)) {
    GO = FALSE
  }

  if (is.null(gseKEGGparam)) {
    KEGG = FALSE
  }

  OrgDb = gseGOparam(object)[["OrgDb"]]

  if(KEGG) {

    usethis::ui_info("GSE KEGG analysis Start. This process will take a few minutes.")

    keggParams <- gseKEGGparam(object)
    gseKEGG_GeneSets = list()

    if("TERM2GENE" %in% names(keggParams)){
      for (i in ok) {
        gseKEGG_GeneSets[[i]] <- suppressWarnings(GSEA_GS(object,which = i,type = "SYMBOL"))
      }
    } else {
      for (i in ok) {
        gseKEGG_GeneSets[[i]] <- suppressWarnings(GSEA_GS(object,which = i,type = "ENTREZID",OrgDb = OrgDb))
      }
    }

    ## 富集分析
    gseKEGG_res <- gse_keggResolve(geneSet_list = gseKEGG_GeneSets,keggParams = keggParams)

  } else {

    gseKEGG_res <- NULL

  }

  if(GO) {

    usethis::ui_info("GSE GO analysis Start. This process will take a few minutes.")

    goParams <- gseGOparam(object)
    gseGO_GeneSets = list()

    for (i in ok) {
      gseGO_GeneSets[[i]] <- suppressWarnings(GSEA_GS(object,which = i,type = "SYMBOL"))
    }

    ## 富集分析
    gseGO_res <- gse_goResolve(geneSet_list = gseGO_GeneSets,goParams = goParams)

  } else {

    gseGO_res = NULL

  }


  ## 保存结果
  tmp <- gseRes(object)
  tmp["gseKEGG_res"] <- list(gseKEGG_res)
  tmp["gseGO_res"] <- list(gseGO_res)
  gseRes(object) <- tmp

  return(object)
}

#' @importFrom clusterProfiler gseKEGG
#' @export
gse_keggResolve <- function(...,geneSet_list,keggParams) {

  keggres <- lapply(seq_along(geneSet_list), function(x){

    gene = geneSet_list[[x]]

    tryCatch(
      expr = {
        gseCore(gene=gene,fparams=keggParams,f = "gseKEGG")
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured. try again.")
        gseCore(gene=gene,fparams=keggParams,f = "gseKEGG")
      },
      finally = {
        usethis::ui_line("GSE KEGG {names(geneSet_list)[x]} analysis done")
      }
    )

  })

  names(keggres) <- names(geneSet_list)

  return(keggres)

}

#' @importFrom clusterProfiler gseGO
#' @export
gse_goResolve <- function(...,geneSet_list,goParams) {

  gores <- lapply(seq_along(geneSet_list), function(x){

    gene = geneSet_list[[x]]

    tryCatch(
      expr = {
        gseCore(gene=gene,fparams=goParams,f = "gseGO")
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured. try again.")
        gseCore(gene=gene,fparams=goParams,f = "gseGO")
      },
      finally = {
        usethis::ui_line("GSE GO {names(geneSet_list)[x]} analysis done")
      }
    )

  })

  names(gores) <- names(geneSet_list)

  return(gores)

}

#' @importFrom clusterProfiler gseKEGG GSEA
gseCore <- function(..., fparams, f = "gseKEGG") {

  params <- list(...)
  fparams <- modifyList(params, fparams)
  f2 = "GSEA"
  f3 = "gseGO2"

  if("ont" %in% names(fparams)&"TERM2GENE" %in% names(fparams)){

    core <- suppressMessages(do.call(f3, modifyList(
      list(),
      fparams)
    ))

  } else if ("TERM2GENE" %in% names(fparams)&!"ont" %in% names(fparams)) {

    core <- suppressMessages(do.call(f2, modifyList(
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
