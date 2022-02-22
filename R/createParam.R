## for hyper module ----
#' @export
Create_hyperParam <- function(goParam = NULL,keggParam = NULL,
                              customGO = FALSE,customKEGG = FALSE) {

  customParams = list(
    pvalueCutoff = 0.99,
    pAdjustMethod = "BH",
    universe = NULL,
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff = 0.99,
    TERM2GENE = NULL, ## if custom* = true ,it must be set
    TERM2NAME = NA)

  ## 设置默认选项 ----
  if (customGO) {
    goParam_default = customParams
  } else {
    goParam_default = list(
      OrgDb = 'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.99,
      pAdjustMethod = "BH",
      universe = NULL,
      qvalueCutoff = 0.99,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = FALSE)
  }

  if(customKEGG) {
    keggParam_default = customParams
  } else {
    keggParam_default = list(
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.99,
      pAdjustMethod = "BH",
      universe = NULL,
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.99,
      use_internal_data = FALSE)
  }

  goParam = paramParse(defaultParams = goParam_default,newParam = goParam)
  keggParam = paramParse(defaultParams = keggParam_default,newParam = keggParam)

  new("hyperParam",
      goParam = goParam,
      keggParam = keggParam)

}
## ----

## for gse module ----
Create_gseParam <- function(goParam = NULL,keggParam = NULL,
                            customGO = FALSE,customKEGG = FALSE) {

  customParams = list(
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps  = 0,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    TERM2GENE = NULL,
    TERM2NAME = NA,
    verbose = TRUE,
    seed = FALSE,
    by = 'fgsea'
  )

  ## 设置默认选项 ----
  if (customGO) {
    goParam_default = customParams
  } else {

    goParam_default = list(
      ont = "ALL",
      OrgDb = NULL,
      keyType = "SYMBOL",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0,
      pvalueCutoff = 0.99,
      pAdjustMethod = "BH",
      verbose = TRUE,
      seed = FALSE,
      by = "fgsea")

  }

  if (customKEGG) {
    keggParam_default = customParams
  } else {

    keggParam_default = list(
      organism = NULL,
      keyType = "kegg",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0,
      pvalueCutoff = 0.99,
      pAdjustMethod = "BH",
      verbose = TRUE,
      use_internal_data = FALSE,
      seed = FALSE,
      by = "fgsea")

  }

  goParam = paramParse(defaultParams = goParam_default,newParam = goParam)
  keggParam = paramParse(defaultParams = keggParam_default,newParam = keggParam)

  new("gseParam",
      goParam = goParam,
      keggParam = keggParam)

}
## ----

## msigdb module
#' Create \code{MSigDB}
#'
#' Create a msigdbParam list. If NULL, default parameters will apply.
#'
#' @param msigdbParam a list contains any paramaters of \code{\link[msigdbr]{msigdbr}}
#'
#' @return a list
#' @export
Create_msigdbParam <- function(msigdbParam = NULL) {

  ## 设置默认选项 ----
  msigdbParam_default = list(
    species = NULL, category = NULL, subcategory = NULL)

  msigdbParam = paramParse(defaultParams = msigdbParam_default, newParam = msigdbParam)

  return(msigdbParam)

}

#' Create \code{MSigDB}
#'
#' Create a GSEAparam list. If NULL, default parameters will apply.
#'
#' @param GSEAparam a list contains any paramaters of \code{\link[clusterProfiler]{GSEA}}
#'
#' @return a list
#' @export
Create_msigdbGSEAparam <- function(GSEAparam = NULL) {

  GSEAparam_default = list(
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 0.99,
    pAdjustMethod = "BH",
    TERM2NAME = NA,
    verbose = TRUE,
    seed = FALSE,
    by = "fgsea"
  )

  GSEAparam = paramParse(defaultParams = GSEAparam_default,newParam = GSEAparam)

  return(GSEAparam)

}

Create_msigdbHyperParam <- function(HyperParam = NULL) {

  HyperParam_default = list(
    pvalueCutoff = 0.99,
    pAdjustMethod = "BH",
    universe = NULL,
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff = 0.99,
    TERM2GENE = NULL, ## if custom* = true ,it must be set
    TERM2NAME = NA)

  HyperParam = paramParse(defaultParams = HyperParam_default,newParam = HyperParam)

  return(HyperParam)

}

## ----

paramParse <- function(defaultParams = list(),newParam = NULL) {

  if (is.null(newParam))
    newParam <- defaultParams

  if (is.list(newParam)&length(newParam) == 0)
    newParam <- defaultParams

  if (is.list(newParam)&length(newParam) > 0) { ## is a list class and include params

    if (!all(names(newParam) %in% names(defaultParams))) { ## is a list class and include params with default param

      usethis::ui_oops("Your {usethis::ui_code('newParam')} will not be applied.\n Please check the parameters of {usethis::ui_code('paramParse')} ")

    } else if (length(newParam) == 1) { ## if only one param include

      defaultParams[[names(newParam)]] <- newParam[names(newParam)][[1]]

      newParam <- defaultParams

    } else if (length(newParam) > 1) { ## if params more than one

      for (x in seq_along(newParam)) {

        defaultParams[[names(newParam)[x]]] <- newParam[names(newParam)[x]][[1]]

      }

      newParam <- defaultParams

    }


  } else { ## if new param is nothing

    usethis::ui_oops("Your {usethis::ui_code('newParam')} will not be applied.\n Please check the usage of {usethis::ui_code('paramParse')} ")

  }

  return(newParam)

}
