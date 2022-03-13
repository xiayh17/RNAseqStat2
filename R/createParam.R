## for hyper module ----
#' Create param for hyper analysis
#'
#' @param goParam arguments in \code{\link[clusterProfiler]{enrichGO}} or \code{\link[clusterProfiler]{enricher}} except \code{gene}.
#' @param keggParam arguments in \code{\link[clusterProfiler]{enrichKEGG}} or \code{\link[clusterProfiler]{enricher}} except \code{gene}.
#' @param customGO logic. FALSE for \code{\link[clusterProfiler]{enrichGO}}, TRUE for \code{\link[clusterProfiler]{enricher}}
#' @param customKEGG logic. FALSE for \code{\link[clusterProfiler]{enrichKEGG}}, TRUE for \code{\link[clusterProfiler]{enricher}}
#' @param skipGO logical, set parameters as NULL for skip this step
#' @param skipKEGG logical, set parameters as NULL for skip this step
#' @return hyperParam
#' @export
#'
#' @examples
#' Create_hyperParam()
Create_hyperParam <- function(goParam = NULL,keggParam = NULL,
                              customGO = FALSE,customKEGG = FALSE,
                              skipGO = FALSE,skipKEGG = FALSE) {

  customKEGGParams = list(
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = NULL,
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff = 1,
    TERM2GENE = NULL, ## if custom* = true ,it must be set
    TERM2NAME = NA)

  customGOParam = list(TERM2GENE = NULL,
                       TERM2NAME = NA,
                       organism = "UNKNOW",
                       keyType = "SYMBOL",
                       ont = "ALL",
                       pvalueCutoff=1,
                       pAdjustMethod="BH",
                       universe = NULL,
                       qvalueCutoff = 1,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pool=FALSE)

  ## 设置默认选项 ----
  if (customGO) {
    goParam_default = customGOParam
  } else {
    goParam_default = list(
      OrgDb = 'org.Hs.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      universe = NULL,
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = FALSE)
  }

  if(customKEGG) {
    keggParam_default = customKEGGParams
  } else {
    keggParam_default = list(
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      universe = NULL,
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 1,
      use_internal_data = FALSE)
  }

  goParam = paramParse(defaultParams = goParam_default,newParam = goParam)
  keggParam = paramParse(defaultParams = keggParam_default,newParam = keggParam)

  if (skipGO) {
    goParam = NULL
  }

  if (skipKEGG) {
    keggParam = NULL
  }

  new("hyperParam",
      goParam = goParam,
      keggParam = keggParam)

}
## ----

## for gse module ----
#' Create param for gse analysis
#'
#' @param goParam arguments in \code{\link[clusterProfiler]{gseGO}} or \code{\link[clusterProfiler]{GSEA}} except \code{geneList}.
#' @param keggParam arguments in \code{\link[clusterProfiler]{gseKEGG}} or \code{\link[clusterProfiler]{GSEA}} except \code{geneList}.
#' @param customGO logic. FALSE for \code{\link[clusterProfiler]{gseGO}}, TRUE for \code{\link[clusterProfiler]{GSEA}}
#' @param customKEGG logic. FALSE for \code{\link[clusterProfiler]{gseKEGG}}, TRUE for \code{\link[clusterProfiler]{GSEA}}
#' @param skipGO logical, set parameters as NULL for skip this step
#' @param skipKEGG logical, set parameters as NULL for skip this step
#'
#' @return gseParam
#' @export
#'
#' @examples
#' Create_gseParam()
Create_gseParam <- function(goParam = NULL,keggParam = NULL,
                            customGO = FALSE,customKEGG = FALSE,
                            skipGO = FALSE,skipKEGG = FALSE) {

  customKEGGParams = list(
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps  = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = NULL,
    TERM2NAME = NA,
    verbose = TRUE,
    seed = FALSE,
    by = 'fgsea'
  )

  customGOParams = list(
    ont           = "ALL",
    TERM2GENE = NULL,
    TERM2NAME = NA,
    organism = "UNKNOW",
    keyType = "SYMBOL",
    exponent      = 1,
    minGSSize     = 10,
    maxGSSize     = 500,
    eps           = 0,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    verbose       = TRUE,
    seed          = FALSE,
    by            = 'fgsea'
  )

  ## 设置默认选项 ----
  if (customGO) {
    goParam_default = customKEGGParams
  } else {

    goParam_default = list(
      ont = "ALL",
      OrgDb = NULL,
      keyType = "SYMBOL",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = TRUE,
      seed = FALSE,
      by = "fgsea")

  }

  if (customKEGG) {
    keggParam_default = customGOParams
  } else {

    keggParam_default = list(
      organism = NULL,
      keyType = "kegg",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 0,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = TRUE,
      use_internal_data = FALSE,
      seed = FALSE,
      by = "fgsea")

  }

  goParam = paramParse(defaultParams = goParam_default,newParam = goParam)
  keggParam = paramParse(defaultParams = keggParam_default,newParam = keggParam)

  if (skipGO) {
    goParam = NULL
  }

  if (skipKEGG) {
    keggParam = NULL
  }

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
#' @param skipMSigDB logical, set parameters as NULL for skip this step
#'
#' @return a list
#' @export
#'
#' @examples
#' Create_msigdbParam()
Create_msigdbParam <- function(msigdbParam = NULL,skipMSigDB = FALSE) {

  ## 设置默认选项 ----
  msigdbParam_default = list(
    species = NULL, category = NULL, subcategory = NULL)

  msigdbParam = paramParse(defaultParams = msigdbParam_default, newParam = msigdbParam)

  if (skipMSigDB) {
    msigdbParam = NULL
  }

  return(msigdbParam)

}

#' Create \code{msigdbGSEAparam}
#'
#' Create a msigdbHyperParam list. If NULL, default parameters will apply.
#'
#' @param GSEAparam a list contains any paramaters of \code{\link[clusterProfiler]{GSEA}} except \code{geneList}
#'
#' @return a list
#' @export
#'
#' @examples
#' Create_msigdbGSEAparam()
Create_msigdbGSEAparam <- function(GSEAparam = NULL) {

  GSEAparam_default = list(
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2NAME = NA,
    verbose = TRUE,
    seed = FALSE,
    by = "fgsea"
  )

  GSEAparam = paramParse(defaultParams = GSEAparam_default,newParam = GSEAparam)

  return(GSEAparam)

}

#' Create \code{msigdbHyperParam}
#'
#' Create a msigdbHyperParam list. If NULL, default parameters will apply.
#'
#' @param GSEAparam a list contains any paramaters of \code{\link[clusterProfiler]{enricher}} except \code{geneList}
#'
#' @return a list
#' @export
#'
#' @examples
#' Create_msigdbHyperParam()
Create_msigdbHyperParam <- function(HyperParam = NULL) {

  HyperParam_default = list(
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = NULL,
    minGSSize=10,
    maxGSSize=500,
    qvalueCutoff = 1,
    TERM2GENE = NULL, ## if custom* = true ,it must be set
    TERM2NAME = NA)

  HyperParam = paramParse(defaultParams = HyperParam_default,newParam = HyperParam)

  return(HyperParam)

}

## ----
## help function
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
