#' Function for treatInfo
#'
#' @param cutFC a numeric or list. Threshold of logFC for finding significantly different genes.
#' @param cutFDR a numeric. Threshold of Pvalue for finding significantly different genes.
#' @param label a character vector. Name for every group
#' @param label_ns a character. Name of no signification group
#' @param sigCol a character vector. sigColor of every group
#' @param sigAlpha a numeric vector. transparency for every group
#'
#' @return
#' @export
Create_treatInfo <- function(cutFC = NULL,
                      cutFDR = 0.05,
                      label = c("Down","Stable","Up"),
                      label_ns = "Stable",
                      sigCol = c("#2874C5", "grey", "#f87669"),
                      sigAlpha = c(1,0.5,1),
                      sigSize = c(0.6,0.5,0.6),
                      sigShape = 16) {

  new("treatInfo",
      cutFC = cutFC,
      cutFDR = cutFDR,
      label = label,
      label_ns = label_ns,
      sigCol = sigCol,
      sigAlpha = sigAlpha,
      sigSize = sigSize,
      sigShape = sigShape
      )

}

Create_vsData <- function(limma,edgeR,DESeq2) {

  limma <- data.frame(matrix(ncol = 0, nrow = 0))
  # colnames(limma) <- c("logFC","AveExpr","t","P.Value","adj.P.Val","B")

  edgeR <- data.frame(matrix(ncol = 0, nrow = 0))
  # colnames(edgeR) <- c("logFC","logCPM","LR","PValue","FDR")

  DESeq2 <- data.frame(matrix(ncol = 0, nrow = 0))
  # colnames(DESeq2) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

  new("vsData",
      limma_res = limma,
      edgeR_res = edgeR,
      DESeq2_res = DESeq2)

}

Create_degResults <- function(vsData,treatInfo = Create_treatInfo()) {

  new("degResults",
      vsData = vsData,
      treatInfo = treatInfo
  )

}

Create_hyperParam <- function(goParam = NULL,keggParam = NULL) {

  ## 设置默认选项 ----
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

  if (is.null(goParam))
    goParam <- goParam_default

  if (is.null(keggParam))
    keggParam <- keggParam_default

  # if (is.null(keggParam[['organism']]))
  #   keggParam[['organism']] = switch(keggParam[["OrgDb"]],
  #                                    'org.Hs.eg.db' = 'hsa',
  #                                    'org.Mm,eg.db' = 'mmu'
  #   )
  #
  # if (is.null(keggParam[['down_label']]))
  #   keggParam[['down_label']] = keggParam[['label']][1]
  ## ----

  ## 读取用户选项 ----
  ## Go
  if (is.list(goParam)&length(goParam) > 0) { ## 判断格式是否正确

    if (!all(names(goParam) %in% names(goParam_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the parameters of {usethis::ui_code('Create_hyperParam')} ")

    } else if (length(goParam) == 1) { ## 判断 长度为1

      goParam_default[[names(goParam)]] <- goParam[names(goParam)][[1]]

      goParam <- goParam_default

    } else if (length(goParam) > 1) { ## 其他长度

      for (x in seq_along(goParam)) {

        goParam_default[[names(goParam)[x]]] <- goParam[names(goParam)[x]][[1]]

      }

      goParam <- goParam_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the usage of {usethis::ui_code('Create_hyperParam')} ")

  }

  ## KEGG
  if (is.list(keggParam)&length(keggParam) > 0) { ## 判断格式是否正确

    if (!all(names(keggParam) %in% names(keggParam_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the parameters of {usethis::ui_code('Create_hyperParam')} ")

    } else if (length(keggParam) == 1) { ## 判断 长度为1

      keggParam_default[[names(keggParam)]] <- keggParam[names(keggParam)][[1]]

      keggParam <- keggParam_default

    } else if (length(keggParam) > 1) { ## 其他长度

      for (x in seq_along(keggParam)) {

        keggParam_default[[names(keggParam)[x]]] <- keggParam[names(keggParam)[x]][[1]]

      }

      keggParam <- keggParam_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the usage of {usethis::ui_code('Create_hyperParam')} ")

  }
  ## ----

  new("hyperParam",
      goParam = goParam,
      keggParam = keggParam)

}

Create_hyperResults <- function(hyperParam=Create_hyperParam()) {

  new("hyperResults",
      hyperRes = list(),
      hyperParam = hyperParam
  )

}

Create_gseParam <- function(goParam = NULL,keggParam = NULL) {

  ## 设置默认选项 ----
  goParam_default = list(
    ont = "ALL",
    OrgDb = NULL,
    keyType = "SYMBOL",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10,
    pvalueCutoff = 0.99,
    pAdjustMethod = "BH",
    verbose = TRUE,
    seed = FALSE,
    by = "fgsea")

  keggParam_default = list(
    organism = NULL,
    keyType = "kegg",
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-10,
    pvalueCutoff = 0.99,
    pAdjustMethod = "BH",
    verbose = TRUE,
    use_internal_data = FALSE,
    seed = FALSE,
    by = "fgsea")

  if (is.null(goParam))
    goParam <- goParam_default

  if (is.null(keggParam))
    keggParam <- keggParam_default

  # if (is.null(keggParam[['organism']]))
  #   keggParam[['organism']] = switch(keggParam[["OrgDb"]],
  #                                    'org.Hs.eg.db' = 'hsa',
  #                                    'org.Mm,eg.db' = 'mmu'
  #   )
  #
  # if (is.null(keggParam[['down_label']]))
  #   keggParam[['down_label']] = keggParam[['label']][1]
  ## ----

  ## 读取用户选项 ----
  ## Go
  if (is.list(goParam)&length(goParam) > 0) { ## 判断格式是否正确

    if (!all(names(goParam) %in% names(goParam_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the parameters of {usethis::ui_code('gseGO')} ")

    } else if (length(goParam) == 1) { ## 判断 长度为1

      goParam_default[[names(goParam)]] <- goParam[names(goParam)][[1]]

      goParam <- goParam_default

    } else if (length(goParam) > 1) { ## 其他长度

      for (x in seq_along(goParam)) {

        goParam_default[[names(goParam)[x]]] <- goParam[names(goParam)[x]][[1]]

      }

      goParam <- goParam_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the usage of {usethis::ui_code('gseGO')} ")

  }

  ## KEGG
  if (is.list(keggParam)&length(keggParam) > 0) { ## 判断格式是否正确

    if (!all(names(keggParam) %in% names(keggParam_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('keggParam')} will not be applied.\n Please check the parameters of {usethis::ui_code('gseKEGG')} ")

    } else if (length(keggParam) == 1) { ## 判断 长度为1

      keggParam_default[[names(keggParam)]] <- keggParam[names(keggParam)][[1]]

      keggParam <- keggParam_default

    } else if (length(keggParam) > 1) { ## 其他长度

      for (x in seq_along(keggParam)) {

        keggParam_default[[names(keggParam)[x]]] <- keggParam[names(keggParam)[x]][[1]]

      }

      keggParam <- keggParam_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('keggParam')} will not be applied.\n Please check the usage of {usethis::ui_code('gseKEGG')} ")

  }
  ## ----

  new("gseParam",
      goParam = goParam,
      keggParam = keggParam)

}

Create_gseResults <- function(gseParam=Create_gseParam()) {

  new("gseResults",
      gseRes = list(),
      gseParam = gseParam
  )

}

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
    species = "Homo sapiens", category = NULL, subcategory = NULL)

  if (is.null(msigdbParam))
    msigdbParam <- msigdbParam_default

  ## 读取用户选项 ----
  if (is.list(msigdbParam)&length(msigdbParam) > 0) { ## 判断格式是否正确

    if (!all(names(msigdbParam) %in% names(msigdbParam_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('msigdbParam')} will not be applied.\n Please check the parameters of {usethis::ui_code('msigdbr')} ")

    } else if (length(msigdbParam) == 1) { ## 判断 长度为1

      msigdbParam_default[[names(msigdbParam)]] <- msigdbParam[names(msigdbParam)][[1]]

      msigdbParam <- msigdbParam_default

    } else if (length(msigdbParam) > 1) { ## 其他长度

      for (x in seq_along(msigdbParam)) {

        msigdbParam_default[[names(msigdbParam)[x]]] <- msigdbParam[names(msigdbParam)[x]][[1]]

      }

      msigdbParam <- msigdbParam_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('msigdbParam')} will not be applied.\n Please check the usage of {usethis::ui_code('msigdbr')} ")

  }

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

  if (is.null(GSEAparam))
    GSEAparam <- GSEAparam_default

  ## 读取用户选项 ----
  if (is.list(GSEAparam)&length(GSEAparam) > 0) { ## 判断格式是否正确

    if (!all(names(GSEAparam) %in% names(GSEAparam_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('GSEAparam')} will not be applied.\n Please check the parameters of {usethis::ui_code('GSEA')} ")

    } else if (length(GSEAparam) == 1) { ## 判断 长度为1

      GSEAparam_default[[names(GSEAparam)]] <- GSEAparam[names(GSEAparam)][[1]]

      GSEAparam <- GSEAparam_default

    } else if (length(GSEAparam) > 1) { ## 其他长度

      for (x in seq_along(GSEAparam)) {

        GSEAparam_default[[names(GSEAparam)[x]]] <- GSEAparam[names(GSEAparam)[x]][[1]]

      }

      GSEAparam <- GSEAparam_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('goParam')} will not be applied.\n Please check the usage of {usethis::ui_code('GSEA')} ")

  }

  return(GSEAparam)

}

#' Create \code{MSigDB}
#'
#' Create a MSigDB object.
#'
#' @param msigdbParam from \code{Create_msigdbParam}
#' @param msigdbGSEAparam from \code{Create_msigdbGSEAparam}
#'
#' @return a \code{MSigDB} object
#' @export
Create_MSigDB <- function(msigdbParam = Create_msigdbParam(),msigdbGSEAparam = Create_msigdbGSEAparam()) {

  new("MSigDB",
      msigdbParam = msigdbParam,
      msigdbData = list(),
      msigdbGSEAparam = msigdbGSEAparam,
      msigdbGSEAresult = list(),
      msigdbGSVAresult = list())

}