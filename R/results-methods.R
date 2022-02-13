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
#'
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

#' @importFrom usethis ui_stop ui_info ui_value
setValidity("treatInfo", function(object) {

  label =  paste0(label(object),sigCollapse = ';')
  if (!label_ns(object) %in% label(object)) {
    usethis::ui_stop("label_ns must be in one of label")
  } else {
    usethis::ui_done("label_ns:{ui_value(label_ns(object))} and label:{ui_value(label)} seems ok")
  }

})

# Methods for treatInfo ---------------------------------------------------
setGeneric(name="treatInfo", def=function(obj) standardGeneric("treatInfo"))
setMethod(f="treatInfo", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo)

setMethod(f="treatInfo", signature="degResults", definition=function(obj) obj@treatInfo)

setGeneric(name="label", def=function(obj) standardGeneric("label"))
setMethod(f="label", signature="treatInfo", definition=function(obj) obj@label)
setMethod(f="label", signature="degResults", definition=function(obj) obj@treatInfo@label)
setMethod(f="label", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@label)

setGeneric(name="label_ns", def=function(obj) standardGeneric("label_ns"))
setMethod(f="label_ns", signature="treatInfo", definition=function(obj) obj@label_ns)
setMethod(f="label_ns", signature="degResults", definition=function(obj) obj@treatInfo@label_ns)
setMethod(f="label_ns", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@label_ns)

setGeneric(name="cutFC", def=function(obj) standardGeneric("cutFC"))
setMethod(f="cutFC", signature="treatInfo", definition=function(obj) obj@cutFC)
setMethod(f="cutFC", signature="degResults", definition=function(obj) obj@treatInfo@cutFC)
setMethod(f="cutFC", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@cutFC)

setGeneric(name="cutFDR", def=function(obj) standardGeneric("cutFDR"))
setMethod(f="cutFDR", signature="treatInfo", definition=function(obj) obj@cutFDR)
setMethod(f="cutFDR", signature="degResults", definition=function(obj) obj@treatInfo@cutFDR)
setMethod(f="cutFDR", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@cutFDR)

setGeneric(name="sigCol", def=function(obj) standardGeneric("sigCol"))
setMethod(f="sigCol", signature="treatInfo", definition=function(obj) obj@sigCol)
setMethod(f="sigCol", signature="degResults", definition=function(obj) obj@treatInfo@sigCol)
setMethod(f="sigCol", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigCol)

setGeneric(name="sigAlpha", def=function(obj) standardGeneric("sigAlpha"))
setMethod(f="sigAlpha", signature="treatInfo", definition=function(obj) obj@sigAlpha)
setMethod(f="sigAlpha", signature="degResults", definition=function(obj) obj@treatInfo@sigAlpha)
setMethod(f="sigAlpha", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigAlpha)

setGeneric(name="sigSize", def=function(obj) standardGeneric("sigSize"))
setMethod(f="sigSize", signature="treatInfo", definition=function(obj) obj@sigSize)
setMethod(f="sigSize", signature="degResults", definition=function(obj) obj@treatInfo@sigSize)
setMethod(f="sigSize", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigSize)

setGeneric(name="sigShape", def=function(obj) standardGeneric("sigShape"))
setMethod(f="sigShape", signature="treatInfo", definition=function(obj) obj@sigShape)
setMethod(f="sigShape", signature="degResults", definition=function(obj) obj@treatInfo@sigShape)
setMethod(f="sigShape", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigShape)

setGeneric("treatInfo<-", function(obj, value) standardGeneric("treatInfo<-"))
setReplaceMethod("treatInfo", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo <- value; validObject(obj); obj})
setReplaceMethod("treatInfo", "degResults",
                 function(obj, value) {obj@treatInfo <- value; validObject(obj); obj})

setGeneric("cutFC<-", function(obj, value) standardGeneric("cutFC<-"))
setReplaceMethod("cutFC", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@cutFC <- value; validObject(obj); obj})
setReplaceMethod("cutFC", "degResults",
                 function(obj, value) {obj@treatInfo@cutFC <- value; validObject(obj); obj})
setReplaceMethod("cutFC", "treatInfo",
                 function(obj, value) {obj@cutFC <- value; validObject(obj); obj})

setGeneric("cutFDR<-", function(obj, value) standardGeneric("cutFDR<-"))
setReplaceMethod("cutFDR", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@cutFDR <- value; validObject(obj); obj})
setReplaceMethod("cutFDR", "degResults",
                 function(obj, value) {obj@treatInfo@cutFDR <- value; validObject(obj); obj})
setReplaceMethod("cutFDR", "treatInfo",
                 function(obj, value) {obj@cutFDR <- value; validObject(obj); obj})

setGeneric("label<-", function(obj, value) standardGeneric("label<-"))
setReplaceMethod("label", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@label <- value; validObject(obj); obj})
setReplaceMethod("label", "degResults",
                 function(obj, value) {obj@treatInfo@label <- value; validObject(obj); obj})
setReplaceMethod("label", "treatInfo",
                 function(obj, value) {obj@label <- value; validObject(obj); obj})

setGeneric("label_ns<-", function(obj, value) standardGeneric("label_ns<-"))
setReplaceMethod("label_ns", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@label_ns <- value; validObject(obj); obj})
setReplaceMethod("label_ns", "degResults",
                 function(obj, value) {obj@treatInfo@label_ns <- value; validObject(obj); obj})
setReplaceMethod("label_ns", "treatInfo",
                 function(obj, value) {obj@label_ns <- value; validObject(obj); obj})

setGeneric("sigCol<-", function(obj, value) standardGeneric("sigCol<-"))
setReplaceMethod("sigCol", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigCol <- value; validObject(obj); obj})
setReplaceMethod("sigCol", "degResults",
                 function(obj, value) {obj@treatInfo@sigCol <- value; validObject(obj); obj})
setReplaceMethod("sigCol", "treatInfo",
                 function(obj, value) {obj@sigCol <- value; validObject(obj); obj})

setGeneric("sigAlpha<-", function(obj, value) standardGeneric("sigAlpha<-"))
setReplaceMethod("sigAlpha", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigAlpha <- value; validObject(obj); obj})
setReplaceMethod("sigAlpha", "degResults",
                 function(obj, value) {obj@treatInfo@sigAlpha <- value; validObject(obj); obj})
setReplaceMethod("sigAlpha", "treatInfo",
                 function(obj, value) {obj@sigAlpha <- value; validObject(obj); obj})

setGeneric("sigSize<-", function(obj, value) standardGeneric("sigSize<-"))
setReplaceMethod("sigSize", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigSize <- value; validObject(obj); obj})
setReplaceMethod("sigSize", "degResults",
                 function(obj, value) {obj@treatInfo@sigSize <- value; validObject(obj); obj})
setReplaceMethod("sigSize", "treatInfo",
                 function(obj, value) {obj@sigSize <- value; validObject(obj); obj})

setGeneric("sigShape<-", function(obj, value) standardGeneric("sigShape<-"))
setReplaceMethod("sigShape", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigShape <- value; validObject(obj); obj})
setReplaceMethod("sigShape", "degResults",
                 function(obj, value) {obj@treatInfo@sigShape <- value; validObject(obj); obj})
setReplaceMethod("sigShape", "treatInfo",
                 function(obj, value) {obj@sigShape <- value; validObject(obj); obj})

# End of treatInfo methods ------------------------------------------------

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

# Methods for vsData ----
setGeneric(name="vsData", def=function(obj) standardGeneric("vsData"))
setMethod(f="vsData", signature="DEGContainer", definition=function(obj) obj@degResults@vsData)
setMethod(f="vsData", signature="degResults", definition=function(obj) obj@vsData)

setGeneric(name="limma_res", def=function(obj) standardGeneric("limma_res"))
setMethod(f="limma_res", signature="vsData", definition=function(obj) obj@limma_res)
setMethod(f="limma_res", signature="degResults", definition=function(obj) obj@vsData@limma_res)
setMethod(f="limma_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@limma_res)

setGeneric(name="edgeR_res", def=function(obj) standardGeneric("edgeR_res"))
setMethod(f="edgeR_res", signature="vsData", definition=function(obj) obj@edgeR_res)
setMethod(f="edgeR_res", signature="degResults", definition=function(obj) obj@vsData@edgeR_res)
setMethod(f="edgeR_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@edgeR_res)

setGeneric(name="DESeq2_res", def=function(obj) standardGeneric("DESeq2_res"))
setMethod(f="DESeq2_res", signature="vsData", definition=function(obj) obj@DESeq2_res)
setMethod(f="DESeq2_res", signature="degResults", definition=function(obj) obj@vsData@DESeq2_res)
setMethod(f="DESeq2_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@DESeq2_res)

setGeneric(name="merge_res", def=function(obj) standardGeneric("merge_res"))
setMethod(f="merge_res", signature="vsData", definition=function(obj) obj@merge_res)
setMethod(f="merge_res", signature="degResults", definition=function(obj) obj@vsData@merge_res)
setMethod(f="merge_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@merge_res)

setGeneric("vsData<-", function(obj, value) standardGeneric("vsData<-"))
setReplaceMethod("vsData", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData <- value; validObject(obj); obj})
setReplaceMethod("vsData", "degResults",
                 function(obj, value) {obj@vsData <- value; validObject(obj); obj})

setGeneric("DESeq2_res<-", function(obj, value) standardGeneric("DESeq2_res<-"))
setReplaceMethod("DESeq2_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@DESeq2_res <- value; validObject(obj); obj})
setReplaceMethod("DESeq2_res", "degResults",
                 function(obj, value) {obj@vsData@DESeq2_res <- value; validObject(obj); obj})
setReplaceMethod("DESeq2_res", "vsData",
                 function(obj, value) {obj@DESeq2_res <- value; validObject(obj); obj})

setGeneric("limma_res<-", function(obj, value) standardGeneric("limma_res<-"))
setReplaceMethod("limma_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@limma_res <- value; validObject(obj); obj})
setReplaceMethod("limma_res", "degResults",
                 function(obj, value) {obj@vsData@limma_res <- value; validObject(obj); obj})
setReplaceMethod("limma_res", "vsData",
                 function(obj, value) {obj@limma_res <- value; validObject(obj); obj})

setGeneric("edgeR_res<-", function(obj, value) standardGeneric("edgeR_res<-"))
setReplaceMethod("edgeR_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@edgeR_res <- value; validObject(obj); obj})
setReplaceMethod("edgeR_res", "degResults",
                 function(obj, value) {obj@vsData@edgeR_res <- value; validObject(obj); obj})
setReplaceMethod("edgeR_res", "vsData",
                 function(obj, value) {obj@edgeR_res <- value; validObject(obj); obj})

setGeneric("merge_res<-", function(obj, value) standardGeneric("merge_res<-"))
setReplaceMethod("merge_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@merge_res <- value; validObject(obj); obj})
setReplaceMethod("merge_res", "degResults",
                 function(obj, value) {obj@vsData@merge_res <- value; validObject(obj); obj})
setReplaceMethod("merge_res", "vsData",
                 function(obj, value) {obj@merge_res <- value; validObject(obj); obj})
# End of vsData methods ----

Create_degResults <- function(vsData,treatInfo = Create_treatInfo()) {

  new("degResults",
      vsData = vsData,
      treatInfo = treatInfo
  )

}

## Methods for degResults ----

## End of Methods for degResults ----

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

## Methods for hyper results ----
setGeneric(name="hyperKEGGparam", def=function(obj) standardGeneric("hyperKEGGparam"))
setMethod(f="hyperKEGGparam", signature="hyperParam", definition=function(obj) obj@keggParam)
setMethod(f="hyperKEGGparam", signature="hyperResults", definition=function(obj) obj@hyperParam@keggParam)
setMethod(f="hyperKEGGparam", signature="DEGContainer", definition=function(obj) obj@hyperResults@hyperParam@keggParam)

setGeneric("hyperKEGGparam<-", function(obj, value) standardGeneric("hyperKEGGparam<-"))
setReplaceMethod("hyperKEGGparam", "DEGContainer",
                 function(obj, value) {obj@hyperResults@hyperParam@keggParam <- value; validObject(obj); obj})
setReplaceMethod("hyperKEGGparam", "hyperResults",
                 function(obj, value) {obj@hyperParam@keggParam <- value; validObject(obj); obj})
setReplaceMethod("hyperKEGGparam", "hyperParam",
                 function(obj, value) {obj@keggParam <- value; validObject(obj); obj})



setGeneric(name="hyperGOparam", def=function(obj) standardGeneric("hyperGOparam"))
setMethod(f="hyperGOparam", signature="hyperParam", definition=function(obj) obj@goParam)
setMethod(f="hyperGOparam", signature="hyperResults", definition=function(obj) obj@hyperParam@goParam)
setMethod(f="hyperGOparam", signature="DEGContainer", definition=function(obj) obj@hyperResults@hyperParam@goParam)

setGeneric("hyperGOparam<-", function(obj, value) standardGeneric("hyperGOparam<-"))
setReplaceMethod("hyperGOparam", "DEGContainer",
                 function(obj, value) {obj@hyperResults@hyperParam@goParam <- value; validObject(obj); obj})
setReplaceMethod("hyperGOparam", "hyperResults",
                 function(obj, value) {obj@hyperParam@goParam <- value; validObject(obj); obj})
setReplaceMethod("hyperGOparam", "hyperParam",
                 function(obj, value) {obj@goParam <- value; validObject(obj); obj})

setGeneric(name="hyperRes", def=function(obj) standardGeneric("hyperRes"))
# setMethod(f="hyperRes", signature="hyperParam", definition=function(obj) obj@keggParam)
setMethod(f="hyperRes", signature="hyperResults", definition=function(obj) obj@hyperRes)
setMethod(f="hyperRes", signature="DEGContainer", definition=function(obj) obj@hyperResults@hyperRes)

setGeneric("hyperRes<-", function(obj, value) standardGeneric("hyperRes<-"))
setReplaceMethod("hyperRes", "DEGContainer",
                 function(obj, value) {obj@hyperResults@hyperRes <- value; validObject(obj); obj})
setReplaceMethod("hyperRes", "hyperResults",
                 function(obj, value) {obj@hyperRes <- value; validObject(obj); obj})
## End of methods for hyper results ----

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

## Methods for gse results ----
setGeneric(name="gseKEGGparam", def=function(obj) standardGeneric("gseKEGGparam"))
setMethod(f="gseKEGGparam", signature="gseParam", definition=function(obj) obj@keggParam)
setMethod(f="gseKEGGparam", signature="gseResults", definition=function(obj) obj@gseParam@keggParam)
setMethod(f="gseKEGGparam", signature="DEGContainer", definition=function(obj) obj@gseResults@gseParam@keggParam)

setGeneric("gseKEGGparam<-", function(obj, value) standardGeneric("gseKEGGparam<-"))
setReplaceMethod("gseKEGGparam", "DEGContainer",
                 function(obj, value) {obj@gseResults@gseParam@keggParam <- value; validObject(obj); obj})
setReplaceMethod("gseKEGGparam", "gseResults",
                 function(obj, value) {obj@gseParam@keggParam <- value; validObject(obj); obj})
setReplaceMethod("gseKEGGparam", "gseParam",
                 function(obj, value) {obj@keggParam <- value; validObject(obj); obj})



setGeneric(name="gseGOparam", def=function(obj) standardGeneric("gseGOparam"))
setMethod(f="gseGOparam", signature="gseParam", definition=function(obj) obj@goParam)
setMethod(f="gseGOparam", signature="gseResults", definition=function(obj) obj@gseParam@goParam)
setMethod(f="gseGOparam", signature="DEGContainer", definition=function(obj) obj@gseResults@gseParam@goParam)

setGeneric("gseGOparam<-", function(obj, value) standardGeneric("gseGOparam<-"))
setReplaceMethod("gseGOparam", "DEGContainer",
                 function(obj, value) {obj@gseResults@gseParam@goParam <- value; validObject(obj); obj})
setReplaceMethod("gseGOparam", "gseResults",
                 function(obj, value) {obj@gseParam@goParam <- value; validObject(obj); obj})
setReplaceMethod("gseGOparam", "gseParam",
                 function(obj, value) {obj@goParam <- value; validObject(obj); obj})

setGeneric(name="gseRes", def=function(obj) standardGeneric("gseRes"))
# setMethod(f="gseRes", signature="gseParam", definition=function(obj) obj@keggParam)
setMethod(f="gseRes", signature="gseResults", definition=function(obj) obj@gseRes)
setMethod(f="gseRes", signature="DEGContainer", definition=function(obj) obj@gseResults@gseRes)

setGeneric("gseRes<-", function(obj, value) standardGeneric("gseRes<-"))
setReplaceMethod("gseRes", "DEGContainer",
                 function(obj, value) {obj@gseResults@gseRes <- value; validObject(obj); obj})
setReplaceMethod("gseRes", "gseResults",
                 function(obj, value) {obj@gseRes <- value; validObject(obj); obj})
## End of methods for gse results ----

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

Create_MSigDB <- function(msigdbParam = Create_msigdbParam(),msigdbGSEAparam = Create_msigdbGSEAparam()) {

  new("MSigDB",
      msigdbParam = msigdbParam,
      msigdbData = list(),
      msigdbGSEAparam = msigdbGSEAparam,
      msigdbGSEAresult = list(),
      msigdbGSVAresult = list())

}

## Methods for MSigDB ----
setGeneric(name="msigdbParam", def=function(obj) standardGeneric("msigdbParam"))
setMethod(f="msigdbParam", signature="MSigDB", definition=function(obj) obj@msigdbParam)
setMethod(f="msigdbParam", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbParam)

setGeneric(name="msigdbGSEAparam", def=function(obj) standardGeneric("msigdbGSEAparam"))
setMethod(f="msigdbGSEAparam", signature="MSigDB", definition=function(obj) obj@msigdbGSEAparam)
setMethod(f="msigdbGSEAparam", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbGSEAparam)

setGeneric("msigdbParam<-", function(obj, value) standardGeneric("msigdbParam<-"))
setReplaceMethod("msigdbParam", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbParam <- value; validObject(obj); obj})
setReplaceMethod("msigdbParam", "MSigDB",
                 function(obj, value) {obj@msigdbParam <- value; validObject(obj); obj})

setGeneric("msigdbGSEAparam<-", function(obj, value) standardGeneric("msigdbGSEAparam<-"))
setReplaceMethod("msigdbGSEAparam", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbGSEAparam <- value; validObject(obj); obj})
setReplaceMethod("msigdbGSEAparam", "MSigDB",
                 function(obj, value) {obj@msigdbGSEAparam <- value; validObject(obj); obj})

setGeneric(name="msigdbData", def=function(obj) standardGeneric("msigdbData"))
setMethod(f="msigdbData", signature="MSigDB", definition=function(obj) obj@msigdbData)
setMethod(f="msigdbData", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbData)

setGeneric("msigdbData<-", function(obj, value) standardGeneric("msigdbData<-"))
setReplaceMethod("msigdbData", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbData <- value; validObject(obj); obj})
setReplaceMethod("msigdbData", "MSigDB",
                 function(obj, value) {obj@msigdbData <- value; validObject(obj); obj})

setGeneric(name="msigdbGSEAresult", def=function(obj) standardGeneric("msigdbGSEAresult"))
setMethod(f="msigdbGSEAresult", signature="MSigDB", definition=function(obj) obj@msigdbGSEAresult)
setMethod(f="msigdbGSEAresult", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbGSEAresult)

setGeneric("msigdbGSEAresult<-", function(obj, value) standardGeneric("msigdbGSEAresult<-"))
setReplaceMethod("msigdbGSEAresult", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbGSEAresult <- value; validObject(obj); obj})
setReplaceMethod("msigdbGSEAresult", "MSigDB",
                 function(obj, value) {obj@msigdbGSEAresult <- value; validObject(obj); obj})

setGeneric(name="msigdbGSVAresult", def=function(obj) standardGeneric("msigdbGSVAresult"))
setMethod(f="msigdbGSVAresult", signature="MSigDB", definition=function(obj) obj@msigdbGSVAresult)
setMethod(f="msigdbGSVAresult", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbGSVAresult)

setGeneric("msigdbGSVAresult<-", function(obj, value) standardGeneric("msigdbGSVAresult<-"))
setReplaceMethod("msigdbGSVAresult", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbGSVAresult <- value; validObject(obj); obj})
setReplaceMethod("msigdbGSVAresult", "MSigDB",
                 function(obj, value) {obj@msigdbGSVAresult <- value; validObject(obj); obj})
## End of methods for gse results ----
