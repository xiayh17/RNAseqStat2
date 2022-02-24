#' create \code{dataInfo}
#'
#' a function to create \code{dataInfo} OB
#'
#' @param species species for your data. `Human` or `Mouse`
#' @param dataType kind of expresses value matrix. `Counts` (Integer) or `Array` (Decimal).
#' @param idType kind of gene id. `ENSEMBL` or `SYMBOL`, If `ENSEMBL`, it will be automatically converted to `SYMBOL`.
#' @param expMatrix expresses value matrix. Should be a data.frame row named by gene ID and column named by Sample
#' @param groupInfo a Character Vectors ordered by samples in matrix.
#' @param caseGroup a Character names of case group.
#' @param filterMethod a function used to filter expresses value matrix.
#' @param matrixFiltered a expresses value matrix after filter by `filterMethod`.
#'
#' @importFrom methods new
#'
#' @return \code{dataInfo}
Create_dataInfo <- function(species,
                     dataType,
                     idType,
                     expMatrix,
                     groupInfo,
                     caseGroup,
                     filterMethod,
                     matrixFiltered) {

  species = as.character(species)
  dataType = as.character(dataType)
  idType = as.character(idType)
  expMatrix <- as.data.frame(expMatrix)
  groupInfo = as.character(groupInfo)
  caseGroup = as.character(caseGroup)
  filterMethod = as.character(filterMethod)
  matrixFiltered <- as.data.frame(matrixFiltered)


  new("dataInfo",
      species = species,
      dataType = dataType,
      idType = idType,
      expMatrix = expMatrix,
      groupInfo = groupInfo,
      caseGroup = caseGroup,
      filterMethod = filterMethod,
      matrixFiltered = matrixFiltered
  )

}

# Low level function to create DEGContainer -------------------------------
#' @title Initialize an object of class \code{DEGContainer}
#'
#' @name newDEGContainer
#' @docType methods
#'
#' @description Constructs a \code{DEGContainer} object. Additional helper
#'   methods for manipulating \code{DEGContainer} objects are  also
#'   described below. We now recommend using
#'   \code{Create_DEGContainer} to create objects, instead.
#'
#' @param dataInfo a dataInfo object
#' @param degResults a degResults object
#' @param hyperResults a hyperResults object
#' @param gseResults a gseResults object
#' @param MSigDB a MSigDB object
#'
#' @return a DEGContainer
#'
#' @export
newDEGContainer <- function(dataInfo=NULL,
                         degResults=NULL,
                         hyperResults=NULL,
                         gseResults=NULL,
                         MSigDB = NULL)
{

  new(
    "DEGContainer",
    dataInfo=dataInfo,
    degResults=degResults,
    hyperResults=hyperResults,
    gseResults=gseResults,
    MSigDB = MSigDB
  )

}

## results OB create
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

## Create zero data frame for start
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

#' @export
Create_degResults <- function(vsData,treatInfo = Create_treatInfo()) {

  new("degResults",
      vsData = vsData,
      treatInfo = treatInfo
  )

}

Create_hyperResults <- function(hyperParam=Create_hyperParam()) {

  new("hyperResults",
      hyperRes = list(),
      hyperParam = hyperParam
  )

}

Create_gseResults <- function(gseParam=Create_gseParam()) {

  new("gseResults",
      gseRes = list(),
      gseParam = gseParam
  )

}

#' Create \code{MSigDB}
#'
#' Create a MSigDB object.
#'
#' @param msigdbParam from \code{Create_msigdbParam}
#' @param msigdbGSEAparam from \code{Create_msigdbGSEAparam}
#' @param msigdbHyperParam from \code{Create_msigdbHyperParam}
#'
#' @return a \code{MSigDB} object
#' @export
Create_MSigDB <- function(msigdbParam = Create_msigdbParam(),
                          msigdbGSEAparam = Create_msigdbGSEAparam(),
                          msigdbHyperParam = Create_msigdbHyperParam()) {

  new("MSigDB",
      msigdbParam = msigdbParam,
      msigdbData = list(),
      msigdbGSEAparam = msigdbGSEAparam,
      msigdbHyperParam = msigdbHyperParam,
      msigdbGSEAresult = list(),
      msigdbHyperResult = list(),
      msigdbGSVAresult = list(),
      msigdbTreat = Create_treatInfo(cutFDR = 1))

}
