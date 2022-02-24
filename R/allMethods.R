# methods for DEGContainer ------------------------------------------------
#' @importFrom usethis ui_line ui_value
setMethod("show", "DEGContainer",
          function(object) {

            ui_line("DEGContainer
                     {ui_value(dataType(object))} Matrix Containing {ui_value(ncol(expMatrix(object)))} samples")
          }
)
#' @export
setMethod(f="dataInfo", signature="DEGContainer", definition=function(obj) obj@dataInfo)
#' @export
setMethod(f="degResults", signature="DEGContainer", definition=function(obj) obj@degResults)
#' @export
setMethod(f="hyperResults", signature="DEGContainer", definition=function(obj) obj@hyperResults)
#' @export
setMethod(f="gseResults", signature="DEGContainer", definition=function(obj) obj@gseResults)
#' @export
setMethod(f="MSigDB", signature="DEGContainer", definition=function(obj) obj@MSigDB)
# methods for DEGContainer ------------------------------------------------

# methods for dataInfo ------------------------------------------------
#' @export
setMethod(f="species", signature="DEGContainer", definition=function(obj) obj@dataInfo@species)

#' @export
setMethod(f="dataType", signature="DEGContainer", definition=function(obj) obj@dataInfo@dataType)

#' @export
setMethod(f="idType", signature="DEGContainer", definition=function(obj) obj@dataInfo@idType)

#' @export
setMethod(f="expMatrix", signature="DEGContainer",definition=function(obj) obj@dataInfo@expMatrix)

#' @export
setMethod("groupInfo", "DEGContainer", function(obj) obj@dataInfo@groupInfo)

#' @export
setMethod(f="caseGroup", signature="DEGContainer",definition=function(obj) obj@dataInfo@caseGroup)

#' @export
setMethod(f="filterMethod", signature="DEGContainer", definition=function(obj) obj@dataInfo@filterMethod)

#' @export
setMethod(f="matrixFiltered", signature="DEGContainer", definition=function(obj) obj@dataInfo@matrixFiltered)

#' @export
setMethod(f="geneNames", signature="DEGContainer", definition=function(obj,filtered = FALSE) {
  if (filtered) {
    genes_name = rownames(obj@dataInfo@matrixFiltered)
  } else {
    genes_name = rownames(obj@dataInfo@expMatrix)
  }
  return(genes_name)})


#' @export
setMethod(f="sampleNames", signature="DEGContainer", definition=function(obj,filtered = FALSE) {
  if (filtered) {
    sample_names = colnames(obj@dataInfo@matrixFiltered)
  } else {
    sample_names = colnames(obj@dataInfo@expMatrix)
  }
  return(sample_names)})

#' @export
setMethod(f="species", signature="dataInfo", definition=function(obj) obj@species)

#' @export
setMethod(f="dataType", signature="dataInfo", definition=function(obj) obj@dataType)

#' @export
setMethod(f="idType", signature="dataInfo", definition=function(obj) obj@idType)

#' @export
setMethod(f="expMatrix", signature="dataInfo", definition=function(obj) obj@expMatrix)

#' @export
setMethod("groupInfo", "dataInfo", function(obj) obj@groupInfo)

#' @export
setMethod("caseGroup", "dataInfo", function(obj) obj@caseGroup)

#' @export
setMethod(f="filterMethod", signature="dataInfo", definition=function(obj) obj@filterMethod)

#' @export
setMethod(f="matrixFiltered", signature="dataInfo", definition=function(obj) obj@matrixFiltered)

#' @export
setMethod(f="sampleNames", signature="dataInfo", definition=function(obj,filtered = FALSE) {

  if (filtered) {
    sample_names = colnames(obj@matrixFiltered)
  } else {
    sample_names = colnames(obj@expMatrix)
  }
  return(sample_names)

})

## setter
#' @export
setReplaceMethod("species", "dataInfo",
                 function(obj, value) {obj@species <- value; validObject(obj); obj})

#' @export
setReplaceMethod("dataType", "dataInfo",
                 function(obj, value) {obj@dataType <- value; validObject(obj); obj})

#' @export
setReplaceMethod("idType", "dataInfo",
                 function(obj, value) {obj@idType <- value; validObject(obj); obj})

#' @export
setReplaceMethod("expMatrix", "dataInfo",
                 function(obj, value) {obj@expMatrix <- value; validObject(obj); obj})

#' @export
setReplaceMethod("groupInfo", "dataInfo",
                 function(obj, value) {obj@groupInfo <- value; validObject(obj); obj})

#' @export
setReplaceMethod("caseGroup", "dataInfo",
                 function(obj, value) {obj@caseGroup <- value; validObject(obj); obj})

#' @export
setReplaceMethod("filterMethod", "dataInfo",
                 function(obj, value) {obj@filterMethod <- value; validObject(obj); obj})

#' @export
setReplaceMethod("matrixFiltered", "dataInfo",
                 function(obj, value) {obj@matrixFiltered <- value; validObject(obj); obj})

#' @export
setReplaceMethod("geneNames", "dataInfo",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     rownames(obj@matrixFiltered) <- value;
                   } else {
                     rownames(obj@expMatrix) <- value;
                   };

                   validObject(obj); obj})

#' @export
setReplaceMethod("sampleNames", "dataInfo",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     colnames(obj@matrixFiltered) <- value;
                   } else {
                     colnames(obj@expMatrix) <- value;
                   };

                   validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("species", "DEGContainer",
                 function(obj, value) {obj@dataInfo@species <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("dataType", "DEGContainer",
                 function(obj, value) {obj@dataInfo@dataType <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("idType", "DEGContainer",
                 function(obj, value) {obj@dataInfo@idType <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("expMatrix", "DEGContainer",
                 function(obj, value) {obj@dataInfo@expMatrix <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("groupInfo", "DEGContainer",
                 function(obj, value) {obj@dataInfo@groupInfo <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("caseGroup", "DEGContainer",
                 function(obj, value) {obj@dataInfo@caseGroup <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("filterMethod", "DEGContainer",
                 function(obj, value) {obj@dataInfo@filterMethod <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("matrixFiltered", "DEGContainer",
                 function(obj, value) {obj@dataInfo@matrixFiltered <- value; validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("geneNames", "DEGContainer",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     rownames(obj@dataInfo@matrixFiltered) <- value;
                   } else {
                     rownames(obj@dataInfo@expMatrix) <- value;
                   };

                   validObject(obj); obj})
#' @rdname DEGContainer
#' @export
setReplaceMethod("sampleNames", "DEGContainer",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     colnames(obj@dataInfo@matrixFiltered) <- value;
                   } else {
                     colnames(obj@dataInfo@expMatrix) <- value;
                   };

                   validObject(obj); obj})
# methods for dataInfo ------------------------------------------------

# Methods for treatInfo ---------------------------------------------------
#' @export
setMethod(f="treatInfo", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo)
#' @export
setMethod(f="treatInfo", signature="degResults", definition=function(obj) obj@treatInfo)

#' @export
setMethod(f="label", signature="treatInfo", definition=function(obj) obj@label)
#' @export
setMethod(f="label", signature="degResults", definition=function(obj) obj@treatInfo@label)
#' @export
setMethod(f="label", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@label)

#' @export
setMethod(f="label_ns", signature="treatInfo", definition=function(obj) obj@label_ns)
#' @export
setMethod(f="label_ns", signature="degResults", definition=function(obj) obj@treatInfo@label_ns)
#' @export
setMethod(f="label_ns", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@label_ns)

#' @export
setMethod(f="cutFC", signature="treatInfo", definition=function(obj) obj@cutFC)
#' @export
setMethod(f="cutFC", signature="degResults", definition=function(obj) obj@treatInfo@cutFC)
#' @export
setMethod(f="cutFC", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@cutFC)

#' @export
setMethod(f="cutFDR", signature="treatInfo", definition=function(obj) obj@cutFDR)
#' @export
setMethod(f="cutFDR", signature="degResults", definition=function(obj) obj@treatInfo@cutFDR)
#' @export
setMethod(f="cutFDR", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@cutFDR)

#' @export
setMethod(f="sigCol", signature="treatInfo", definition=function(obj) obj@sigCol)
#' @export
setMethod(f="sigCol", signature="degResults", definition=function(obj) obj@treatInfo@sigCol)
#' @export
setMethod(f="sigCol", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigCol)

#' @export
setMethod(f="sigAlpha", signature="treatInfo", definition=function(obj) obj@sigAlpha)
#' @export
setMethod(f="sigAlpha", signature="degResults", definition=function(obj) obj@treatInfo@sigAlpha)
#' @export
setMethod(f="sigAlpha", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigAlpha)

#' @export
setMethod(f="sigSize", signature="treatInfo", definition=function(obj) obj@sigSize)
#' @export
setMethod(f="sigSize", signature="degResults", definition=function(obj) obj@treatInfo@sigSize)
#' @export
setMethod(f="sigSize", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigSize)

#' @export
setMethod(f="sigShape", signature="treatInfo", definition=function(obj) obj@sigShape)
#' @export
setMethod(f="sigShape", signature="degResults", definition=function(obj) obj@treatInfo@sigShape)
#' @export
setMethod(f="sigShape", signature="DEGContainer", definition=function(obj) obj@degResults@treatInfo@sigShape)

#' @export
setReplaceMethod("treatInfo", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo <- value; validObject(obj); obj})
#' @export
setReplaceMethod("treatInfo", "degResults",
                 function(obj, value) {obj@treatInfo <- value; validObject(obj); obj})

#' @export
setReplaceMethod("cutFC", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@cutFC <- value; validObject(obj); obj})
#' @export
setReplaceMethod("cutFC", "degResults",
                 function(obj, value) {obj@treatInfo@cutFC <- value; validObject(obj); obj})
#' @export
setReplaceMethod("cutFC", "treatInfo",
                 function(obj, value) {obj@cutFC <- value; validObject(obj); obj})

#' @export
setReplaceMethod("cutFDR", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@cutFDR <- value; validObject(obj); obj})
#' @export
setReplaceMethod("cutFDR", "degResults",
                 function(obj, value) {obj@treatInfo@cutFDR <- value; validObject(obj); obj})
#' @export
setReplaceMethod("cutFDR", "treatInfo",
                 function(obj, value) {obj@cutFDR <- value; validObject(obj); obj})

#' @export
setReplaceMethod("label", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@label <- value; validObject(obj); obj})
#' @export
setReplaceMethod("label", "degResults",
                 function(obj, value) {obj@treatInfo@label <- value; validObject(obj); obj})
#' @export
setReplaceMethod("label", "treatInfo",
                 function(obj, value) {obj@label <- value; validObject(obj); obj})

#' @export
setReplaceMethod("label_ns", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@label_ns <- value; validObject(obj); obj})
#' @export
setReplaceMethod("label_ns", "degResults",
                 function(obj, value) {obj@treatInfo@label_ns <- value; validObject(obj); obj})
#' @export
setReplaceMethod("label_ns", "treatInfo",
                 function(obj, value) {obj@label_ns <- value; validObject(obj); obj})

#' @export
setReplaceMethod("sigCol", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigCol <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigCol", "degResults",
                 function(obj, value) {obj@treatInfo@sigCol <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigCol", "treatInfo",
                 function(obj, value) {obj@sigCol <- value; validObject(obj); obj})

#' @export
setReplaceMethod("sigAlpha", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigAlpha <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigAlpha", "degResults",
                 function(obj, value) {obj@treatInfo@sigAlpha <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigAlpha", "treatInfo",
                 function(obj, value) {obj@sigAlpha <- value; validObject(obj); obj})

#' @export
setReplaceMethod("sigSize", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigSize <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigSize", "degResults",
                 function(obj, value) {obj@treatInfo@sigSize <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigSize", "treatInfo",
                 function(obj, value) {obj@sigSize <- value; validObject(obj); obj})

#' @export
setReplaceMethod("sigShape", "DEGContainer",
                 function(obj, value) {obj@degResults@treatInfo@sigShape <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigShape", "degResults",
                 function(obj, value) {obj@treatInfo@sigShape <- value; validObject(obj); obj})
#' @export
setReplaceMethod("sigShape", "treatInfo",
                 function(obj, value) {obj@sigShape <- value; validObject(obj); obj})
# Methods for treatInfo ---------------------------------------------------

# Methods for vsData ----
#' @export
setMethod(f="vsData", signature="DEGContainer", definition=function(obj) obj@degResults@vsData)
#' @export
setMethod(f="vsData", signature="degResults", definition=function(obj) obj@vsData)


#' @export
setMethod(f="limma_res", signature="vsData", definition=function(obj) obj@limma_res)
#' @export
setMethod(f="limma_res", signature="degResults", definition=function(obj) obj@vsData@limma_res)
#' @export
setMethod(f="limma_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@limma_res)

#' @export
setMethod(f="edgeR_res", signature="vsData", definition=function(obj) obj@edgeR_res)
#' @export
setMethod(f="edgeR_res", signature="degResults", definition=function(obj) obj@vsData@edgeR_res)
#' @export
setMethod(f="edgeR_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@edgeR_res)


#' @export
setMethod(f="DESeq2_res", signature="vsData", definition=function(obj) obj@DESeq2_res)
#' @export
setMethod(f="DESeq2_res", signature="degResults", definition=function(obj) obj@vsData@DESeq2_res)
#' @export
setMethod(f="DESeq2_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@DESeq2_res)


#' @export
setMethod(f="merge_res", signature="vsData", definition=function(obj) obj@merge_res)
#' @export
setMethod(f="merge_res", signature="degResults", definition=function(obj) obj@vsData@merge_res)
#' @export
setMethod(f="merge_res", signature="DEGContainer", definition=function(obj) obj@degResults@vsData@merge_res)


#' @export
setReplaceMethod("vsData", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData <- value; validObject(obj); obj})
#' @export
setReplaceMethod("vsData", "degResults",
                 function(obj, value) {obj@vsData <- value; validObject(obj); obj})


#' @export
setReplaceMethod("DESeq2_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@DESeq2_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("DESeq2_res", "degResults",
                 function(obj, value) {obj@vsData@DESeq2_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("DESeq2_res", "vsData",
                 function(obj, value) {obj@DESeq2_res <- value; validObject(obj); obj})


#' @export
setReplaceMethod("limma_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@limma_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("limma_res", "degResults",
                 function(obj, value) {obj@vsData@limma_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("limma_res", "vsData",
                 function(obj, value) {obj@limma_res <- value; validObject(obj); obj})


#' @export
setReplaceMethod("edgeR_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@edgeR_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("edgeR_res", "degResults",
                 function(obj, value) {obj@vsData@edgeR_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("edgeR_res", "vsData",
                 function(obj, value) {obj@edgeR_res <- value; validObject(obj); obj})


#' @export
setReplaceMethod("merge_res", "DEGContainer",
                 function(obj, value) {obj@degResults@vsData@merge_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("merge_res", "degResults",
                 function(obj, value) {obj@vsData@merge_res <- value; validObject(obj); obj})
#' @export
setReplaceMethod("merge_res", "vsData",
                 function(obj, value) {obj@merge_res <- value; validObject(obj); obj})
# Methods for vsData ----

## Methods for hyper results ----
#' @export
setMethod(f="hyperKEGGparam", signature="hyperParam", definition=function(obj) obj@keggParam)
#' @export
setMethod(f="hyperKEGGparam", signature="hyperResults", definition=function(obj) obj@hyperParam@keggParam)
#' @export
setMethod(f="hyperKEGGparam", signature="DEGContainer", definition=function(obj) obj@hyperResults@hyperParam@keggParam)

#' @export
setReplaceMethod("hyperKEGGparam", "DEGContainer",
                 function(obj, value) {obj@hyperResults@hyperParam@keggParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("hyperKEGGparam", "hyperResults",
                 function(obj, value) {obj@hyperParam@keggParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("hyperKEGGparam", "hyperParam",
                 function(obj, value) {obj@keggParam <- value; validObject(obj); obj})

#' @export
setMethod(f="hyperGOparam", signature="hyperParam", definition=function(obj) obj@goParam)
#' @export
setMethod(f="hyperGOparam", signature="hyperResults", definition=function(obj) obj@hyperParam@goParam)
#' @export
setMethod(f="hyperGOparam", signature="DEGContainer", definition=function(obj) obj@hyperResults@hyperParam@goParam)

#' @export
setReplaceMethod("hyperGOparam", "DEGContainer",
                 function(obj, value) {obj@hyperResults@hyperParam@goParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("hyperGOparam", "hyperResults",
                 function(obj, value) {obj@hyperParam@goParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("hyperGOparam", "hyperParam",
                 function(obj, value) {obj@goParam <- value; validObject(obj); obj})

#' @export
setMethod(f="hyperRes", signature="hyperResults", definition=function(obj) obj@hyperRes)
#' @export
setMethod(f="hyperRes", signature="DEGContainer", definition=function(obj) obj@hyperResults@hyperRes)

#' @export
setReplaceMethod("hyperRes", "DEGContainer",
                 function(obj, value) {obj@hyperResults@hyperRes <- value; validObject(obj); obj})
#' @export
setReplaceMethod("hyperRes", "hyperResults",
                 function(obj, value) {obj@hyperRes <- value; validObject(obj); obj})
## Methods for hyper results ----

## Methods for gse results ----
#' @export
setMethod(f="gseKEGGparam", signature="gseParam", definition=function(obj) obj@keggParam)
#' @export
setMethod(f="gseKEGGparam", signature="gseResults", definition=function(obj) obj@gseParam@keggParam)
#' @export
setMethod(f="gseKEGGparam", signature="DEGContainer", definition=function(obj) obj@gseResults@gseParam@keggParam)

#' @export
setReplaceMethod("gseKEGGparam", "DEGContainer",
                 function(obj, value) {obj@gseResults@gseParam@keggParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("gseKEGGparam", "gseResults",
                 function(obj, value) {obj@gseParam@keggParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("gseKEGGparam", "gseParam",
                 function(obj, value) {obj@keggParam <- value; validObject(obj); obj})

#' @export
setMethod(f="gseGOparam", signature="gseParam", definition=function(obj) obj@goParam)
#' @export
setMethod(f="gseGOparam", signature="gseResults", definition=function(obj) obj@gseParam@goParam)
#' @export
setMethod(f="gseGOparam", signature="DEGContainer", definition=function(obj) obj@gseResults@gseParam@goParam)

#' @export
setReplaceMethod("gseGOparam", "DEGContainer",
                 function(obj, value) {obj@gseResults@gseParam@goParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("gseGOparam", "gseResults",
                 function(obj, value) {obj@gseParam@goParam <- value; validObject(obj); obj})
#' @export
setReplaceMethod("gseGOparam", "gseParam",
                 function(obj, value) {obj@goParam <- value; validObject(obj); obj})

#' @export
setMethod(f="gseRes", signature="gseResults", definition=function(obj) obj@gseRes)
#' @export
setMethod(f="gseRes", signature="DEGContainer", definition=function(obj) obj@gseResults@gseRes)

#' @export
setReplaceMethod("gseRes", "DEGContainer",
                 function(obj, value) {obj@gseResults@gseRes <- value; validObject(obj); obj})
#' @export
setReplaceMethod("gseRes", "gseResults",
                 function(obj, value) {obj@gseRes <- value; validObject(obj); obj})
## Methods for gse results ----

## Methods for MSigDB ----
#' @rdname MSigDB
#' @export
setMethod(f="msigdbParam", signature="MSigDB", definition=function(obj) obj@msigdbParam)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbParam", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbParam)

#' @rdname MSigDB
#' @export
setMethod(f="msigdbGSEAparam", signature="MSigDB", definition=function(obj) obj@msigdbGSEAparam)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbGSEAparam", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbGSEAparam)

#' @rdname MSigDB
#' @export
setMethod(f="msigdbHyperParam", signature="MSigDB", definition=function(obj) obj@msigdbHyperParam)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbHyperParam", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbHyperParam)

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbParam", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbParam <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbParam", "MSigDB",
                 function(obj, value) {obj@msigdbParam <- value; validObject(obj); obj})

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbGSEAparam", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbGSEAparam <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbGSEAparam", "MSigDB",
                 function(obj, value) {obj@msigdbGSEAparam <- value; validObject(obj); obj})

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbHyperParam", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbHyperParam <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbHyperParam", "MSigDB",
                 function(obj, value) {obj@msigdbHyperParam <- value; validObject(obj); obj})

#' @rdname MSigDB
#' @export
setMethod(f="msigdbData", signature="MSigDB", definition=function(obj) obj@msigdbData)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbData", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbData)

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbData", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbData <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbData", "MSigDB",
                 function(obj, value) {obj@msigdbData <- value; validObject(obj); obj})

#' @rdname MSigDB
#' @export
setMethod(f="msigdbGSEAresult", signature="MSigDB", definition=function(obj) obj@msigdbGSEAresult)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbGSEAresult", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbGSEAresult)

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbGSEAresult", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbGSEAresult <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbGSEAresult", "MSigDB",
                 function(obj, value) {obj@msigdbGSEAresult <- value; validObject(obj); obj})

#' @rdname MSigDB
#' @export
setMethod(f="msigdbHyperResult", signature="MSigDB", definition=function(obj) obj@msigdbHyperResult)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbHyperResult", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbHyperResult)

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbHyperResult", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbHyperResult <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbHyperResult", "MSigDB",
                 function(obj, value) {obj@msigdbHyperResult <- value; validObject(obj); obj})

#' @rdname MSigDB
#' @export
setMethod(f="msigdbGSVAresult", signature="MSigDB", definition=function(obj) obj@msigdbGSVAresult)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbGSVAresult", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbGSVAresult)

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbGSVAresult", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbGSVAresult <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbGSVAresult", "MSigDB",
                 function(obj, value) {obj@msigdbGSVAresult <- value; validObject(obj); obj})

#' @rdname MSigDB
#' @export
setMethod(f="msigdbTreat", signature="MSigDB", definition=function(obj) obj@msigdbTreat)
#' @rdname DEGContainer
#' @export
setMethod(f="msigdbTreat", signature="DEGContainer", definition=function(obj) obj@MSigDB@msigdbTreat)

#' @rdname DEGContainer
#' @export
setReplaceMethod("msigdbTreat", "DEGContainer",
                 function(obj, value) {obj@MSigDB@msigdbTreat <- value; validObject(obj); obj})
#' @rdname MSigDB
#' @export
setReplaceMethod("msigdbTreat", "MSigDB",
                 function(obj, value) {obj@msigdbTreat <- value; validObject(obj); obj})

#' @export
setMethod(f="cutFC", signature="MSigDB", definition=function(obj) obj@msigdbTreat@cutFC)

#' @export
setReplaceMethod("cutFC", "MSigDB",
                 function(obj, value) {obj@msigdbTreat@cutFC <- value; validObject(obj); obj})

## Methods for MSigDB ----