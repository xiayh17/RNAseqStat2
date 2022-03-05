# Generic for DEGContainer ------------------------------------------------
#' @export
setGeneric(name="dataInfo", def=function(obj) standardGeneric("dataInfo"))

#' @export
setGeneric(name="degResults", def=function(obj) standardGeneric("degResults"))

#' @export
setGeneric(name="hyperResults", def=function(obj) standardGeneric("hyperResults"))

#' @export
setGeneric(name="gseResults", def=function(obj) standardGeneric("gseResults"))

#' @export
setGeneric(name="MSigDB", def=function(obj) standardGeneric("MSigDB"))
# Generic for DEGContainer ------------------------------------------------

# Generic for every slot in dataInfo --------------------------------------------------
#' @rdname dataInfo
#' @export
setGeneric(name="species", def=function(obj) standardGeneric("species"))

#' @rdname dataInfo
#' @export
setGeneric(name="dataType", def=function(obj) standardGeneric("dataType"))

#' @rdname dataInfo
#' @export
setGeneric(name="idType", def=function(obj) standardGeneric("idType"))

#' @rdname dataInfo
#' @export
setGeneric(name="expMatrix", def=function(obj) standardGeneric("expMatrix"))

#' @rdname dataInfo
#' @export
setGeneric(name="groupInfo", def=function(obj) standardGeneric("groupInfo"))

#' @rdname dataInfo
#' @export
setGeneric(name="caseGroup", def=function(obj) standardGeneric("caseGroup"))

#' @rdname dataInfo
#' @export
setGeneric(name="filterMethod", def=function(obj) standardGeneric("filterMethod"))

#' @rdname dataInfo
#' @export
setGeneric(name="matrixFiltered", def=function(obj) standardGeneric("matrixFiltered"))

#' @rdname dataInfo
#' @export
setGeneric(name="geneNames", def=function(obj,filtered) standardGeneric("geneNames"))

#' @rdname dataInfo
#' @export
setGeneric(name="sampleNames", def=function(obj,filtered) standardGeneric("sampleNames"))

## setter
#' @rdname dataInfo
#' @export
setGeneric("species<-", function(obj, value) standardGeneric("species<-"))

#' @rdname dataInfo
#' @export
setGeneric("dataType<-", function(obj, value) standardGeneric("dataType<-"))

#' @rdname dataInfo
#' @export
setGeneric("idType<-", function(obj, value) standardGeneric("idType<-"))

#' @rdname dataInfo
#' @export
setGeneric("expMatrix<-", function(obj, value) standardGeneric("expMatrix<-"))

#' @rdname dataInfo
#' @export
setGeneric("groupInfo<-", function(obj, value) standardGeneric("groupInfo<-"))

#' @rdname dataInfo
#' @export
setGeneric("caseGroup<-", function(obj, value) standardGeneric("caseGroup<-"))

#' @rdname dataInfo
#' @export
setGeneric("filterMethod<-", function(obj, value) standardGeneric("filterMethod<-"))

#' @rdname dataInfo
#' @export
setGeneric("matrixFiltered<-", function(obj, value) standardGeneric("matrixFiltered<-"))

#' @rdname dataInfo
#' @export
setGeneric("geneNames<-", function(obj, value, filtered) standardGeneric("geneNames<-"))

#' @rdname dataInfo
#' @export
setGeneric("sampleNames<-", function(obj, value, filtered) standardGeneric("sampleNames<-"))
# Generic for every slot in dataInfo --------------------------------------------------

# Generic for treatInfo ---------------------------------------------------

#' @export
setGeneric(name="treatInfo", def=function(obj) standardGeneric("treatInfo"))

#' @export
setGeneric(name="label", def=function(obj) standardGeneric("label"))

#' @export
setGeneric(name="label_ns", def=function(obj) standardGeneric("label_ns"))

#' @export
setGeneric(name="cutFC", def=function(obj) standardGeneric("cutFC"))

#' @export
setGeneric(name="cutFDR", def=function(obj) standardGeneric("cutFDR"))

#' @export
setGeneric(name="sigCol", def=function(obj) standardGeneric("sigCol"))

#' @export
setGeneric(name="sigAlpha", def=function(obj) standardGeneric("sigAlpha"))

#' @export
setGeneric(name="sigSize", def=function(obj) standardGeneric("sigSize"))

#' @export
setGeneric(name="sigShape", def=function(obj) standardGeneric("sigShape"))

#' @export
setGeneric("treatInfo<-", function(obj, value) standardGeneric("treatInfo<-"))

#' @export
setGeneric("cutFC<-", function(obj, value) standardGeneric("cutFC<-"))

#' @export
setGeneric("cutFDR<-", function(obj, value) standardGeneric("cutFDR<-"))

#' @export
setGeneric("label<-", function(obj, value) standardGeneric("label<-"))

#' @export
setGeneric("label_ns<-", function(obj, value) standardGeneric("label_ns<-"))

#' @export
setGeneric("sigCol<-", function(obj, value) standardGeneric("sigCol<-"))

#' @export
setGeneric("sigAlpha<-", function(obj, value) standardGeneric("sigAlpha<-"))

#' @export
setGeneric("sigSize<-", function(obj, value) standardGeneric("sigSize<-"))

#' @export
setGeneric("sigShape<-", function(obj, value) standardGeneric("sigShape<-"))
# Generic for treatInfo ---------------------------------------------------

# Generic for vsData ----

#' @export
setGeneric(name="vsData", def=function(obj) standardGeneric("vsData"))

#' @export
setGeneric(name="limma_res", def=function(obj) standardGeneric("limma_res"))

#' @export
setGeneric(name="edgeR_res", def=function(obj) standardGeneric("edgeR_res"))

#' @export
setGeneric(name="DESeq2_res", def=function(obj) standardGeneric("DESeq2_res"))

#' @export
setGeneric(name="merge_res", def=function(obj) standardGeneric("merge_res"))

#' @export
setGeneric("vsData<-", function(obj, value) standardGeneric("vsData<-"))

#' @export
setGeneric("DESeq2_res<-", function(obj, value) standardGeneric("DESeq2_res<-"))

#' @export
setGeneric("limma_res<-", function(obj, value) standardGeneric("limma_res<-"))

#' @export
setGeneric("edgeR_res<-", function(obj, value) standardGeneric("edgeR_res<-"))

#' @export
setGeneric("merge_res<-", function(obj, value) standardGeneric("merge_res<-"))

# Generic for vsData ----

## Generic for hyper results ----
#' @export
setGeneric(name="hyperKEGGparam", def=function(obj) standardGeneric("hyperKEGGparam"))

#' @export
setGeneric("hyperKEGGparam<-", function(obj, value) standardGeneric("hyperKEGGparam<-"))

#' @export
setGeneric(name="hyperGOparam", def=function(obj) standardGeneric("hyperGOparam"))

#' @export
setGeneric("hyperGOparam<-", function(obj, value) standardGeneric("hyperGOparam<-"))

#' @export
setGeneric(name="hyperRes", def=function(obj) standardGeneric("hyperRes"))

#' @export
setGeneric("hyperRes<-", function(obj, value) standardGeneric("hyperRes<-"))
## Generic for hyper results ----

## Generic for gse results ----
#' @export
setGeneric(name="gseKEGGparam", def=function(obj) standardGeneric("gseKEGGparam"))

#' @export
setGeneric("gseKEGGparam<-", function(obj, value) standardGeneric("gseKEGGparam<-"))

#' @export
setGeneric(name="gseGOparam", def=function(obj) standardGeneric("gseGOparam"))

#' @export
setGeneric("gseGOparam<-", function(obj, value) standardGeneric("gseGOparam<-"))

#' @export
setGeneric(name="gseRes", def=function(obj) standardGeneric("gseRes"))

#' @export
setGeneric("gseRes<-", function(obj, value) standardGeneric("gseRes<-"))
## Generic for gse results ----

## Generic for MSigDB ----
#' @export
setGeneric(name="msigdbParam", def=function(obj) standardGeneric("msigdbParam"))

#' @export
setGeneric("msigdbParam<-", function(obj, value) standardGeneric("msigdbParam<-"))

#' @export
setGeneric(name="msigdbGSEAparam", def=function(obj) standardGeneric("msigdbGSEAparam"))

#' @export
setGeneric("msigdbGSEAparam<-", function(obj, value) standardGeneric("msigdbGSEAparam<-"))

#' @export
setGeneric(name="msigdbHyperParam", def=function(obj) standardGeneric("msigdbHyperParam"))

#' @export
setGeneric("msigdbHyperParam<-", function(obj, value) standardGeneric("msigdbHyperParam<-"))

#' @export
setGeneric(name="msigdbData", def=function(obj) standardGeneric("msigdbData"))

#' @export
setGeneric("msigdbData<-", function(obj, value) standardGeneric("msigdbData<-"))

#' @export
setGeneric(name="msigdbGSEAresult", def=function(obj) standardGeneric("msigdbGSEAresult"))

#' @export
setGeneric("msigdbGSEAresult<-", function(obj, value) standardGeneric("msigdbGSEAresult<-"))

#' @export
setGeneric(name="msigdbHyperResult", def=function(obj) standardGeneric("msigdbHyperResult"))

#' @export
setGeneric("msigdbHyperResult<-", function(obj, value) standardGeneric("msigdbHyperResult<-"))


#' @export
setGeneric(name="msigdbGSVAresult", def=function(obj) standardGeneric("msigdbGSVAresult"))

#' @export
setGeneric("msigdbGSVAresult<-", function(obj, value) standardGeneric("msigdbGSVAresult<-"))

#' @export
setGeneric(name="msigdbTreat", def=function(obj) standardGeneric("msigdbTreat"))

#' @export
setGeneric("msigdbTreat<-", function(obj, value) standardGeneric("msigdbTreat<-"))

## Generic for MSigDB ----
