setClassUnion("data.frame_OR_NULL", c("data.frame", "NULL"))
setClassUnion("character_OR_NULL", c("character", "NULL"))
setClassUnion("numeric_OR_NULL", c("numeric", "NULL"))
setClassUnion("numeric_OR_list_OR_NULL", c("list","numeric", "NULL"))

#' a S4 class contains inputs info for workflow
#'
#' contains begin info of RNAseqStat2
#'
#' @importFrom methods setClass new
#'
#' @rdname dataInfo
#'
#' @family DEGContainer
#'
#' @export
setClass(Class="dataInfo",
         slots = c(
           species = "character", ## 物种
           dataType = "character", ## 数据类型，counts 或者 array
           idType = "character", ## id 类型
           expMatrix = "data.frame", ## 矩阵数据框
           groupInfo = "character",  ## 分组信息字符串
           caseGroup = "character", ## 实验组名称
           filterMethod = "character_OR_NULL", ## 过滤函数
           matrixFiltered = "data.frame_OR_NULL" ## 过滤后的矩阵
         )
)

#' a S4 class contains config info for workflow
#'
#' contains configurations info of RNAseqStat2
#'
#' @importFrom methods setClass new
#'
#' @rdname dataInfo
#'
#' @family DEGContainer
#'
setClass(Class = "treatInfo", ##
         slots = c(
           cutFC = "numeric_OR_list_OR_NULL", ## FC 阈值
           cutFDR = "numeric_OR_list_OR_NULL", ## FDR 阈值
           label = "character", ## 分组名字
           label_ns = "character", ## 不显著分组名字
           sigCol = "character", ## 分组颜色
           sigAlpha = "numeric", ## 分组透明度
           sigSize = "numeric", ## 分组点的大小
           sigShape = "numeric" ## 分组点的形状
         ))
setClassUnion("treatInfo_OR_NULL", c("treatInfo", "NULL"))

setClass(Class = "vsData",
         slots = c(
           limma_res = "data.frame_OR_NULL",
           edgeR_res = "data.frame_OR_NULL",
           DESeq2_res = "data.frame_OR_NULL",
           merge_res = "data.frame_OR_NULL"
         ))
setClassUnion("vsData_OR_NULL", c("vsData", "NULL"))

#' a S4 class contains degResults from deg analysis
#'
#' contains deg analysis results in RNAseqStat2
#'
#' @importFrom methods setClass new
#'
#' @rdname degResults
#'
#' @family DEGContainer
#'
setClass(Class="degResults",
         slots = c(
           vsData = "vsData_OR_NULL", ## 差异分析表 含分组
           treatInfo = "treatInfo_OR_NULL"
         )
)


# hyper -------------------------------------------------------------------
setClassUnion("list_OR_NULL", c("list", "NULL"))
setClass(Class = "hyperParam",
         slots = c(goParam = "list",
                   keggParam = "list"))
setClassUnion("hyperParam_OR_NULL", c("hyperParam", "NULL"))

setClass(Class = "hyperResults",
         slots = c(
           hyperRes = "list_OR_NULL", ## 分析结果
           hyperParam = "hyperParam_OR_NULL"
         ))
setClassUnion("hyperResults_OR_NULL", c("hyperResults", "NULL"))


# gse ---------------------------------------------------------------------
setClass(Class = "gseParam",
         slots = c(goParam = "list",
                   keggParam = "list"))
setClassUnion("gseParam_OR_NULL", c("gseParam", "NULL"))

setClass(Class = "gseResults",
         slots = c(
           gseRes = "list_OR_NULL", ## 分析结果
           gseParam = "gseParam_OR_NULL"
         ))
setClassUnion("gseResults_OR_NULL", c("gseResults", "NULL"))

setClassUnion("dataInfo_OR_NULL", c("dataInfo", "NULL"))
setClassUnion("degResults_OR_NULL", c("degResults", "NULL"))
setClassUnion("gseResults_OR_NULL", c("gseResults", "NULL"))

# msigdb ------------------------------------------------------------------
#' @title Class \code{MSigDB}
#' @aliases MSigDB-class
#'
#' @description This was the class for storing data from MSigDB and analysis
#'   results. We now generally reommend using the
#'   \code{\link{Create_DEGContainer}} to create it.
#'
#' @slot msigdbParam list_OR_NULL.
#' @slot msigdbData list.
#' @slot msigdbGSEAparam list_OR_NULL.
#' @slot msigdbGSEAresult list.
#' @slot msigdbGSVAresult list.
#'
#' @return The accessor functions \code{msigdbParam}, \code{msigdbData},
#'   \code{msigdbGSEAparam}, \code{msigdbGSEAresult}, \code{msigdbGSVAresult}
#'   return the corresponding elements of a
#'   \code{DEGContainer} or \code{MSigDB}.
#'
#' @export
setClass(Class = "MSigDB",
         slots = c(msigdbParam = "list_OR_NULL",
                   msigdbData = "list",
                   msigdbGSEAparam = "list_OR_NULL",
                   msigdbGSEAresult = "list",
                   msigdbGSVAresult = "list"
         ))
setClassUnion("MSigDB_OR_NULL", c("MSigDB", "NULL"))

#' @title Class \code{DEGContainer}
#' @aliases DEGContainer-class
#'
#' @description This was the universal class for storing data and
#'   results. We now generally reommend using the
#'   \code{\link{Create_DEGContainer}} to create it.
#'
#' @slot dataInfo store input data of workflow
#' @slot degResults store parameters and results in \code{runDEG} module
#' @slot hyperResults store parameters and results in \code{runHyper} module
#' @slot gseResults store parameters and results in \code{runGSEA} module
#' @slot MSigDB store parameters and results in \code{runMSigDB} module
#'
#' @return The accessor functions \code{dataInfo}, \code{degResults},
#'   \code{hyperResults}, \code{gseResults}, \code{MSigDB} return
#'   the corresponding elements of a
#'   \code{DEGContainer}.
#'
#' @export
setClass(Class = "DEGContainer",
          slots = c(
            dataInfo = "dataInfo_OR_NULL",
            degResults = "degResults_OR_NULL",
            hyperResults = "hyperResults_OR_NULL",
            gseResults = "gseResults_OR_NULL",
            MSigDB = "MSigDB_OR_NULL"
          )
)

# methods for DEGContainer ------------------------------------------------
#' @importFrom usethis ui_line ui_value
setMethod("show", "DEGContainer",
          function(object) {

            ui_line("DEGContainer
                     {ui_value(dataType(object))} Matrix Containing {ui_value(ncol(expMatrix(object)))} samples")
          }
)

#' @rdname DEGContainer
#' @export
setGeneric(name="dataInfo", def=function(obj) standardGeneric("dataInfo"))
setMethod(f="dataInfo", signature="DEGContainer", definition=function(obj) obj@dataInfo)

#' @rdname DEGContainer
#' @export
setGeneric(name="degResults", def=function(obj) standardGeneric("degResults"))
setMethod(f="degResults", signature="DEGContainer", definition=function(obj) obj@degResults)

#' @rdname DEGContainer
#' @export
setGeneric(name="hyperResults", def=function(obj) standardGeneric("hyperResults"))
setMethod(f="hyperResults", signature="DEGContainer", definition=function(obj) obj@hyperResults)

#' @rdname DEGContainer
#' @export
setGeneric(name="gseResults", def=function(obj) standardGeneric("gseResults"))
setMethod(f="gseResults", signature="DEGContainer", definition=function(obj) obj@gseResults)

#' @rdname DEGContainer
#' @export
setGeneric(name="MSigDB", def=function(obj) standardGeneric("MSigDB"))
setMethod(f="MSigDB", signature="DEGContainer", definition=function(obj) obj@MSigDB)
