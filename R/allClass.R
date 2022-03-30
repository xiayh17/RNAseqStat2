setClassUnion("data.frame_OR_NULL", c("data.frame", "NULL"))
setClassUnion("character_OR_NULL", c("character", "NULL"))
setClassUnion("numeric_OR_NULL", c("numeric", "NULL"))
setClassUnion("numeric_OR_list_OR_NULL", c("list","numeric", "NULL"))
setClassUnion("list_OR_NULL", c("list", "NULL"))

#' @title Class \code{dataInfo}
#' @aliases dataInfo-class
#'
#' @family DEGContainer
#' @name dataInfo
#' @docType methods
#'
#' @description This was the class for storing input data for workflow.
#'   We now generally recommend using the \code{\link{Create_DEGContainer}}
#'   to create it in \code{DEGContainer} obj.
#'
#' @slot species species for your data. `Human` or `Mouse`.
#' @slot dataType kind of expresses value matrix. `Counts` (Integer) or `Array` (Decimal).
#' @slot idType kind of gene id. `ENSEMBL` or `SYMBOL`, If `ENSEMBL`, it will be automatically converted to `SYMBOL`.
#' @slot expMatrix expresses value matrix. Should be a data.frame row named by gene ID and column named by Sample
#' @slot groupInfo a Character Vectors ordered by samples in matrix.
#' @slot caseGroup a Character names of case group.
#' @slot filterMethod a function used to filter expresses value matrix. Or disable filter by set as `NULL`.
#' @slot matrixFiltered expresses value matrix after apply \code{filterMethod}. Should be a data.frame row named by gene ID and column named by Sample
#'
#' @importFrom methods setClass new
#'
#' @return The accessor functions \code{species}, \code{dataType}, \cr
#'   \code{idType}, \code{expMatrix}, \code{groupInfo}, \cr
#'   \code{caseGroup}, \code{filterMethod}, \code{matrixFiltered} \cr
#'   return the corresponding elements of a \cr
#'   \code{DEGContainer} or \code{dataInfo}. \cr
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
setClassUnion("dataInfo_OR_NULL", c("dataInfo", "NULL"))
setClassUnion("degResults_OR_NULL", c("degResults", "NULL"))
setClassUnion("gseResults_OR_NULL", c("gseResults", "NULL"))

# msigdb ------------------------------------------------------------------
#' @title Class \code{MSigDB}
#' @aliases MSigDB-class
#'
#' @name MSigDB
#' @docType methods
#' @family DEGContainer
#'
#' @description This was the class for storing data from MSigDB and analysis
#'   results.
#'   We now generally recommend using the \code{\link{Create_DEGContainer}}
#'   to create it in \code{DEGContainer} obj.
#'
#' @slot msigdbParam list_OR_NULL. Store Param of \code{\link[msigdbr]{msigdbr}}
#' @slot msigdbData list. Save data download from MSigDB by \code{\link[msigdbr]{msigdbr}}.
#' @slot msigdbGSEAparam list_OR_NULL. Store Param of \code{\link[clusterProfiler]{GSEA}}
#' @slot msigdbHyperParam list_OR_NULL. Store Param of \code{\link[clusterProfiler]{enricher}}
#' @slot msigdbGSEAresult list. Store results of GSEA.
#' @slot msigdbHyperResult list. Store results of Hyper
#' @slot msigdbGSVAresult list. Store results of GSVA.
#'
#' @return The accessor functions \code{msigdbParam}, \code{msigdbData},
#'   \code{msigdbGSEAparam}, \code{msigdbGSEAresult}, \code{msigdbGSVAresult}
#'   return the corresponding elements of a
#'   \code{DEGContainer} or \code{MSigDB}.
#'
#' @export
setClass(Class = "MSigDB",
         slots = c(msigdbParam = "list_OR_NULL",
                   msigdbData = "list_OR_NULL",
                   msigdbGSEAparam = "list_OR_NULL",
                   msigdbHyperParam = "list_OR_NULL",
                   msigdbGSEAresult = "list_OR_NULL",
                   msigdbHyperResult = "list_OR_NULL",
                   msigdbGSVAresult = "list_OR_NULL",
                   msigdbTreat = "treatInfo_OR_NULL"
         ))
setClassUnion("MSigDB_OR_NULL", c("MSigDB", "NULL"))

# setClass(Class = "IDstat",
# slots = c(
#   OrgDb = "character_OR_NULL",
#   OriginalIDtype = "character_OR_NULL",
#   CurrentIDtype = "character_OR_NULL",
#   LoseOriginal2Current = "character_OR_NULL",
#   Lose2ENTREZID = "character_OR_NULL"
# ))

#' @title Class \code{DEGContainer}
#' @aliases DEGContainer-class
#'
#' @name DEGContainer
#' @docType methods
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

# check Input data --------------------------------------------------------
#' @importFrom usethis ui_info ui_stop ui_value ui_code
setValidity("dataInfo", function(object) {

  species = species(object)
  dataType = dataType(object)
  idType = idType(object)
  group_list = groupInfo(object)
  case_group = caseGroup(object)
  control_group = setdiff(group_list,case_group)
  counts_data = expMatrix(object)
  ## species
  if (length(species)==1) {
    usethis::ui_done("species: {ui_value(species)}")
  } else {
    usethis::ui_stop("Please make sure your {ui_code('species')} is a {ui_value('Character')}.")
  }
  ## dataType
  if (dataType %in% c("Counts","Array")&length(dataType)==1) {
    usethis::ui_done("dataType: {ui_value(dataType)}")
  } else {
    usethis::ui_stop("Please make sure your {ui_code('dataType')} is one of {ui_value('Counts')} or {ui_value('Array')}")
  }
  ## idType
  if (idType %in% c("SYMBOL","ENSEMBL")&length(idType)==1) {
    usethis::ui_done("idType: {ui_value(idType)}")
  } else {
    usethis::ui_stop("Please make sure your {ui_code('idType')} is one of {ui_value('SYMBOL')} or {ui_value('ENSEMBL')}")
  }
  ## counts
  if (identical(dataType,"Counts") & any(class(counts_data) == "data.frame") & all(apply(counts_data, 2, is.integer))) {
    usethis::ui_done("Counts data frame seems ok")
  } else if (identical(dataType,"Array") & any(class(counts_data) == "data.frame") & all(apply(counts_data, 2, is.numeric))) {
    usethis::ui_done("Array data frame seems ok")
  } else {
    usethis::ui_stop("Please check your data frame!
                     Is it an integer data frame for {ui_value('Counts')} or numeric for {ui_value('Array')}?")
  }
  ## group
  if (all(c(case_group, control_group) %in% group_list) & ncol(counts_data) == length(group_list)) {
    usethis::ui_done("{ui_value(case_group)}_VS_{ui_value(control_group)} seems ok")
  } else {
    usethis::ui_stop("Please check {ui_code('groupInfo')} and {ui_code('caseGroup')}")
  }
  ## names
  usethis::ui_info("Please make sure your data frame is rownamed by Gene Symbol")
  ## ----

})

#' @importFrom usethis ui_stop ui_info ui_value
setValidity("treatInfo", function(object) {

  label =  paste0(label(object),sigCollapse = ';')
  if (!label_ns(object) %in% label(object)) {
    usethis::ui_stop("label_ns must be in one of label")
  } else {
    usethis::ui_done("label_ns:{ui_value(label_ns(object))} and label:{ui_value(label)} seems ok")
  }

})
