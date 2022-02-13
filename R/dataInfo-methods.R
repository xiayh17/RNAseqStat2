#' dataInfo create function
#'
#' a function to create dataInfo OB
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
#' @return a DEG_container ob
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


# Methods for every slot in dataInfo --------------------------------------------------
##
setGeneric(name="species", def=function(obj) standardGeneric("species"))
setMethod(f="species", signature="DEGContainer", definition=function(obj) obj@dataInfo@species)

setGeneric(name="dataType", def=function(obj) standardGeneric("dataType"))
setMethod(f="dataType", signature="DEGContainer", definition=function(obj) obj@dataInfo@dataType)

setGeneric(name="idType", def=function(obj) standardGeneric("idType"))
setMethod(f="idType", signature="DEGContainer", definition=function(obj) obj@dataInfo@idType)

setGeneric(name="expMatrix", def=function(obj) standardGeneric("expMatrix"))
setMethod(f="expMatrix", signature="DEGContainer",definition=function(obj) obj@dataInfo@expMatrix)

setGeneric(name="groupInfo", def=function(obj) standardGeneric("groupInfo"))
setMethod("groupInfo", "DEGContainer", function(obj) obj@dataInfo@groupInfo)

setGeneric(name="caseGroup", def=function(obj) standardGeneric("caseGroup"))
setMethod(f="caseGroup", signature="DEGContainer",definition=function(obj) obj@dataInfo@caseGroup)

setGeneric(name="filterMethod", def=function(obj) standardGeneric("filterMethod"))
setMethod(f="filterMethod", signature="DEGContainer", definition=function(obj) obj@dataInfo@filterMethod)

setGeneric(name="matrixFiltered", def=function(obj) standardGeneric("matrixFiltered"))
setMethod(f="matrixFiltered", signature="DEGContainer", definition=function(obj) obj@dataInfo@matrixFiltered)

setGeneric(name="geneNames", def=function(obj,filtered) standardGeneric("geneNames"))
setMethod(f="geneNames", signature="DEGContainer", definition=function(obj,filtered = FALSE) {
  if (filtered) {
    genes_name = rownames(obj@dataInfo@matrixFiltered)
  } else {
    genes_name = rownames(obj@dataInfo@expMatrix)
  }
  return(genes_name)})

setGeneric(name="sampleNames", def=function(obj,filtered) standardGeneric("sampleNames"))
setMethod(f="sampleNames", signature="DEGContainer", definition=function(obj,filtered = FALSE) {
  if (filtered) {
    sample_names = colnames(obj@dataInfo@matrixFiltered)
  } else {
    sample_names = colnames(obj@dataInfo@expMatrix)
  }
  return(sample_names)})

setMethod(f="species", signature="dataInfo", definition=function(obj) obj@species)
setMethod(f="dataType", signature="dataInfo", definition=function(obj) obj@dataType)
setMethod(f="idType", signature="dataInfo", definition=function(obj) obj@idType)
setMethod(f="expMatrix", signature="dataInfo", definition=function(obj) obj@expMatrix)
setMethod("groupInfo", "dataInfo", function(obj) obj@groupInfo)
setMethod("caseGroup", "dataInfo", function(obj) obj@caseGroup)
setMethod(f="filterMethod", signature="dataInfo", definition=function(obj) obj@filterMethod)
setMethod(f="matrixFiltered", signature="dataInfo", definition=function(obj) obj@matrixFiltered)
setMethod(f="sampleNames", signature="dataInfo", definition=function(obj,filtered = FALSE) {
  if (filtered) {
    sample_names = colnames(obj@matrixFiltered)
  } else {
    sample_names = colnames(obj@expMatrix)
  }
  return(sample_names)})

## setter
setGeneric("species<-", function(obj, value) standardGeneric("species<-"))
setReplaceMethod("species", "dataInfo",
                 function(obj, value) {obj@species <- value; validObject(obj); obj})

setGeneric("dataType<-", function(obj, value) standardGeneric("dataType<-"))
setReplaceMethod("dataType", "dataInfo",
                 function(obj, value) {obj@dataType <- value; validObject(obj); obj})

setGeneric("idType<-", function(obj, value) standardGeneric("idType<-"))
setReplaceMethod("idType", "dataInfo",
                 function(obj, value) {obj@idType <- value; validObject(obj); obj})

setGeneric("expMatrix<-", function(obj, value) standardGeneric("expMatrix<-"))
setReplaceMethod("expMatrix", "dataInfo",
                 function(obj, value) {obj@expMatrix <- value; validObject(obj); obj})

setGeneric("groupInfo<-", function(obj, value) standardGeneric("groupInfo<-"))
setReplaceMethod("groupInfo", "dataInfo",
                 function(obj, value) {obj@groupInfo <- value; validObject(obj); obj})

setGeneric("caseGroup<-", function(obj, value) standardGeneric("caseGroup<-"))
setReplaceMethod("caseGroup", "dataInfo",
                 function(obj, value) {obj@caseGroup <- value; validObject(obj); obj})

setGeneric("filterMethod<-", function(obj, value) standardGeneric("filterMethod<-"))
setReplaceMethod("filterMethod", "dataInfo",
                 function(obj, value) {obj@filterMethod <- value; validObject(obj); obj})

setGeneric("matrixFiltered<-", function(obj, value) standardGeneric("matrixFiltered<-"))
setReplaceMethod("matrixFiltered", "dataInfo",
                 function(obj, value) {obj@matrixFiltered <- value; validObject(obj); obj})

setGeneric("geneNames<-", function(obj, value, filtered) standardGeneric("geneNames<-"))
setReplaceMethod("geneNames", "dataInfo",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     rownames(obj@matrixFiltered) <- value;
                   } else {
                     rownames(obj@expMatrix) <- value;
                   };

                   validObject(obj); obj})
setGeneric("sampleNames<-", function(obj, value, filtered) standardGeneric("sampleNames<-"))
setReplaceMethod("sampleNames", "dataInfo",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     colnames(obj@matrixFiltered) <- value;
                   } else {
                     colnames(obj@expMatrix) <- value;
                   };

                   validObject(obj); obj})

setReplaceMethod("species", "DEGContainer",
                 function(obj, value) {obj@dataInfo@species <- value; validObject(obj); obj})

setReplaceMethod("dataType", "DEGContainer",
                 function(obj, value) {obj@dataInfo@dataType <- value; validObject(obj); obj})

setReplaceMethod("idType", "DEGContainer",
                 function(obj, value) {obj@dataInfo@idType <- value; validObject(obj); obj})

setReplaceMethod("expMatrix", "DEGContainer",
                 function(obj, value) {obj@dataInfo@expMatrix <- value; validObject(obj); obj})

setReplaceMethod("groupInfo", "DEGContainer",
                 function(obj, value) {obj@dataInfo@groupInfo <- value; validObject(obj); obj})

setReplaceMethod("caseGroup", "DEGContainer",
                 function(obj, value) {obj@dataInfo@caseGroup <- value; validObject(obj); obj})

setReplaceMethod("filterMethod", "DEGContainer",
                 function(obj, value) {obj@dataInfo@filterMethod <- value; validObject(obj); obj})

setReplaceMethod("matrixFiltered", "DEGContainer",
                 function(obj, value) {obj@dataInfo@matrixFiltered <- value; validObject(obj); obj})

setReplaceMethod("geneNames", "DEGContainer",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     rownames(obj@dataInfo@matrixFiltered) <- value;
                   } else {
                     rownames(obj@dataInfo@expMatrix) <- value;
                   };

                   validObject(obj); obj})

setReplaceMethod("sampleNames", "DEGContainer",
                 function(obj, value, filtered = FALSE) {

                   if (filtered) {
                     colnames(obj@dataInfo@matrixFiltered) <- value;
                   } else {
                     colnames(obj@dataInfo@expMatrix) <- value;
                   };

                   validObject(obj); obj})

# End of Methods for every slot in dataInfo --------------------------------------------------

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
  if (species %in% c("Human","Mouse")&length(species)==1) {
    usethis::ui_done("species: {ui_value(species)}")
  } else {
    usethis::ui_stop("Please make sure your {ui_code('species')} is one of {ui_value('Human')} or {ui_value('Mouse')}")
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
  if (any(class(counts_data) == "data.frame") & all(apply(counts_data, 2, is.integer))) {
    usethis::ui_done("Counts data frame seems ok")
  } else {
    usethis::ui_stop("Please check your data frame! Is it an integer data frame?")
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
