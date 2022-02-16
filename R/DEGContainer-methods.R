# Low level function to create DEGContainer -------------------------------
#' @title Initialize an object of class \code{DEGContainer}
#'
#' @name DEGContainer
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
DEGContainer <- function(dataInfo=NULL,
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

#' Create \code{DEGContainer}
#'
#' Begin in DEGContainer.
#'
#' @param species species for your data. `Human` or `Mouse`.
#' @param dataType kind of expresses value matrix. `Counts` (Integer) or `Array` (Decimal).
#' @param idType kind of gene id. `ENSEMBL` or `SYMBOL`, If `ENSEMBL`, it will be automatically converted to `SYMBOL`.
#' @param expMatrix expresses value matrix. Should be a data.frame row named by gene ID and column named by Sample
#' @param groupInfo a Character Vectors ordered by samples in matrix.
#' @param caseGroup a Character names of case group.
#' @param filterMethod a function used to filter expresses value matrix. Or disable filter by set as `NULL`.
#'
#' @importFrom usethis ui_info ui_done ui_code
#'
#' @return a \code{DEGContainer} object
#' @export
#'
#' @examples
#' Create_DEGContainer(expMatrix = counts_input,groupInfo = c("C","C","C","T","T","T"),caseGroup = "T")
Create_DEGContainer <- function(species = "Human",
                                dataType = "Counts",
                                idType = "SYMBOL",
                                expMatrix,
                                groupInfo,
                                caseGroup,
                                filterMethod = "rowSums(expMatrix > 0) >= ncol(expMatrix)/2"){

  expMatrix <- as.data.frame(expMatrix)

  ## convert ids
  if (idType == "ENSEMBL") {
    expMatrix <- suppressWarnings(toSYMBOL(row_counts = expMatrix))
    idType = "SYMBOL"
    ui_done("ENSEMBL named expression matrix renamed by SYMBOL")
  }

  ## check filter
  if(is.null(filterMethod)){

    ui_info("Nothing filtered in your data.")

  } else {

    keep_feature <- eval(parse(text = filterMethod))
    matrixFiltered <- as.data.frame(expMatrix[keep_feature,])
    left = nrow(matrixFiltered)


    if (left == 0) {

      matrixFiltered = NULL
      ui_stop("Your {ui_code('filterMethod')} have some problem OR too Strict.
              It Should be something like {ui_code('rowSums(expMatrix > 2) >= ncol(expMatrix)')} in a character.")

    } else {

      ui_done("Your data have filtered by {ui_code(filterMethod)}")
      ui_info("{left} Keeped, {nrow(expMatrix) - left} removed. ")

    }
  }

  ## Create basic data info
  dataInfo <- Create_dataInfo(
    species = species,
    dataType = dataType,
    idType = idType,
    expMatrix = expMatrix,
    groupInfo = groupInfo,
    caseGroup = caseGroup,
    filterMethod = filterMethod,
    matrixFiltered = matrixFiltered
  )

  ## Create zero for DEG
  degResults <- suppressWarnings(Create_degResults(
    vsData = Create_vsData(),
    treatInfo = Create_treatInfo()
  ))

  ## Create zero for enrich
  OrgDb = switch (species,
                  'Human' = "org.Hs.eg.db",
                  'Mouse' = "org.Mm.eg.db"
  )

  organism = switch (species,
                     'Human' = "hsa",
                     'Mouse' = "mmu"
  )
  hyperResults <- suppressWarnings(Create_hyperResults(hyperParam=Create_hyperParam(goParam = list(OrgDb = OrgDb),keggParam = list(organism = organism))))
  gseResults <- suppressWarnings(Create_gseResults(gseParam=Create_gseParam(goParam = list(OrgDb = OrgDb),keggParam = list(organism = organism))))

  msi_species <- switch (species,
                         'Human' = "Homo sapiens",
                         'Mouse' = "Mus musculus"
  )

  MSigDB = suppressWarnings(Create_MSigDB(msigdbParam = Create_msigdbParam(msigdbParam = list(species = msi_species,category = "H")),msigdbGSEAparam = Create_msigdbGSEAparam()))

  DEGContainer(dataInfo = dataInfo,
               degResults = degResults,
               hyperResults = hyperResults,
               gseResults = gseResults,
               MSigDB = MSigDB)

}
