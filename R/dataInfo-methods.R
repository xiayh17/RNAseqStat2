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
