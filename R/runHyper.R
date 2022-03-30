#' runHyper
#'
#' @param obj a DEGContainer
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param top top for enrich result filter
#' @param GO logic, run GO in enrich step
#' @param KEGG logic, run KEGG in enrich step
#'
#' @import enrichplot
#' @export
#'
#' @return a DEGContainer
#' @export
#'
#' @examples
#' runHyper(DEGContainer)
setGeneric(name="runHyper", def=function(obj, dir = ".", prefix = "3-runHyper",top = 10,GO = FALSE, KEGG = TRUE) standardGeneric("runHyper"))

setMethod(f="runHyper", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "3-runHyper",top = 10,GO = FALSE, KEGG = TRUE) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## hyper resolve of limma edgeR DESeq2
  if (length(hyperRes(obj)) == 0) {
    obj <- hyperResolve(obj = obj,GO = GO,KEGG=KEGG)
  }

  tryCatch(
    expr = {
      hyperSummary(obj = obj, dir = dir, prefix = prefix,top = top)
    },
    error = function(e){
      usethis::ui_oops("Something wrong occured in Hyper Summary. try again later by {ui_code(hyperSummary)}.")
    },
    finally = {
      usethis::ui_line("Hyper analysis step done")
    }
  )

  return(obj)

})
