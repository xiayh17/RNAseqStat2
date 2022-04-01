#' runGSEA Module
#'
#' runGSEA Module control GO and KEGG
#'
#' @param obj a DEGContainer
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param top top rows of up and down
#' @param GO run GO
#' @param KEGG run KEGG
#'
#' @import enrichplot
#' @import ggplot2
#' @importFrom DOSE theme_dose
#' @export
#'
#' @return
#' @export
#'
#' @examples
#' runGSEA(DEGContainer)
setGeneric(name="runGSEA", def=function(obj, dir = ".", prefix = "4-runGSEA",top =10, GO = FALSE, KEGG = TRUE) standardGeneric("runGSEA"))

setMethod(f="runGSEA", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "4-runGSEA",top = 10, GO = FALSE, KEGG = TRUE) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## hyper resolve of limma edgeR DESeq2
  if (length(gseRes(obj)) == 0) {
    obj <- gseResolve(obj = obj, GO = GO, KEGG = KEGG)
  }

  tryCatch(
    expr = {
      gseSummary(obj = obj, dir = dir, prefix = prefix,top = top)
    },
    error = function(e){
      usethis::ui_oops("Something wrong occured in GSEA Summary. try again later by {ui_code('gseSummary')}.")
    },
    finally = {
      usethis::ui_line("GSEA analysis step done")
    }
  )

  return(obj)

})
