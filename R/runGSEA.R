#' @import enrichplot
#' @import ggplot2
#' @importFrom DOSE theme_dose
#' @export
setGeneric(name="runGSEA", def=function(obj, dir = ".", prefix = "4-runGSEA",top =10) standardGeneric("runGSEA"))

setMethod(f="runGSEA", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "4-runGSEA",top = 10) {

  ## hyper resolve of limma edgeR DESeq2
  if (length(gseRes(obj)) == 0) {
    obj <- gseResolve(obj = obj)
  }

  gesSummary(obj = obj, dir = dir, prefix = prefix,top = top)

  return(obj)

})
