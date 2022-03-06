#' @import enrichplot
#' @export
setGeneric(name="runHyper", def=function(obj, dir = ".", prefix = "3-runHyper",top = 10,GO = FALSE, KEGG = TRUE) standardGeneric("runHyper"))

setMethod(f="runHyper", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "3-runHyper",top = 10,GO = FALSE, KEGG = TRUE) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## hyper resolve of limma edgeR DESeq2
  if (length(hyperRes(obj)) == 0) {
    obj <- hyperResolve(obj = obj,GO = GO,KEGG=KEGG)
  }

  hyperSummary(obj = obj, dir = dir, prefix = prefix,top = top)

  return(obj)

})
