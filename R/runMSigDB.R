#' @export
setGeneric(name="runMSigDB", def=function(obj, dir = ".", prefix = "5-runMSigDB", top =10) standardGeneric("runMSigDB"))
setMethod(f="runMSigDB", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "5-runMSigDB",top = 10) {

  ## Download  data
  if (length(msigdbData(obj)) == 0) {
    obj <- msigdbGet(obj)
  }

  ## Do gse
  if (length(msigdbGSEAresult(obj)) == 0) {
    obj <- gseMSigDB(obj)
  }

  ## Do Hyper
  if (length(msigdbGSEAresult(obj)) == 0) {
    obj <- hyperMSigDB(obj)
  }

  ## gsva
  if (length(msigdbGSVAresult(obj)) == 0) {
    obj <- gsvaResolve(obj)
  }



  return(obj)

})
