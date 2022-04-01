#' run MSigDB module
#'
#' @param obj a DEGContainer after MSigDB analysis
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param top top
#'
#' @return
#' @export
#'
#' @examples
#' runMSigDB(data_gse)
setGeneric(name="runMSigDB", def=function(obj, dir = ".", prefix = "5-runMSigDB", top =10) standardGeneric("runMSigDB"))
setMethod(f="runMSigDB", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "5-runMSigDB",top = 10) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  if (!is.null(msigdbParam(obj))) {

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

    obj@MSigDB <- degGroup(obj@MSigDB)

    tryCatch(
      expr = {
        MSigDBSummary(obj, dir = dir, prefix = prefix,top =top)
      },
      error = function(e){
        usethis::ui_oops("Something wrong occured in Hyper Summary. try again later by {ui_code('MSigDBSummary')}.")
      },
      finally = {
        usethis::ui_line("Hyper analysis step done")
      }
    )

  } else {
    ui_info("MSigDB step skiped.")
  }

  return(obj)

})
