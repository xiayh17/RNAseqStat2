#' Group degResults
#'
#' group results of DEG analysis
#'
#' @param obj a DEGContainer
#'
#' @return a DEGContainer
#' @export
#'
#' @examples
#' degGroup(DEGContainer)
setGeneric(name="degGroup", def=function(obj) standardGeneric("degGroup"))
setMethod(f="degGroup", signature="DEGContainer", definition=function(obj) {

  test <- deg_here(obj)

  ok <- names(test)[which(test == TRUE)]

  main <- setdiff(ok,"merge")
  merge_data <- intersect(ok,"merge")

  if(length(main)>=1){

    for (i in main) {

      obj <- cutMuch(obj,which = i)

    }

  } else {
    ui_info("None available results of DEG.")
  }

  if (length(main)>=2) {

    merge_res(obj) <- merge_deg(obj)

    ui_done("Common Genes updated")

  }

  return(obj)

})

setMethod(f="degGroup", signature="MSigDB", definition=function(obj) {

  main = names(msigdbGSVAresult(obj)[['GSVA_diff']])

  if(length(main)>=1){

    for (i in main) {

      obj <- cutMuch(obj,which = "MSigDB",category = i)

    }

  } else {
    ui_info("None available results of MSigDB DEG analysis.")
  }

  return(obj)

})
