#' run DEG analysis
#'
#' run DEG analysis and summary results
#'
#' @param obj a DEGContainer
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param parallel use parallel in DESeq2
#' @param qc qc in DESeq2
#' @param PointVolcanoParam for volcano plot
#'
#' @return
#' @export
#'
#' @examples
#' runDEG(DEGContainer)
setGeneric(name="runDEG", def=function(obj, dir = ".", prefix = "2-runDEG",
                                            parallel = TRUE,qc=TRUE,
                                       PointVolcanoParam = list(gene = 10,light = NULL,
                                                                light_color = "#24ac56",
                                                                light_label_color = "#24ac56",
                                                                expend = c(0.12, 0.12))) standardGeneric("runDEG"))

setMethod(f="runDEG", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "2-runDEG",
                                                                    parallel = TRUE,qc=TRUE,
                                      PointVolcanoParam = list(gene = 10,light = NULL,
                                                               light_color = "#24ac56",
                                                               light_label_color = "#24ac56",
                                                               expend = c(0.12, 0.12))) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## DEG resolve of limma edgeR DESeq2
  if (dataType(obj) == "Array") {

    if (all(deg_here(obj)[c(1,2,4)])) {
      usethis::ui_info("DEG analysis seems done. ignored")
    } else {
      obj <- degResolveArray(obj = obj, dir = dir, prefix = prefix)
    }

  } else {

    if (all(deg_here(obj))) {
      usethis::ui_info("DEG analysis seems done. ignored")
    } else {
      obj <- degResolve(obj = obj, dir = dir, prefix = prefix,parallel = parallel, qc =qc)
    }


  }

  ## group DEG Results
  obj <- degGroup(obj = obj)

  degSummary(obj = obj, dir = dir,prefix = prefix,PointVolcanoParam = PointVolcanoParam)

  return(obj)

})
