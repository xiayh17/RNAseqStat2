#' @export
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

  ## DEG resolve of limma edgeR DESeq2
  if (dataType(obj) == "Array") {
    obj <- degResolveArray(obj = obj, dir = dir, prefix = prefix)
  } else {

    obj <- degResolve(obj = obj, dir = dir, prefix = prefix,parallel = parallel, qc =qc)

  }

  ## group DEG Results
  obj <- degGroup(obj = obj)

  degSummary(obj = obj, dir = dir,prefix = prefix,PointVolcanoParam = PointVolcanoParam)

  return(obj)

})
