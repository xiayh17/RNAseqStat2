setGeneric(name="runDEG", def=function(obj, dir = ".", prefix = "2-runDEG",
                                            parallel = TRUE,qc=TRUE,
                                             gene = 10,light = NULL,label_light = TRUE,
                                             light_color = "#24ac56",
                                             light_label_color = "#24ac56",
                                             expend = c(0.12, 0.12)) standardGeneric("runDEG"))

setMethod(f="runDEG", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "2-runDEG",
                                                                    parallel = TRUE,qc=TRUE,
                                                                    gene = 10,light = NULL,label_light = TRUE,
                                                                    light_color = "#24ac56",
                                                                    light_label_color = "#24ac56",
                                                                    expend = c(0.12, 0.12)) {

  ## DEG resolve of limma edgeR DESeq2

  obj <- degResolve(obj = obj, dir = dir, prefix = prefix,parallel = parallel, qc =qc)

  ## group DEG Results
  obj <- degGroup(obj = obj)

  test <- deg_here(obj)

  ok <- names(test)[which(test == TRUE)]

  main <- setdiff(ok,"merge")
  merge_data <- intersect(ok,"merge")

  ## plot and export volcano
  if(length(main)>=1){

    for (i in main) {

      which = i
      volcano_plot <-  PointVolcano(obj = obj,
                                    which = which,
                                    gene = gene,
                                    light = light,
                                    label_light = label_light,
                                    light_color = light_color,
                                    light_label_color = light_label_color,
                                    expend = expend)

      volcano_file = glue("{dir}/{prefix}_{which}_volcano.pdf")
      ggsave(volcano_plot,filename = volcano_file, width = 1600,height = 1600,units = "px",limitsize = FALSE,device = cairo_pdf)
      ui_done(glue("{which} volcano results were store in {ui_path(volcano_file)}."))

    }

  } else {
    ui_info("None available results of DEG.")
  }

  return(obj)

})
