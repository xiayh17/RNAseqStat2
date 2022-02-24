#' @export
degSummary <- function(obj, dir = ".", prefix = "2-runDEG",
                       PointVolcanoParam = list(gene = 10,light = NULL,label_light = TRUE,
                                                light_color = "#24ac56",
                                                light_label_color = "#24ac56",
                                                expend = c(0.12, 0.12))) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  test <- deg_here(obj)
  ok <- names(test)[which(test == TRUE)]
  main <- setdiff(ok,"merge")
  # merge_data <- intersect(ok,"merge")

  if(test["merge"]) {

    dat <- dataDEG(obj,which = "merge")
    csv_name = glue('{dir}/{prefix}_merge_results.csv')
    write.csv(dat,file = csv_name)
    ui_done(glue("merge DEG results in csv format is stored in {usethis::ui_path(csv_name)}"))

    geneSymbol_list <- hyper_GS(obj,"merge",type = "SYMBOL")
    lapply(seq_along(geneSymbol_list), function(x){

      file_name <- glue("{dir}/{prefix}_{names(geneSymbol_list)[x]}_merge_Gene.txt")
      writeLines(geneSymbol_list[[x]],con = file_name)
      ui_done(glue("merge {names(geneSymbol_list)[x]} Gene in txt format is stored in {usethis::ui_path(file_name)}"))

    })

  }

  ## plot and export data
  if(length(main)>=1){

    for (i in main) {

      which = i

      geneSymbol_list <- hyper_GS(object = obj,which = which,type = "SYMBOL")
      lapply(seq_along(geneSymbol_list), function(x){

        file_name <- glue("{dir}/{prefix}_{names(geneSymbol_list)[x]}_{which}_Gene.txt")
        writeLines(geneSymbol_list[[x]],con = file_name)
        ui_done(glue("{which} {names(geneSymbol_list)[x]} Gene in txt format is stored in {usethis::ui_path(file_name)}"))

      })

      dat <- dataDEG(obj,which = which)
      csv_name = glue('{dir}/{prefix}_{which}_results.csv')
      write.csv(dat,file = csv_name)
      ui_done(glue("{which} DEG results in csv format is stored in {usethis::ui_path(csv_name)}"))

      volcano_plot <-  do.call("PointVolcano",c(list(obj = obj,which = i), PointVolcanoParam))
      volcano_file = glue("{dir}/{prefix}_{which}_volcano.pdf")
      ggsave(volcano_plot,filename = volcano_file, width = 1600,height = 1600,units = "px",limitsize = FALSE,device = cairo_pdf)
      ui_done(glue("{which} volcano results were store in {ui_path(volcano_file)}."))

      topHeatmap_file = glue("{dir}/{prefix}_{which}_volcano.pdf")
      topHeatmap <- do.call("DEGtopHeatmap",list(obj = obj,which = i,filename = topHeatmap_file))
      ui_done(glue("{which} top 50 heatmap were store in {ui_path(topHeatmap_file)}."))

    }

  } else {
    ui_info("None available results of DEG.")
  }

}

dataDEG <- function(obj,which) {

  if (which == "limma") {
    deg_data = limma_res(obj)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(obj)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(obj)
  } else if (which == "merge") {
    deg_data = merge_res(obj)
  } else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2')}")
  }

  return(deg_data)

}
