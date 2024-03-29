#' Summary results of DEG analysis
#'
#' @param obj a DEGContainer object
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param PointVolcanoParam more parameters for \code{\link{PointVolcano}}
#'
#' @return a batch of files
#' @export
#'
#' @examples
#' degSummary(DEGContainer)
degSummary <- function(obj, dir = ".", prefix = "2-runDEG",
                       PointVolcanoParam = list(gene = 10,light = NULL,
                                                light_color = "#24ac56",
                                                light_label_color = "#24ac56",
                                                expend = c(0.12, 0.12))) {

  ## create dir if not exists
  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## check deg results
  test <- deg_here(obj)
  ok <- names(test)[which(test == TRUE)]
  ## except merge
  main <- setdiff(ok,"merge")

  ## when merge exists
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

  ## other dataset
  if(length(main)>=1){

    for (i in main) {

      which = i

      geneSymbol_list <- hyper_GS(object = obj,which = which,type = "SYMBOL")
      ## gene names of every  regulate group
      lapply(seq_along(geneSymbol_list), function(x){

        file_name <- glue("{dir}/{prefix}_{names(geneSymbol_list)[x]}_{which}_Gene.txt")
        writeLines(geneSymbol_list[[x]],con = file_name)
        ui_done(glue("{which} {names(geneSymbol_list)[x]} Gene in txt format is stored in {usethis::ui_path(file_name)}"))

      })

      ## deg results
      dat <- dataDEG(obj,which = which)
      csv_name = glue('{dir}/{prefix}_{which}_results.csv')
      write.csv(dat,file = csv_name)
      ui_done(glue("{which} DEG results in csv format is stored in {usethis::ui_path(csv_name)}"))

      ## volcano plot
      volcano_plot <-  do.call("PointVolcano",c(list(obj = obj,which = i), PointVolcanoParam))
      volcano_file = glue("{dir}/{prefix}_{which}_volcano.pdf")
      ggsave(volcano_plot,filename = volcano_file, width = 1600,height = 1600,units = "px",limitsize = FALSE,device = cairo_pdf)
      ui_done(glue("{which} volcano results were store in {ui_path(volcano_file)}."))

      topHeatmap_file = glue("{dir}/{prefix}_{which}_top500_heatmap.pdf")
      topHeatmap <- do.call("DEGtopHeatmap",list(obj = obj,which = i,filename = topHeatmap_file))
      ui_done(glue("{which} top 500 heatmap were store in {ui_path(topHeatmap_file)}."))

    }

    ## top heatmap all
    heatmap_ls <- lapply(main, function(x){

      do.call("DEGtopHeatmap",list(obj = obj, which = x,filename = NA,show_gene = FALSE,
                                   legend = FALSE,annotation_legend = FALSE,annotation_names_col = F))

    })
    heatmap_all <- aplot::plot_list(gglist = heatmap_ls,labels = main,guides = 'collect')
    heatmap_all_file = glue("{dir}/{prefix}_ALL_top500heatmap.pdf")
    ggsave(filename = heatmap_all_file,plot = heatmap_all,device = cairo_pdf,width = 4.5*length(main),height = 4.5)
    ui_done(glue("DEG top heatmap plot were store in {ui_path(heatmap_all_file)}."))

    ## venn plot
    if (length(main) > 1) {
      index <- c(setdiff(label(obj),label_ns(obj)),"diff")

      geneSets <- lapply(main, function(x){
        geneSymbol_list <- hyper_GS(object = obj,which = x,type = "SYMBOL")
      })
      names(geneSets) <- main

      geneSets_ls <- list()
      for (i in index){
        tmp <- lapply(main, function(x){geneSets[[x]][[i]]})
        names(tmp) <- main
        geneSets_ls[[i]] <- tmp
      }

      p_list <- lapply(seq_along(geneSets_ls), function(x){
        DEGvenn(geneSets = geneSets_ls[[x]])
      })
      p_l <- aplot::plot_list(gglist = p_list,labels = index)

      venn_file = glue("{dir}/{prefix}_vennplot.pdf")
      ggsave(filename = venn_file,plot = p_l,device = cairo_pdf,width = 4.5*length(index),height = 4.5)
      ui_done(glue("DEG venn plot were store in {ui_path(venn_file)}."))
    }

  } else {
    ui_info("None available results of DEG.")
  }

}
