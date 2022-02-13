#' @import enrichplot
#' @export
setGeneric(name="runHyper", def=function(obj, dir = ".", prefix = "3-runHyper",top = 10) standardGeneric("runHyper"))

setMethod(f="runHyper", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "3-runHyper",top = 10) {

  ## hyper resolve of limma edgeR DESeq2

  if (length(hyperRes(obj)) == 0) {
    obj <- hyperResolve(obj = obj)
  }

  ## get data
  ui_info("Start plot")
  data <- hyperRes(obj = obj)

  ## KEGG data
  index <- c(setdiff(label(obj),label_ns(obj)),"diff")
  lapply(index, function(i){

    plot_list <- lapply(seq_along(data$hyperKEGG_res), function(j){


      x = data$hyperKEGG_res[[j]]
      eob <- x[[i]]

      if(!is.null(eob)){

        dat <- enrichplot:::fortify.enrichResult(model = eob, showCategory = top, by = "Count")
        if(nrow(dat) != 0) {

          barplot(eob)+ theme(legend.position="none") + ggtitle(names(data$hyperKEGG_res)[j])

        }

      }



    })

    legend <- cowplot::get_legend(
      # create some space to the left of the legend
      plot_list[[1]] + ggplot2::theme(legend.position="right",legend.box.margin = ggplot2::margin(0, 0, 0, 12))
    )

    p <- cowplot::plot_grid(plotlist = plot_list,legend,ncol = 5, rel_widths = c(.4,3,3,3,3))

    ggplot2::ggsave(p,filename = glue::glue("{dir}/{prefix}_{i}.pdf"), width = 6400,height = 1200*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf)

  })

  # ## group DEG Results
  # obj <- degGroup(obj = obj)
  #
  # test <- deg_here(obj)
  #
  # ok <- names(test)[which(test == TRUE)]
  #
  # main <- setdiff(ok,"merge")
  # merge_data <- intersect(ok,"merge")
  #
  # ## plot and export volcano
  # if(length(main)>=1){
  #
  #   for (i in main) {
  #
  #     which = i
  #     volcano_plot <-  PointVolcano(obj = obj,
  #                                   which = which,
  #                                   gene = gene,
  #                                   light = light,
  #                                   label_light = label_light,
  #                                   light_color = light_color,
  #                                   light_label_color = light_label_color,
  #                                   expend = expend)
  #
  #     volcano_file = glue("{dir}/{prefix}_{which}_volcano.pdf")
  #     ggsave(volcano_plot,filename = volcano_file, width = 1600,height = 1600,units = "px",limitsize = FALSE,device = cairo_pdf)
  #     ui_done(glue("{which} volcano results were store in {ui_path(volcano_file)}."))
  #
  #   }
  #
  # } else {
  #   ui_info("None available results of DEG.")
  # }

  return(obj)

})
