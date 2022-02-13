#' @import enrichplot
#' @import ggplot2
#' @importFrom DOSE theme_dose
#' @export
setGeneric(name="runGSEA", def=function(obj, dir = ".", prefix = "4-runGSEA",top =10) standardGeneric("runGSEA"))

setMethod(f="runGSEA", signature="DEGContainer", definition=function(obj, dir = ".", prefix = "4-runGSEA",top = 10) {

  ## hyper resolve of limma edgeR DESeq2

  if (length(gseRes(obj)) == 0) {
    obj <- gseResolve(obj = obj)
  }

  ## get data
  ui_info("Start plot")
  res_list <- gseRes(obj = obj)$gseKEGG_res

  # ## KEGG data
  plot_list <- lapply(seq_along(res_list), function(i){
    x = res_list[[i]]
    GSEAbar(x,top = top) +
      theme(legend.position="none") +
      ggtitle(names(res_list)[i])

  })

  legend <- cowplot::get_legend(
    # create some space to the left of the legend
    plot_list[[1]] + ggplot2::theme(legend.position="right",legend.box.margin = ggplot2::margin(0, 0, 0, 12))
  )

  p <- cowplot::plot_grid(plotlist = plot_list,legend,ncol = 5, rel_widths = c(.4,3,3,3,3))

  ggplot2::ggsave(p,filename = glue::glue("{dir}/{prefix}_KEGG_bar.pdf"), width = 6400,height = 1200*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf)

  tmp <- lapply(seq_along(res_list), function(i){

    x = res_list[[i]]

    gseaplots_l <- GSEAplot(x,top =top)

    ## 保存gse富集趋势图到多页pdf
    mod <- names(res_list)[i]
    down_plots_pdf <- glue("{dir}/{prefix}_KEGG_{mod}_Down_gseplot.pdf")
    pdf(down_plots_pdf,height = 3,width = 4)
    invisible(lapply(gseaplots_l[["down_plots"]], print))
    dev.off()
    ui_done(glue("Enrichplot of KEGG head {top} {mod} in Down ploted in {ui_path({down_plots_pdf})}"))

    up_plots_pdf <- glue("{dir}/{prefix}_KEGG_{mod}_Up_gseplot.pdf")
    pdf(up_plots_pdf,height = 3,width = 4)
    invisible(lapply(gseaplots_l[["up_plots"]], print))
    dev.off()
    ui_done(glue("Enrichplot of KEGG head {top} {mod} in Up ploted in {ui_path({up_plots_pdf})}"))

  })

  return(obj)

})
