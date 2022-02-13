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

  ## gsva
  if (length(msigdbGSVAresult(obj)) == 0) {
    obj <- gsvaResolve(obj)
  }

  ## plot
  ## gse
  gse_list <- msigdbGSEAresult(obj)
  tmp <- lapply(msigdbParam(obj)[["category"]], function(j){

    plot_list <- lapply(seq_along(gse_list), function(i){

      y = gse_list[[i]]
      x = y[[j]]
      GSEAbar(x,top = top) +
        theme(legend.position="none") +
        ggtitle(names(gse_list)[i])

    })

    legend <- cowplot::get_legend(
      # create some space to the left of the legend
      plot_list[[1]] + ggplot2::theme(legend.position="right",legend.box.margin = ggplot2::margin(0, 0, 0, 12))
    )

    p <- cowplot::plot_grid(plotlist = plot_list,legend,ncol = 5, rel_widths = c(.4,3,3,3,3))

    cat <- switch (j,
      "H" = "HALLMARK",
      "C1" = "positional",
      "C2" = "curated",
      "C3" = "regulatory_target",
      "C4" = "computational",
      "C5" = "ontology",
      "C6" = "oncogenic_signature",
      "C7" = "immunologic_signature",
      "C8" = "cell_type_signature"
    )
    ggplot2::ggsave(p,filename = glue::glue("{dir}/{prefix}_GSEA_bar_{cat}.pdf"), width = 6400,height = 1200*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf)


    tmp <- lapply(seq_along(gse_list), function(i){

      y = gse_list[[i]]
      x = y[[j]]

      gseaplots_l <- GSEAplot(x,top =top)

      ## 保存gse富集趋势图到多页pdf
      mod <- names(gse_list)[i]
      down_plots_pdf <- glue("{dir}/{prefix}_GSEA_{mod}_{cat}_Down_gseplot.pdf")
      pdf(down_plots_pdf,height = 3,width = 4)
      invisible(lapply(gseaplots_l[["down_plots"]], print))
      dev.off()
      ui_done(glue("Enrichplot of GSEA head {top} {mod} {cat} in Down ploted in {ui_path({down_plots_pdf})}"))

      up_plots_pdf <- glue("{dir}/{prefix}_GSEA_{mod}_{cat}_Up_gseplot.pdf")
      pdf(up_plots_pdf,height = 3,width = 4)
      invisible(lapply(gseaplots_l[["up_plots"]], print))
      dev.off()
      ui_done(glue("Enrichplot of GSEA head {top} {mod} {cat} in Up ploted in {ui_path({up_plots_pdf})}"))

    })


  })

  ## gsva
  gsvares_list <- msigdbGSVAresult(obj)[["GSVA_matrix"]]
  lapply(seq_along(gsvares_list), function(x){

    pdf_file = glue("{dir}/{prefix}_GSVA_{names(gsvares_list)[x]}.pdf")
    ac=data.frame(Groups=groupInfo(obj))
    rownames(ac)=sampleNames(obj,filtered = T)
    pheatmap::pheatmap(gsvares_list[[x]],
                       annotation_col = ac,
                       filename = pdf_file)

  })


  return(obj)

})
