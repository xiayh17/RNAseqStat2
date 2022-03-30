#' Summary MSigDB modules results
#'
#' @param obj a DEGContainer
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param top top of hyper or gsea
#'
#' @return
#' @export
#'
#' @examples
#' MSigDBSummary(data_msigdb)
MSigDBSummary <- function(obj, dir = ".", prefix = "5-runMSigDB",top =10) {

  ## plot
  ## gse ----
  gse_list <- msigdbGSEAresult(obj)

  if(length(gse_list) != 0){

    tmp <- lapply(msigdbParam(obj)[["category"]], function(j){

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

      tmp <- lapply(seq_along(gse_list), function(i){

        y = gse_list[[i]]
        x = y[[j]]
        dat = x@result

        res_name <- glue("{dir}/{prefix}_GSEA_{names(gse_list)[i]}_{j}.csv")
        write.csv(dat,file = res_name)
        ui_done(glue("{names(gse_list)[i]} {j} GSEA result in csv format is stored in {usethis::ui_path(res_name)}"))

      })

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

  }
  ## ----

  ## hyper ----
  hyper_list <- msigdbHyperResult(obj)

  index <- c(setdiff(label(obj),label_ns(obj)),"diff")

  if (length(hyper_list) != 0) {

    tmp <- lapply(index, function(i){

      lapply(msigdbParam(obj)[["category"]], function(j){

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

        plot_list <- lapply(seq_along(hyper_list), function(k){

          x = hyper_list[[k]]
          set <- x[[j]]
          eob <- set[[i]]

          res_name <- glue("{dir}/{prefix}_Hyper_{names(hyper_list)[k]}_{j}_{i}.csv")
          write.csv(eob@result,file = res_name)
          ui_done(glue("{names(hyper_list)[k]} {j} {i} Hyper result in csv format is stored in {usethis::ui_path(res_name)}"))

          if(!is.null(eob)){
            hyperBar(eob,top = top)+ theme(legend.position="none") + ggtitle(names(hyper_list)[j])
          }

        })

        legend <- cowplot::get_legend(
          # create some space to the left of the legend
          plot_list[[1]] + ggplot2::theme(legend.position="right",legend.box.margin = ggplot2::margin(0, 0, 0, 12))
        )

        p <- cowplot::plot_grid(plotlist = plot_list,legend,ncol = 5, rel_widths = c(.4,3,3,3,3))

        plot_path = glue::glue("{dir}/{prefix}_Hyper_{j}_{i}.pdf")
        ggplot2::ggsave(p,filename = plot_path, width = 6400,height = 1200*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf)

        ui_done(glue("{j} {i} Genes Hyper result is stored in {usethis::ui_path(plot_path)}"))

      })

    })

  }


  ## ----

  ## gsva
  gsvares_list <- msigdbGSVAresult(obj)[["GSVA_matrix"]]
  invisible(lapply(seq_along(gsvares_list), function(x){

    pdf_file = glue("{dir}/{prefix}_GSVA_{names(gsvares_list)[x]}.pdf")
    ac=data.frame(Groups=groupInfo(obj))
    rownames(ac)=sampleNames(obj,filtered = T)
    pheatmap::pheatmap(gsvares_list[[x]],
                       annotation_col = ac,
                       filename = pdf_file)
    ui_done("{names(gsvares_list)[x]} GSVA heatmap is stored in {usethis::ui_path(pdf_file)}")

    csv_file = glue("{dir}/{prefix}_GSVA_{names(gsvares_list)[x]}_matrix.csv")
    write.csv(gsvares_list[[x]],file = csv_file)
    ui_done("{names(gsvares_list)[x]} GSVA matrix is stored in {usethis::ui_path(csv_file)}")

  }))

  gsvadiff_list <- msigdbGSVAresult(obj)[["GSVA_diff"]]
  invisible(lapply(seq_along(gsvadiff_list), function(x){

    csv_file = glue("{dir}/{prefix}_GSVA_{names(gsvadiff_list)[x]}_limmaDEG.csv")
    write.csv(gsvadiff_list[[x]],file = csv_file)
    ui_done("{names(gsvadiff_list)[x]} GSVA analysis by limma is stored in {usethis::ui_path(csv_file)}")

    volcano_file = glue("{dir}/{prefix}_GSVA_{names(gsvadiff_list)[x]}_volcano.pdf")
    p <- PointVolcano(object = obj,which = "MSigDB",category = names(gsvadiff_list)[x],gene = 5,expend = c(0.4,0.4))
    ggsave(p,filename = volcano_file, width = 1600,height = 1600,units = "px",limitsize = FALSE,device = cairo_pdf)
    ui_done("Volcano of {names(gsvadiff_list)[x]} GSVA analysis by limma is plot in {usethis::ui_path(volcano_file)}")

  }))

}
