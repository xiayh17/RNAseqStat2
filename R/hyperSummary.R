#' @export
hyperSummary <- function(obj, dir = ".", prefix = "3-runHyper",top = 10) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## get data
  # ui_info("Start plot")
  data <- hyperRes(obj = obj)

  ## KEGG data
  index <- c(setdiff(label(obj),label_ns(obj)),"diff")

  subRes <- data$hyperKEGG_res

  if (!is.null(subRes)) {
    tmp <- lapply(index, function(i){

      plot_list <- lapply(seq_along(subRes), function(j){

        x = subRes[[j]]
        eob <- x[[i]]

        if(!is.null(eob)){

          dat <- enrichplot:::fortify.enrichResult(model = eob, showCategory = top, by = "Count")
          if(nrow(dat) != 0) {

            barplot(eob)+ theme(legend.position="none") + ggtitle(names(subRes)[j])

          }

        }

      })

      legend <- cowplot::get_legend(
        # create some space to the left of the legend
        plot_list[[1]] + ggplot2::theme(legend.position="right",legend.box.margin = ggplot2::margin(0, 0, 0, 12))
      )

      p <- cowplot::plot_grid(plotlist = plot_list,legend,ncol = 5, rel_widths = c(.4,3,3,3,3))

      plot_path = glue::glue("{dir}/{prefix}_KEGG_{i}.pdf")
      ggplot2::ggsave(p,filename = plot_path, width = 6400,height = 1200*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf)

      ui_done(glue("{i} Genes Hyper KEGG result is stored in {usethis::ui_path(plot_path)}"))

    })
  }

  tmp <- lapply(index, function(i){ ## Up Down diff

    lapply(seq_along(subRes), function(j){ ## limma edgeR DESeq2

      x = subRes[[j]]
      eob <- x[[i]]

      res_name <- glue("{dir}/{prefix}_KEGG_{names(subRes)[j]}_{i}_Gene.csv")
      write.csv(eob@result,file = res_name)
      ui_done(glue("{names(subRes)[j]} {i} Hyper KEGG result in csv format is stored in {usethis::ui_path(res_name)}"))

    })

  })

  ## GO data
  subRes <- data$hyperGO_res

  if (!is.null(subRes)) {
    tmp <- lapply(index, function(i){

      plot_list <- lapply(seq_along(subRes), function(j){


        x = subRes[[j]]
        eob <- x[[i]]

        if(!is.null(eob)){

          dat <- enrichplot:::fortify.enrichResult(model = eob, showCategory = top, by = "Count")
          if(nrow(dat) != 0) {

            barplot(eob)+ theme(legend.position="none") + ggtitle(names(subRes)[j])

          }

        }

      })

      legend <- cowplot::get_legend(
        # create some space to the left of the legend
        plot_list[[1]] + ggplot2::theme(legend.position="right",legend.box.margin = ggplot2::margin(0, 0, 0, 12))
      )

      p <- cowplot::plot_grid(plotlist = plot_list,legend,ncol = 5, rel_widths = c(.4,3,3,3,3))

      plot_path = glue::glue("{dir}/{prefix}_GO_{i}.pdf")
      ggplot2::ggsave(p,filename = plot_path, width = 6400,height = 1200*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf)

      ui_done(glue("{i} Genes Hyper GO result is stored in {usethis::ui_path(plot_path)}"))

    })
  }

  tmp <- lapply(index, function(i){ ## Up Down diff

    lapply(seq_along(subRes), function(j){ ## limma edgeR DESeq2

      x = subRes[[j]]
      eob <- x[[i]]

      res_name <- glue("{dir}/{prefix}_GO_{names(subRes)[j]}_{i}_Gene.csv")
      write.csv(eob@result,file = res_name)
      ui_done(glue("{names(subRes)[j]} {i} Hyper GO result in csv format is stored in {usethis::ui_path(res_name)}"))

    })

  })

}