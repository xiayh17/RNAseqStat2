#' Summary Hyper results
#'
#' @param obj a DEGContainer
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param top top rows of up and down
#'
#' @importFrom usethis ui_info
#' @importFrom fs dir_create
#' @importFrom cowplot plot_grid get_legend
#' @importFrom ggplot2 ggsave
#'
#' @return files
#' @export
#'
#' @examples
#' hyperSummary(DEGContainer)
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

res_l <- modelEnrich(obj,dataBase = "KEGG",orderBy = "pvalue",head = top)

  tmp <- lapply(index, function(x){

    dat_d <- res_l[[x]]
    cc_file_name = glue("{dir}/{prefix}_KEGG_compareEnrichCircle_{x}.pdf")
    compareEnrichCircle(result_g = dat_d,filename = cc_file_name,mar = c(8,0,0,19))

  })

    tmp <- lapply(index, function(i){

      plot_list <- lapply(seq_along(subRes), function(j){

        x = subRes[[j]]
        eob <- x[[i]]

        if(!is.null(eob)){

            hyperBar(eob,top = top)+ theme(legend.position="none") + ggtitle(names(subRes)[j])


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

    res_l <- modelEnrich(obj,dataBase = "GO",orderBy = "pvalue",head = top)

    tmp <- lapply(index, function(x){

      dat_d <- res_l[[x]]
      cc_file_name = glue("{dir}/{prefix}_GO_compareEnrichCircle_{x}.pdf")
      compareEnrichCircle(result_g = dat_d,filename = cc_file_name,mar = c(8,0,0,19))

    })

    tmp <- lapply(index, function(i){

      plot_list <- lapply(seq_along(subRes), function(j){

        x = subRes[[j]]
        eob <- x[[i]]

        if(!is.null(eob)){

          enrichBar(eob,top = top,plot_title = names(subRes)[j])

        }

      })

      p <- cowplot::plot_grid(plotlist = plot_list, ncol = 4)

      plot_path = glue::glue("{dir}/{prefix}_GO_bar_{i}.pdf")
      ggplot2::ggsave(p,filename = plot_path, width = 7700,height = 1800*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf,dpi = 300)
      ui_done(glue("{i} Genes Hyper GO result is stored in {usethis::ui_path(plot_path)}"))

    })
  }

  tmp <- lapply(index, function(i){ ## Up Down diff

    lapply(seq_along(subRes), function(j){ ## limma edgeR DESeq2

      x = subRes[[j]]
      eob <- x[[i]]

      if (is.null(eob)) {

        ui_info("Hyper GO of {j} {i} Data not avalible")

      } else {

        res_name <- glue("{dir}/{prefix}_GO_{names(subRes)[j]}_{i}_Gene.csv")
        write.csv(eob@result,file = res_name)
        ui_done(glue("{names(subRes)[j]} {i} Hyper GO result in csv format is stored in {usethis::ui_path(res_name)}"))

      }

    })

  })

}
