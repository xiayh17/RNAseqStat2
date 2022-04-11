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

  pro = NULL

  if (length(hyperRes(obj)[["hyperKEGG_res"]]) > 0) {
    pro = "KEGG"
  }

  if (length(hyperRes(obj)[["hyperGO_res"]]) > 0) {
    pro = c(pro,"GO")
  }

  if(is.null(pro)){

    ui_info("NO KEGG or GO data found!")

  } else {

    lapply(pro, function(x){
      ui_info("Start to Summary Hyper {x}")

      ## Bar plot of Count, order by pvalue, fill by qvalue
      tryCatch(
        expr = {
          hyperSummaryHyperBar(obj = obj,dataBase = x,top = top,dir = dir,prefix = prefix)
        },
        error = function(e){
          usethis::ui_oops("Something wrong occured in hyperSummaryHyperBar {x}.
                       try again later by {ui_code('hyperSummaryHyperBar')}.")
        }
      )

      ## Circlize plot of KEGG ids in Different methods
      tryCatch(
        expr = {
          hyperSummaryHyperCompareCircle(obj = obj,dataBase = x,orderBy = "pvalue",top = top,dir = dir,prefix = prefix)
        },
        error = function(e){
          usethis::ui_oops("Something wrong occured in hyperSummaryHyperCompareCircle {x}.
                       try again later by {ui_code('hyperSummaryHyperCompareCircle')}.")
        }
      )

      ## CSV files for results
      tryCatch(
        expr = {
          hyperSummaryCSV(obj = obj, dataBase = x, dir = dir, prefix = prefix)
        },
        error = function(e){
          usethis::ui_oops("Something wrong occured in hyperSummaryCSV {x}.
                       try again later by {ui_code('hyperSummaryCSV')}.")
        }
      )

      if(identical(x,"GO")){
        tryCatch(
          expr = {
            ## GO enrich barplot Split by Ontology
            hyperSummaryGOBar(obj = obj, top = top, dir = dir, prefix = prefix)
          },
          error = function(e){
            usethis::ui_oops("Something wrong occured in hyperSummaryGOBar.
                       try again later by {ui_code('hyperSummaryGOBar')}.")
          }
        )
      }

    })

  }

}

## Hyper Summary  Enrich bar plot
hyperSummaryGOBar <- function(obj,top = 10,dir = ".",prefix = "3-runHyper"){

  index <- c(setdiff(label(obj),label_ns(obj)),"diff")
  ## get data
  # ui_info("Start plot")
  data <- hyperRes(obj = obj)

  subRes <- data$hyperGO_res

  if (is.null(subRes)) {

    ui_oops("NO data avalible of {ui_value('GO')}")

  } else {

    tmp <- lapply(index, function(i) {

      plot_list <- lapply(seq_along(subRes), function(j){

        x = subRes[[j]]
        eob <- x[[i]]

        if(!is.null(eob)){

          enrichBar(eob,top = top,plot_title = names(subRes)[j])

        } else {

          ui_info("Hyper GO of {j} {i} Data not avalible")

        }

      })

      p <- patchwork::wrap_plots(plot_list,nrow = 1,guides = "collect")

      plot_path = glue::glue("{dir}/{prefix}_GO_{i}_OntologySplitBar.pdf")
      ggplot2::ggsave(p,filename = plot_path, width = 2566.66*seq_along(subRes),height = 1800*(top*2/10),units = "px",limitsize = FALSE,device = cairo_pdf,dpi = 300)
      ui_done(glue("{i} Genes Hyper GO enrich barplot Split by Ontology is stored in {usethis::ui_path(plot_path)}"))

    })

  }

}

## CSV files for results
hyperSummaryCSV <- function(obj,dataBase,dir = ".",prefix = "3-runHyper"){

  if (missing(dataBase)) {

    ui_oops("Please select from {ui_value('KEGG')} or {ui_value('GO')} for ui_code('dataBase')")

  }

  index <- c(setdiff(label(obj),label_ns(obj)),"diff")
  ## get data
  # ui_info("Start plot")
  data <- hyperRes(obj = obj)

  if (dataBase == "KEGG") {
    subRes <- data$hyperKEGG_res
  } else if (dataBase == "GO") {
    subRes <- data$hyperGO_res
  } else {

    ui_oops("Please select from {ui_value('KEGG')} or {ui_value('GO')} for ui_code('dataBase')")

  }

  if (is.null(subRes)) {

    ui_oops("NO data avalible of {ui_value(dataBase)}")

  } else {

    tmp <- lapply(index, function(i){ ## Up Down diff

      lapply(seq_along(subRes), function(j){ ## limma edgeR DESeq2

        x = subRes[[j]]
        eob <- x[[i]]

        if (is.null(eob)) {

          ui_info("Hyper {dataBase} of {j} {i} Data not avalible")

        } else {

          res_name <- glue("{dir}/{prefix}_{dataBase}_{names(subRes)[j]}_{i}_Gene.csv")
          write.csv(eob@result,file = res_name)
          ui_done(glue("{names(subRes)[j]} {i} Hyper {dataBase} result in csv format is stored in {usethis::ui_path(res_name)}"))

        }

      })

    })

  }

}

## Circlize plot of KEGG and GO in Different methods ----
hyperSummaryHyperCompareCircle <- function(obj,dataBase,orderBy = "pvalue",top = 10,dir = ".",prefix = "3-runHyper"){

  index <- c(setdiff(label(obj),label_ns(obj)),"diff")
  ## check deg results
  test <- deg_here(obj)
  ok <- names(test)[which(test == TRUE)]
  ## except merge
  main <- setdiff(ok,"merge")

  ## convert data format for compareEnrichCircle
  res_l <- modelEnrich(obj,dataBase = dataBase,orderBy = orderBy,head = top)

  if (is.null(res_l)) {
    ui_oops("NO data avaliable for {ui_code('compareEnrichCircle')}")
  } else {
    tmp <- lapply(index, function(x){
      tryCatch(
        expr = {
          dat_d <- res_l[[x]]
          cc_file_name = glue("{dir}/{prefix}_{dataBase}_compareEnrichCircle_{x}.pdf")
          compareEnrichCircle(result_g = dat_d,filename = cc_file_name,mar = c(8,0,0,19))
        },
        error = function(e){
          usethis::ui_oops("Something wrong occured in {ui_code('compareEnrichCircle')} of {x}")
        }
      )
    })
  }

}

## barplot of KEGG and GO ----
hyperSummaryHyperBar <- function(obj,dataBase,top = 10,dir=".",prefix = "3-runHyper"){

  if (missing(dataBase)) {

    ui_oops("Please select from {ui_value('KEGG')} or {ui_value('GO')} for ui_code('dataBase')")

  }

  index <- c(setdiff(label(obj),label_ns(obj)),"diff")
  ## get data
  # ui_info("Start plot")
  data <- hyperRes(obj = obj)

  if (dataBase == "KEGG") {
    subRes <- data$hyperKEGG_res
  } else if (dataBase == "GO") {
    subRes <- data$hyperGO_res
  } else {

    ui_oops("Please select from {ui_value('KEGG')} or {ui_value('GO')} for ui_code('dataBase')")

  }

  if (is.null(subRes)) {

    ui_oops("NO data avalible of {ui_value(dataBase)}")

  } else {

    ## Down up diff split
    tmp <- lapply(index, function(i){

      ## plot for every results
      plot_list <- lapply(seq_along(subRes), function(j){

        tryCatch(
          expr = {

            ## get limma edgeR DESeq2
            x = subRes[[j]]
            ## get up down diff
            eob <- x[[i]]
            ## plot
            if(!is.null(eob)) {
              hyperBar(eob,top = top) + ggtitle(names(subRes)[j])
            } else {
              ui_oops("NO data avalible of {ui_value(dataBase)} {ui_value(names(subRes)[j])} {ui_value(i)}")
            }

          },
          error = function(e){
            usethis::ui_oops("Something wrong occured in {ui_code('hyperBar')} {ui_value(names(subRes)[j])} {ui_value(i)}.")
          }
        )

      })

      ## collect plots
      p <- patchwork::wrap_plots(plot_list,nrow = 1,guides = "collect")

      ## save plot
      plot_path = glue::glue("{dir}/{prefix}_{dataBase}_{i}.pdf")
      ggplot2::ggsave(p,filename = plot_path, width = 6400,height = 1200*(top*2.5/10),units = "px",limitsize = FALSE,device = cairo_pdf)

      ui_done(glue("{i} Genes Hyper {dataBase} barplot is stored in {usethis::ui_path(plot_path)}"))

    })

  }

}
