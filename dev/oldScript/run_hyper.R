## GO and KEGG at the same time
run_hyper <- function(deg_data, x, y, optionsGO, optionsKEGG) {

  ## 设置默认选项 ----
  optionsGO_default = list(cut_FC = 1, cut_FDR = 0.05, showCategory = 10, dir = ".", prefix = "3-EnrichGO",
                OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "ALL", simplify = TRUE,
                pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,minGSSize = 10,
                maxGSSize = 500, readable = FALSE, pool = FALSE,
                label = c("Down", "Stable", "Up"),
                label_ns = "Stable",
                mc.cores = 1L)

  optionsKEGG_default = list(cut_FC = 1, cut_FDR = 0.05, top = 10,
                  dir= ".", prefix = "4-EnrichKEGG",
                  organism = NULL,
                  keyType = "kegg",
                  OrgDb = 'org.Hs.eg.db',
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  minGSSize = 10,
                  maxGSSize = 500,
                  qvalueCutoff = 0.2,
                  use_internal_data = FALSE,
                  label = c("Down", "Stable", "Up"),
                  label_ns = "Stable",down_label = NULL,
                  mc.cores = 1L)

  if (missing(optionsGO))
    optionsGO <- optionsGO_default

  if (missing(optionsKEGG))
    optionsKEGG <- optionsKEGG_default

  if (is.null(optionsKEGG[['organism']]))
    optionsKEGG[['organism']] = switch(optionsKEGG[["OrgDb"]],
                                       'org.Hs.eg.db' = 'hsa',
                                       'org.Mm,eg.db' = 'mmu'
    )

  if (is.null(optionsKEGG[['down_label']]))
    optionsKEGG[['down_label']] = optionsKEGG[['label']][1]
  ## ----

  ## 读取用户选项 ----
  ## Go
  if (is.list(optionsGO)&length(optionsGO) > 0) { ## 判断格式是否正确

    if (!all(names(optionsGO) %in% names(optionsGO_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('optionsGO')} will not be applied.\n Please check the parameters of {usethis::ui_code('run_hyperGO')} ")

    } else if (length(optionsGO) == 1) { ## 判断 长度为1

      optionsGO_default[[names(optionsGO)]] <- optionsGO[names(optionsGO)][[1]]

    } else if (length(optionsGO) > 1) { ## 其他长度

       for (x in seq_along(optionsGO)) {

        optionsGO_default[[names(optionsGO)[x]]] <- optionsGO[names(optionsGO)[x]][[1]]

       }

      optionsGO <- optionsGO_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('optionsGO')} will not be applied.\n Please check the usage of {usethis::ui_code('run_hyper')} ")

  }

  ## KEGG
  if (is.list(optionsKEGG)&length(optionsKEGG) > 0) { ## 判断格式是否正确

    if (!all(names(optionsKEGG) %in% names(optionsKEGG_default))) { ## 判断有效性

      usethis::ui_oops("Your {usethis::ui_code('optionsGO')} will not be applied.\n Please check the parameters of {usethis::ui_code('run_hyperGO')} ")

    } else if (length(optionsKEGG) == 1) { ## 判断 长度为1

      optionsKEGG_default[[names(optionsKEGG)]] <- optionsKEGG[names(optionsKEGG)][[1]]

    } else if (length(optionsKEGG) > 1) { ## 其他长度

      for (x in seq_along(optionsKEGG)) {

        optionsKEGG_default[[names(optionsKEGG)[x]]] <- optionsKEGG[names(optionsKEGG)[x]][[1]]

      }

      optionsKEGG <- optionsKEGG_default

    }


  } else {

    usethis::ui_oops("Your {usethis::ui_code('optionsGO')} will not be applied.\n Please check the usage of {usethis::ui_code('run_hyper')} ")

  }
  ## ----

  ## 打印最后的关键选项
  ## 打印差异分析分组选项 FC P
  usethis::ui_info(glue::glue(" {optionsGO[['cut_FC']]}"))

  run_hyperGO(deg_data = deg_data, x = x, y = y, !!!optionsGO)
  run_hyperKEGG(deg_data = deg_data, x = x, y = y, !!!optionsKEGG)

}

#' Enrichment analysis of go
#'
#' run hyper_go and plot results
#'
#' @inheritParams hyper_go
#' @param showCategory Category numbers to show
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#'
#' @importFrom utils write.table
#' @importFrom graphics strwidth
#'
#' @return a list result files of GO
#' @export
#'
#' @examples
#' \dontrun{
#' run_hyperGO(DEG_df, x = "log2FoldChange", y = "pvalue", dir = tempdir())
#' }
run_hyperGO <- function(deg_data, x, y, cut_FC = 1, cut_FDR = 0.05, showCategory = 10, dir = ".", prefix = "3-EnrichGO",
                        OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL", ont = "ALL", simplify = TRUE,
                        pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,minGSSize = 10,
                        maxGSSize = 500, readable = FALSE, pool = FALSE,
                        label = c("Down", "Stable", "Up"),
                        label_ns = "Stable",
                        mc.cores = 1L) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  usethis::ui_info(glue("Enrich GO analysis Start. This process will take a few minutes."))

  go_resl <- hyper_go(deg_data = deg_data, x = x, y = y, cut_FC = cut_FC, cut_FDR = cut_FDR,
                      OrgDb = OrgDb, keyType = keyType, ont = ont, simplify = simplify,
                      pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,minGSSize = minGSSize,
                      maxGSSize = maxGSSize, readable = readable, pool = pool,
                      label = label,
                      label_ns = label_ns,
                      mc.cores = 1L)

  go_res_file <- glue("{dir}/{prefix}_go_result.Rdata")
  save(go_resl,file = go_res_file)
  usethis::ui_done(glue::glue("Result of enrichGO stored in {go_res_file}"))

  go_csv_file <- glue::glue('{dir}/{prefix}_gene_{names(go_resl)[[x]]}_GO_enrichment.csv')
  tmp <- lapply(seq_along(go_resl), function(x)
    write.table(go_resl[[x]]@result,file = go_csv_file)
  )
  usethis::ui_done(glue::glue("Result of enrichGO in csv format stored in {go_csv_file}"))

  if (ont == "ALL") {
    plots <- lapply(go_resl, function(x)
      enhance_barplot(x@result,showCategory=showCategory,split = 'ONTOLOGY')
    )
    text_w_l <- lapply(go_resl, function(x)
      max(strwidth(x@result$Description,units = "inch"))
    )
    lapply(seq_along(plots), function(x) {

      text_w <- text_w_l[[x]]
      ggsave(plot = plots[[x]],
             filename = glue::glue("{dir}/{prefix}_barplot-gene_{names(plots)[[x]]}_GO_enrichment.pdf"),
             height = 880*2.5/300, width = text_w+1, dpi = 300)

    }
    )
  } else {
    plots <- lapply(go_resl, function(x)
      enhance_barplot(x@result,showCategory=showCategory)
    )
    text_w_l <- lapply(go_resl, function(x)
      max(strwidth(x@result$Description,units = "inch"))
    )
    lapply(seq_along(plots), function(x) {

      text_w <- text_w_l[[x]]
      ggsave(plot = plots[[x]], filename = glue::glue("{dir}/{prefix}_barplot-gene_{names(plots)[[x]]}_GO_enrichment.pdf"),
             height = 880*2.5/300/3, width = text_w+1, dpi = 300)

    }

    )
  }

  usethis::ui_done(glue::glue("Result of enrichGO ploted in {dir}/{prefix}.*_GO_enrichment.pdf"))

}

#' Enrichment analysis of kegg
#'
#' run hyper_kegg and plot results
#'
#' @inheritParams hyper_kegg
#' @inheritParams kegg_barplot
#' @param top top rows for up and down
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#'
#' @return a list result files of KEGG
#' @export
#'
#' @examples
#' \dontrun{
#' run_hyperKEGG(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue", cut_FC = 1,
#' cut_FDR = 0.05, top = 10)
#' }
run_hyperKEGG <- function(deg_data, x, y, cut_FC = 1, cut_FDR = 0.05, top = 10,
                          dir= ".", prefix = "4-EnrichKEGG",
                          organism = switch (OrgDb,
                                             'org.Hs.eg.db' = 'hsa',
                                             'org.Mm,eg.db' = 'mmu'
                          ),
                          keyType = "kegg",
                          OrgDb = 'org.Hs.eg.db',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 500,
                          qvalueCutoff = 0.2,
                          use_internal_data = FALSE,
                          label = c("Down", "Stable", "Up"),
                          label_ns = "Stable",
                          mc.cores = 1L,down_label = "Down"
) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }
  # hyper_kegg
  kegg_resl <- hyper_kegg(deg_data = deg_data, x = x, y = y, cut_FC = cut_FC,
                          cut_FDR = cut_FDR,
                          organism = organism,
                          keyType = keyType,
                          OrgDb = OrgDb,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          qvalueCutoff = qvalueCutoff,
                          use_internal_data = use_internal_data,
                          label = label,
                          label_ns = label_ns,
                          mc.cores = mc.cores)

  kegg_res_file <- glue("{dir}/{prefix}_kegg_result.Rdata")
  save(kegg_resl,file = kegg_res_file)
  usethis::ui_done(glue::glue("Result of hyper_kegg stored in {kegg_res_file}"))

  text_w <- max(strwidth(kegg_resl$Down@result$Description,units = "inch"),
                strwidth(kegg_resl$Up@result$Description,units = "inch"))

  p <- kegg_barplot(kegg_resl, top = top, down_label = down_label)
  kegg_barplot_filename = glue::glue("{dir}/{prefix}_up_and_down_KEGG.pdf")
  ggsave(plot = p, filename = kegg_barplot_filename,height = 5, width = text_w + 3)
  usethis::ui_done(glue::glue("Result of enrichKEGG ploted in {kegg_barplot_filename}"))

}


