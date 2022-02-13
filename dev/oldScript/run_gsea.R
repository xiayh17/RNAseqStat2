#' run enrichgesKEGG
#'
#' run enrichgesKEGG and make output files
#'
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#' @inheritParams gsea_kegg
#' @inheritParams gesa_barplot
#'
#' @importFrom fs dir_exists
#' @importFrom usethis ui_info
#' @importFrom glue glue
#'
#' @return a list of file
#' @export
#'
#' @examples
#' \dontrun{
#' run_gseKEGG(deg_data = DEG_df, x = "log2FoldChange", dir = tempdir(),eps = 0)
#' }
run_gseKEGG <- function(deg_data,x,dir= ".", prefix = "5-GSEA",
                              pvalue_cut = 0.1, enrichmentScore_cut = 0.5, top = 10,
                              organism = "hsa",
                              keyType = "kegg",
                              OrgDb = 'org.Hs.eg.db',
                              exponent = 1,
                              minGSSize = 10,
                              maxGSSize = 500,
                              eps = 1e-10,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              verbose = TRUE,
                              use_internal_data = FALSE,
                              seed = FALSE,
                              by = "fgsea") {

  if (!dir_exists(dir)) {
    dir_create(dir)
  }

  usethis::ui_info(glue("Start GSEA of KEGG. This process will take a few minutes."))

  gsekegg_res <- gsea_kegg(deg_data = deg_data, x = x,
                           organism = organism,
                           keyType = keyType,
                           OrgDb = OrgDb,
                           exponent = exponent,
                           minGSSize = minGSSize,
                           maxGSSize = maxGSSize,
                           eps = eps,
                           pvalueCutoff = pvalueCutoff,
                           pAdjustMethod = pAdjustMethod,
                           verbose = verbose,
                           use_internal_data = use_internal_data,
                           seed = seed,
                           by = by)

  gsea_result_summary(gseaResult = gsekegg_res,
                      type = "KEGG",
                      prefix = prefix,
                      dir = dir,
                      top = top,
                      pvalue_cut = pvalue_cut,
                      enrichmentScore_cut = enrichmentScore_cut)

}

#' run enrich gesOG
#'
#' run enrich gesGO and make output files
#'
#' @param dir where to save results files
#' @param prefix a prefix of file names in this step
#' @inheritParams gsea_go
#' @inheritParams gesa_barplot
#'
#' @importFrom fs dir_exists
#' @importFrom usethis ui_info
#' @importFrom glue glue
#'
#' @return a list of file
#' @export
#'
#' @examples
#' \dontrun{
#' run_gseGO(deg_data = DEG_df, x = "log2FoldChange", dir = tempdir())
#' }
run_gseGO <- function(deg_data,x,dir= ".", prefix = "5-GSEA",
                      pvalue_cut = 0.1, enrichmentScore_cut = 0.5, top = 10,
                        ont = "BP",
                        keyType = "SYMBOL",
                        OrgDb = 'org.Hs.eg.db',
                        exponent = 1,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 1e-10,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        verbose = TRUE,
                        seed = FALSE,
                        by = "fgsea") {

  if (!dir_exists(dir)) {
    dir_create(dir)
  }

  usethis::ui_info(glue("Start GSEA of GO. This process will take a few minutes."))

  gsego_res <- gsea_go(deg_data = deg_data,
                        x = x,
                        ont = ont,
                        keyType = keyType,
                        OrgDb = OrgDb,
                        exponent = exponent,
                        minGSSize = minGSSize,
                        maxGSSize = maxGSSize,
                        eps = eps,
                        pvalueCutoff = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        verbose = verbose,
                        seed = seed,
                        by = by)

  gsea_result_summary(gseaResult = gsego_res,
                      type = "GO",
                      prefix = prefix,
                      dir = dir,
                      top = top,
                      pvalue_cut = pvalue_cut,
                      enrichmentScore_cut = enrichmentScore_cut)

}

#' @inheritParams gesa_barplot
#' @inheritParams enhance_gseplot
#' @importFrom usethis ui_done
#' @importFrom glue glue
#' @importFrom utils write.csv
gsea_result_summary <- function(gseaResult,
                                type,
                                prefix = "5-GSEA",
                                dir = ".",
                                top = 10,
                                pvalue_cut = 0.1,
                                enrichmentScore_cut = 0.5) {

  # 保存Rdata
  gsea_res_file <- glue("{dir}/{prefix}_{type}_result.Rdata")
  save(gseaResult,file = gsea_res_file)
  ui_done(glue("Result of GSEA_{type} stored in {ui_path({gsea_res_file})}"))

  ## 保存csv
  gsea_csv <- glue('{dir}/{prefix}_{type}_result.csv')
  write.csv(gseaResult@result,file = gsea_csv)
  ui_done(glue("Data in CSV format of GSEA_{type} writed in {ui_path({gsea_res_file})}"))

  ## 绘制和保存柱状图
  gsea_bar <- glue("{dir}/{prefix}_{type}_barplot.pdf")
  text_w <- max(strwidth(gseaResult@result$Description,units = "inch"))
  p = gesa_barplot(gseaResult,pvalue_cut = pvalue_cut, enrichmentScore_cut = enrichmentScore_cut, top = top)
  ggsave(plot = p,filename = gsea_bar,height = 5, width = text_w+3)
  ui_done(glue("Barplot of head {top} GSEA_{type} ploted in {ui_path({gsea_bar})}"))

  ## 绘制gse富集趋势图
  plots_l <- enhance_gseplot(gseaResult, top = top,
                             pvalue_cut = pvalue_cut, enrichmentScore_cut = enrichmentScore_cut)

  # 试图拼图
  # li_down = structure(plots_l[["down_plots"]], class = c("gglist", "ggplot"))
  # print.gglist = function(x, ...) {plyr::l_ply(x, print, ...)}
  # ggsave(li_down, file = glue("{dir}/{prefix}_gseKEGG_down_gseplot.pdf"))
  #
  # li_up = structure(plots_l[["up_plots"]], class = c("gglist", "ggplot"))
  # ggsave(li_up, file = glue("{dir}/{prefix}_gseKEGG_up_gseplot.pdf"))

  ## 保存gse富集趋势图到多页pdf
  down_plots_pdf <- glue("{dir}/{prefix}_{type}_Down_gseplot.pdf")
  pdf(down_plots_pdf,height = 3,width = 4)
  invisible(lapply(plots_l[["down_plots"]], print))
  dev.off()
  ui_done(glue("Enrichplot of head {top} GSEA_{type} in Down ploted in {ui_path({down_plots_pdf})}"))

  up_plots_pdf <- glue("{dir}/{prefix}_{type}_Up_gseplot.pdf")
  pdf(up_plots_pdf,height = 3,width = 4)
  invisible(lapply(plots_l[["up_plots"]], print))
  dev.off()
  ui_done(glue("Enrichplot of head {top} GSEA_{type} in Up ploted in {ui_path({up_plots_pdf})}"))

}
