#' hyper_go
#'
#' a funtion group deg data frame and make enrichment analysis in one step
#'
#' @inheritParams cut_much
#' @inheritParams enhance_enrichGO
#' @param mc.cores param for mclapply, choose cores to use in non-Windows machine
#'
#' @importFrom parallel mclapply
#' @importFrom glue glue
#'
#' @return enrichGO result
#' @export
#'
#' @examples
#' \dontrun{
#' hyper_go(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue")
#' }
hyper_go <- function(deg_data, x, y, cut_FC = 1, cut_FDR = 0.05,
                      OrgDb = 'org.Hs.eg.db',
                      keyType = "SYMBOL",
                      ont = "ALL",
                      simplify = TRUE,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.2,
                      minGSSize = 10,
                      maxGSSize = 500,
                      readable = FALSE,
                      pool = FALSE,
                      label = c("Down", "Stable", "Up"),
                      label_ns = "Stable",
                      mc.cores = 1L) {

  deg_df_g <- cut_much(deg_data, x = x, y = y,cut_FC = cut_FC,cut_FDR = cut_FDR)

  g <- setdiff(label,label_ns)

  gene_list <- list()
  for (i in g) {
    gene_list[[i]] = row.names(deg_df_g[which(deg_df_g$group == i),])
  }
  rm(list = "i")

  gene_list[["diff"]] <- unique(unlist(gene_list))

  test <- mclapply(gene_list, function(x)
    suppressMessages(enhance_enrichGO(gene = x,OrgDb = OrgDb,keyType = keyType,ont = ont,simplify = simplify,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          qvalueCutoff = qvalueCutoff,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          readable = FALSE,
                          pool = FALSE)),mc.cores = getOption("mc.cores", mc.cores))

  return(test)

}


