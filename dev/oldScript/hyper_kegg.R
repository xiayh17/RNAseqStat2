#' hyper_go
#'
#' a funtion group deg data frame and make enrichment analysis in one step
#'
#' @inheritParams cut_much
#' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
#' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
#' @param OrgDb OrgDb
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param use_internal_data logical, use KEGG.db or latest online KEGG data
#' @param mc.cores param for mclapply, choose cores to use in non-Windows machine
#'
#' @importFrom parallel mclapply
#' @importFrom clusterProfiler enrichKEGG bitr
#' @importFrom glue glue
#'
#' @return enrichKEGG result
#' @export
#'
#' @examples
#' \dontrun{
#' hyper_kegg(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue")
#' }
hyper_kegg <- function(deg_data, x, y, cut_FC = 1, cut_FDR = 0.05,
                      organism = switch(OrgDb,
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
                      mc.cores = 1L) {

  deg_df_g <- cut_much(deg_data, x = x, y = y,cut_FC = cut_FC,cut_FDR = cut_FDR)

  g <- setdiff(label,label_ns)

  gene_list <- list()
  for (i in g) {
    SYMBOLS_id = row.names(deg_df_g[which(deg_df_g$group == i),])
    ENTREZ_id = clusterProfiler::bitr(SYMBOLS_id,
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = OrgDb)
    gene_list[[i]] = ENTREZ_id[,"ENTREZID"]
  }
  rm(list = "i")

  usethis::ui_info("Enrich KEGG analysis Start. This process will take a few minutes.")

  test <- mclapply(gene_list, function(x)
    suppressMessages(enrichKEGG(
      gene = x,
      organism = organism,
                                keyType = keyType,
                                pvalueCutoff = pvalueCutoff,
                                pAdjustMethod = pAdjustMethod,
                                minGSSize = minGSSize,
                                maxGSSize = maxGSSize,
                                qvalueCutoff = qvalueCutoff,
                                use_internal_data = use_internal_data)),mc.cores = getOption("mc.cores", mc.cores))

  return(test)

}
