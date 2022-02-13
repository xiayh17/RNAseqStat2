#' gseGO from DEG data
#'
#' run \code{\link[clusterProfiler]{gseGO}}. of DEG data
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param keyType keytype of gene
#' @param OrgDb OrgDb
#' @param exponent weight of each step
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of genes annotated for testing
#' @param eps This parameter sets the boundary for calculating the p value.
#' @param pvalueCutoff pvalue Cutoff
#' @param pAdjustMethod pvalue adjustment method
#' @param verbose print message or not
#' @param seed logical
#' @param by one of 'fgsea' or 'DOSE'
#' @param ... other parameter
#'
#' @importFrom clusterProfiler gseGO
#'
#' @return a gsego result
#' @export
#'
#' @examples
#' \dontrun{
#' gsea_go(DEG_df,x = "log2FoldChange")
#' }
gsea_go <- function(deg_data,x,
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
                      by = "fgsea",...) {

  # gene <- bitr(rownames(deg_data), fromType = "SYMBOL",
  #              toType =  "ENTREZID",
  #              OrgDb = OrgDb)
  # gene$logfc <- deg_data[,x][match(gene$SYMBOL,rownames(deg_data))]

  geneList = deg_data[,x]
  names(geneList) = rownames(deg_data)
  geneList=sort(geneList,decreasing = TRUE)

  kk_gse <- gseGO(
    geneList = geneList,
    ont = ont,
    OrgDb = OrgDb,
    keyType = keyType,
    exponent = exponent,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    eps = eps,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    verbose = verbose,
    seed = seed,
    by = by,
    ...
  )

  # kk=DOSE::setReadable(kk_gse, OrgDb=OrgDb,keyType='ENTREZID')

  return(kk_gse)
}

