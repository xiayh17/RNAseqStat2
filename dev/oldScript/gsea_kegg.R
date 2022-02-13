#' gseKEGG from DEG data
#'
#' run \code{\link[clusterProfiler]{gseKEGG}}. of DEG data
#'
#' @param deg_data a DEG data frame contains logFC and p value
#' @param x which column is log FC
#' @param organism supported organism listed in 'http://www.genome.jp/kegg/catalog/org_list.html'
#' @param keyType one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
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
#' @param use_internal_data logical, use KEGG.db or latest online KEGG data
#' @param ... other parameter
#'
#' @importFrom clusterProfiler bitr gseKEGG
#' @importFrom DOSE setReadable
#'
#' @return a gsekegg result
#' @export
#'
#' @examples
#' \dontrun{
#' gsea_kegg(DEG_df,x = "log2FoldChange")
#' }
gsea_kegg <- function(deg_data,x,
                           organism = switch (OrgDb,
                                              'org.Hs.eg.db' = 'hsa',
                                              'org.Mm,eg.db' = 'mmu'
                           ),
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
                           by = "fgsea",...) {



  gene <- bitr(rownames(deg_data), fromType = "SYMBOL",
               toType =  "ENTREZID",
               OrgDb = OrgDb)
  gene$logfc <- deg_data[,x][match(gene$SYMBOL,rownames(deg_data))]

  geneList=gene$logfc
  names(geneList)=gene$ENTREZID
  geneList=sort(geneList,decreasing = TRUE)

  kk_gse <- gseKEGG(geneList = geneList,
    organism = organism,
                    keyType = keyType,
                    exponent = exponent,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    eps = eps,
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod = pAdjustMethod,
                    verbose = verbose,
                    use_internal_data = use_internal_data,
                    seed = seed,
                    by = by,...)

  kk=DOSE::setReadable(kk_gse, OrgDb=OrgDb,keyType='ENTREZID')

  return(kk)
}

