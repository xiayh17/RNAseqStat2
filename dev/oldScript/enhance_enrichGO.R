##' GO Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment GO
##' categories after FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param OrgDb OrgDb
##' @param keyType keytype of input gene
##' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
##' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
##' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
##' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant.  Tests must pass i) \code{pvalueCutoff} on unadjusted pvalues, ii) \code{pvalueCutoff} on adjusted pvalues and iii) \code{qvalueCutoff} on qvalues to be reported.
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param readable whether mapping gene ID to gene Name
##' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
##' @param simplify whether simplify
##' @return An \code{enrichResult} instance.
##' @importClassesFrom DOSE enrichResult
##' @importFrom DOSE setReadable
##' @importFrom furrr future_map furrr_options
##' @importFrom clusterProfiler simplify
##' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @examples
##' \dontrun{
##'   data(geneList, package = "DOSE")
##' 	de <- names(geneList)[1:100]
##' 	yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
##' 	head(yy)
##' }
enhance_enrichGO <- function(gene,
                     OrgDb,
                     keyType = "ENTREZID",
                     ont="MF",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500,
                     readable=FALSE, pool=FALSE, simplify = FALSE) {

  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))
  GO_DATA <- enhance_get_GO_data(OrgDb, ont, keyType)

  if (missing(universe))
    universe <- NULL

  if (ont == "ALL" && !pool) {
    lres <- lapply(c("BP", "CC", "MF"), function(ont)
      suppressMessages(enhance_enrichGO(gene, OrgDb, keyType, ont,
                                pvalueCutoff, pAdjustMethod, universe,
                                qvalueCutoff, minGSSize, maxGSSize
      ))
    )

    lres <- lres[!vapply(lres, is.null, logical(1))]
    if (length(lres) == 0)
      return(NULL)

    if (simplify) {
      lres <- lapply(lres, function(x) clusterProfiler::simplify(x))
    }

    df <- do.call('rbind', lapply(lres, as.data.frame))
    geneSets <- lres[[1]]@geneSets
    if (length(lres) > 1) {
      for (i in 2:length(lres)) {
        geneSets <- append(geneSets, lres[[i]]@geneSets)
      }
    }
    res <- lres[[1]]
    res@result <- df
    res@geneSets <- geneSets
  } else {
    res <- enhance_enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             USER_DATA = GO_DATA
    )

    if (is.null(res))
      return(res)
  }

  res@keytype <- keyType
  res@organism <- get_organism(OrgDb)
  if(readable) {
    res <- setReadable(res, OrgDb)
  }
  res@ontology <- ont


  if (ont != "ALL" & simplify) {
    res <- clusterProfiler::simplify(res)
  }

  if (ont == "ALL") {
    res <- add_GO_Ontology(res, GO_DATA)
  }

  return(res)
}

enricher <- clusterProfiler::enricher

add_GO_Ontology <- clusterProfiler:::add_GO_Ontology
