##' Gene Set Enrichment Analysis of Gene Ontology
##'
##'
##' @title gseGO2
##' @param geneList order ranked geneList
##' @param ont one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
##' @param keyType keytype of gene
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param eps This parameter sets the boundary for calculating the p value.
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @param ... other parameter
##' @importClassesFrom DOSE gseaResult
##' @export
##' @return gseaResult object
gseGO2 <- function(geneList,
                  ont           = "ALL",
                  TERM2GENE,
                  TERM2NAME = NA,
                  organism = "UNKNOW",
                  keyType = "SYMBOL",
                  exponent      = 1,
                  minGSSize     = 10,
                  maxGSSize     = 500,
                  eps           = 1e-10,
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  verbose       = TRUE,
                  seed          = FALSE,
                  by            = 'fgsea',
                  ...) {

  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

  GO_DATA <- get_GO_data2(TERM2GENE, TERM2NAME, keyType, ont, organism)

  res <-  DOSE:::GSEA_internal(geneList      = geneList,
                        exponent      = exponent,
                        minGSSize     = minGSSize,
                        maxGSSize     = maxGSSize,
                        eps           = eps,
                        pvalueCutoff  = pvalueCutoff,
                        pAdjustMethod = pAdjustMethod,
                        verbose       = verbose,
                        USER_DATA     = GO_DATA,
                        seed          = seed,
                        by            = by,
                        ...)



  if (is.null(res))
    return(res)

  res@organism <- organism
  res@setType <- ont
  res@keytype <- keyType

  if (ont == "ALL") {
    res <- clusterProfiler:::add_GO_Ontology(res, GO_DATA)
  }
  return(res)
}

#' GO Enrichment Analysis of a gene set.
#'
#' Given a vector of genes, this function will return the enrichment GO
#' categories after FDR control. \code{\link[clusterProfiler]{enrichGO}} only support OrgDb and
#' \code{\link[clusterProfiler]{enricher}} is too simple for GO.
#' In this case, a modified version of enrichGO here
#'
#' @param gene a vector of entrez gene id.
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
#' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column with term and name
#' @param keyType keytype of input gene
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param pvalueCutoff adjusted pvalue cutoff on enrichment tests to report
#' @param pAdjustMethod  one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant.  Tests must pass i) \code{pvalueCutoff} on unadjusted pvalues, ii) \code{pvalueCutoff} on adjusted pvalues and iii) \code{qvalueCutoff} on qvalues to be reported.
#' @param minGSSize minimal size of genes annotated by Ontology term for testing.
#' @param maxGSSize maximal size of genes annotated for testing
#' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
#' @return An \code{enrichResult} instance.
#' @importClassesFrom DOSE enrichResult
#' @importFrom DOSE setReadable
#' @keywords manip
#' @export
#' @examples
#' \dontrun{
#'  enrichGO2(gene,TERM2GENE,TERM2NAME)
#' }
enrichGO2 <- function(gene,
                      TERM2GENE,
                      TERM2NAME = NA,
                      organism = "UNKNOW",
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pvalueCutoff=0.05,
                      pAdjustMethod="BH",
                      universe,
                      qvalueCutoff = 0.2,
                      minGSSize = 10,
                      maxGSSize = 500,
                      pool=FALSE) {

  ont %<>% toupper
  ont <- match.arg(ont, c("BP", "MF", "CC", "ALL"))

  GO_DATA <- get_GO_data2(TERM2GENE, TERM2NAME, keyType, ont, organism)

  # GO_DATA <- get_GO_data(OrgDb, ont, keyType)

  if (missing(universe))
    universe <- NULL

  if (ont == "ALL" && !pool) {
    lres <- lapply(c("BP", "CC", "MF"), function(ont)
      suppressMessages(enrichGO2(gene, TERM2GENE,TERM2NAME,organism, keyType, ont,
                                pvalueCutoff, pAdjustMethod, universe,
                                qvalueCutoff, minGSSize, maxGSSize
      ))
    )

    lres <- lres[!vapply(lres, is.null, logical(1))]
    if (length(lres) == 0)
      return(NULL)

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
    res <- DOSE:::enricher_internal(gene,
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
  res@organism <- organism
  # if(readable) {
  #   res <- clusterProfiler::setReadable(res, OrgDb)
  # }
  res@ontology <- ont

  if (ont == "ALL") {
    res <- clusterProfiler:::add_GO_Ontology(res, GO_DATA)
  }
  return(res)

}

get_GO_data2 <- function(TERM2GENE, TERM2NAME = NA, keyType, ont, organism = "UNKNOW") {

  GO_Env <- clusterProfiler:::get_GO_Env()
  use_cached <- FALSE
  keytype = keyType

  ont2 <- NULL
  if (exists("ont", envir = GO_Env, inherits = FALSE))
    ont2 <- get("ont", envir = GO_Env)

  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE) &&
      !is.null(ont2)) {

    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)

    if (org == organism &&
        keytype == kt &&
        (ont == ont2 || ont2 == "ALL") &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){

      use_cached <- TRUE
    }
  }

  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
    if (!is.null(ont2) && ont2 != ont) { ## ont2 == "ALL"
      goAnno <- goAnno[goAnno$ONTOLOGYALL == ont,]
    }
  } else {

    goOnt.df <- clusterProfiler::go2ont(as.character(TERM2GENE[,1]))
    goterms <- goOnt.df[,2]
    names(goterms) <- goOnt.df[,1]

    if (ont != "ALL") {
      goterms <- goterms[goterms == ont]
    }

    goAnno <- TERM2GENE[,c(2,1)]

    colnames(goAnno) <- c(keytype, "GOALL")
    goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
    goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]

    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("ont", ont, envir = GO_Env)
    assign("organism", organism, envir=GO_Env)
  }

  GO2GENE <- unique(goAnno[, c(2,1)])

  GO_DATA <- DOSE:::build_Anno(GO2GENE, TERM2NAME)

  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique

  if (!is.null(ont2) && ont2 == "ALL") {
    return(GO_DATA)
  }

  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)

  return(GO_DATA)
}

