#' run all workflow
#'
#' runALL include runCheck runDEG runHyper runGSEA runMSigDB
#'
#' @param object a DEGContainer
#' @param dir a directory to store results
#' @param top top for enrich result filter
#' @param parallel use parallel in DESeq2
#' @param GO logic, run GO in enrich step
#' @param KEGG logic, run KEGG in enrich step
#'
#' @return a DEGContainer
#'
#' @export
runALL <- function(object,dir = "output",top = 10,parallel = T,GO = FALSE,KEGG = TRUE) {

  runCheck(object = object,dir = dir)

  object_g <- runDEG(obj = object,dir =dir, parallel = parallel)

  object_h <- runHyper(obj = object_g,dir = dir,top = top,GO = GO,KEGG = KEGG)

  object_gs <- runGSEA(obj = object_h,dir = dir,top = top,GO = GO,KEGG= KEGG)

  object_va <- runMSigDB(obj = object_gs,dir = dir,top = top)

  return(object_va)

}
