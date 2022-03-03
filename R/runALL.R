#' @export
runALL <- function(object,dir = "output",top = 10,parallel = T,GO = FALSE) {

  runCheck(object = object,dir = dir)

  object_g <- runDEG(obj = object,dir =dir, parallel = parallel)

  object_h <- runHyper(obj = object_g,dir = dir,top = top,GO = GO)

  object_gs <- runGSEA(obj = object_h,dir = dir,top = top)

  object_va <- runMSigDB(obj = object_gs,dir = dir,top = top)

  return(object_va)

}
