#' run all step in default setting
#'
#' Only need count data, group, OrgDb, dir
#'
#' @param count_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param OrgDb OrgDb
#' @param dir a directory to store results
#' @param case_group which one is test group in your group_list
#' @param control_group which one is control group in your group_list
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @importFrom usethis ui_info ui_done
#'
#' @return a dir contains all results
#' @export
#'
#' @examples
#' \dontrun{
#' runAll(count_data = counts_input, group_list = group_list, OrgDb = 'org.Hs.eg.db', dir = tempdir())
#' }
# runAll <- function(count_data, group_list, OrgDb = 'org.Hs.eg.db', dir = ".",case_group = "T", control_group = "C",parallel = FALSE) {
#
#   ui_info("Step1: Check you data")
#   run_check(counts_data = count_data, group_list = group_list, dir = dir)
#   ui_done("All Quick check of your data have done!")
#
#   ui_info("Step2: DEG analysis")
#   deg_results <- run_deg(count_data, group_list, case_group = case_group, control_group = control_group,dir = dir,parallel = parallel)
#   ui_done("All DEG analysis have done!")
#
#   ui_info("Step3: DEG Vision")
#   deg_group <- run_degVision(deg_results)
#   ui_done("All DEG Vision have done")
#
#   ui_info("Step3: Hyper analysis")
#
#
#   message(glue("Step4: EnrichKEGG analysis"))
#
#
#   message(glue("Step5: GSEA_KEGG analysis"))
#
#
# }
