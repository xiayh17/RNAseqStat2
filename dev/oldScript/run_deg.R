#' run DEG
#'
#' run DEG module in limma, DESeq2 and edgeR
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param case_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @importFrom fs dir_exists dir_create
#' @importFrom glue glue
#' @importFrom utils write.csv
#' @importFrom usethis ui_done
#'
#' @return a csv file and a DEG_container ob
#' @export
#'
#' @examples
#' \dontrun{
#' run_deg(counts_input,group_list,case_group = "T", control_group = "C",dir = tempdir())
#' }
run_deg <- function(counts_data,group_list,case_group = "T", control_group = "C",dir,prefix = "2-DEG",parallel = FALSE) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  ## 检查数据 ----
  ## group
  if (all(c(case_group, control_group) %in% group_list) & ncol(counts_data) == length(group_list)) {
    usethis::ui_info(glue("{case_group}_VS_{control_group} info seems ok"))
  } else {
    usethis::ui_stop(glue("Please check group_list = {group_list},case_group = {case_group}, control_group = {control_group}"))
  }

  ## counts
  if (any(class(counts_data) == "data.frame") & all(apply(counts_data, 2, is.integer))) {
    usethis::ui_info(glue("Counts data frame seems ok"))
  } else {
    usethis::ui_stop(glue("Please check your data frame! Is it an integer data frame?"))
  }
  ## names
  usethis::ui_info("Please make sure your data frame is rownamed by Gene Symbol")
  ## ----

  usethis::ui_info(glue("Start DESeq2 analysis."))
  deg_df_DESeq2 <- run_DESeq2(counts_data,group_list,parallel = parallel,
            case_group = case_group, control_group = control_group, qc = TRUE)

  usethis::ui_info(glue("Start edgeR analysis."))
  deg_df_edgeR <- run_edgeR(counts_data, group_list,
            control_group = control_group)

  usethis::ui_info(glue("Start limma analysis."))
  deg_df_limma <- run_limma(counts_data,group_list,
            case_group = case_group, control_group = control_group)

  usethis::ui_info(glue("Merge data above."))
  allg <- Reduce(intersect, list(rownames(deg_df_DESeq2),
                         rownames(deg_df_edgeR),
                         rownames(deg_df_limma)))

  deg_df_intersect=cbind(deg_df_limma[allg,c(1,4)],
              deg_df_edgeR[allg,c(1,5)],
              deg_df_DESeq2[allg,c(2,6)])

  colnames(deg_df_intersect) <- paste0(colnames(deg_df_intersect),"_",rep(c("limma","edgeR","DESeq2"),each=2))

  csv_name = glue('{dir}/{prefix}_intersect_results.csv')
  write.csv(deg_df_intersect,file = csv_name)
  ui_done(glue("DEG results in csv format is stored in {usethis::ui_path(csv_name)}"))

  deg_results <- create_DEG_container(
    deg_df_limma = deg_df_limma,
    deg_df_edgeR = deg_df_edgeR,
    deg_df_DESeq2 = deg_df_DESeq2,
    deg_df_intersect = deg_df_intersect,
    # counts_data_filtered = counts_data_filtered,
    group_list = group_list
  )

  rdata_name = glue('{dir}/{prefix}_results.Rdata')
  save(deg_results,file = rdata_name)
  ui_done(glue("DEG results in Rdata is stored in {usethis::ui_path(rdata_name)}"))

  return(deg_results)

}

#' a S4 class
#'
#' contains results of DEG module
#'
#' @importFrom methods setClass
#'
#' @rdname DEG_container
#'
#' @export
setClass(Class="DEG_container",
         slots = representation(
           deg_df_limma = "data.frame",
           deg_df_edgeR = "data.frame",
           deg_df_DESeq2 = "data.frame",
           deg_df_intersect = "data.frame",
           # counts_data_filtered = "data.frame",
           group_list = "character"
         )
)


#' DEG_container create function
#'
#' a function to create DEG_container ob
#'
#' @param deg_df_limma a deg data frame was produced by limma
#' @param deg_df_edgeR a deg data frame was produced by edgeR
#' @param deg_df_DESeq2 a deg data frame was produced by DESeq2
#' @param deg_df_intersect a deg data frame was combined by results of limma edgeR DESeq2
#' @param group_list a list ordered by samples in counts_data
#'
#' @importFrom methods new
#'
#' @return a DEG_container ob
#' @noRd
create_DEG_container <- function(
  deg_df_limma = NA,
  deg_df_edgeR = NA,
  deg_df_DESeq2 = NA,
  deg_df_intersect = NA,
  # counts_data_filtered = NA,
  group_list = NA
) {

  new("DEG_container",
      deg_df_limma = as.data.frame(deg_df_limma),
      deg_df_edgeR = as.data.frame(deg_df_edgeR),
      deg_df_DESeq2 = as.data.frame(deg_df_DESeq2),
      deg_df_intersect = as.data.frame(deg_df_intersect),
      # counts_data_filtered = as.data.frame(counts_data_filtered),
      group_list = as.character(group_list)
  )
}
