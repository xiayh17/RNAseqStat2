#' run DEG minimal
#'
#' run DEG module in limma, DESeq2 and edgeR
#'
#' @param object a DEGContainer contains at least dataInfo
#' @param dir a directory to store results
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
#' degResolve(counts_input,group_list,case_group = "T", control_group = "C",dir = tempdir())
#' }
degResolveArray <- function(object,dir = ".", prefix = "2-DEG") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  if (is.null(matrixFiltered(object))) {
    expr_data = expMatrix(object)
  } else {
    expr_data = matrixFiltered(object)
  }

  group_list = groupInfo(object)
  case_group = caseGroup(object)
  control_group = setdiff(group_list,case_group)

  # usethis::ui_info(glue("Start DESeq2 analysis."))
  # deg_df_DESeq2 <- DESeq2_resolve(expr_data,group_list,parallel = parallel,dir=dir,prefix = prefix,
  #                                 case_group = case_group, control_group = control_group, qc = qc)

  # usethis::ui_info(glue("Start edgeR analysis."))
  # deg_df_edgeR <- edgeR_resolveArray(expr_data, group_list,
  #                               control_group = control_group)

  usethis::ui_info(glue("Start limma analysis."))
  deg_df_limma <- limma_resolveArray(expr_data,group_list,
                                case_group = case_group, control_group = control_group)

  # usethis::ui_info(glue("Merge data above."))
  # allg <- Reduce(intersect, list(
  #                                rownames(deg_df_edgeR),
  #                                rownames(deg_df_limma)))
  #
  # deg_df_intersect=cbind(deg_df_limma[allg,c("logFC","P.Value","adj.P.Val")],
  #                        deg_df_edgeR[allg,c("logFC","PValue","FDR")])
  #
  # colnames(deg_df_intersect) <- paste0(colnames(deg_df_intersect),"_",rep(c("limma","edgeR"),each=2))

  # csv_name = glue('{dir}/{prefix}_intersect_results.csv')
  # write.csv(deg_df_intersect,file = csv_name)
  # ui_done(glue("DEG results in csv format is stored in {usethis::ui_path(csv_name)}"))

  limma_res(object) = deg_df_limma
  # edgeR_res(object) = deg_df_edgeR
  # DESeq2_res(object) = deg_df_DESeq2
  # merge_res(object) = deg_df_intersect
  ui_done(glue("DEG results is updated in {ui_code('DEGContainer')}"))
  # rdata_name = glue('{dir}/{prefix}_results.Rdata')
  # save(deg_results,file = rdata_name)
  # ui_done(glue("DEG results in Rdata is stored in {usethis::ui_path(rdata_name)}"))



  return(object)

}

#' Basic produce of limma
#'
#' A basic function to get data and produce results of limma
#'
#' @param expr_data a counts data frame of rows in genes and columns in samples
#' @param group_list a character vector ordered by samples in counts_data
#' @param case_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#'
#' @importFrom stats model.matrix na.omit
#' @importFrom edgeR DGEList cpm calcNormFactors
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#'
#' @return a DEG data frame
#' @export
#'
#' @examples
#' limma_resolve(counts_input, group_list, control_group= "C", case_group = "T")
limma_resolveArray <- function(expr_data, group_list, control_group, case_group) {

  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(expr_data)

  fit <- lmFit(expr_data, design)

  con=paste0(case_group,'-',control_group)

  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)

  tempOutput = topTable(fit2, coef=con, number=Inf)
  DEG_limma_voom = na.omit(tempOutput)

  return(DEG_limma_voom)

}
