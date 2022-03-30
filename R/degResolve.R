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
degResolve <- function(object,dir = ".", prefix = "2-DEG",parallel = FALSE,qc=TRUE) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  if (is.null(matrixFiltered(object))) {
    counts_data = expMatrix(object)
  } else {
    counts_data = matrixFiltered(object)
  }

  group_list = groupInfo(object)
  case_group = caseGroup(object)
  control_group = setdiff(group_list,case_group)

  usethis::ui_info(glue("Start DESeq2 analysis."))
  deg_df_DESeq2 <- DESeq2_resolve(counts_data,group_list,parallel = parallel,dir=dir,prefix = prefix,
                              case_group = case_group, control_group = control_group, qc = qc)

  usethis::ui_info(glue("Start edgeR analysis."))
  deg_df_edgeR <- edgeR_resolve(counts_data, group_list,
                            control_group = control_group)

  usethis::ui_info(glue("Start limma analysis."))
  deg_df_limma <- limma_resolve(counts_data,group_list,
                            case_group = case_group, control_group = control_group)

  usethis::ui_info(glue("Merge data above."))
  allg <- Reduce(intersect, list(rownames(deg_df_DESeq2),
                                 rownames(deg_df_edgeR),
                                 rownames(deg_df_limma)))

  deg_df_intersect=cbind(deg_df_limma[allg,c("logFC","P.Value","adj.P.Val")],
                         deg_df_edgeR[allg,c("logFC","PValue","FDR")],
                         deg_df_DESeq2[allg,c("log2FoldChange","pvalue","padj")])

  colnames(deg_df_intersect) <- paste0(colnames(deg_df_intersect),"_",rep(c("limma","edgeR","DESeq2"),each=3))

  # csv_name = glue('{dir}/{prefix}_intersect_results.csv')
  # write.csv(deg_df_intersect,file = csv_name)
  # ui_done(glue("DEG results in csv format is stored in {usethis::ui_path(csv_name)}"))

  limma_res(object) = deg_df_limma
  edgeR_res(object) = deg_df_edgeR
  DESeq2_res(object) = deg_df_DESeq2
  merge_res(object) = deg_df_intersect
  ui_done(glue("DEG results is updated in {ui_code('DEGContainer')}"))
  # rdata_name = glue('{dir}/{prefix}_results.Rdata')
  # save(deg_results,file = rdata_name)
  # ui_done(glue("DEG results in Rdata is stored in {usethis::ui_path(rdata_name)}"))



  return(object)

}

#' Basic produce of edgeR
#'
#' A basic function to get data and produce results of edgeR
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a character vector ordered by samples in counts_data
#' @param control_group the name of the denominator level for the fold change (Control group)
#'
#' @importFrom stats relevel model.matrix
#' @importFrom edgeR DGEList cpm calcNormFactors estimateGLMCommonDisp estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit glmLRT topTags
#'
#' @return a DEG data frame
#' @export
#'
#' @examples
#' edgeR_resolve_(counts_input, group_list, control_group= "C")
edgeR_resolve <- function(counts_data, group_list, control_group) {
  g=factor(group_list)
  g=relevel(g,control_group)

  d <- DGEList(counts=counts_data,group=g)
  keep <- rowSums(cpm(d)>1) >= 2

  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)

  dge=d
  design <- model.matrix(~0+g)
  # rownames(design)<-colnames(dge)
  # colnames(design)<-levels(g)

  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)

  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit, contrast=c(-1,1))
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)

  return(nrDEG)

}

#' Basic produce of limma
#'
#' A basic function to get data and produce results of limma
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
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
limma_resolve <- function(counts_data, group_list, control_group, case_group) {

  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(counts_data)

  dge <- DGEList(counts=counts_data)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)

  v <- voom(dge,design,plot=TRUE, normalize.method="quantile")
  fit <- lmFit(v, design)

  con=paste0(case_group,'-',control_group)

  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)

  tempOutput = topTable(fit2, coef=con, number=Inf)
  DEG_limma_voom = na.omit(tempOutput)

  return(DEG_limma_voom)

}

#' Basic produce of DESeq2
#'
#' A basic function to get data and produce results of DESeq2
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a character vector ordered by samples in counts_data
#' @param case_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param qc qc plots
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom stats na.omit
#'
#' @return a DEG data frame
#' @export
#'
#' @examples
#' \dontrun{
#' DESeq2_resolve(counts_input, group_list,case_group = "T", control_group = "C", dir = tempdir())
#' }
DESeq2_resolve <- function(counts_data,group_list,case_group,control_group,qc = TRUE,dir = ".",prefix = "2-DEG_DEseq2",parallel = FALSE) {

  colData <- data.frame(row.names=colnames(counts_data),
                        group_list=group_list)

  colData$group_list <- factor(colData$group_list)

  dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)

  if (qc) {
    DESeq2_qc(counts_data,dds,dir = dir,prefix = prefix)
  }

  res <- results(dds,parallel = parallel,
                 contrast=c("group_list",case_group,control_group))
  resOrdered <- res[order(res$padj),]

  DEG =as.data.frame(resOrdered)
  DEG = na.omit(DEG)
  return(DEG)
}

#' QC for DESeq2
#'
#' plot dispersions and RAWvsNORM
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param dds a DESeqDataSet class data set
#'
#' @importFrom glue glue
#' @importFrom usethis ui_done
#' @importFrom grDevices pdf dev.off rainbow
#' @importFrom graphics par boxplot hist
#' @importFrom DESeq2 plotDispEsts rlogTransformation varianceStabilizingTransformation
#' @importMethodsFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @return dispersions and RAWvsNORM firures
#'
#' @noRd
#' @examples
#' \dontrun{
#' DESeq2_qc(counts_input,dds, dir = tempdir())
#' }
DESeq2_qc <- function(counts_data, dds, dir = ".", prefix = "2-DEG") {

  pdf(glue("{dir}/{prefix}_DESeq2qc_dispersions.pdf"), 18, 18, pointsize=35)
  plotDispEsts(dds, main="Dispersion plot",cex = 0.45)
  dev.off()

  if (length(rownames(dds@colData)) >=50 ) {
    rld <- varianceStabilizingTransformation(dds)
  } else {

    rld <- rlogTransformation(dds)

  }

  exprMatrix_rlog=assay(rld)

  pdf(glue("{dir}/{prefix}_DESeq2qc_RAWvsNORM.pdf"),height = 8,width = 8)
  par(cex = 0.7)
  n.sample=ncol(counts_data)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(counts_data, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(counts_data))
  hist(exprMatrix_rlog)
  dev.off()
  ui_done("QC for DESeq2 is done")
}
