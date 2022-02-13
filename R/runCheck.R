#' Check you data by PCA and Heatmaps
#'
#' A previous check function for overview the data sets.
#'
#' @param object a DEGContainer contains at least dataInfo
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#' @importFrom usethis ui_done ui_path ui_info
#' @importFrom glue glue
#' @importFrom edgeR cpm
#' @importFrom fs dir_exists dir_create
#'
#' @return
#' @export
#'
#' @examples
#' runCheck(object, tempdir())
runCheck <- function(object,
                      dir = ".",
                      prefix = "1-run_check",
                      palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  exprSet= expMatrix(object)
  group_list =  groupInfo(object)

  if(dataType(object) == "Counts") {

    ui_info("Your Counts will be Convert to CPM and make a check!")
    # dat=cpm(exprSet, prior.count = 2)
    dat = log2(cpm(exprSet)+1)

    pca_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("PCA checking have done, a plot was store in {ui_path(dir)}."))
    corall_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("Correlation checking have done, a plot was store in {ui_path(dir)}."))
    cor500_check(exprSet,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("Correlation to top 500 genes checking have done, a plot was store in {ui_path(dir)}."))
    top1000_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("Standard Deviation top 1000 genes checking have done, a plot was store in {ui_path(dir)}."))

  } else {

    ui_info("Your Array will be make a check!")

    pca_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("PCA checking have done, a plot was store in {ui_path(dir)}."))
    corall_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("Correlation checking have done, a plot was store in {ui_path(dir)}."))
    cor500_check(exprSet,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("Correlation to top 500 genes checking have done, a plot was store in {ui_path(dir)}."))
    top1000_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("Standard Deviation top 1000 genes checking have done, a plot was store in {ui_path(dir)}."))

  }

}

#' PCA QC check
#'
#' check data quality by PCA plot
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param list a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#'
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggplot2 theme element_text labs ggsave
#' @importFrom RColorBrewer brewer.pal
#' @importFrom glue glue
#'
#' @return a figure of PCA
#'
#' @examples
#' pca_check(data, list)
#' @noRd
pca_check <- function(data, list, dir = ".", prefix = "1-run_check", palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_all_samples_PCA_by_type.pdf')
  dat=t(data)
  dat.pca <- PCA(dat , graph = FALSE)
  p <- fviz_pca_ind(dat.pca,
                    geom.ind = "point", # show points only (nbut not "text")
                    col.ind =  list, # color by groups
                    palette = palette,
                    addEllipses = TRUE, # Concentration ellipses
                    legend.title = "Groups"
  ) + theme(plot.title = element_text(face = "bold")) +
    labs(title = "all samples - PCA")
  ggsave(filename,p, width = 400/100, height = 350/100, dpi = 300, units = "in", limitsize = FALSE)
}

#' Correlation of all samples and all genes
#'
#' Check data quality by calculating correlation of all samples
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param list a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#'
#' @importFrom glue glue
#' @importFrom pheatmap pheatmap
#'
#' @return a Heatmap shows correlation of all samples
#'
#' @examples
#' corall_check(data, list)
#' @noRd
corall_check <- function(data, list, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_cor_all.pdf')
  corExprHeatmap(expr = data,group_list = list,palette = palette,filename = filename,main = "Correlation by all genes")
}

#' Correlation of all samples and top 500 genes
#'
#' Check data quality by calculating correlation of all samples but top 500 genes
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param list a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#'
#' @importFrom glue glue
#' @importFrom pheatmap pheatmap
#' @importFrom stats cor mad
#'
#' @return a Heatmap shows correlation of all samples to 500 genes
#'
#' @examples
#' cor500_check(counts_input, list)
#' @noRd
cor500_check <- function(data, list, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_cor_top500.pdf')
  corExprHeatmap(expr = data,group_list = list,palette = palette,top=500,filename = filename,main = "Correlation by 500 genes")
}

#' Heatmap of all samples and Top1000 genes
#'
#' Check data quality by plotting heatmap of top 1000 genes
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param list a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#'
#' @importFrom glue glue
#' @importFrom pheatmap pheatmap
#' @importFrom stats sd
#' @importFrom utils tail
#'
#' @return a Heatmap shows top100 genes of all samples
#'
#' @examples
#' top1000_check(data, list)
#' @noRd
top1000_check <- function(data, list, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_heatmap_top1000_sd.pdf')
  topExprHeatmap(expr = data,group_list = list,filename = filename,palette = palette,top = 1000,main = "SD Top 1000 genes")
}
