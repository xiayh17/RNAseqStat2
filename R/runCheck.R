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

  exprSet= matrixFiltered(object)
  group_list =  groupInfo(object)

  if(dataType(object) == "Counts") {

    ui_info("Your Counts will be Convert to CPM and make a check!")
    # dat=cpm(exprSet, prior.count = 2, log = TRUE)
    dat = log2(cpm(exprSet)+1)

    box_check(dat,group_list,dir = dir,prefix = prefix,palette = palette,y = "log2(cpm(count)+1)")
    ui_done(glue("Box checking have done, a plot was store in {ui_path(dir)}."))
    suppressMessages(density_check(dat,group_list,dir = dir,prefix = prefix,palette = palette,y = "log2(cpm(count)+1)"))
    ui_done(glue("Density checking have done, a plot was store in {ui_path(dir)}."))

  } else if (dataType(object) == "Array") {

    ui_info("Your Array will be make a check!")
    dat = exprSet

    box_check(dat,group_list,dir = dir,prefix = prefix,palette = palette,y = "expression value")
    ui_done(glue("Box checking have done, a plot was store in {ui_path(dir)}."))
    suppressMessages(density_check(dat,group_list,dir = dir,prefix = prefix,palette = palette,y = "expression value"))
    ui_done(glue("Density checking have done, a plot was store in {ui_path(dir)}."))


  }

  if(species(object) == "Human"){
    HKG_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
    ui_done(glue("HKG checking have done, a plot was store in {ui_path(dir)}."))
  }

  pca_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
  ui_done(glue("PCA checking have done, a plot was store in {ui_path(dir)}."))
  corall_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
  ui_done(glue("Correlation checking have done, a plot was store in {ui_path(dir)}."))
  cor500_check(exprSet,group_list,dir = dir,prefix = prefix,palette = palette)
  ui_done(glue("Correlation to top 500 genes checking have done, a plot was store in {ui_path(dir)}."))
  top1000_check(dat,group_list,dir = dir,prefix = prefix,palette = palette)
  ui_done(glue("Standard Deviation top 1000 genes checking have done, a plot was store in {ui_path(dir)}."))

}

#' PCA QC check
#'
#' check data quality by PCA plot
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#' @importFrom glue glue
#'
#' @return a figure of PCA
#' @export
#'
#' @examples
#' pca_check(data, group)
pca_check <- function(data, group, dir = ".", prefix = "1-run_check", palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_all_samples_PCA_by_type.pdf')
  exprPCA(expr = data, group_list = group,
          palette = palette,
          filename = filename, main = "all samples - PCA",
          width=4, height = 4)
}

#' Correlation of all samples and all genes
#'
#' Check data quality by calculating correlation of all samples
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#' @importFrom glue glue
#'
#' @return a Heatmap shows correlation of all samples
#' @export
#'
#' @examples
#' corall_check(data, group)
corall_check <- function(data, group, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_cor_all.pdf')
  exprCorHeatmap(expr = data,group_list = group,palette = palette,filename = filename,main = "Correlation by all genes")
}

#' Correlation of all samples and top 500 genes
#'
#' Check data quality by calculating correlation of all samples but top 500 genes
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#' @importFrom glue glue
#'
#' @return a Heatmap shows correlation of all samples to 500 genes
#' @export
#'
#' @examples
#' cor500_check(counts_input, group)
cor500_check <- function(data, group, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_cor_top500.pdf')
  exprCorHeatmap(expr = data,group_list = group,palette = palette,top=500,filename = filename,main = "Correlation by 500 genes")
}

#' Heatmap of all samples and Top1000 genes
#'
#' Check data quality by plotting heatmap of top 1000 genes
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#' @importFrom glue glue
#' @importFrom pheatmap pheatmap
#' @importFrom stats sd
#' @importFrom utils tail
#'
#' @return a Heatmap shows top100 genes of all samples
#' @export
#'
#' @examples
#' top1000_check(data, group)
top1000_check <- function(data, group, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_heatmap_top1000_sd.pdf')
  exprTopHeatmap(expr = data,group_list = group,filename = filename,palette = palette,top = 1000,main = "SD Top 1000 genes")
}


#' Boxplot of expression matrix
#'
#' Check data by box plot
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#' @param ... more parameters in \code{\link{exprBox}}
#'
#' @return a ridges plot file
#' @export
#'
#' @examples
#' density_check(data,group)
box_check <- function(data, group, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],...) {
  filename = glue('{dir}/{prefix}_boxplot.pdf')
  exprBox(expr = data,group_list = group,filename = filename,palette = palette,main = "Boxplot of gene expression",...)
}
#' Density of expression matrix
#'
#' Check data density by ridges plot
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#' @param ... more parameters in \code{\link{exprRidges}}
#'
#' @return a ridges plot file
#' @export
#'
#' @examples
#' density_check(data,group)
density_check <- function(data, group, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],...) {
  filename = glue('{dir}/{prefix}_density.pdf')
  exprRidges(expr = data,group_list = group,filename = filename,palette = palette,main = "Density of gene expression",...)
}

#' Human housekeeping genes expression
#'
#' Check Human housekeeping genes expression
#'
#' @param data a cpm data frame of rows in genes and columns in samples
#' @param group a list ordered by samples in data
#' @param dir a directory to store results
#' @param prefix a prefix of file names in this step
#' @param palette a color palette for plots
#'
#' @return a heatmap plot file
#' @export
#'
#' @examples
#' HKG_check(data,group)
HKG_check <- function(data, group, dir = ".", prefix = "1-run_check",palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {
  filename = glue('{dir}/{prefix}_HKG.pdf')
  exprHKGheatmap(expr = data,group_list = group,filename = filename,palette = palette,main = "Human housekeeping genes expression")
}
