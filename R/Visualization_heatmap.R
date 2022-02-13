#' Correlation Heatmap between samples
#'
#' plot a heatmap of correlation between samples
#'
#' @param expr a expression matrix
#' @param group_list a character vector order by samples
#' @param filename file name you want to save plot
#' @param top if not NA, only keep top of genes by mad
#' @param main plot title
#' @param palette annotation palette
#' @param anno_title annotation title
#' @param annotation_col annotation dataframe
#' @param width width of plot
#' @param height height of plot
#' @param heatmapParam more parameters in \code{\link[pheatmap]{pheatmap}}
#' @param ... more parameters in \code{\link[pheatmap]{pheatmap}}
#' @importFrom utils modifyList
#' @importFrom stats cor mad
#'
#' @details
#' Use \code{\link[stats]{cor}} to get correlation of samples. Can filter genes by `top`.
#' `top` means the top mad (by \code{\link[stats]{mad}}) of expression matrix.
#'
#' @returns a heatmap plot or a file
#' @export
#'
#' @examples
#' corExprHeatmap(counts_input,group_list=group_list)
corExprHeatmap <- function(expr,group_list,filename = NA,top = NA,main = "Correlation by all genes",
                        palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],anno_title = "Group",
                        annotation_col = ac_(expr,group_list),
                        width = ncol(expr)*0.3+2.2,height = ncol(expr)*0.3+2.2,
                        ...,heatmapParam = list(show_rownames = F)) {


  names(palette) <- unique(group_list)
  palette = split(palette,anno_title)

  if (is.na(top)) {
    m=cor(expr)
  } else if (is.numeric(top)){
    expr2=expr[names(sort(apply(expr, 1,mad),decreasing = T)[1:top]),]
    m=cor(expr2)
  }


  width = width
  height = height


  param = list()
  heatmapParam = modifyList(param,heatmapParam)
  heatmap <- do.call("pheatmap",modifyList(
    list(mat = m,width = width, height = height, main = main,
         annotation_col = annotation_col, annotation_colors = palette,filename = filename),
    heatmapParam))

  return(heatmap)

}

#' Expression Heatmap
#'
#' plot a heatmap of expression matrix
#'
#' @param expr a expression matrix
#' @param group_list a character vector order by samples
#' @param filename file name you want to save plot
#' @param top if not NA, only keep top of genes by mad
#' @param main plot title
#' @param palette annotation palette
#' @param anno_title annotation title
#' @param annotation_col annotation dataframe
#' @param width width of plot
#' @param height height of plot
#' @param heatmapParam more parameters in \code{\link[pheatmap]{pheatmap}}
#' @param ... more parameters in \code{\link[pheatmap]{pheatmap}}
#'
#' @importFrom stats cor mad
#'
#' @return a heatmap or a file
#' @export
#'
#' @details
#' Can filter genes by `top`.
#' `top` means the top sd (by \code{\link[stats]{sd}}) of expression matrix.
#'
#' @examples
#' topExprHeatmap(counts_input,group_list=group_list)
topExprHeatmap <- function(expr,group_list,filename = NA,top = 1000,main = "SD Top 1000 genes",
                       palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],anno_title = "Group",
                       annotation_col = ac_(expr,group_list),
                       width = ncol(expr)*0.3+2.2,height = ncol(expr)*0.3+2.2,
                       ...,heatmapParam = list(show_rownames = F)) {

  names(palette) <- unique(group_list)
  palette = split(palette,anno_title)

  if (is.na(top)) {

    m=t(scale(t(expr)))
    m[m>2]=2
    m[m< -2]= -2

  } else if (is.numeric(top)){

    expr2=expr[names(sort(apply(expr, 1,sd),decreasing = T)[1:top]),]
    m=t(scale(t(expr2)))
    m[m>2]=2
    m[m< -2]= -2

  }

  param = list()
  heatmapParam = modifyList(param,heatmapParam)
  heatmap <- do.call("pheatmap",modifyList(
    list(mat = m,width = width, height = height, main = main,
         annotation_col = annotation_col, annotation_colors = palette,filename = filename),
    heatmapParam))

  return(heatmap)

}

## 构建分组信息
ac_ <- function(expr,group_list) {
  colD=data.frame(Groups=group_list)
  rownames(colD)=colnames(expr)
  return(colD)
}
