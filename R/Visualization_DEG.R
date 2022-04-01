#' Heatmap for DEG data frame
#'
#' default will return a top 100 deg heatmap in p value = 0.05
#'
#' @param object a counts data frame of rows in genes and columns in samples
#' @param top a single number or a length of 2 numeric vector, if 2 numeric vector, first one is top max logFC.
#' @param which which model of deg analysis. kinds of DEG; can be "limma", "edgeR", "DESeq2" or "MSigDB"
#' @param palette a color palette for plots
#' @param filename NA or a file path
#' @param show_gene logical,show gene name
#' @param category MSigDB collection abbreviation, such as H or C1.
#' @param ... More \code{\link[pheatmap]{pheatmap}} parameters.
#'
#' @importFrom edgeR cpm
#' @importFrom pheatmap pheatmap
#' @importFrom glue glue
#'
#' @return a heatmap plot file
#' @export
#'
#' @examples
#' DEGtopHeatmap(object,which = "limma")
DEGtopHeatmap <- function(object, which, top = 50, filename = NA, show_gene = TRUE,category = "H",
                            palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],...) {

  if (is.null(matrixFiltered(object))) {
    counts_data = expMatrix(object)
  } else {
    counts_data = matrixFiltered(object)
  }

  deg_data <- dataDEG(obj = object,which = which,category = category)

  group_list = groupInfo(object)
  x = FC_Identify(res = deg_data)
  y = pvalue_Identify(res = deg_data)
  top = top

  choose_gene <- topGene(object = object, topSig = top, which = which)

  filename = filename

  # exprSet=log2(edgeR::cpm(counts_data)+1)
  exprSet=cpm(counts_data, prior.count = 2, log = TRUE)
  choose_matrix=exprSet[choose_gene,]
  choose_matrix=t(scale(t(choose_matrix)))
  choose_matrix[choose_matrix>2]=2
  choose_matrix[choose_matrix< -2]= -2
  colD=data.frame(Groups=group_list)
  rownames(colD)=colnames(exprSet)
  names(palette) <- unique(group_list)
  pheatmap(choose_matrix,...,
           annotation_col = colD,fontsize = 12,
           width = (ncol(choose_matrix)*0.3+2.2) *2,
           height = 550/100*3*nrow(choose_matrix)/100,
           annotation_colors = list(Groups = palette),show_rownames = show_gene,border_color = NA,
           filename = filename)
}

#' Plot venn for a set of gene
#'
#' @param geneSets a list of gene sets
#' @param palette palette
#'
#' @importFrom VennDiagram venn.diagram
#' @importFrom ggplotify as.ggplot
#' @importFrom cowplot as_grob
#'
#' @return a plot of venn
#' @export
#'
#' @examples
#' DEGvenn(geneSets)
DEGvenn <- function(geneSets, palette = c('#1f78b4','#33a02c','#ff7f00')) {

  p =  venn.diagram(x = geneSets,
                                 filename= NULL,
                                 disable.logging = FALSE,
                                 col=NA,
                                 fill=palette,
                                 cat.col=palette,
                                 cat.cex = 1,
                                 cat.dist = -0.15,
                                 rotation.degree = 0,
                                 main.cex = 1,
                                 cex=1,
                                 alpha = 0.5,
                                 reverse=TRUE)
  file.remove(dir(pattern = ("^VennDiagram.*log$")))
  p = as.ggplot(as_grob(p))

  return(p)

}
