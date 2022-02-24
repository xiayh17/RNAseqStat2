# PCA plot ----------------------------------------------------------------
#' PCA plot for Grouped Samples
#'
#' plot a pca for expression matrix
#'
#' @param expr a expression matrix
#' @param group_list a character vector order by samples
#' @param palette group palette
#' @param filename file name you want to save plot
#' @param main plot title
#' @param width width of plot
#' @param height height of plot
#' @param ... more parameters in \code{\link[factoextra]{fviz_pca_ind}}
#'
#' @importFrom factoextra fviz_pca_ind
#' @importFrom FactoMineR PCA
#' @import ggplot2
#'
#' @details
#' Use \code{\link[FactoMineR]{PCA}} to get PCA of grouped samples.
#'
#' @returns a PCA plot or a file
#' @export
#'
#' @examples
#' exprPCA(counts_input,group_list=group_list)
exprPCA <- function(expr, group_list,
                    palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],
                    filename = NA,main = "all samples - PCA",
                    width=3.5, height = 4,...) {
  dat=t(expr)
  dat.pca <- PCA(dat , graph = FALSE)
  p <- fviz_pca_ind(dat.pca,
                    geom.ind = "point", # show points only (nbut not "text")
                    col.ind =  group_list, # color by groups
                    palette = palette,
                    addEllipses = TRUE, # Concentration ellipses
                    legend.title = "Groups",...
  ) + theme(plot.title = element_text(face = "bold")) +
    labs(title = main)
  if(is.na(filename)){
    return(p)
  } else (
    ggsave(filename,p, width = width, height = height, dpi = 300, units = "in", limitsize = FALSE)
  )
}
# PCA plot ----------------------------------------------------------------

# heatmap plot ------------------------------------------------------------
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
#'
#' @importFrom pheatmap pheatmap
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
exprCorHeatmap <- function(expr,group_list,filename = NA,top = NA,main = "Correlation by all genes",
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
#' @importFrom utils modifyList
#' @importFrom pheatmap pheatmap
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
exprTopHeatmap <- function(expr,group_list,filename = NA,top = 1000,main = "SD Top 1000 genes",
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
# heatmap plot ------------------------------------------------------------

# boxplot plot ------------------------------------------------------------
#' Boxplot of all samples
#'
#' plot a boxplot of all samples
#'
#' @param expr a expression matrix
#' @param group_list a character vector order by samples
#' @param palette group palette
#' @param filename file name you want to save plot
#' @param main plot title
#' @param width width of plot
#' @param height height of plot
#' @param x x axis title
#' @param y y axis title
#'
#' @return
#' @export
#'
#' @examples
#' exprBox(counts_input,group_list=group_list)
exprBox <- function(expr,group_list,palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],
                    main = "Boxplot of all samples",
                    filename = NA,x= "Samples",y = "log2(cpm(count)+1)",
                    height = 4.46,width = 8.3) {
  # format data
  boxplot_data <- exprLong(expr = expr,group_list = group_list,palette = palette)

  n = length(unique(boxplot_data$sample))

  # plot
  if (n <= 9) {
    myboxplot <- plot_boxplot_base(boxplot_data) + theme_normal() +
      coord_cartesian(clip = "off") +
      labs(x = x, y = y,
           subtitle = main)
  } else {
    myboxplot <- plot_boxplot_base(boxplot_data) +
      scale_y_continuous(position = "right") +
      coord_flip(clip = "off") +
      theme_wide() +
      labs(x = y, y = x,
           subtitle = main)
  }

  ## save file
  if (is.na(filename)) {
    return(myboxplot)
  } else {
    if (n > 9) {
      ggplot2::ggsave(filename, myboxplot,height = (height-1.032066)*n/12+1.032066,width = width,dpi = 300,limitsize = FALSE)
    } else {
      ggplot2::ggsave(filename,myboxplot,height = height,width = width,dpi = 300)
    }
  }
}

# boxplot base function
#' @importFrom ggdist stat_halfeye median_qi stat_slabinterval
plot_boxplot_base <- function(boxplot_data) {
  p <- ggplot(boxplot_data, aes(x = label, y = value), show.legend = F) +
    stat_halfeye(
      adjust = .5,
      width = .6,
      .width = 0, show.legend = F,
      justification = -.3,
      point_colour = NA) +
    geom_boxplot(
      lwd=0.2,
      fatten = T,
      show.legend = F,
      width = .25,
      outlier.shape = NA
    )
    # gghalves::geom_half_point(
    #   ## draw jitter on the left
    #   side = "l",
    #   ## control range of jitter
    #   range_scale = .4,
    #   ## add some transparency
    #   alpha = .3
    # )

  return(p)
}

## for normal  boxplot
#' @importFrom ggtext element_markdown
theme_normal <- function(...) theme(...,
                                    # plot
                                    plot.margin = margin(t = .3, r = 1, b = .3, l = .3, unit = "cm"),
                                    panel.grid.major = element_line(colour = "snow4", size = 0.1),
                                    panel.grid.minor = element_line(linetype = "blank"),
                                    panel.background = element_rect(fill = NA),
                                    # strip
                                    strip.background = element_blank(),
                                    strip.text = element_blank(),
                                    # title
                                    plot.subtitle = element_text(colour = "olivedrab"),
                                    # legend
                                    legend.position=c(0.5,0.95),
                                    legend.direction = "horizontal",
                                    legend.background = element_rect(colour = "white", fill = "white"),
                                    legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
                                    legend.key = element_rect(colour = "transparent", fill = "transparent"),
                                    legend.text = element_markdown(lineheight = .8),
                                    legend.key.height = unit(1, "cm"),
                                    # axis
                                    axis.ticks = element_line(linetype = "blank"),
                                    axis.text.x = element_markdown()
)

## for wide boxplot
theme_wide <- function(...)  theme(...,
                                   # plot
                                   plot.margin = margin(t = .3, r = 1, b = .3, l = .3, unit = "cm"),
                                   panel.grid.major = element_line(colour = "snow4", size = 0.1),
                                   panel.grid.minor = element_line(linetype = "blank"),
                                   panel.background = element_rect(fill = NA),
                                   # strip
                                   strip.background = element_blank(),
                                   strip.text = element_blank(),
                                   # title
                                   plot.subtitle = element_text(colour = "olivedrab"),
                                   # legend
                                   legend.position= "top",
                                   legend.direction = "horizontal",
                                   legend.background = element_rect(colour = "white", fill = "white"),
                                   legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
                                   legend.key = element_rect(colour = "transparent", fill = "transparent"),
                                   legend.text = element_markdown(lineheight = .8),
                                   legend.key.height = unit(1, "cm"),
                                   # axis
                                   axis.title.y = element_text(hjust = 1, vjust = 1),
                                   axis.title.x = element_text(hjust = 0, vjust = 1),
                                   axis.text.x = element_markdown(),
                                   axis.ticks = element_line(linetype = "blank"))
# boxplot plot ------------------------------------------------------------

# ridges plot -------------------------------------------------------------
#' Title
#'
#' @param expr a expression matrix
#' @param group_list a character vector order by samples
#' @param palette group palette
#' @param filename file name you want to save plot
#' @param main plot title
#' @param width width of plot
#' @param height height of plot
#' @param x x axis title
#' @param y y axis title
#'
#' @importFrom ggridges geom_density_ridges theme_ridges
#'
#' @return a Ridges plot or a file
#' @export
#'
#' @examples
#' exprRidges(expr,group_list)
exprRidges <- function(expr,group_list,palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],
                       main = "Density of gene expression",
                       filename = NA,x= "Samples",y = "log2(cpm(count)+1)",
                       height = 5,width = 8.3) {

  expr_long = exprLong(expr,group_list = group_list,palette = palette)
  names(palette) <- unique(group_list)
  n = length(unique(expr_long$sample))

  p <- ggplot(expr_long, aes(x = value, y = label)) +
    ggridges::geom_density_ridges(
      aes(color = group, fill = group),
      jittered_points = F, scale = .95, rel_min_height = .01,
      size = 0.25
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = palette, labels =names(palette), guide = "none") +
    scale_color_manual(values = palette, labels =names(palette), guide = "none") +
    coord_cartesian(clip = "off") +
    guides(fill = guide_legend(
      override.aes = list(
        fill = palette,
        color = palette)
    )) +
    # ggtitle("Density of gene expression") +
    ggridges::theme_ridges(center = TRUE) +
    theme(axis.text.y = element_markdown(),legend.position = "none")+
    labs(x = y, y = x,
         subtitle = main)

  ## save file
  if (is.na(filename)) {
    return(p)
  } else {

  ggplot2::ggsave(filename, p,height = (height-1.032066)*n/12+1.032066,width = width,dpi = 300,limitsize = FALSE)

  }

}
# ridges plot -------------------------------------------------------------


# hallmark heatmap --------------------------------------------------------
#' Human housekeeping genes Heatmap
#'
#' plot a heatmap of highly uniform and strongly expressed Human housekeeping genes
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
#' @details
#' Human housekeeping genes revisited in
#' https://www.tau.ac.il/~elieis/HKG/
#' E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013)
#'
#' @return a heatmap or a file
#' @export
#'
#' @examples
#' exprHKGheatmap(expr,group_list)
exprHKGheatmap <- function(expr,group_list,filename = NA, main = "Human housekeeping genes",
                           palette = RColorBrewer::brewer.pal(3,"Set2")[1:2],anno_title = "Group",
                           annotation_col = ac_(expr,group_list),
                           width = ncol(expr)*0.3+2.2,height = 10*0.3+2.2,
                           ...,heatmapParam = list(show_rownames = T,cluster_cols = F)){

  # highly uniform and strongly expressed genes
  hk=c('C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')
  hk = hk[hk %in% rownames(expr)]
  # https://www.tau.ac.il/~elieis/HKG/
  # "Human housekeeping genes revisited"
  # E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013)
  # pheatmap::pheatmap(expr[hk,])

  m = expr[hk,]
  names(palette) <- unique(group_list)
  palette = split(palette,anno_title)

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
# hallmark heatmap --------------------------------------------------------

# help function -----------------------------------------------------------
## for heatmap构建分组信息
ac_ <- function(expr,group_list) {
  colD=data.frame(Groups=group_list)
  rownames(colD)=colnames(expr)
  return(colD)
}

## convert matrix to long data
#' @importFrom data.table setDT melt
#' @importFrom edgeR cpm
#' @importFrom glue glue
exprLong <- function(expr,group_list,palette = RColorBrewer::brewer.pal(3,"Set2")[1:2]) {

  expr <- data.table::setDT(as.data.frame(expr),keep.rownames = T)
  expr_long <- data.table::melt(expr, measure.vars = setdiff(colnames(expr),"rn"),
                                variable.name = "sample", value.name = "value")

  names(palette) <- unique(group_list)

  col_data <- data.table(
    sample = setdiff(colnames(expr),"rn"),
    group = group_list
  )

  col_data$color = palette[col_data$group]
  col_data$label = glue("<i style='color:{col_data$color}'>{col_data$sample}</i><br>(**{col_data$group}**)")

  expr_long2 <- merge(expr_long,col_data)

  return(expr_long2)

}
# help function -----------------------------------------------------------
