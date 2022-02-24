# basical volcano ---------------------------------------------------------
#' Basic DEG volcano basic plot
#'
#' more beautiful and clear DEG volcano plot
#'
#' @param object a grouped DEGContainer
#'
#' @import ggplot2
#'
#' @return a ggplot object
#'
#' @noRd
#' @examples
#' BaseVolcanoPlot(object)
BaseVolcanoPlot <- function(object,which,category = "H") {

  if (which == "limma") {
    deg_data = limma_res(object)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(object)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(object)
  } else if (which == "MSigDB") {
    if(!is.null(category)) {

      deg_data = msigdbGSVAresult(object)[["GSVA_diff"]][[category]]

    } else {

      ui_stop("when {ui_code('which')} set as {ui_value('MSigDB')}, {ui_code('category')} should be one of category in {ui_value('MSigDB')}")

    }
  }else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2, MSigDB')}")
  }

  x = FC_Identify(deg_data)
  y = pvalue_Identify(deg_data)

  if (which == "MSigDB") {

    treat = msigdbTreat(object)
    label = treat@label
    label_ns = treat@label_ns
    palette = treat@sigCol
    cut_FDR = treat@cutFDR
    cut_FC = treat@cutFC

    if(is.list(cut_FC)){
      cut_FC = cut_FC[[category]]
    } else {
      cut_FC = cut_FC
    }

  } else {

    label = label(object)
    label_ns = label_ns(object)


    if(is.list(cutFC(object))){
      cut_FC = cutFC(object)[[which]]
    } else {
      cut_FC = cutFC(object)
    }
    cut_FDR = cutFDR(object)
    palette =  sigCol(object)

  }

  names(palette) <- label

  if (length(cut_FC) == 1) {
    cut_FC <- c(-cut_FC, cut_FC)
  }

  FC_data <- data.frame(
    cut_FC = cut_FC
  )

  volcano <- create_volcano(object= object,which = which, category = category)

  group_count <- table(deg_data$group)
  updown_label <-  setdiff(label,label_ns)

  p <- ggplot(volcano,aes(x = get(x), y = log10Pvalue)) +
    geom_point(aes(
      color = point.color,
      size = point.size,
      alpha = point.alpha,
      shape = point.shape))  +
    scale_colour_identity(name = "Group",labels = levels(volcano$group),breaks = get_breaks(volcano,"point.color"),guide = "legend") +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_shape_identity() +
    geom_vline(data=FC_data, mapping=aes(xintercept=cut_FC),
               color = palette[updown_label], linetype="longdash",size = 0.1) +
    geom_hline(aes(yintercept = -log10(cut_FDR)),
               color = palette[label_ns], linetype="longdash", size = 0.1)+
    theme_volcano() +
    coord_cartesian(expand = T,clip = 'off') +
    labs(colour = "Group",x = "log2FoldChange", y = "-log10(PValue)",
         subtitle = NULL,
         caption = paste(sprintf('%s  %.3f;',y, cut_FDR),
                         sprintf('FC %.3f;',cut_FC),
                         paste(updown_label,group_count[updown_label],collapse = "; "),
                         sprintf('; Total: %1.0f',nrow(deg_data))
         ))

  return(p)

}
# basical volcano ---------------------------------------------------------

# theme of volcano --------------------------------------------------------
#' theme_volcano
#'
#' a nice theme for DEG volcano
#'
#' @param ... param passed from theme
#'
#' @importFrom ggplot2 theme element_line element_text element_rect
#'
#' @return a ggplot theme
#'
#' @noRd
theme_volcano <- function(...) {theme(...,
                                      axis.line = element_line(size = 0.2, linetype = "solid"),
                                      axis.ticks = element_line(size = 0.2),
                                      panel.grid.major = element_line(colour = "gray", size = 0.05,linetype = "dotted"),
                                      panel.grid.minor = element_line(linetype = "blank"),
                                      panel.border = element_rect(colour = "black",fill = NA,size = 0.1),
                                      axis.title = element_text(family = "Times"),
                                      axis.text = element_text(family = "Times",color = "black"),
                                      axis.line.y = element_line(linetype = "blank"),
                                      axis.line.x = element_line(linetype = "blank"),
                                      axis.text.x = element_text(family = "Times"),
                                      axis.text.y = element_text(family = "Times"),
                                      axis.title.x = element_text(family = "Times",face = "bold"),
                                      axis.title.y = element_text(family = "Times",face = "bold"),
                                      legend.text = element_text(family = "Times"),
                                      legend.title = element_text(family = "Times"),
                                      panel.background = element_rect(fill = NA),
                                      legend.key = element_rect(fill = NA),
                                      legend.background = element_rect(fill = NA),
                                      legend.position = "top", legend.direction = "horizontal",
                                      plot.margin = margin(10, 25, 10, 10),
                                      plot.caption = element_text(hjust = 1,family = "Times", size = 6, face = "italic", colour = "black"))}
# theme of volcano --------------------------------------------------------

# point volcano -----------------------------------------------------------
#' Showing gene labels and highlighting gene points in a volcano plot in ggplot
#'
#' specific some genes to label
#'
#' @param object an DEGContainer
#' @param which limma edgeR or DESeq2
#' @param gene number or gene name
#' @param light number or gene name
#' @param light_color character
#' @param light_label_color character
#' @param expend c(0.12, 0.12)
#'
#' @importFrom data.table data.table
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @return a ggplot ob
#'
#' @noRd
PointVolcano <- function(object,which, category = "H",
                         gene,light = NULL,
                         light_color = "#24ac56",
                         light_label_color = "#24ac56",
                         expend = c(0.12, 0.12)) {


  ## basic plot
  p <- BaseVolcanoPlot(object,which = which,category = category)

  if (which == "limma") {
    deg_data = limma_res(object)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(object)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(object)
  } else if (which == "MSigDB") {
    if(!is.null(category)) {

      deg_data = msigdbGSVAresult(object)[["GSVA_diff"]][[category]]

    } else {

      ui_stop("when {ui_code('which')} set as {ui_value('MSigDB')}, {ui_code('category')} should be one of category in {ui_value('MSigDB')}")

    }
  }else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2, MSigDB')}")
  }

  x = FC_Identify(deg_data)

  if (!is.null(gene)) {
    ## gene label data
    gene_label <- genePoint(object = object,which = which,gene = gene, category = category)

    ## up and down
    label_up <- gene_label[which(gene_label[,x] >= 0),]
    label_down <- gene_label[which(gene_label[,x] < 0),]

    ## add gene label
    p <- p +
      scale_x_continuous(expand = c(0.12, 0.12)) +
      geom_volcano_text(data = label_up,
                        mapping = aes(label = rn,color = point.color),
                        nudge_y      = -6.5,
                        hjust        = 0,
                        min.segment.length = 0,
                        nudge_x =  volcano_nudge_x_up(object = object, label_data = label_up,which = which, category = category)) +
      geom_volcano_text(data = label_down,
                        mapping = aes(label = rn,color = point.color),
                        nudge_y      = -4,
                        hjust        = 1,
                        min.segment.length = 0,
                        nudge_x =  volcano_nudge_x_down(object = object, label_data = label_down,which = which, category = category))
  }

  if (!is.null(light)) {
    ## highlight data
    light_label <- genePoint(object = object,which = which,gene = light, category = category)

    ## add highlight point
    p <- p + geom_volcano_point(data = light_label,stroke = 0.1,color = light_color)

    ## highlight up and down
    light_up <- light_label[which(light_label[,x] >= 0),]
    light_down <- light_label[which(light_label[,x] < 0),]

    ## add highlight label
    p <- p +
      geom_volcano_text(data = light_up,
                        mapping = aes(label = rn),color = light_label_color,
                        nudge_y      = -6.5,
                        hjust        = 0,
                        min.segment.length = 0,
                        nudge_x =  volcano_nudge_x_up(object = object, label_data = light_up,which = which)) +
      geom_volcano_text(data = light_down,
                        mapping = aes(label = rn),color = light_label_color,
                        nudge_y      = -4,
                        hjust        = 1,
                        min.segment.length = 0,
                        nudge_x =  volcano_nudge_x_down(object = object, label_data = light_down,which = which))
  }

  return(p)

}
# point volcano -----------------------------------------------------------

# help function -----------------------------------------------------------
volcano_nudge_x_up <- function(object,label_data,which,just = 0, category = "H") {

  if (which == "limma") {
    deg_data = limma_res(object)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(object)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(object)
  } else if (which == "MSigDB") {
    if(!is.null(category)) {

      deg_data = msigdbGSVAresult(object)[["GSVA_diff"]][[category]]

    } else {

      ui_stop("when {ui_code('which')} set as {ui_value('MSigDB')}, {ui_code('category')} should be one of category in {ui_value('MSigDB')}")

    }
  }else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2, MSigDB')}")
  }

  x = FC_Identify(deg_data)

  max_lfc <- max(deg_data[,x],na.rm = T)

  ## up
  value <- label_data[which(label_data[,x] > 0),][,x] + max_lfc+max_lfc/4 + just

  return(value)

}

volcano_nudge_x_down <- function(object,label_data,which,just = 0, category = "H") {

  if (which == "limma") {
    deg_data = limma_res(object)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(object)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(object)
  } else if (which == "MSigDB") {
    if(!is.null(category)) {

      deg_data = msigdbGSVAresult(object)[["GSVA_diff"]][[category]]

    } else {

      ui_stop("when {ui_code('which')} set as {ui_value('MSigDB')}, {ui_code('category')} should be one of category in {ui_value('MSigDB')}")

    }
  }else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2, MSigDB')}")
  }

  x = FC_Identify(deg_data)

  min_lfc <- min(deg_data[,x],na.rm = T)

  ## down
  value <- label_data[which(label_data[,x] < 0),][,x] + min_lfc+min_lfc/4 + just

  return(value)

}

geom_volcano_text <- function(data,
                              mapping = NULL,
                              nudge_x = NULL,
                              nudge_y      = -1,
                              hjust        = 0,
                              size = 1.8,
                              direction    = "y",
                              segment.size = 0.1,
                              segment.linetype = 6,
                              max.overlaps = 10,
                              max.iter = 1000000,
                              max.time = 10,
                              min.segment.length = 0,
                              fontface = "bold",
                              family = "Times", ...)  {
  geom_text_repel(data = data,
                  size = size,
                  nudge_y = nudge_y,
                  hjust = hjust,
                  nudge_x = nudge_x,
                  direction    = direction,
                  segment.size = segment.size,
                  segment.linetype = segment.linetype,
                  max.overlaps = max.overlaps,
                  max.iter = max.iter,
                  max.time = max.time,
                  min.segment.length = min.segment.length,
                  fontface = fontface,
                  family = family, mapping = mapping, ...)
}

geom_volcano_point <- function(data,
                               shape = 1,
                               stroke = 0.15,
                               color = "#24ac56",
                               fill = NA,
                               ...,mapping = NULL) {

  geom_point(data = data,
             shape = shape,
             stroke = stroke,
             color = color,
             fill = fill,
             ...,mapping = mapping)

}

get_breaks <- function(volcano, column_name) {
  pat <- unique(volcano[,c("group",column_name)])
  pat <- pat[match(levels(volcano$group),pat[,"group"]),]
  res <- pat[,column_name]
  return(res)
}

genePoint <- function(object,which, gene, category = "H") {

  volcano <- create_volcano(object,which = which,category = category)

  if (is.character(gene)&any(gene %in% volcano[,"rn"])) {
    gene <- data.frame(
      rn = gene
    )
    # look up piont for labels
    label_data <- merge(gene,volcano)
  } else if (is.numeric(gene)&length(gene)==1) {

    gene = topGene(object = object,topSig = gene,which = which,category = category)

    gene <- data.frame(
      rn = gene
    )
    # look up piont for labels
    label_data <- merge(gene,volcano)

    ui_info("top {nrow(gene)} gene of Up and Down were Choosed.")

  } else {

    usethis::ui_oops("Make Sure You Gene Name match on you data! ")

  }

  if (nrow(label_data)>0) {

    ui_done(glue("{nrow(label_data)/nrow(gene)*100}% label match on your data."))

  } else {

    usethis::ui_oops("Make Sure You Gene Name match on you data! ")

  }

  return(label_data)

}

create_volcano <- function(object, which = 'MSigDB', category = "H") {

  deg_data_list <- msigdbGSVAresult(test)[["GSVA_diff"]]

  if (which == "limma") {
    deg_data = limma_res(object)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(object)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(object)
  } else if (which == "MSigDB") {
    if(!is.null(category)) {

      deg_data = msigdbGSVAresult(object)[["GSVA_diff"]][[category]]

    } else {

      ui_stop("when {ui_code('which')} set as {ui_value('MSigDB')}, {ui_code('category')} should be one of category in {ui_value('MSigDB')}")

    }
  }else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2, MSigDB')}")
  }

  if (which == "MSigDB") {

    x = FC_Identify(deg_data)
    y = pvalue_Identify(deg_data)
    treat = msigdbTreat(object)
    sigGroup = treat@label
    sigCol = treat@sigCol
    sigAlpha = treat@sigAlpha
    sigSize = treat@sigSize
    sigShape = treat@sigShape

  } else {

    x = FC_Identify(deg_data)
    y = pvalue_Identify(deg_data)
    sigGroup = label(object)
    sigCol = sigCol(object)
    sigAlpha = sigAlpha(object)
    sigSize = sigSize(object)
    sigShape = sigShape(object)

  }

  if(length(sigCol) == length(sigGroup)) {
    names(sigCol) = sigGroup
    # deg_data[,"point.color"] = factor(sigCol[deg_data$group], levels = sigCol)
  } else if (length(sigCol) == 1) {
    sigCol = rep(sigCol, length(sigGroup))
    names(sigCol) = sigGroup
  }

  if(length(sigAlpha) == length(sigGroup)) {
    names(sigAlpha) = sigGroup
    # deg_data[,"point.alpha"] = factor(sigAlpha[deg_data$group], levels = sigAlpha)
  } else if (length(sigAlpha) == 1) {
    sigAlpha = rep(sigAlpha, length(sigGroup))
    names(sigAlpha) = sigGroup
  }

  if(length(sigSize) == length(sigGroup)) {
    names(sigSize) = sigGroup
    # deg_data[,"point.size"] = factor(sigSize[deg_data$group], levels = sigSize)
  } else if (length(sigSize) == 1) {
    sigSize = rep(sigSize, length(sigGroup))
    names(sigSize) = sigGroup
  }

  if(length(sigShape) == length(sigGroup)) {
    names(sigShape) = sigGroup
    # deg_data[,"point.shape"] = factor(sigShape[deg_data$group], levels = sigShape)
  } else if (length(sigShape) == 1) {
    sigShape = rep(sigShape, length(sigGroup))
    names(sigShape) = sigGroup
  }


  deg_data[,"point.color"] = sigCol[deg_data$group]
  deg_data[,"point.alpha"] = sigAlpha[deg_data$group]
  deg_data[,"point.size"] = sigSize[deg_data$group]
  deg_data[,"point.shape"] = sigShape[deg_data$group]
  deg_data[,"rn"] = rownames(deg_data)

  y = -log10(deg_data[,y])

  deg_data[,"log10Pvalue"] = y

  return(deg_data)
}
# help function -----------------------------------------------------------
