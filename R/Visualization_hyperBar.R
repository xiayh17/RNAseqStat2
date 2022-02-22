#' Enhanced barplot
#'
#' plot enrich result in bar plot
#'
#' @param data 'enrichResult' object, enrichGO result
#' @param fillstrip fill strip rect or not
#' @param split_color color for CC, MF, and BP
#' @param bar_color color for bar
#' @param print logic for print plot
#' @param showCategory Category numbers to show
#' @param by one of Count and GeneRatio
#' @param order logical
#' @param drop logical
#' @param split separate result by 'split' variable
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggfun element_roundrect
#' @importFrom grid grid.draw
#'
#' @return ggplot or grid plot
#' @export
#'
#' @examples
#' \dontrun{
#' test <- enrich_go(deg_data = DEG_df, x = "log2FoldChange", y = "pvalue")
#' enhance_barplot(test$Down,showCategory=10,split = "ONTOLOGY")
#' enhance_barplot(test$Down,showCategory=30)
#' }
enhance_barplot <- function(data, showCategory = 10,by = "Count",split=NULL,
                            order = FALSE,
                            drop = FALSE,
                            fillstrip = TRUE,
                            print = FALSE,
                            split_color = RColorBrewer::brewer.pal(3,"Dark2"),
                            bar_color = rev(RColorBrewer::brewer.pal(5,"GnBu")[3:5])
) {

  if (missing(split)){
    dat <- fortify.enrichResult(model = data, showCategory = showCategory, by = by, order = order, drop = drop) %>%
      sort_goTerms(by=by,split = split)

    p <- barplot_base2(dat, bar_color = bar_color, text_color = split_color[[1]])

    return(p)

  } else {
    dat <- fortify.enrichResult(model = data, showCategory = showCategory, by = by, order = order, drop = drop, split = split) %>%
      sort_goTerms(by=by,split = split)

    p <- barplot_base(dat, bar_color = bar_color, split_color = split_color)
    p <- p +
      facet_grid(ONTOLOGY~., scales="free", space="free_y", switch = "both") +
      theme(strip.background=element_roundrect(fill=NA, color=NA, r=0.31415,size = 0.5,linetype = "dotted")
            # ,legend.background = element_roundrect(fill=NA, color="grey80", r=0.31415,size = 0.5,linetype = "solid")
      )

    g <- ggplot_gtable(ggplot_build(p))
    strip_both <- which(grepl('strip-', g$layout$name))

    k <- 1

    if (fillstrip) {
      for (i in strip_both) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        m <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- split_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- split_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lty <- "solid"
        g$grobs[[i]]$grobs[[1]]$children[[m]]$children[[1]]$gp$col <- "white"
        k <- k+1
      }
    } else {
      for (i in strip_both) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        m <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- split_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- NA
        g$grobs[[i]]$grobs[[1]]$children[[m]]$children[[1]]$gp$col <- split_color[k]
        k <- k+1
      }
    }

    if (print) {
      grid.draw(g)
    }

    return(g)

  }

}

#' barplot_base
#'
#' plot basic barplot for enhance barplot
#'
#' @param data a data frame from enrichGO result
#' @param bar_color color for bar
#' @param split_color color for ont
#'
#' @import ggplot2
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return ggplot ob
barplot_base <- function(data, bar_color, split_color) {
  p <- ggplot(data,aes(Count, myY, fill = p.adjust)) +
    geom_shadowtext(aes(color = ONTOLOGY, label = Description),
                    x = 0,
                    nudge_y = -0.5,
                    hjust = -0.01, show.legend = F,
                    # position=position_dodge2(width=0.9),
                    size = 4,
                    bg.colour='#ffffff',bg.r = NA) +
    scale_color_manual(values = split_color)+
    new_scale_color() +
    geom_segment(aes(x=0,y=myY-1,xend = Count,yend = myY-1,color = p.adjust),size = 1, lineend = 'round') +
    scale_y_continuous(name = "Description", breaks = data$myY, labels = data$Description, expand = c(0, 0.5)) +
    # geom_col(width = .15,alpha = 1,aes(color = p.adjust),
    #          position=position_dodge2(padding = 0.9)) +
    scale_color_gradientn(colours = bar_color) +
    # geom_vline(aes(xintercept = 0,colour = ONTOLOGY),size = 0.5,show.legend = F)+
    # geom_text(aes(label = Description,color = ONTOLOGY),
    #           x = 0,
    #           hjust = 0,
    #           size = 4) +
    scale_x_continuous(expand = c(0, 0))+
    theme_minimal()+
    theme(legend.position = "right",
          # plot.margin=margin(t= 100,b=2,l=-2,r= 2,unit="pt"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dotted"),
          panel.grid.minor = element_blank(),
          # axis.line = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.ticks.y = element_blank())

  return(p)
}

#' barplot_base2
#'
#' plot basic barplot for enhance barplot not split
#'
#' @param data a data frame from enrichGO result
#' @param bar_color color for bar
#' @param text_color color for text
#'
#' @import ggplot2
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return ggplot ob
barplot_base2 <- function(data, bar_color, text_color) {
  p <- ggplot(data,aes(Count, myY, fill = p.adjust)) +
    geom_vline(aes(xintercept = 0,colour = text_color),size = 0.5,show.legend = F)+
    geom_shadowtext(aes(color = text_color, label = Description),
                    x = 0,
                    nudge_y = -0.5,
                    hjust = -0.01, show.legend = F,
                    # position=position_dodge2(width=0.9),
                    size = 4,
                    bg.colour='#ffffff',bg.r = NA) +
    scale_color_manual(values = text_color)+
    new_scale_color() +
    geom_segment(aes(x=0,y=myY-1,xend = Count,yend = myY-1,color = p.adjust),size = 1, lineend = 'round') +
    scale_y_continuous(name = "Description", breaks = data$myY, labels = data$Description, expand = c(0, 0.5)) +
    # geom_col(width = .15,alpha = 1,aes(color = p.adjust),
    #          position=position_dodge2(padding = 0.9)) +
    scale_color_gradientn(colours = bar_color) +
    # geom_text(aes(label = Description,color = ONTOLOGY),
    #           x = 0,
    #           hjust = 0,
    #           size = 4) +
    scale_x_continuous(expand = c(0, 0))+
    theme_minimal()+
    theme(legend.position = "right",
          legend.background = element_roundrect(fill=NA, color="grey80", r=0.31415,size = 0.5,linetype = "solid"),
          # plot.margin=margin(t= 100,b=2,l=-2,r= 2,unit="pt"),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dotted"),
          panel.grid.minor = element_blank(),
          # axis.line = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.ticks.y = element_blank())

  return(p)
}
