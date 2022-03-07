hyperBar <- function(res, top = 10) {

  dat <- res@result

  dat <- dat[order(dat[,"pvalue"])[1:top],]

  p <- ggplot(dat,
              aes(Count, fct_reorder(stringr::str_wrap(Description,25), Count),
                  fill=qvalue)) +
    geom_col() +
    scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                         guide=guide_colorbar(reverse=TRUE)) +
    theme_dose(12) +
    xlab("Count") +
    ylab(NULL)

  return(p)

}

#' @importFrom ggfittext geom_fit_text
#' @importFrom ggnewscale new_scale_color
#' @importFrom ggfun element_roundrect
#' @import ggplot2
enrichBar <- function(enrichResult,
                      top = 10, group = "ONTOLOGY", space = 0.9,
                      order_by = "pvalue", FDR = "qvalue", bar = "logP", bar_size = 2,
                      point = "gr", point_shape = 21, point_size = 2.5,
                      point_color = c('#de2d26','#fc9272','#fee0d2'),
                      bar_color = rev(RColorBrewer::brewer.pal(5,"GnBu")[3:5]),
                      group_color = RColorBrewer::brewer.pal(3,"Dark2"),
                      group_names = c("Biological Process", "Cellular Component", "Molecular Function"),
                      group_color_name = c("BP","CC","MF"),
                      group_title = "ONTOLOGY",
                      bar_title = "Qvalue",
                      point_title = "GeneRatio",
                      x_title = "-log10(pvalue)",
                      y_title = "Description",
                      plot_title = NULL,
                      legend_text_size = 12,
                      FDR_color = c('#de2d26','#fc9272','#fee0d2')) {

  dat <- textBarData(enrichResult = enrichResult,top = top,group = group,order_by = order_by)

  max_x <- max(dat[,bar],na.rm = T)

  if(is.null(group)) {

    p <- ggplot(dat,aes_(x = as.name(bar),y = ~helpY)) +
      geom_rect(
        xmin = 0, xmax = max_x,
        aes(ymin = helpY - space,
            ymax = helpY + space),
        fill = NA
      )+
      geom_fit_text(
        xmin = 0, xmax = max_x,
        aes_(ymin = ~helpY - space, ymax = ~helpY + space,
             label = ~Description),
        color = group_color[1],
        grow = TRUE, reflow = TRUE,fullheight = TRUE,
        place = "left", show.legend = F
      )

  } else {

    if(group %in% colnames(dat)) {

      p <- ggplot(dat,aes_(x = as.name(bar),y = ~helpY)) +
        geom_rect(
          xmin = 0, xmax = max_x,
          aes(ymin = helpY - space,
              ymax = helpY + space),
          fill = NA
        )+
        geom_fit_text(
          xmin = 0, xmax = max_x,
          aes_(ymin = ~helpY - space, ymax = ~helpY + space,
               label = ~Description,color = as.name(group)),
          grow = TRUE, reflow = TRUE,fullheight = TRUE,
          place = "left", show.legend = TRUE
        )+
        scale_color_manual(values = group_color,
                           labels = group_names,
                           breaks = group_color_name,
                           guide = guide_enrichLegend(order = 1, title = group_title,
                                                      label.theme = element_text(size = legend_text_size)
                           ))

    } else {

      p <- ggplot(dat,aes_(x = as.name(bar),y = ~helpY)) +
        geom_rect(
          xmin = 0, xmax = max_x,
          aes(ymin = helpY - space,
              ymax = helpY + space),
          fill = NA
        )+
        geom_fit_text(
          xmin = 0, xmax = max_x,
          aes_(ymin = ~helpY - space, ymax = ~helpY + space,
               label = ~Description),
          color = group_color[1],
          grow = TRUE, reflow = TRUE,fullheight = TRUE,
          place = "left", show.legend = F
        )

    }

  }

  p <- p + scale_x_continuous(name = x_title, expand = c(0, 0),limits = c(0,max_x*1.1))+
    new_scale_color() +
    geom_segment(
      aes_(xend = as.name(bar),color = as.name(FDR),
           y = ~helpY - space, yend = ~helpY - space),
      x = 0,
      size = bar_size,lineend = 'round',
    ) +
    scale_color_gradientn(colours = bar_color,
                          guide = guide_enrichBar(order = 2,title = bar_title,
                                                  label.theme = element_text(size = legend_text_size)
                          ))+ ## bar 的颜色
    new_scale_color() +
    geom_point(
      aes_(y = ~helpY-space, x = as.name(bar),
           fill = as.name(point),color = as.name(point)),
      shape = point_shape,
      size = point_size
    ) +
    scale_fill_gradientn(colours = point_color,guide = "none")+ # 点的填充色
    scale_color_gradientn(colours = point_color,
                          guide = guide_enrichBar(order = 3,title = point_title,
                                                  label.theme = element_text(size = legend_text_size)
                          ))+ # 点的边界色
    scale_y_continuous(name = y_title, breaks = dat$helpY, expand = c(0, 0.5)) +
    theme_enrichBar() +
    ggtitle(plot_title)

  ## 如果有分组，添加
  if(is.null(group)) {

    g <- p

  } else {

    if(group %in% colnames(dat)) {

      p <- p +
        facet_grid(get(group)~., scale="free", switch="both") +
        theme(strip.background=element_roundrect(fill=NA, color=NA, r=0.31415,size = 0.5,linetype = "dotted"))

      g <- ggplot_gtable(ggplot_build(p))
      strip_both <- which(grepl('strip-', g$layout$name))

      k <- 1

      for (i in strip_both) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        m <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- group_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lty <- "solid"
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- group_color[k]
        g$grobs[[i]]$grobs[[1]]$children[[m]]$children[[1]]$gp$col <- "white"
        k <- k+1
      }

    } else {

      g <- p

    }

  }

  return(g)
}

guide_enrichLegend <- function(...) {

  guide_legend(...,
               title.position = 'right',
               title.vjust = 0,
               direction = "horizontal",
               reverse = T,
               keywidth = unit(8, 'lines'),
               keyheight = unit(.1, 'lines'),
               override.aes = list(size = 1),
               label.position = "top")

}

guide_enrichBar <- function(...) {

  guide_colorbar(...,
                 title.position = 'right',
                 title.vjust = 1,
                 reverse = T,
                 barwidth = unit(25, 'lines'),
                 barheight = unit(.5, 'lines'))

}

theme_enrichBar <- function(...) {
  theme_minimal()+
    theme(..., legend.position = "bottom",legend.box = "vertical",
          legend.box.just = "left",
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

#' @importFrom tibble add_column
#' @importFrom dplyr pull arrange slice_head mutate group_by across
#' @importFrom forcats fct_reorder
#' @importFrom usethis ui_info ui_oops
textBarData <- function(enrichResult,
                        top = 10, group = "ONTOLOGY",
                        order_by = "pvalue"
) {
  ## 判断数据类型
  if(identical(class(enrichResult)[1],"enrichResult")) {

    ## 添加richFactor
    result <- mutate(enrichResult@result, richFactor = Count / as.numeric
                     (sub("/\\d+", "", BgRatio)))

  } else {

    ## 添加richFactor
    result <- mutate(enrichResult, richFactor = Count / as.numeric
                     (sub("/\\d+", "", BgRatio)))

  }

  if (!is.null(group)) {

    if(group %in% colnames(result)) {

      usethis::ui_info("grouped by {group} and choose top {top} of every group")
      ## 分组排序取top
      dat_go <- result %>%
        group_by(across({{ group }})) %>%
        arrange({{ order_by }}, .by_group = TRUE) %>% # 排序 pvalue 升序
        slice_head(n = top) # 排序后数据的top

    } else {

      usethis::ui_oops("group cant found in enrichResult, group is ignored...")
      ## 排序取top
      dat_go <- result %>%
        arrange({{ order_by }}, .by_group = TRUE) %>% # 排序 pvalue 升序
        slice_head(n = top) # 排序后数据的top

    }

  } else {

    usethis::ui_info("top {top} choosed")
    ## 排序取top
    dat_go <- result %>%
      arrange({{ order_by }}, .by_group = TRUE) %>% # 排序 pvalue 升序
      slice_head(n = top) # 排序后数据的top

  }

  ## 辅助绘图数据处理
  dat_go$Description <- factor(pull(dat_go, "Description"))
  dat_go$Description <- fct_reorder(pull(dat_go, "Description"), pull(dat_go, order_by), .desc = FALSE)
  dat_go <- dat_go %>% tibble::add_column(helpY = seq(2,length(.$Description)*2,2))
  dat_go <- dat_go %>% tibble::add_column("logP" = -log10(pull(., order_by)))
  dat <- dat_go %>% tibble::add_column(gr =
                                         as.numeric(sub("/\\d+", "", .$GeneRatio))/as.numeric(sub("\\d+/", "", .$GeneRatio)))

  ##
  dat <- as.data.frame(dat)
  return(dat)

}

