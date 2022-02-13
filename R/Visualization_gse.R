#' gseKEGG bar plot
#'
#' plot for enrich_gsekegg result
#'
#' @param data a enrich_gsekegg result
#' @param pvalue_cut filter by pvalue
#' @param enrichmentScore_cut filter by enrichmentScore
#' @param top top rows of up and down
#'
#' @import ggplot2
#' @importFrom dplyr top_n mutate
#' @importFrom forcats fct_reorder
#' @importFrom RColorBrewer brewer.pal
#' @importFrom shadowtext geom_shadowtext
#' @importFrom ggnewscale new_scale_color
#'
#' @return ggplot ob
#' @export
#'
#' @examples
#' \dontrun{
#' test <- enrich_gsekegg(DEG_df,x = "log2FoldChange")
#' geskegg_barplot(test)
#' }
GSEAbar <- function(res, top = 10) {

  dat <- res@result

  down_kegg <- dat[dat$enrichmentScore < 0,];
  up_kegg <- dat[dat$enrichmentScore >= 0,];


  down_kegg <- down_kegg[order(down_kegg[,"pvalue"])[1:top],]
  up_kegg <- up_kegg[order(up_kegg[,"pvalue"])[1:top],]

  dat=rbind(up_kegg,down_kegg)

  p <- ggplot(dat,
         aes(NES, fct_reorder(stringr::str_wrap(Description,25), NES),
             fill=qvalues)) +
    geom_col() +
    scale_fill_gradientn(colours=c("#b3eebe",
                                   "#46bac2", "#371ea3"),
                         guide=guide_colorbar(reverse=TRUE)) +
    theme_dose(12) +
    xlab("Normalized Enrichment Score") +
    ylab(NULL)

  return(p)

}

#' plot gseKEGG plot
#'
#' plot gseKEGG plots by  up and down
#'
#' @param res output from gse
#' @param top filter top by pvalue
#' @param pvalue_cut filter cut of pvalue
#' @param enrichmentScore_cut filter cut of enrichmentScore
#'
#' @importFrom enrichplot gseaplot2
#'
#' @return a list contains up and down
#' @export
#'
#' @examples
#' \dontrun{
#' gsekegg_res <- enrich_gsekegg(DEG_df,x = "log2FoldChange")
#' plots_l <- enhance_gseplot(gsekegg_res)
#' }
GSEAplot <- function(res, top = 10) {

  data = res
  down_kegg<-data[data$enrichmentScore < 0,]
  up_kegg<-data[data$enrichmentScore >= 0,]

  down_kegg <- down_kegg[order(down_kegg[,"pvalue"])[1:top],]
  up_kegg <- up_kegg[order(up_kegg[,"pvalue"])[1:top],]
  # down_kegg <- res %>% top_n(-top,wt = pvalue)
  # up_kegg <- res %>% top_n(top,wt = pvalue)

  plot_gseplot <- function(data,data_ud,x) {
    gseaplot2(res,data_ud$ID[x],
              title=stringr::str_wrap(data_ud$Description[x],30),pvalue_table = FALSE,subplots = 1:3,color = "blue") +
      labs(caption = paste(sprintf('Pvalue:  %.3f _ P.adjust: %.3f', data_ud$pvalue[x], data_ud$p.adjust[x]))) +
      theme(plot.caption = element_text(family = "Times", size = 8, face = "italic", colour = "dodgerblue"))
  }

  down_plots <- lapply(1:nrow(down_kegg), function(x)
    plot_gseplot(data,down_kegg,x)
  )

  up_plots <- lapply(1:nrow(up_kegg), function(x)
    plot_gseplot(data,up_kegg,x)
  )

  resl <- list(up_plots, down_plots)

  names(resl) <- c("up_plots", "down_plots")

  return(resl)

}
