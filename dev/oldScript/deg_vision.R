#' Title
#'
#' @param counts_data
#' @param deg_container
#' @param group_list
#' @param x
#' @param y
#' @param group
#' @param label
#' @param label_ns
#' @param volcano_palette
#' @param cut_FC
#' @param cut_FDR
#' @param volcano_top
#' @param size
#' @param expand
#' @param genes_list
#' @param highlight
#' @param heatmap_top
#' @param dir
#' @param prefix
#' @param heatmap_palette
#'
#'
#' @return
#' @export
#'
#' @examples
deg_vision <- function(counts_data,
                       deg_container,
                       group_list,
                       x = NULL,
                       y = NULL,
                       group = "group",
                       label = c("Down", "Stable", "Up"),
                       label_ns = "Stable",
                       volcano_palette = c("#2874C5", "grey", "#f87669"),
                       cut_FC = "auto",
                       cut_FDR = 0.05,
                       volcano_top = 10,
                       size = 2,
                       expand = c(0.25, 0.25),
                       genes_list = "top",
                       highlight = NULL,
                       heatmap_top = 50,
                       dir = ".",
                       prefix = "3-DEG_Vision",
                       heatmap_palette = RColorBrewer::brewer.pal(3, "Set2")[1:2]) {

  ana_type <- c("limma","edgeR","DESeq2")
  deg_group_l <- lapply(ana_type, function(x){

    index <- grep(x,slotNames(deg_container))
    index_name <- slotNames(deg_container)[index]
    deg_data <- slot(deg_container,index_name)

    ## 分组
    deg_group <- deg_group(deg_data)

    ## 火山图
    res <- deg_volcano(
      deg_data,
      x = x,
      y = y,
      group = group,
      label = label,
      label_ns = label_ns,
      palette = volcano_palette,
      cut_FC = cut_FC,
      cut_FDR = cut_FDR,
      top = volcano_top,
      size = size,
      expand = expand,
      genes_list = genes_list,
      highlight = highlight
    )

    volcano_file = glue("{dir}/{prefix}_{x}_volcano.pdf")
    ggsave(res,filename = volcano_file, width = 1600,height = 1600,units = "px",limitsize = FALSE)
    ui_done(glue("{x} volcano results were store in {volcano_file}."))

    ## 热图
    heatmap_prefix = glue("{prefix}_{x}")
    top_heatmap(
      counts_data = counts_data,
      deg_data = deg_data,
      group_list= group_list,
      x = x,
      y = y,
      top = heatmap_top,
      cut_FDR = cut_FDR,
      dir = dir,
      prefix = heatmap_prefix,
      palette = heatmap_palette
    )
    ui_done(glue("{x} top heatmap results were store in {heatmap_prefix}_top{top*2}_heatmap.pdf."))

    return(deg_group)
  })

  names(deg_data_l) <- ana_type

  deg_container@deg_df_limma <- deg_group_l[[]]

}

## This module will generate a grouped expression matrices and make some visualizations.

top_deg <- degTool::top_deg
deg_group <- degTool::deg_group
deg_volcano <- degPoint::deg_volcano
top_heatmap <- degPoint::top_heatmap

#' Run DEseq2
#'
#' A integrated function for run DEseq2 in a counts data and return results files.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param case_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param qc qc plots
#' @param x which column is log FC
#' @param y which column is P value
#' @param prefix a prefix of file names in this step
#' @param parallel if FALSE, no parallelization. if TRUE, parallel execution using BiocParallel
#'
#' @importFrom glue glue
#' @importFrom ggplot2 ggsave
#' @importFrom fs dir_exists dir_create
#'
#' @return a directory contains figures and csv files and a deg data frame
#' @export
#'
#' @examples
#' \dontrun{
#' deg_DESeq2(counts_input,group_list,
#'           case_group = "T", control_group = "C", qc = TRUE,
#'            x = "log2FoldChange", y = "pvalue",
#'            dir = tempdir(), prefix = "2-DEG_DEseq2")
#' }
deg_DESeq2 <- function(counts_data,group_list, parallel = F,
                       case_group,control_group,qc = TRUE,x,y,
                       dir = ".",prefix = "2-DEG_DEseq2") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  # deg_data <- run_DESeq2(counts_data = counts_data,group_list = group_list, parallel = parallel,
  #                        case_group = case_group,control_group=control_group,qc = qc,dir = dir,prefix = prefix)

  enhance_heatmap(counts_data, deg_data, group_list, x = x, y = y, dir = dir, prefix = prefix)
  message(glue("{emoji('deciduous_tree')} DESeq2 heatmap results were store in {dir}."))

  res <- enhance_volcano(deg_data,x = x, y = y,
                         label = c("Down","Stable","Up"), label_ns = "Stable",
                         palette =  c("#2874C5", "grey", "#f87669"),
                         cut_FC = "auto",cut_FDR = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
                         genes_list = "top", highlight = NULL)
  ggsave(res,filename = glue("{dir}/{prefix}_volcano.pdf"), width = 1600,height = 1600,units = "px",limitsize = FALSE)
  message(glue("{emoji('volcano')} DESeq2 volcano results were store in {dir}."))

  return(deg_data)

}

#' Run degeR
#'
#' A integrated function for run edgeR in a counts data and return results files.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param x which column is log FC
#' @param y which column is P value
#' @param prefix a prefix of file names in this step
#'
#' @importFrom glue glue
#' @importFrom ggplot2 ggsave
#' @importFrom fs dir_exists dir_create
#'
#' @return a directory contains figures and csv files  and a deg data frame
#' @export
#'
#' @examples
#' deg_edgeR(counts_input,group_list,
#'           control_group = "C",
#'            x = "logFC", y = "PValue",
#'            dir = tempdir(), prefix = "2-DEG_edgeR")
deg_edgeR <- function(counts_data,group_list,
                      control_group,x = "logFC", y = "PValue",
                      dir = ".",prefix = "2-DEG_edgeR") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  # deg_data <- run_edgeR(counts_data, group_list, control_group= control_group)

  enhance_heatmap(counts_data, deg_data, group_list, x = x, y = y, dir = dir, prefix = prefix)
  message(glue("{emoji('deciduous_tree')} edgeR heatmap results were store in {dir}."))

  res <- enhance_volcano(deg_data,x = x, y = y,
                         label = c("Down","Stable","Up"), label_ns = "Stable",
                         palette =  c("#2874C5", "grey", "#f87669"),
                         cut_FC = "auto",cut_FDR = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
                         genes_list = "top", highlight = NULL)
  ggsave(res,filename = glue("{dir}/{prefix}_volcano.pdf"), width = 1600,height = 1600,units = "px",limitsize = FALSE)
  message(glue("{emoji('volcano')} edgeR volcano results were store in {dir}."))

  return(deg_data)

}

#' Run limma
#'
#' A integrated function for run limma in a counts data and return results files.
#'
#' @param counts_data a counts data frame of rows in genes and columns in samples
#' @param group_list a list ordered by samples in counts_data
#' @param dir a directory to store results
#' @param case_group the name of the numerator level for the fold change (Test group)
#' @param control_group the name of the denominator level for the fold change (Control group)
#' @param x which column is log FC
#' @param y which column is P value
#' @param prefix a prefix of file names in this step
#'
#' @importFrom glue glue
#' @importFrom ggplot2 ggsave
#' @importFrom fs dir_exists dir_create
#'
#' @return a directory contains figures and csv files  and a deg data frame
#' @export
#'
#' @examples
#' deg_limma(counts_input,group_list,
#'           case_group = "T", control_group = "C",
#'            x = "logFC", y = "P.Value",
#'            dir = tempdir(), prefix = "2-DEG_limma")
deg_limma <- function(counts_data,group_list,
                      case_group,control_group,x,y,
                      dir = ".",prefix = "2-DEG_limma") {

  if (!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }

  # deg_data <- run_limma(counts_data = counts_data,group_list = group_list,
  #                       case_group = case_group,control_group=control_group)

  enhance_heatmap(counts_data, deg_data, group_list, x = x, y = y, dir = dir, prefix = prefix)
  message(glue("{emoji('deciduous_tree')} limma heatmap results were store in {dir}."))

  res <- enhance_volcano(deg_data,x = x, y = y,
                         label = c("Down","Stable","Up"), label_ns = "Stable",
                         palette =  c("#2874C5", "grey", "#f87669"),
                         cut_FC = "auto",cut_FDR = 0.05,top = 10, size = 2.0,expand = c(0.25,0.25),
                         genes_list = "top", highlight = NULL)
  ggsave(res,filename = glue("{dir}/{prefix}_volcano.pdf"), width = 1600,height = 1600,units = "px",limitsize = FALSE)
  message(glue("{emoji('volcano')} limma volcano results were store in {dir}."))

  return(deg_data)

}
