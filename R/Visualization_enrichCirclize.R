#' Plot circos of common ID
#'
#' After a enrich analysis across limma, edgeR, DESeq2 results. How about the distribution among
#' different datset.
#'
#' @param result_g results of \code{\link{modelEnrich}} or similar format data
#' @param group which column used to group data
#' @param height height of plot
#' @param width width of plot
#' @param mar margin of bottom, left, top, right around circos
#' @param filename filename
#' @param IDpalette color for unique ID
#' @param heatCircle1 which column used to first outer circle
#' @param heatCircle2 which column used to second outer circle
#' @param heatCircle1Pal color of first outer circle
#' @param heatCircle2Pal color of second outer circle
#' @param groupPal color for every group
#' @param group_name group name
#'
#' @importFrom circlize circos.clear chordDiagram circos.track
#' @importFrom Cairo Cairo
#' @importFrom grDevices dev.off pdf
#' @importFrom grid unit pushViewport viewport grid.draw upViewport
#'
#' @return a file store plot
#' @export
#'
#' @examples
#' compareEnrichCircle(result_g)
compareEnrichCircle <- function(result_g,group = "model",
                                height = 12000,width = 10000,
                                mar = c(8,0,0,17),filename = "circlize_cross.pdf",
                                IDpalette = as.character(MetBrewer::met.brewer("Signac",length(unique(result_g$ID)))),
                                heatCircle1 = "Count",
                                heatCircle2 = "qvalue",
                                heatCircle1Pal = c('#e5f5e0','#a1d99b','#31a354'),
                                heatCircle2Pal = c('#efedf5','#bcbddc','#756bb1'),
                                groupPal = c('#1f78b4','#33a02c','#ff7f00','#b15928'),
                                group_name = c("edgeR","limma","DESeq2","merge")){

  # 构建link 矩阵
  mat_df <- linkCir(df = result_g,group = group)

  # 分组
  group <- groupCir(matrix = mat_df)

  # 调色
  pal <- paletteCir(df = result_g,group = group,palette = IDpalette)

  # 初步绘图
  pdf(file = tempfile(fileext = "pdf"))
  cnm_o = chordDiagram(mat_df, grid.col = pal,
                       directional = FALSE, annotationTrack = "grid",
                       big.gap = 10, small.gap = 1,
                       preAllocateTracks = list(track.height = 0.05),
                       link.target.prop = FALSE)
  circos.clear()
  dev.off()

  # 隐藏孤儿
  col2 = hideCir(matrix = mat_df,preData = cnm_o)

  # 注释数据
  heatCount <- heatCir(matrix = mat_df, df = result_g,look = heatCircle1)
  heatQvalue <- heatCir(matrix = mat_df, df = result_g,look = heatCircle2)

  # 注释配色
  countPal <- heatCircle1Pal
  qvaluePal <- heatCircle2Pal
  groupPal = groupPal
  names(groupPal) <- group_name

  Cairo::Cairo(height = height, width = width, file=filename, type="pdf", bg="white",dpi = 300,units = "px")
  par(mar = mar)
  chordDiagram(mat_df,group = group,
               grid.col = pal ,col=col2, directional = F,transparency = 0,
               big.gap = 10, small.gap = 1,
               link.sort = TRUE,
               link.decreasing = TRUE,
               annotationTrack = "grid", preAllocateTracks = 1)

  circos.track(ylim = c(0, 1), panel.fun =

                 rectCir(preData = cnm_o, annoData = heatCount, palatte = countPal, track.index = 1,
                         ylim = c(0,1), y1Scale = 0.1, y2Scale = 0.5)

               , track.height = 0.0,bg.border = NA, track.index = 1,track.margin = c(0,0))

  circos.track(ylim = c(0, 1), panel.fun =

                 rectCir(preData = cnm_o, annoData = heatQvalue, palatte = qvaluePal, track.index = 1,
                         ylim = c(0,1), y1Scale = 0.5, y2Scale = 0.9)

               , track.height = 0,bg.border = NA,track.index = 1,track.margin = c(0,0))

  highlightGroup(df = result_g,track.index = 1,palatte = groupPal, lwd = 1,
                 text.vjust = "10mm", cex = 0.9,text.col = "black",
                 padding = c(0, 0, 0, 0))

  ## 构建注释
  countLegend <- heatLegend(annoData = heatCount,palatte = countPal,title = heatCircle1,int = TRUE,legend_width = unit(6, 'cm'),
                            direction = "horizontal",title_position = "lefttop")
  qvalueLegend <- heatLegend(annoData = heatQvalue,palatte = qvaluePal,title = heatCircle2,int = F,legend_width = unit(6, 'cm'),
                             direction = "horizontal",title_position = "lefttop")
  idLegend <- idLegend(df = result_g,palatte = pal,labels_gp = gpar(fontsize = 10))

  ## 绘制注释

  ## 热图图例在圆的正下方
  y_coord <- 0.1
  x_coord <- 0.3
  pushViewport(viewport(x = x_coord, y = y_coord))
  grid.draw(countLegend)
  upViewport()

  pushViewport(viewport(x = x_coord, y = y_coord + 0.1))
  grid.draw(qvalueLegend)
  upViewport()

  ## ID 图例在正右方
  pushViewport(viewport(x = 0.75, y = 0.75))
  grid.draw(idLegend)
  upViewport()

  ## 清理和关闭
  circos.clear()
  dev.off()

}


# help functions -----------------------------------------------------------
# 用于生成基本link矩阵的函数---------------------------------------------------------------------
#' Create matrix for circlize chordDiagram
#'
#' @param df a data frame contains multigroup results of enrich
#' @param group which column is group data
#'
#' @importFrom reshape2 dcast
#' @importFrom textshape column_to_rownames
#'
#' @return a matrix contains relations of IDs in group
#' @export
#'
#' @examples
#' linkCir(df)
linkCir <- function(df, group = "model") {

  # 构建一个分组与ID合并的列
  dat_d = df
  dat_d[,"u_id"] <- paste(dat_d[,get(group)],dat_d$ID,sep = "_")

  # 形成一个ID之间的组合, 每个组合的初始值为零
  mat_o <- merge(dat_d$u_id,dat_d$u_id)
  mat_o$value <- 0

  # 计算ID的重复数量，重复大于一的为有link ，值为1
  count_id <- as.data.frame(table(dat_d$ID))
  count_ids <- count_id$Var1[which(count_id$Freq > 1)] %>% as.character()
  mat_o$value[gsub(".*_","",mat_o$x) %in% count_ids] <- 1

  # 去除 相同分组相同ID的 组合
  z_1 <- grepl(TRUE,gsub("_.*","",mat_o$x) == gsub("_.*","",mat_o$y))
  mat_o$value[z_1] <- 0

  # 去除 不同分组 并且 不同ID的 组合
  z_2 <- grepl(TRUE,(gsub("_.*","",mat_o$x) != gsub("_.*","",mat_o$y))&(gsub(".*_","",mat_o$x) != gsub(".*_","",mat_o$y)))
  mat_o$value[z_2] <- 0

  ## 将长数据转为宽数据
  mat <- reshape2::dcast(mat_o,x~y)
  ## 第一列作为行名
  mat <- textshape::column_to_rownames(mat,1)
  ## 去掉对称link
  mat[lower.tri(mat)] <- 0
  ## 去掉矩阵对角线link
  diag(mat) <- 0

  ## keep poor 有一些ID没有任何link 但是需要暂时给一个link 才能在图中显示
  count_poor <- count_id$Var1[which(count_id$Freq == 1)] %>% as.character()
  count_poor_all <- data.frame(
    id = gsub(".*_","",mat_o$x)[gsub(".*_","",mat_o$x) %in% count_poor],
    index =  grep(TRUE,gsub(".*_","",mat_o$x) %in% count_poor)
  )
  count_poor_unique <- count_poor_all[!duplicated(count_poor_all$id),]
  count_poor_unique_ids <- count_poor_unique$id

  ispoor <- gsub(".*_","",rownames(mat)) %in% count_poor_unique_ids

  ispoor_m <- matrix(rep(ispoor,ncol(mat)),ncol(mat),ncol(mat))

  poor_unique <- upper.tri(mat)&ispoor_m
  poor_unique <- apply(poor_unique, 1, function(x){

    grep(T,x)[1]

  })

  if (any(ispoor_m[,nrow(ispoor_m)])) {

    poor_unique[length(poor_unique)] <- length(poor_unique)

  }

  pu_index_row <- grep("[0-9]",poor_unique)
  pu_index_col <- na.omit(poor_unique)
  # print(length(pu_index_row))
  # print(length(pu_index_col))

  for (i in seq_along(pu_index_row)) {

    mat[pu_index_row[i],pu_index_col[i]] <- 2

  }

  # mat_o$value[count_poor_unique$index] <- 2

  ## 转化为matrix格式
  mat <- as.matrix(mat)

  return(mat)

}

# 用于调色的函数 -----------------------------------------------------------------
#' Create a palette of ID across different group
#'
#' @param df a data frame contains multigroup results of enrich
#' @param group which column is group data
#' @param palette a color vecter, should be length with unique ID
#'
#' @importFrom MetBrewer met.brewer
#'
#' @return a named palette
#' @export
#'
#' @examples
#' paletteCir(df)
paletteCir <- function(df, group = "model",
                       palette = as.character(MetBrewer::met.brewer("Signac",length(unique(df$ID))))) {

  if (length(unique(df$ID)) != length(palette)) {
    usethis::ui_stop("Colors number should be {length(unique(df$ID))}")
  }

  dat_d = df
  dat_d[,"u_id"] <- paste(dat_d$model,dat_d$ID,sep = "_")
  ## 所有需要的色块
  u_mod <- unique(dat_d$u_id)
  uu_id <- gsub(".*_","",u_mod)

  ## 有哪些色块是相同的ID 他们的颜色也要一样
  p_ids <- unique(uu_id)
  ## 初始颜色的数量应该和独立ID的数量一致
  palattes <- palette
  names(palattes) <- p_ids

  ## 将颜色放回色块并且命名
  pal <- palattes[uu_id]
  names(pal) <- u_mod
  return(pal)

}

# 构建分组信息 ----------------------------------------------------------------
#' Group circlize matrix
#'
#' @param matrix row names is from, column names is to
#'
#' @return a named character vector
#' @export
#'
#' @examples
#' groupCir(matrix)
groupCir <- function(matrix) {

  nm = unique(unlist(dimnames(matrix)))
  group = structure(gsub("_.*", "", nm), names = nm)
  return(group)

}

# 隐藏多余信息 ------------------------------------------------------------------
hideCir <- function(matrix,preData) {

  count_poor_unique_id_ls <- apply(matrix, 1, function(x) which(x == 2))
  count_poor_unique_id <- gsub("\\..*","",names(unlist(count_poor_unique_id_ls)))
  idx <- which(gsub(".*_","",preData$rn) %in% gsub(".*_","",count_poor_unique_id), arr.ind = TRUE)

  col2<- preData$col
  col2[idx] <- 'transparent'

  return(col2)

}

# 热图注释数据 ------------------------------------------------------------------
#' heatmap values for annotation circlize
#'
#' @param matrix row names is from, column names is to
#' @param df a data frame contains multigroup results of enrich
#' @param look which column as heatmap value
#'
#' @return a list for heatmap in circlize
#' @export
#'
#' @examples
#' heatCir(matrix,df)
heatCir <- function(matrix,df,look = "Count") {

  ## logical for na
  if (any(is.na(df[,get(look)]))) {

    ui_info("{ui_value(look)} contains NA, {ui_value('p.adjust')} will be choosed.")
    fill = "p.adjust"

  }

  if (any(is.na(df[,get(look)]))) {

    ui_info("{ui_value(look)} contains NA, {ui_value('pvalue')} will be choosed.")
    fill = "pvalue"

  }

  # list 的names是from from 下面的子元素是to
  f_col <- apply(matrix, 1, function(x) which(x > 0))
  # list 的names是to to 下面的子元素是from
  t_col <- apply(matrix, 2, function(x) which(x > 0))

  values_f <- heatValues(f_col, df, look)
  values_t <- heatValues(t_col, df, look)

  # 将查询到的值替换到矩阵中
  # from
  mat_f <- matrix
  for(i in seq_along(f_col)) {

    mat_f[names(f_col)[i],f_col[[i]]] <- values_f[[i]]

  }
  # to
  mat_t <- matrix
  for(i in seq_along(t_col)) {

    ## 需要原始位置的对称位置
    mat_t[t_col[[i]], names(t_col)[i]] <- values_t[[i]]

  }

  res_l <- list(
    "mat_from" = mat_f,
    "mat_to" = mat_t
  )
  return(res_l)

}

heatValues <- function(index, df, look) {

  df[,"u_id"] <- paste(df$model,df$ID,sep = "_")
  # from
  values <- lapply(seq_along(index), function(x){

    ## from 是哪个就查哪个
    values = df[which(df$u_id == names(index)[x]),][,get(look)]
    ## to 有几个就重复几次
    values = rep(values, length(index[[x]]))

    return(values)

  })

  return(values)

}

# 画外圈注释热图 -----------------------------------------------------------------
#' Plot rect in continuous value
#'
#' @param preData data after pre plot
#' @param annoData annotation data
#' @param palatte palatte
#' @param track.index index number
#' @param ylim y range
#' @param y1Scale adjust y start
#' @param y2Scale adjust y end
#'
#' @return a rect plot of circlize
#' @export
#'
#' @examples
#' rectCir(preData, annoData, palatte)
rectCir <- function(preData, annoData, palatte, track.index = 1,
                    ylim = c(0,1), y1Scale = 0, y2Scale = 0.5) {

  cnm_o = preData
  heatCount = annoData

  for(i in seq_len(nrow(cnm_o))) {
    mat_f <- heatCount[["mat_from"]]
    mat_t <- heatCount[["mat_to"]]
    abs_max = quantile(abs(range(c(mat_f,mat_t)[c(mat_f,mat_t) != 0])), c(0,0.5,1), na.rm = TRUE)
    col_fun = colorRamp2(abs_max, palatte)
    y1 = ylim[1] + (ylim[2] - ylim[1])*y1Scale
    y2 = ylim[2]*y2Scale
    # from
    circos.rect(cnm_o[i, "x1"], y1, (cnm_o[i, "x1"] - abs(cnm_o[i, "value1"])), y2,
                col = col_fun(mat_f[cnm_o$rn[i], cnm_o$cn[i]]),
                border = col_fun(mat_f[cnm_o$rn[i], cnm_o$cn[i]]),
                sector.index = cnm_o$rn[i], track.index = track.index)
    # # to
    circos.rect(cnm_o[i, "x2"], y1, (cnm_o[i, "x2"] - abs(cnm_o[i, "value1"])), y2,
                col = col_fun(mat_t[cnm_o$rn[i], cnm_o$cn[i]]),
                border = col_fun(mat_t[cnm_o$rn[i], cnm_o$cn[i]]),
                sector.index = cnm_o$cn[i], track.index = track.index)
  }

}

#' plot rect in discrete value
#'
#' @param preData preData
#' @param annoData annoData
#' @param palatte palatte
#' @param track.index track.index
#' @param ylim ylim
#' @param y1Scale y1Scale
#' @param y2Scale y2Scale
#'
#' @importFrom circlize circos.rect
#'
#' @return a rect plot of circlize
#' @export
#'
#' @examples
#' rectCirDiscrete(preData, annoData, palatte)
rectCirDiscrete <- function(preData, annoData, palatte, track.index = 1,
                         ylim = c(0,1), y1Scale = 0, y2Scale = 0.5){

  cnm_o = preData
  heatCount = annoData
  # palatte_f = palattes[gsub("","",)]

  for(i in seq_len(nrow(cnm_o))) {
    mat_f <- heatCount[["mat_from"]]
    mat_t <- heatCount[["mat_to"]]
    y1 = ylim[1] + (ylim[2] - ylim[1])*y1Scale
    y2 = ylim[2]*y2Scale
    circos.rect(cnm_o[i, "x1"], y1, cnm_o[i, "x1"] - abs(cnm_o[i, "value1"]), y2,
                col = palatte[gsub("_.*","",cnm_o$rn[i])], border = palatte[gsub("_.*","",cnm_o$rn[i])],
                sector.index = cnm_o$rn[i], track.index = track.index)
    circos.rect(cnm_o[i, "x2"], y1, cnm_o[i, "x2"] - abs(cnm_o[i, "value1"]), y2,
                col = palatte[gsub("_.*","",cnm_o$cn[i])], border = palatte[gsub("_.*","",cnm_o$cn[i])],
                sector.index = cnm_o$cn[i], track.index = track.index)

  }

}


# 构建图例 --------------------------------------------------------------------

# 不同通路的注释图例
#' Legend for ID
#'
#' @param df df
#' @param palatte p
#' @param type t
#' @param title t
#' @param title_position t
#' @param text_width t
#'
#' @importFrom grid gpar
#'
#' @return a legend of circlize
#' @export
#'
#' @examples
#' idLegend(df, palatte)
idLegend <- function(df, palatte, type = "points",
                     title = 'Pathways', title_position = "topleft",...) {



  p_n <- gsub(".*_","",names(palatte))
  p_n <- unique(p_n)
  palatte = unique(palatte)

  ## 查找通路名字
  dat_u <- unique(df[,c("ID","Description")])
  dat_o <- dat_u[match(p_n, dat_u$ID),]
  ids <- unlist(dat_o[,"Description"])
  Legend(at =ids, type = type,
         grid_width = unit(0.4, 'cm'),grid_height = unit(0.4, 'cm'),
         pch = NA,...,
         background = palatte,
         legend_gp = gpar(fill = palatte,color = palatte), title_position = title_position,
         title = title)

}

# 分组图例
groupLegend <- function(palatte, type = "points", title = "Group", title_position = "topleft") {

  Legend(at = names(palatte), type = type,
         grid_width = unit(0.4, 'cm'),grid_height = unit(0.4, 'cm'),
         pch = NA,
         background = palatte,
         legend_gp = gpar(fill = palatte,color = palatte),
         title_position = title_position,
         title = title)

}

# 热图图例
#' Legend for heatmap
#'
#' @param annoData a
#' @param palatte p
#' @param ... more Legend parameters
#' @param title t
#' @param title_position t
#' @param int integer value?
#'
#' @importFrom stats quantile
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Legend
#'
#' @return a legend of circlize
#' @export
#'
#' @examples
#' heatLegend(annoData, palatte)
heatLegend <- function(annoData, palatte,...,
                       title = "Values", title_position = "topleft", int = T) {

  mat_f <- annoData[["mat_from"]]
  mat_t <- annoData[["mat_to"]]
  abs_max = quantile(abs(range(c(mat_f,mat_t)[c(mat_f,mat_t) != 0])), c(0,0.5,1), na.rm = TRUE)
  col_fun = colorRamp2(abs_max, palatte)

  if (!int) {

    abs_max_label = format(abs_max,scientific = T,digits = 3)
    Legend(at = abs_max, col_fun = col_fun, labels = abs_max_label,
           title_position = title_position, title = title,...)

  } else {

    Legend(at = abs_max, col_fun = col_fun,
           title_position = title_position, title = title,...)

  }

}

# highlight Group ---------------------------------------------------------
#' Group highlight
#'
#' @param df d
#' @param track.index t
#' @param palatte p
#' @param cex c
#' @param text.col t
#' @param niceFacing f
#' @param text.vjust t
#' @param ... more highlight.sector parameters
#'
#' @importFrom circlize highlight.sector
#'
#' @return a highlight of circlize
#' @export
#'
#' @examples
#' highlightGroup(df,track.index = 1,palatte)
highlightGroup <- function(df,track.index = 1,palatte,
                           cex = 1, text.col = "white", niceFacing = TRUE,text.vjust = "4mm",...) {

  df[,"u_id"] <- paste(df$model,df$ID,sep = "_")
  sector_names <- df$u_id
  for (i in unique(df$model)) {

    highlight.sector(sector_names[grep(i,sector_names)], track.index = track.index,
                     border = palatte[i],
                     text = i, cex = cex, text.col = text.col,col = NA,niceFacing = T,
                     text.vjust = text.vjust,...)

  }


}

#' Get circlize data from obj after hyper analysis
#'
#' @param obj obj
#' @param dataBase KEGG or GO
#' @param orderBy pvalue or other column
#' @param head head after order
#'
#' @importFrom usethis ui_oops
#' @importFrom data.table rbindlist
#'
#' @return a list
#' @export
#'
#' @examples
#' modelEnrich(obj)
modelEnrich <- function(obj,dataBase,orderBy = "pvalue", head = 3) {

  if (missing(dataBase)) {

    ui_oops("Please select from {ui_value('KEGG')} or {ui_value('GO')} for ui_code('dataBase')")

  }

  # 验证方法得到的结果
  test <- deg_here(obj)
  ok <- names(test)[which(test == TRUE)]
  main <- setdiff(ok,"merge")

  # 验证上下调分组
  index <- c(setdiff(label(obj),label_ns(obj)),"diff")

  if(length(main) >= 2) { ## 基本条件

    if (dataBase == "KEGG") {
      res_l <- hyperRes(obj = obj)[['hyperKEGG_res']]
    } else if (dataBase == "GO") {
      res_l <- hyperRes(obj = obj)[['hyperGO_res']]
    }

    if(is.null(res_l)) {
      usethis::ui_oops("No available Data in hyper {dataBase}")
      subres_ls <- NULL

    } else {

      subres_ls <- list()
      subres_ls <- lapply(index, function(i){

        ## 提取数据框 完成筛选
        tmp <- lapply(ok, function(x){
          dat <- res_l[[x]][[i]]
          dat_df <- dat@result
          dat_or <- dat_df[order(dat_df[,orderBy]),]
          dat_top <- head(dat_or,head)
          dat_top[,"model"] <- x
          return(dat_top)
        })

        ## 合并每个上下调中的多种方法的数据
        tmp <- data.table::rbindlist(tmp)
        return(tmp)

      })
      names(subres_ls) <- index

    }
  } else {

    ui_oops("NO data avaliable for {ui_code('modelEnrich')}")
    subres_ls <- NULL

  }

  return(subres_ls)

}
