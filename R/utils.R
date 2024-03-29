## res is one of result from limma or DESeq2 or edgeR
## cutFC is from treatInfo
## cutFDR is from treatInfo

#' Identify the column where the FC is located
#'
#' @param res one of result from limma or DESeq2 or edgeR
#'
#' @return a column name
#' @export
#'
#' @examples
#' FC_Identify(deg_data)
FC_Identify <- function(res) {
  x = NULL
  deg_data = res
  ## 自动检测 FC 所在的列 ----
  fc_pattern <- "logFC|logFC|log2FoldChange"
  x = colnames(deg_data)[grep(fc_pattern,colnames(deg_data))]
  return(x)
}

#' Identify the column where the pvalue is located
#'
#' @param res one of result from limma or DESeq2 or edgeR
#'
#' @return a column name
#' @export
#'
#' @examples
#' pvalue_Identify(deg_data)
pvalue_Identify <- function(res) {
  y = NULL
  deg_data = res
  ## 自动检测 Pvalue 所在的列 ----
  fdr_pattern <- "P.Value|PValue|pvalue"
  y = colnames(deg_data)[grep(fdr_pattern,colnames(deg_data))]
  return(y)
}

#' Verify the rationality of the FC threshold
#'
#' @param res one of result from limma or DESeq2 or edgeR
#' @param cut_FC if is null, will give a recommend value.
#' @param scale adjust recommend value
#'
#' @return a value
#' @export
#'
#' @examples
#' cutFC_Verify(deg_data)
cutFC_Verify <- function(res, cut_FC, scale = 1) {

  deg_data = res
  x <- FC_Identify(res = deg_data)

  ## FC阈值 ----
  if (is.null(cut_FC)) { ## 如果是自动
    cut_FC = (mean(abs(deg_data[,x])) + 2*sd(abs(deg_data[,x]))) * scale
    ui_info(glue("{ui_code('|cutFC|')} = {ui_value(cut_FC)} calculated automatically by {ui_code(glue('(mean(abs({x})) + 2*sd(abs({x})))*{scale}'))}
                 Your can also set by {ui_code('cutFC')} argument."))
  } else if (is.numeric(cut_FC)) { ## 如果是数字

    if (max(cut_FC) > max(deg_data[,x])) { ## 如果数字不合理

      usethis::ui_oops(glue("{ui_code('|cutFC|')} = {ui_value(cut_FC)} specified by you is larger than Maximum({max(deg_data[,x])}) of (log)FC column"))
      cut_FC = (mean(abs(deg_data[,x])) + 2*sd(abs(deg_data[,x]))) * scale
      ui_info(glue("{ui_code('|cutFC|')} = {ui_value(cut_FC)} calculated automatically by {ui_code(glue('(mean(abs({x})) + 2*sd(abs({x})))*{scale}'))}"))

    } else {

      ui_info(glue("{ui_code('|cutFC|')} = {ui_value(cut_FC)} specified by you seems ok."))

    }


  } else { ## 如果是其他情况

    usethis::ui_oops(glue("{ui_code('|cutFC|')} = {ui_value(cut_FC)} specified by you will not applied"))
    cut_FC = (mean(abs(deg_data[,x])) + 2*sd(abs(deg_data[,x]))) * scale
    ui_info(glue("{ui_code('|cutFC|')} = {ui_value(cut_FC)} calculated automatically by {ui_code(glue('(mean(abs({x})) + 2*sd(abs({x})))*{scale}'))}"))

  }

  return(cut_FC)

}

upDate_cutFC_NULL <- function(object,which,new_cutFC) {

  if (is.null(cutFC(object))) {

    cutFC = list()
    cutFC[[which]] = new_cutFC
    cutFC(object) = cutFC

  } else if (is.list(cutFC(object))) {

    cutFC = cutFC(object)
    cutFC[[which]] = new_cutFC
    cutFC(object) = cutFC

  } else if (length(cutFC(object)) == 1&is.numeric(cutFC(object))) {

    cutFC = list()
    cutFC[[which]] = new_cutFC
    cutFC(object) = cutFC

  }

  ui_info("cutFC updated for object!")

  return(object)

}

#' get deg data of DEGContainer
#'
#' @param obj a DEGContainer object
#' @param which kinds of DEG; can be "limma", "edgeR", "DESeq2" or "MSigDB"
#' @param category MSigDB collection abbreviation, such as H or C1.
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' dataDEG(DEGContainer)
dataDEG <- function(obj,which,category = "H") {

  if (which == "limma") {
    deg_data = limma_res(obj)
  } else if (which == "edgeR") {
    deg_data = edgeR_res(obj)
  } else if (which == "DESeq2") {
    deg_data = DESeq2_res(obj)
  } else if (which == "merge") {
    deg_data = merge_res(obj)
    deg_data <- commonGroup(merge_data = deg_data)
  } else if (which == "MSigDB") {
    if(!is.null(category)) {

      deg_data = msigdbGSVAresult(obj)[["GSVA_diff"]][[category]]

    } else {

      ui_stop("when {ui_code('which')} set as {ui_value('MSigDB')}, {ui_code('category')} should be one of category in {ui_value('MSigDB')}")

    }
  } else {
    ui_stop("{ui_code('which')} should be one of {ui_value('limma, edgeR, DESeq2, merge, MSigDB')}")
  }

  return(deg_data)

}

#' Group result of DEG analysis results
#'
#' @param object a DEGContainer
#' @param which kinds of DEG; can be "limma", "edgeR", "DESeq2" or "MSigDB"
#' @param category MSigDB collection abbreviation, such as H or C1.
#'
#' @return a DEGContainer
#' @export
#'
#' @examples
#' cutMuch(DEGContainer,"limma","H")
cutMuch <- function(object,which,category = 'H') {

  deg_data <- dataDEG(obj = object,which = which,category = category)

  x = FC_Identify(deg_data)
  y = pvalue_Identify(deg_data)

  if (which == "MSigDB") {

    treat = msigdbTreat(object)
    label = treat@label
    label_ns = treat@label_ns
    cut_FDR = treat@cutFDR
    cut_FC = treat@cutFC

  } else {

    label = label(object)
    label_ns = label_ns(object)
    cut_FDR = cutFDR(object)
    cut_FC = cutFC(object)

  }

  if (is.null(cut_FC)) {
    if (which == 'MSigDB') {
      cut_FC <- cutFC_Verify(res = deg_data,cut_FC,scale = 0.5)
    } else {
      cut_FC <- cutFC_Verify(res = deg_data,cut_FC)
    }
    object <- upDate_cutFC_NULL(object = object,which = which,new_cutFC = cut_FC)
    cut_FC <- c(-cut_FC, cut_FC)
  } else {

    if (length(cut_FC) == 1&is.numeric(cut_FC)) {

      if (which == 'MSigDB') {
        cut_FC <- cutFC_Verify(res = deg_data,cut_FC,scale = 0.5)
      } else {
        cut_FC <- cutFC_Verify(res = deg_data,cut_FC)
      }
      if (cut_FC != cutFC(object)) {
        object <- upDate_cutFC_NULL(object = object,which = which,new_cutFC = cut_FC)
      }
      cut_FC <- c(-cut_FC, cut_FC)

    } else if(is.list(cut_FC)) {

      cut_FC <- cut_FC[[which]]
      if (is.null(cut_FC)) {
        if (which == 'MSigDB') {
          cut_FC <- cutFC_Verify(res = deg_data,cut_FC,scale = 0.5)
        } else {
          cut_FC <- cutFC_Verify(res = deg_data,cut_FC)
        }
        object <- upDate_cutFC_NULL(object = object,which = which,new_cutFC = cut_FC)
      } else {

        if (which == 'MSigDB') {
          cut_FC_new <- cutFC_Verify(res = deg_data,cut_FC,scale = 0.5)
        } else {
          cut_FC_new <- cutFC_Verify(res = deg_data,cut_FC)
        }
        if (cut_FC != cut_FC_new) {
          object <- upDate_cutFC_NULL(object = object,which = which,new_cutFC = cut_FC_new)
          cut_FC = cut_FC_new
        }
      }



      if (length(cut_FC) == 1&is.numeric(cut_FC)) {
        cut_FC <- c(-cut_FC, cut_FC)
      }
    }

  }

  if (length(cut_FDR) == 1&is.numeric(cut_FDR)) {
    cut_FDR <- rep(cut_FDR,length(cut_FC))
  } else {
    ui_stop("{ui_code('cut_FDR')} should be a single number.")
  }

  ## group data
  label_cg <- setdiff(label,label_ns)
  names(cut_FDR) <- label_cg
  deg_data$group <- cut(deg_data[,x],
                        breaks = c(-Inf,cut_FC,Inf),
                        labels = label)
  index = list()
  for (i in label_cg) {
    index[[i]] <- setdiff(which(deg_data$group == i), which(deg_data$group == i & deg_data[,y] < cut_FDR[i]))
    deg_data$group[index[[i]]] <- label_ns
  }

  group_count <- table(deg_data$group)
  group_count <- as.character(group_count)
  # ## 报告分组情况
  usethis::ui_info(glue("The threshold of the FC is {ui_value(unique(abs(cut_FC)))} for {ui_value(which)}"))
  usethis::ui_info(glue("The threshold of the FDR is {ui_value(unique(abs(cut_FDR)))} for {ui_value(which)}"))
  usethis::ui_line(glue('{paste0(label," ","{ui_value(",group_count,")}",collapse = "; ")}'))
  usethis::ui_done(glue("Group done of {ui_value(which)}!"))

  # 存到对象
  if (which == "limma") {
    limma_res(object) = deg_data
  } else if (which == "edgeR") {
    edgeR_res(object) = deg_data
  } else if (which == "DESeq2") {
    DESeq2_res(object) = deg_data
  } else if (which == "MSigDB") {
    msigdbGSVAresult(object)[["GSVA_diff"]][[category]] = deg_data
  }

  return(object)

}

#' Choose top genes from DEG data frame of limma DESeq2 edgeR
#'
#' For example, if topSig = 50, head 50 and tail 50 will return from a
#' ordered DEG data frame after filtering by threshold value of P value
#'
#' @param object a DEG data frame contains logFC and p value
#' @param topSig a single number or a length of 2 numeric vector, if 2 numeric vector, first one is top max logFC.
#' @param which a single number for threshold value of P value
#'
#' @importFrom utils head tail
#' @importFrom data.table data.table .SD
#'
#' @return a character vector of top genes
#' @export
#'
#' @examples
#' topGene(object,  topSig = 50, which = "limma")
topGene <- function(object, topSig = 50, which, category = 'H') {

  deg_data <- dataDEG(obj = object,which = which,category = category)

  if (which == "merge") {

    x = FC_Identify(deg_data)
    y = pvalue_Identify(deg_data)
    cut_FDR = cutFDR(object)
    top = topSig

    data_f <- deg_data[which(apply((deg_data[,y] <= cut_FDR), 1, all)),]

    gene_named <- data.frame(
      SYMBOL = row.names(data_f),
      value = data_f[,x]
    )

    gene_named$value <- apply(gene_named[,setdiff(colnames(gene_named),c("SYMBOL","ENTREZID"))],1,mean)

    d <- data.table(gene_named, key='value', keep.rownames = TRUE)

    if (length(top) == 1) {
      td <- d[, head(.SD, top)] ## most mini
      hd <-  d[, tail(.SD, top)] ## most large
    } else if (length(top) == 2) {
      td <- d[, head(.SD, top[2])] # most mini
      hd <-  d[, tail(.SD, top[1])] # the descending order
    } else {
      stop("top should be a single number or a length of 2 numeric vector")
    }

  } else if (which == "MSigDB") {

    x = FC_Identify(deg_data)
    y = pvalue_Identify(deg_data)
    treat = msigdbTreat(object)
    cut_FDR = cut_FDR = treat@cutFDR
    top = topSig

    data_f <- deg_data[which(deg_data[,y] <= cut_FDR),]

    d <- data.table(data_f, key=x, keep.rownames = TRUE)

    if (length(top) == 1) {
      td <- d[, head(.SD, top)]
      hd <-  d[, tail(.SD, top)]
    } else if (length(top) == 2) {
      td <- d[, head(.SD, top[2])]
      hd <-  d[, tail(.SD, top[1])] # the descending order
    } else {
      stop("top should be a single number or a length of 2 numeric vector")
    }

  } else {

    x = FC_Identify(deg_data)
    y = pvalue_Identify(deg_data)
    cut_FDR = cutFDR(object)
    top = topSig

    data_f <- deg_data[which(deg_data[,y] <= cut_FDR),]

    d <- data.table(data_f, key=x, keep.rownames = TRUE)

    if (length(top) == 1) {
      td <- d[, head(.SD, top)]
      hd <-  d[, tail(.SD, top)]
    } else if (length(top) == 2) {
      td <- d[, head(.SD, top[2])]
      hd <-  d[, tail(.SD, top[1])] # the descending order
    } else {
      stop("top should be a single number or a length of 2 numeric vector")
    }

  }



  hd_names <- hd$rn
  names(hd_names) <- rep("head",length(hd_names))
  td_names <- td$rn
  names(td_names) <- rep("tail",length(hd_names))

  choose_names <- c(hd_names, td_names)

  return(choose_names)
}


#' Test for exists of DEG results
#'
#' @param object a DEGContainer
#'
#' @return logic vector
#' @export
#'
#' @examples
#' deg_here(DEGContainer)
deg_here <- function(object) {

  obj = object
  test <- list()
  if (!is.null(limma_res(obj))) {
    if (nrow(limma_res(obj))> 3) {
      test[["limma"]] <- TRUE
    } else {
      test[["limma"]] <- FALSE
    }
  } else {
    test[["limma"]] <- FALSE
  }

  if (!is.null(edgeR_res(obj))) {
    if (nrow(edgeR_res(obj))> 3) {
      test[["edgeR"]] <- TRUE
    } else {
      test[["edgeR"]] <- FALSE
    }
  } else {
    test[["edgeR"]] <- FALSE
  }

  if (!is.null(DESeq2_res(obj))) {
    if (nrow(DESeq2_res(obj))> 3) {
      test[["DESeq2"]] <- TRUE
    } else {
      test[["DESeq2"]] <- FALSE
    }
  } else {
    test[["DESeq2"]] <- FALSE
  }

  if (!is.null(merge_res(obj))) {
    if (nrow(merge_res(obj))> 3) {
      test[["merge"]] <- TRUE
    } else {
      test[["merge"]] <- FALSE
    }
  } else {
    test[["merge"]] <- FALSE
  }

  test = unlist(test)
  return(test)

}

merge_deg <- function(object) {

    deg_df_limma = limma_res(object)

    deg_df_edgeR = edgeR_res(object)

    deg_df_DESeq2 = DESeq2_res(object)

    data_list_o <- list(deg_df_limma,deg_df_edgeR,deg_df_DESeq2)

    names(data_list_o) <- c("limma","edgeR","DESeq2")

    data_list <- list()
    for (x in c("limma","edgeR","DESeq2")) {
      if(is.null(data_list_o[[x]])){

      } else {
        data_list[[x]] <- data_list_o[[x]]
        names(data_list[x]) <- x
      }
    }

    rm(ls="x")

    rownames_list <- lapply(data_list, rownames)

  allg <- Reduce(intersect, rownames_list)

  limma_col <- c("logFC","P.Value","adj.P.Val","group")
  edgeR_col <- c("logFC","PValue","FDR","group")
  DESeq2_col <- c("log2FoldChange","pvalue","padj","group")

  for(x in names(data_list)) {

    if(x == "limma") {
      get_col = intersect(colnames(data_list[[x]]),limma_col)
      data_list[[x]] = data_list[[x]][allg,get_col]
      # colnames(data_list[[x]]) <- paste0(colnames(data_list[[x]]),"_",x)
    }

    if(x == "edgeR") {
      get_col = intersect(colnames(data_list[[x]]),edgeR_col)
      data_list[[x]] = data_list[[x]][allg,get_col]
      # colnames(data_list[[x]]) <- paste0(colnames(data_list[[x]]),"_",x)
    }

    if(x == "DESeq2") {
      get_col = intersect(colnames(data_list[[x]]),DESeq2_col)
      data_list[[x]] = data_list[[x]][allg,get_col]
      # colnames(data_list[[x]]) <- paste0(colnames(data_list[[x]]),"_",x)
    }


  }
  rm(ls="x")

  deg_df_intersect <- do.call(cbind,data_list)

  return(deg_df_intersect)

}

#' Title
#'
#' @param merge_data
#'
#' @return
#' @export
#'
#' @examples
#' commonGroup(merge_data)
commonGroup <- function(merge_data) {
  Test <- merge_data[,grepl("group",colnames(merge_data))]

  com_index <- apply(Test, 1, function(x){

    length(unique(x)) == 1

  })

  merge_data2 <- merge_data[com_index,]
  merge_data3 <- merge_data2[,!grepl("group",colnames(merge_data2))]
  merge_data3$group <- merge_data2[,grep("group",colnames(merge_data2))[1]]

  return(merge_data3)
}

#' Get a significant difference gene list
#'
#' @param object a DEGContainer
#' @param which kinds of DEG; can be "limma", "edgeR", "DESeq2" or "MSigDB"
#' @param type "ENTREZID" or "SYMBOL"
#' @param OrgDb OrgDb
#' @param category MSigDB collection abbreviation, such as H or C1.
#'
#' @return
#' @export
#'
#' @examples
#' hyper_GS(object, which, type)
hyper_GS <- function(object,which,type,OrgDb = NULL,category = "H") {

  deg_data <- dataDEG(obj = object,which = which,category = category)

  if(is.null(OrgDb)){

    OrgDb = switch (species(object),
                    'Human' = "org.Hs.eg.db",
                    'Mouse' = "org.Mm.eg.db",
                    'Rat' = "org.Rn.eg.db"
    )

  }

  g <- setdiff(label(object),label_ns(object))
  deg_df_g <- deg_data
  gene_list <- list()

  if (!("group" %in% colnames(deg_df_g))) {

    usethis::ui_stop("Try to group data with{ui_code(degGroup)}")

  }

  if (type == "ENTREZID") {

    for (i in g) {
      SYMBOLS_id = row.names(deg_df_g[which(deg_df_g$group == i),])
      ENTREZ_id = suppressMessages(clusterProfiler::bitr(SYMBOLS_id,
                                                         fromType = "SYMBOL",
                                                         toType = "ENTREZID",
                                                         OrgDb = OrgDb))
      gene_list[[i]] = ENTREZ_id[,"ENTREZID"]
    }
    rm(list = "i")

  } else if (type == "SYMBOL") {

    for (i in g) {
      gene_list[[i]] = row.names(deg_df_g[which(deg_df_g$group == i),])
    }
    rm(list = "i")

  } else {

    ui_stop("{ui_code(type)} should be one of {ui_value('SYMBOL')} OR {ui_value('ENTREZID')}")

  }

  gene_list[["diff"]] <- unique(unlist(gene_list))

  return(gene_list)

}

GSEA_GS <- function(object,which,type,OrgDb = NULL) {

  deg_data <- dataDEG(obj = object,which = which,category = category)

  x = FC_Identify(deg_data)

  if(is.null(OrgDb)){

    OrgDb = switch (species(object),
                  'Human' = "org.Hs.eg.db",
                  'Mouse' = "org.Mm.eg.db",
                  'Rat' = "org.Rn.eg.db"
  )

  }

  if (which == "merge") {

    if(type == "ENTREZID") {

      gene_data <- data.frame(
        SYMBOL = row.names(deg_data),
        value = deg_data[,x]
      )

      ENTREZ_id = suppressMessages(clusterProfiler::bitr(gene_data$SYMBOL,
                                                         fromType = "SYMBOL",
                                                         toType = "ENTREZID",
                                                         OrgDb = OrgDb))

      gene_named <- merge(gene_data,ENTREZ_id)
      gene_named$value <- apply(gene_named[,setdiff(colnames(gene_named),c("SYMBOL","ENTREZID"))],1,mean)
      geneList <- gene_named$value
      names(geneList) <- gene_named$ENTREZID

    } else if (type == "SYMBOL") {

      gene_data <- data.frame(
        SYMBOL = row.names(deg_data),
        value = deg_data[,x]
      )
      gene_data$value <- apply(gene_data[,setdiff(colnames(gene_data),c("SYMBOL","ENTREZID"))],1,mean)
      geneList = gene_data$value
      names(geneList) = gene_data$SYMBOL

    } else {

      ui_stop("{ui_code(type)} should be one of {ui_value('SYMBOL')} OR {ui_value('ENTREZID')}")

    }



  } else {

    if(type == "ENTREZID") {
      gene_data <- data.frame(
        SYMBOL = row.names(deg_data),
        value = deg_data[,x]
      )

      ENTREZ_id = suppressMessages(clusterProfiler::bitr(gene_data$SYMBOL,
                                                         fromType = "SYMBOL",
                                                         toType = "ENTREZID",
                                                         OrgDb = OrgDb))

      gene_named <- merge(gene_data,ENTREZ_id)
      geneList <- gene_named$value
      names(geneList) <- gene_named$ENTREZID

    } else if (type == "SYMBOL") {

      geneList = deg_data[,x]
      names(geneList) = rownames(deg_data)

    }

  }

  geneList=sort(geneList,decreasing = TRUE)
  return(geneList)

}


#' Convert a ENSEMBL named data to SYMBOL named
#'
#' larger median will be retain for duplicated
#'
#' @param row_counts count data frame
#' @param species support Human Mouse Rat
#'
#' @return a renamed dataframe
#' @export
#'
#' @examples
#' toSYMBOL(count,b2g)
toSYMBOL <- function(row_counts,species) {

  row_counts <- as.data.frame(row_counts)
  row_counts[,"ENSEMBL"] <- rownames(row_counts)

  symbols <- AnnoProbe::annoGene(IDs = rownames(row_counts),ID_type = "ENSEMBL",species = tolower(species))
  symbols <- symbols[,c("SYMBOL","ENSEMBL")]

  counts_ids <- merge(symbols,row_counts,by="ENSEMBL")

  data = counts_ids

  # symbols$SYMBOL[grep(TRUE,duplicated(symbols$SYMBOL))]

  if (anyDuplicated(data[,"SYMBOL"]) == 0) {
    data <- tibble::column_to_rownames(data,"SYMBOL")
  } else if (anyDuplicated(data[,"SYMBOL"]) != 0) {
    data$median=apply(data[,2:ncol(data)],1,median)
    data <- data[order(data[,"SYMBOL"],data$median,decreasing = T),]
    data <- data[!duplicated(data[,"SYMBOL"]),][,-ncol(data)]
  }

  rownames(data) <- NULL
  data[,"ENSEMBL"] <-NULL
  data <- tibble::column_to_rownames(data,"SYMBOL")

  return(data)

}

#' Convert a  data to SYMBOL named
#'
#' larger median will be retain for duplicated
#'
#' @param row_counts count data frame
#' @param bad2good a dataframe, first column should be row name of row_counts, second should be SYMBOL
#'
#' @return a renamed dataframe
#' @export
#'
#' @examples
#' cleanSYMBOL(count,b2g)
cleanSYMBOL <- function(row_counts,bad2good) {

  row_counts <- as.data.frame(row_counts)
  row_counts[,"badSYMBOL"] <- rownames(row_counts)

  symbols <- bad2good
  colnames(symbols) <- c("badSYMBOL","SYMBOL")
  data <- merge(symbols,row_counts,by="badSYMBOL")

  # symbols$SYMBOL[grep(TRUE,duplicated(symbols$SYMBOL))]

  if (anyDuplicated(data[,"SYMBOL"]) == 0) {
    data <- tibble::column_to_rownames(data,"SYMBOL")
  } else if (anyDuplicated(data[,"SYMBOL"]) != 0) {
    data$median=apply(data[,2:ncol(data)],1,median)
    data <- data[order(data[,"SYMBOL"],data$median,decreasing = T),]
    data <- data[!duplicated(data[,"SYMBOL"]),][,-ncol(data)]
  }

  rownames(data) <- NULL
  data[,"badSYMBOL"] <-NULL
  data <- tibble::column_to_rownames(data,"SYMBOL")

  return(data)

}


#' Convert EGGnog to GENE2TERM
#'
#'
#'
#' @param data data from eggnog
#' @param SYMBOL which column is SYMBOL
#' @param TERM which column is Pathway
#'
#' @return a data frame
#' @export
#'
#' @examples
#' eggnogG2T(data,SYMBOL,TERM)
eggnogG2T <- function(data,SYMBOL,TERM) {

  list_res <- lapply(seq_along(data[,get(SYMBOL)]), function(x){

    gos <- strsplit(data[,get(TERM)][x],",",perl = T)
    query <- rep(data[,get(SYMBOL)][x],length(gos[[1]]))

    df <- data.frame(
      Term = gos[[1]],
      SYMBOL = query
    )
    return(df)

  })

  res_df <- data.table::rbindlist(list_res)
  res_df <- as.data.frame(res_df)

  return(res_df)

}
