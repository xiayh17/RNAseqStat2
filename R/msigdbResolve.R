#' @importFrom kit fduplicated
msigdbGet <- function(object) {

  ## 下载数据
  msigdbParam = msigdbParam(object)
  gene_sets = msigdbGetCore(msigdbParam = msigdbParam)

  ## 整理数据
  ## 根据设置，得到要获取的子集
  cat_param <- msigdbParam(object)[["category"]]
  if (is.null(cat_param)) {

    gs_cat = c("H",paste0("C",1:8))

  } else {

    gs_cat = cat_param

  }

  ## 分割成子集,并做精简
  MSigDB_l <- lapply(gs_cat, function(x){

    dat <- gene_sets[gene_sets$gs_cat==x,]
    ## 只留 gs_name 和 gene_symbol 并去重
    dat <- dat[,c("gs_name","gene_symbol")]
    # dat <- unique(dat)
    dat <- as.data.frame(dat[!fduplicated(dat),])

    return(dat)

  })
  ## 给子集命名
  names(MSigDB_l) <- gs_cat

  ## 前缀去除,只有 H 好换
  if ("H" %in% names(MSigDB_l)) {
    tmp <- MSigDB_l[["H"]]
    tmp$gs_name <- gsub("HALLMARK_","",tmp$gs_name)
    MSigDB_l[["H"]] <- tmp
  }

  ## 存储到对象
  msigdbData(object) <- MSigDB_l

  return(object)

}

#' @importFrom msigdbr msigdbr
msigdbGetCore <- function(...,msigdbParam){

  params <- list(...)
  msigdbParam <- modifyList(params, msigdbParam)
  msigdb_core <- suppressMessages(do.call("msigdbr", modifyList(
    list(),
    msigdbParam)
  ))

  return(msigdb_core)

}

gseMSigDB <- function(object) {

  msigdbGSEAparam <- msigdbGSEAparam(object)

  ## 获取GeneList
  test <- deg_here(object)
  ok <- names(test)[which(test == TRUE)] ## 取有效数据

  gseKEGG_GeneSets = list()
  for (i in ok) {
    gseKEGG_GeneSets[[i]] <- suppressWarnings(GSEA_GS(object,which = i,type = "SYMBOL"))
  }

  ## 对每个子集进行富集分析, 每个子集又可能存在多个geneList 做富集分析
  t2g_l <- msigdbData(object)
  res_l <- lapply(gseKEGG_GeneSets, function(geneList){

    gse_res_l <- lapply(seq_along(t2g_l), function(x){

     gseMSigDB_Core(geneList = geneList,TERM2GENE = t2g_l[[x]],msigdbGSEAparam =msigdbGSEAparam)

    })
    names(gse_res_l) <- names(t2g_l)
    return(gse_res_l)

  })
  names(res_l) <- names(gseKEGG_GeneSets)
  ## 结果存储
  msigdbGSEAresult(object) <- res_l
  return(object)

}

#' @importFrom clusterProfiler GSEA
gseMSigDB_Core <- function(...,msigdbGSEAparam){

  params <- list(...)
  msigdbGSEAparam <- modifyList(params, msigdbGSEAparam)
  msigdb_core <- suppressMessages(do.call("GSEA", modifyList(
    list(),
    msigdbGSEAparam)
  ))

  return(msigdb_core)

}

#' @importFrom GSEABase GeneSet GeneSetCollection EntrezIdentifier KEGGCollection
#' @importFrom GSVA gsva
gsvaResolve <- function(object) {

  # 重新提取数据
  msi_data_list <- msigdbData(object)

  # 将数据转换格式
  ids_list <- lapply(msi_data_list, function(x){

    res <- split(x[,"gene_symbol"],x[,"gs_name"])
    return(res)

  })
  names(ids_list) <- names(msi_data_list)

  # 构建 GeneSetCollection
  GSC_list <- lapply(ids_list, function(x){

    GeneSetCollection(mapply(function(geneIds, keggId) {
      GeneSet(geneIds, geneIdType=EntrezIdentifier(),
              collectionType=KEGGCollection(keggId),
              setName=keggId)
    }, x,  names(x)))

  })
  ## 运行gsva
  gsva_matrix <- as.matrix(matrixFiltered(object))

  gsvaRes_list <- lapply(GSC_list, function(geneset){

    es.max <- gsva(gsva_matrix, geneset,
                    mx.diff=FALSE,
                    verbose=FALSE,
                    kcdf="Poisson",
                    parallel.sz= parallel::detectCores()/2)

    return(es.max)

  })

  gsvaDiff_list <- lapply(gsvaRes_list, function(es.max){

    gsva_limma_resolve(es.max,group_list = groupInfo(object),case_group = caseGroup(object))

  })

  ## 存储结果
  msigdbGSVAresult(object) <- list(
    "GSVA_matrix" = gsvaRes_list,
    "GSVA_diff" = gsvaDiff_list
  )

  return(object)

}

# 对GSVA 差异结果
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
gsva_limma_resolve <- function(gsva_data, group_list, case_group) {
  control_group = setdiff(group_list,case_group)
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(gsva_data)

  # dge <- DGEList(counts=gsva_data)
  # dge <- calcNormFactors(dge)
  # logCPM <- cpm(dge, log=TRUE, prior.count=3)

  # v <- voom(dge,design,plot=TRUE, normalize.method="quantile")
  fit <- lmFit(gsva_data, design)

  con=paste0(case_group,'-',control_group)

  cont.matrix=makeContrasts(contrasts=c(con),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)

  tempOutput = topTable(fit2, ,coef=con,adjust='BH', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)

  return(DEG_limma_voom)

}
