#' Create \code{DEGContainer}
#'
#' Begin in DEGContainer.
#'
#' @param species species for your data.
#' @param dataType kind of expresses value matrix. `Counts` (Integer) or `Array` (Decimal).
#' @param idType kind of gene id. `ENSEMBL` or `SYMBOL`, If `ENSEMBL`,\code{species} supposed be one of (`Human`, `Mouse` and `Rat`), it will be automatically converted to `SYMBOL`.
#' @param expMatrix expresses value matrix. Should be a data.frame row named by gene ID and column named by Sample
#' @param groupInfo a Character Vectors ordered by samples in matrix.
#' @param caseGroup a Character names of case group.
#' @param filterMethod a function used to filter expresses value matrix. Or disable filter by set as `NULL`.
#' @param OrgDb Select an OrgDb.
#' @param organism Select an organism.
#' @param GOTERM2GENE custom enrich data of GO.
#' @param GOTERM2NAME custom enrich data of GO.
#' @param KEGGTERM2GENE custom enrich data of GSE.
#' @param KEGGTERM2NAME custom enrich data of GSE.
#' @param msi_species one of \code{msigdbr::msigdbr_species}
#'
#' @details `Human`, `Mouse` and `Rat` is full supported. For other species, more setting is needed.
#'     If data can't provided, related steps in workflow will ignore.
#'     OrgDb can be find from here \link{https://bioconductor.org/packages/release/BiocViews.html#___OrgDb}
#'     organism can be find from here \link{https://www.genome.jp/kegg/catalog/org_list.html}
#'     msi_species can be get from \code{msigdbr::msigdbr_species}
#'
#' @importFrom usethis ui_info ui_done ui_code ui_oops ui_stop
#'
#' @return a \code{DEGContainer} object
#' @export
#'
#' @examples
#' Create_DEGContainer(expMatrix = counts_input,groupInfo = c("C","C","C","T","T","T"),caseGroup = "T")
Create_DEGContainer <- function(species = "Human",
                                dataType = "Counts",
                                idType = "SYMBOL",
                                expMatrix,
                                groupInfo,
                                caseGroup,
                                filterMethod = "rowSums(expMatrix > 0) >= ncol(expMatrix)/2",
                                OrgDb = NULL, msi_species = NULL, organism = NULL,
                                GOTERM2GENE = NULL, GOTERM2NAME = NA,
                                KEGGTERM2GENE = NULL, KEGGTERM2NAME = NA){

  expMatrix <- as.data.frame(expMatrix)

  ## convert ids
  if (idType == "ENSEMBL"&presetSpecies(species)) {

    expMatrix <- suppressWarnings(toSYMBOL(row_counts = expMatrix,species = species))
    idType = "SYMBOL"
    ui_done("ENSEMBL named expression matrix renamed by SYMBOL")
  } else if (idType != "SYMBOL"&!presetSpecies(species)) {
    ui_stop("{ui_code('idType')} should be 'SYMBOL'")
  }

  ## check filter
  if(is.null(filterMethod)){

    ui_info("Nothing filtered in your data.")

  } else {

    keep_feature <- eval(parse(text = filterMethod))
    matrixFiltered <- as.data.frame(expMatrix[keep_feature,])
    left = nrow(matrixFiltered)

    if (left == 0) {

      matrixFiltered = NULL
      ui_stop("Your {ui_code('filterMethod')} have some problem OR too Strict.
              It Should be something like {ui_code('rowSums(expMatrix > 0) >= ncol(expMatrix)/2')} in a character.")

    } else {

      ui_done("Your data have filtered by {ui_code(filterMethod)}")
      ui_info("{left} Keeped, {nrow(expMatrix) - left} removed. ")

    }
  }

  ## Create basic data info
  dataInfo <- Create_dataInfo(
    species = species,
    dataType = dataType,
    idType = idType,
    expMatrix = expMatrix,
    groupInfo = groupInfo,
    caseGroup = caseGroup,
    filterMethod = filterMethod,
    matrixFiltered = matrixFiltered
  )

  ## Create zero for DEG
  degResults <- suppressWarnings(Create_degResults(
    vsData = Create_vsData(),
    treatInfo = Create_treatInfo()
  ))

  ## Create zero for enrich
  if(presetSpecies(species)&is.null(GOTERM2GENE)&is.null(KEGGTERM2GENE)) {

    OrgDb = switch (species,
                    'Human' = "org.Hs.eg.db",
                    'Mouse' = "org.Mm.eg.db",
                    "Rat" = "org.Rn.eg.db"
    )

    organism = switch (species,
                       'Human' = "hsa",
                       'Mouse' = "mmu",
                       'Rat' = "rno"
    )

    msi_species <- switch (species,
                           'Human' = "Homo sapiens",
                           'Mouse' = "Mus musculus",
                           'Rat' = "Rattus norvegicus"
    )

    hyperParam = Create_hyperParam(goParam = list(OrgDb = OrgDb),
                      keggParam = list(organism = organism))
    gseParam = Create_gseParam(goParam = list(OrgDb = OrgDb),
                            keggParam = list(organism = organism))

    msigdbParam = Create_msigdbParam(msigdbParam = list(species = msi_species,category = "H"))

    shotSpecies(species)

  } else {

    if(!is.null(GOTERM2GENE)) {
      GOcustom = TRUE
      goParam = list(TERM2GENE = GOTERM2GENE, TERM2NAME = GOTERM2NAME)
    } else if (!is.null(OrgDb)) {
      GOcustom = FALSE
      goParam = list(OrgDb = OrgDb)
    } else {
      GOcustom = FALSE
      goParam = NULL
      ui_oops("hyperGO and gseGO step will skip for lack data")
    }

    if(!is.null(KEGGTERM2GENE)) {
      KEGGcustom = TRUE
      keggParam = list(TERM2GENE = KEGGTERM2GENE, TERM2NAME = KEGGTERM2NAME)
    } else if (!is.null(organism)) {
      KEGGcustom = FALSE
      keggParam = list(organism = organism)
    } else {
      KEGGcustom = FALSE
      keggParam = NULL
      ui_oops("hyperKEGG and gseKEGG step will skip for lack data")
    }

    if(!is.null(msi_species)) {
      msigdbParam = Create_msigdbParam(msigdbParam = list(species = msi_species,category = "H"))
    } else {
      msigdbParam = Create_msigdbParam(msigdbParam = list(species = msi_species,category = "H"))
      ui_oops("MSigDB step will skip for lack data")
    }

    if(all(is.null(GOTERM2GENE),is.null(OrgDb),is.null(KEGGTERM2GENE),is.null(organism),is.null(msi_species))){

      shotSpecies(species)

    }

    hyperParam = Create_hyperParam(goParam = goParam,keggParam = keggParam,
               customGO = GOcustom,customKEGG = KEGGcustom)

    gseParam = Create_gseParam(goParam = goParam,keggParam = keggParam,
                               customGO = GOcustom,customKEGG = KEGGcustom)

    msigdbParam = Create_msigdbParam(msigdbParam = list(species = msi_species,category = "H"))

  }

  hyperResults <- suppressWarnings(Create_hyperResults(hyperParam=hyperParam))
  gseResults <- suppressWarnings(Create_gseResults(gseParam=gseParam))

  MSigDB = suppressWarnings(Create_MSigDB(msigdbParam = msigdbParam,
                                          msigdbGSEAparam = Create_msigdbGSEAparam(),
                                          msigdbHyperParam = Create_msigdbHyperParam()))

  newDEGContainer(dataInfo = dataInfo,
                  degResults = degResults,
                  hyperResults = hyperResults,
                  gseResults = gseResults,
                  MSigDB = MSigDB)

}

# https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
# https://www.genome.jp/kegg/catalog/org_list.html
shotSpecies <- function(species) {
  if(presetSpecies(species)) {
    usethis::ui_done("Your species will be resolved all automaticly")
  } else {
    usethis::ui_oops("The species {ui_value(species)} is out one of preset {ui_value('Human,Mouse,Rat')}");
    usethis::ui_todo("We need you provide more parameters for your species.");
    usethis::ui_todo("{ui_code('OrgDb')} parameters should be set one of below page");
    usethis::ui_todo("{ui_path('https://bioconductor.org/packages/release/BiocViews.html#___OrgDb')}");
    usethis::ui_todo("{ui_code('organism')} parameters should be set one of below page");
    usethis::ui_todo("{ui_path('https://www.genome.jp/kegg/catalog/org_list.html')}");
    usethis::ui_todo("Or you can use custom species by {ui_code('GOTERM2GENE')} and {ui_code('GOTERM2GENE')}
Here is a help page to more details. {ui_path('www.xiayh17.com')}");
  }
}

presetSpecies <- function(species) {

  species %in% c("Human","Mouse","Rat")

}
