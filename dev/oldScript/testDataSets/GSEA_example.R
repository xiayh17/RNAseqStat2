library(msigdbr)
library(clusterProfiler)

data(geneList, package="DOSE")
head(geneList)

gene <- names(geneList)[abs(geneList) > 2]
head(gene)

msigdbr_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em <- enricher(gene, TERM2GENE=m_t2g)
head(em)

C3_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>%
  dplyr::select(gs_name, entrez_gene)
head(C3_t2g)

em2 <- GSEA(geneList, TERM2GENE = C3_t2g)
head(em2)
