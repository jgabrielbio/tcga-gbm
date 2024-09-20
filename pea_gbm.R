# Pathway Expression Analysis on TCGA-GBM cohort
# Based on Cluster Profile vignette - https://yulab-smu.top/biomedical-knowledge-mining-book/
# João Gabriel - 01/03/2024
packages <- c('clusterProfiler', 'ggplot2', 'forcats', 'DOSE', 'org.Hs.eg.db', 'dplyr', 'tidyr', "stringr", "UpSetR", "pathview")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
rm(packages)
load("~/gbm/dea_dataframe_idh_status_alpha_0.01_lfc_1.5.rda")
gbm_dea <- data.frame(ensembls = substr(rownames(deaIDHLimited), 1 ,15), lfc = abs(deaIDHLimited$log2FoldChange))
gbm_dea$entrez_id <- mapIds(org.Hs.eg.db, keys = gbm_dea$ensembls, keytype="ENSEMBL", column = "ENTREZID")
gbm_dea <- na.omit(gbm_dea)
#Symbols genes ID
head(gbm_dea$entrez_id)
#Creating geneList
## assume 1st column is ID
## 2nd column is FC

#obs: deletar coluna type

## feature 1: numeric vector
gbmList = gbm_dea[,2]

## feature 2: named vector
names(gbmList) = as.character(gbm_dea[,3])

## feature 3: decreasing order
gbmList = sort(gbmList, decreasing = TRUE)

#saving list
save(gbmList, file = "gbmList.rda")
load("~/gbm/gbmList.rda")

############ Several analysis to put in Upset Plot ##########
#### Gene ontology ---- 
ego3 <- gseGO(geneList     = gbmList,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
#upsetplot(ego3, 8)
# #Over-representation analysis for disease ontology
# x <- enrichDO(gene          = gbm_dea$entrez_id,
#               ont           = "DO",
#               pvalueCutoff  = 0.05,
#               pAdjustMethod = "BH",
#               universe      = gbm_dea$entrez_id,
#               minGSSize     = 5,
#               maxGSSize     = 500,
#               qvalueCutoff  = 0.05,
#               readable      = TRUE)
# head(x)
# 
# gene2 <- names(gbmList)[abs(gbmList) > 2.5]
# ncg <- enrichNCG(gene2)
# ncg@geneSets[["glioblastoma"]]
# 
# dgn <- enrichDGN(gene2)
# head(dgn)
### Disease gene set enrichment analysis ----
#Disease ontology
#library(enrichplot)
#### Disease Ontology ---- 
y <- gseDO(gbmList,
           minGSSize     = 120,
           pvalueCutoff  = 0.05,
           pAdjustMethod = "BH",
           verbose       = FALSE)
y_geneset <- setReadable(y, 'org.Hs.eg.db')
head(y_geneset, 3)
upsetplot(y_geneset)
gseDO_df <- as.data.frame(y_geneset@result)
gbm_gseDO_entrez <- gseDO_df %>% filter(ID == "DOID:863")
gbm_gseDO_entrez <- str_split(gbm_gseDO_entrez, "/")
#gbm_gseDO_entrez <- gbm_gseDO_entrez[1]
others_entrez <- gseDO_df %>% 
  filter(ID != "DOID:863") 
#  select(core_enrichment) 
others_entrez <- str_split(others_entrez$core_enrichment, "/")
cell_type_cancer_entrez <- (others_entrez[1])
listInput <- list(one = gbm_gseDO_entrez, two = cell_type_cancer_entrez)
#comparing <- data.frame(Col1 = cell_type_cancer_entrez, Col2 = c('1', "2"))
upset(fromList(listInput), sets = gbm_gseDO_entrez, cell_type_cancer_entrez)
#upset(others_entrez, sets = others_entrez[1])
##### Network cancer genes ---- 
ncg_geneset <- gseNCG(gbmList,
              pvalueCutoff  = 1,
              pAdjustMethod = "BH",
              verbose       = FALSE)
ncg_geneset <- setReadable(ncg, 'org.Hs.eg.db')
gbm_entrez_ngc <- as.data.frame(ncg_geneset@geneSets[["glioblastoma"]])
upsetplot(ncg_geneset@result) 
##### Disease gene network ----
dgn <- gseDGN(gbmList,
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              verbose       = FALSE)
dgn <- setReadable(dgn, 'org.Hs.eg.db')
head(dgn, 3) 
results_df <- as.data.frame(dgn@result)
selection <- results_df[grepl("Glioblastoma", results_df$Description),]
selection_gene_ids <- selection %>% select(core_enrichment)

# running enrichment results
#kegg_gsea <- enrichKEGG(data(geneList, package="DOSE")
#gene <- names(geneList)[abs(geneList) > 2]

# kk <- enrichKEGG(gene         = names(gbmList),
#                  organism     = 'hsa',
#                  pvalueCutoff = 0.05)
# head(kk)


############ Analysis for Biological Process visualization ##########
#### GSEA KEGG GBM ----
kk2 <- gseKEGG(geneList     = gbmList,
               organism     = 'hsa',
               #minGSSize    = 120,
               pvalueCutoff = 0.05
               #verbose      = FALSE)
)
head(kk2)
save(kk2, file = 'gsea_kegg_gbm.rda')
load("~/gbm/gsea_kegg_gbm.rda")

upsetplot(kk2)
#Tidying
kk2_tidy <- kk2[, c("ID", "setSize", "Description", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue")]
kk2_arrange <- arrange(kk2_tidy, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:10)
save(kk2_tidy, file = "kk2_tidy.rda")
save(kk2_arrange, file = "kk2_arrange.rda")
load("~/gbm/kk2_tidy.rda")
load("~/gbm/kk2_arrange.rda")

#Biological Process visualization
ggplot(kk2_tidy, showCategory = 10,
       aes(setSize,
           fct_reorder(Description, setSize))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = setSize)) +
  scale_color_gradientn (colours=c("#f7ca64", "#46bac2","#7e62a3"),
                         trans = "log10",
                         guide=guide_colorbar(reverse=TRUE,
                                              order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Set Size") +
  #lab(NULL) +
  ggtitle("Biological Processes")

# KEGGpathways Visualization
ggplot(kk2_arrange, showCategory=10,
       aes(NES, fct_reorder(Description, NES),
           fill=qvalue)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe",
                                 "#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("KEGGPathways")
#### GSEA WIKI GBM ----
geneSet_wiki <- gseWP(gbmList, organism = "Homo sapiens", pvalueCutoff = 0.05) 
save(geneSet_wiki, file = "gsea_wiki_gbm.rda")
load("~/gbm/gsea_wiki_gbm.rda")
head(geneSet_wiki)
#Tidying
wiki2 <- geneSet_wiki[, c("ID", "setSize", "Description", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue")]
ewp2 <- arrange(wiki2, desc(abs(NES))) %>%
  group_by(sign(NES)) %>%
  slice(1:5)
save(wiki2, file = "wiki2.rda")
load("~/gbm/wiki2.rda")
save(ewp2, file = "ewp2.rda")
load("~/gbm/ewp2.rda")
#Biological Process visualization
ggplot(wiki2, showCategory = 10,
       aes(setSize,
           fct_reorder(Description, setSize))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = setSize)) +
  scale_color_gradientn (colours=c("#f7ca64", "#46bac2","#7e62a3"),
                         trans = "log10",
                         guide=guide_colorbar(reverse=TRUE,
                         order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Set Size") +
  #lab(NULL) +
  ggtitle("Biological Processes")

# Wikipathways Visualization
# Obs: arg fill in aes with problem.
ggplot(ewp2, showCategory=10,
       aes(NES, fct_reorder(Description, NES)#,
           #fill=qvalue)) 
           ))+
  geom_col() +
  scale_fill_gradientn(colours=c("#b3eebe",
                                 "#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("WikiPathways")
### Enrichment GO ----
enrich_go_gbm <- gseGO(geneList     = gbmList,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              #minGSSize    = 100,
              #maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
upsetplot(enrich_go_gbm, n = 8)
#### Enrichments maps ----
#KEGG Brite of glioma: 05214
kegg_map <- browseKEGG(kk2, 'hsa05214')
hsa05214 <- pathview(gene.data  = gbmList,
                     pathway.id = "hsa05214",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gbmList)), cpd=1))

#Wikipathways
if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways")
}
library(rWikiPathways)

if(!"RCy3" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)

# Use install.packages() for the following, if necessary:
library(magrittr)
