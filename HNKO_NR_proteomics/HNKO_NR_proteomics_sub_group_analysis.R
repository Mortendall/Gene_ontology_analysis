#####Continued analysis of data from NR Proteomics analysis no.  1 #####

library(tidyverse)
library(openxlsx)
library(clusterProfiler)

limma_results <- vector(mode = "list", length = 4)
temporary_result <- data.frame(NA, NA, NA, NA)

for (i in 1:4){
temporary_result <- read.xlsx("C:/Users/tvb217/Documents/R/Lars Proteomic Data/limma_results.xlsx", i)
limma_results[[i]]<- temporary_result
}

names(limma_results)[[1]]<- "Treatment in WT"
names(limma_results)[[2]]<- "Treatment in KO"
names(limma_results)[[3]]<- "Genotype in control"
names(limma_results)[[4]]<- "Genotype in NR"

View(limma_results)

NR_effect <- limma_results[[2]]
HNKO_effect <- limma_results[[3]]

NR_effect_sig <- NR_effect %>%
  dplyr::filter(adj.P.Val < 0.05)
View(NR_effect_sig)

HNKO_effect_sig <- HNKO_effect %>%
  dplyr::filter(adj.P.Val < 0.05)

overlap <- HNKO_effect_sig %>%
  filter(HNKO_effect_sig$Gene %in% NR_effect_sig$Gene)

View(overlap)
#Analyze terms found in both 
Background<- bitr(HNKO_effect$Gene, 
                              fromType = "SYMBOL", 
                              toType = "ENTREZID", 
                              OrgDb = "org.Mm.eg.db",
                              drop = T)
overlap_entrez <- bitr(overlap$Gene,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = "org.Mm.eg.db",
                       drop = T)
go_Results_overlap <- enrichGO(gene = overlap_entrez$ENTREZID,
                             universe = Background$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")
clusterProfiler::dotplot(go_Results_overlap, showCategory = 10, title = "Overlap between NR treatment and HNKO control")

goResults_interaction <- enrichGO(gene = Interaction_enztrez$ENTREZID,
                                  universe = Interaction_enztrez_bg$ENTREZID,
                                  OrgDb = org.Mm.eg.db,
                                  ont = "BP")
#Analyze unique terms in NR group
Unique_NR <- NR_effect_sig %>%
  dplyr::filter(!NR_effect_sig$Gene %in% HNKO_effect_sig$Gene)
View(Unique_NR)
Unique_NR_entrez <- overlap_entrez <- bitr(Unique_NR$Gene,
                                           fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = "org.Mm.eg.db",
                                           drop = T)
View(Unique_NR_entrez)
go_Results_NR_unique <- enrichGO(gene = Unique_NR_entrez$ENTREZID,
                               universe = Background$ENTREZID,
                               OrgDb = org.Mm.eg.db,
                               ont = "BP")

clusterProfiler::dotplot(go_Results_NR_unique,  title = "Unique Terms of NR in HNKO")

View(go_Results_NR_unique@result)
#No significant data in the unique NR Response
#Analysis of unique HNKO control terms
Unique_HNKO <- HNKO_effect_sig %>%
  dplyr::filter(!HNKO_effect_sig$Gene %in% NR_effect_sig$Gene)

Unique_HNKO_entrez <- overlap_entrez <- bitr(Unique_HNKO$Gene,
                                           fromType = "SYMBOL",
                                           toType = "ENTREZID",
                                           OrgDb = "org.Mm.eg.db",
                                           drop = T)

go_Results_HNKO_unique <- enrichGO(gene = Unique_HNKO_entrez$ENTREZID,
                                 universe = Background$ENTREZID,
                                 OrgDb = org.Mm.eg.db,
                                 ont = "BP")

clusterProfiler::dotplot(go_Results_HNKO_unique,  title = "Unique Terms in HNKO in control")

#analyse effect of NR in HNKO
HNKO_entrez <- bitr(HNKO_effect_sig$Gene,
                                             fromType = "SYMBOL",
                                             toType = "ENTREZID",
                                             OrgDb = "org.Mm.eg.db",
                                             drop = T)
go_Results_HNKO <- enrichGO(gene = HNKO_entrez$ENTREZID,
                                   universe = Background$ENTREZID,
                                   OrgDb = org.Mm.eg.db,
                                   ont = "BP")
NR_entrez <- bitr(NR_effect_sig$Gene,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Mm.eg.db",
                    drop = T)
NR_entrez_bg <- bitr(NR_effect$Gene,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Mm.eg.db",
                  drop = T)

go_Results_NR <- enrichGO(gene = NR_entrez$ENTREZID,
                            universe = NR_entrez_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP")
enrichplot::dotplot(go_Results_NR)+ggtitle("HNKO Con vs NR" )

con_entrez <- bitr(HNKO_effect_sig$Gene,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = "org.Mm.eg.db",
                  drop = T)


go_Results_con <- enrichGO(gene = con_entrez$ENTREZID,
                          universe = NR_entrez_bg$ENTREZID,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP")
enrichplot::dotplot(go_Results_con)+ggtitle("Control WT vs HNKO" )

#####Upsetplot#####
limma_results_sig <- limma_results
  for (i in 1:4){
    limma_results_sig[[i]]<-limma_results_sig[[i]] %>% 
      dplyr::filter(adj.P.Val < 0.05) 
    limma_results_sig[[i]]<-limma_results_sig[[i]]$Gene
  } 

order_upset <- c("Genotype in NR", "Genotype in control", "Treatment in KO","Treatment in WT")
UpSetR::upset(fromList(limma_results_sig),
              sets = order_upset,
              order.by = "freq", 
              keep.order = T,
              text.scale = 2,
)
grid::grid.text("Significantly altered proteins", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 20))

#####heatmap of ox-red prod####
oxphos_proteins <- go_Results_NR@result$geneID[[1]]
oxphos_proteins <- unlist(str_split(oxphos_proteins, "/"))

cpm_key <-   clusterProfiler::bitr(
  oxphos_proteins,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = "org.Mm.eg.db"
)


expressions <- data.table::fread(here::here("HNKO_NR_proteomics/expressions.csv"), header = TRUE)
setup <- data.table::fread(here::here("HNKO_NR_proteomics/setup.csv"))
data.table::setnames(setup, c("sample", "Genotype", "Treatment"))
setup[, group:=paste(Genotype, Treatment, sep = "_")]
View(setup)
setup <- setup[sample != "330"]
#330 is the steatotic mouse
expressions <- expressions %>%
  dplyr::select(!"330")

res <- limma::normalizeBetweenArrays(log(as.matrix(expressions[,-c(1:2)])), method = "quantile")
res <- as.data.frame(res)
res <- res %>% 
  dplyr::mutate(Gene = expressions$Gene)

go_results <- list("Treatment in WT" = NA,
                   "Treatment in KO" = NA, 
                   "Genotype in control" = NA, 
                   "Genotype in NR" = NA)
for (i in 1:4){
  go_results[[i]]<- openxlsx::read.xlsx(here::here("HNKO_NR_proteomics/goData_NR_prot.xlsx"),i)
}
                   
  

oxphos_proteins <- go_results[[2]]$geneID[1]
oxphos_proteins <- unlist(str_split(oxphos_proteins, "/"))

cpm_key <-   clusterProfiler::bitr(
  oxphos_proteins,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = "org.Mm.eg.db"
)

res_ox <- res %>% 
  dplyr::filter(Gene %in% cpm_key$SYMBOL) %>% 
  dplyr::distinct(Gene, .keep_all = T)
rownames(res_ox) <- res_ox$Gene 
res_ox <- res_ox %>%
  dplyr::select(-Gene)

setup_ordered <- setup
order <- c("WT_Control", "KO_Control", "WT_NR", "KO_NR")

setup_ordered <- setup_ordered %>% 
  dplyr::arrange(Treatment,desc(Genotype))

class(colnames(res_ox))
setup_ordered$sample <- as.character(setup_ordered$sample)
res_ox <- res_ox %>% 
  dplyr::select(setup_ordered$sample)

key <- as.data.frame(setup_ordered)

key <- key %>% 
  dplyr::select(group)
rownames(key) <- setup_ordered$sample



pheatmap::pheatmap(res_ox,
         treeheight_col = 0,
         treeheight_row = 0,
         scale = "row",
         legend = T,
         na_col = "white",
         Colv = NA,
         na.rm = T,
         cluster_cols = F,
         fontsize_row = 5,
         fontsize_col = 8,
         cellwidth = 3,
         cellheight = 2,
         annotation_col = key,
         show_colnames = F,
         show_rownames = F,
         main = "Oxidation-reduction process"
)

#try to run the GO-analysis function from prim hep
#' Gene ontology enrichment analysis of genes generated from a results file
#'
#' @param result_list list of data.tables generated from edgeR. Must be data.table and contain a SYMBOL annotation column
#'
#' @return a list containing enrichresults for each element in the results file list

goAnalysis <- function(result_list){
  bg <- result_list[[1]]
  bg_list <- clusterProfiler::bitr(
    bg$Gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db",
    drop = T
  )
  
  goResult_list <- vector(mode = "list", length = length(result_list))
  for(i in 1:length(result_list)){
    sig_list<- result_list[[i]] %>%
      dplyr::filter(adj.P.Val<0.05)
    
    eg <- clusterProfiler::bitr(
      sig_list$Gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db",
      drop = T
    )
    goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                           universe = bg_list$ENTREZID,
                                           OrgDb = org.Mm.eg.db,
                                           ont = "BP")
    goResult_list[[i]]<- goResults
  }
  for (i in 1:length(goResult_list)){
    names(goResult_list)[i]<-names(result_list)[i]
  }
  return(goResult_list)
  
}
goTest_MF <- goAnalysis(limma_results)
write.xlsx(goTest, here::here("HNKO_NR_proteomics/goData_NR_prot.xlsx"))
