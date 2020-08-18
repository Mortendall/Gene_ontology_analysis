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
?dotplot
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
