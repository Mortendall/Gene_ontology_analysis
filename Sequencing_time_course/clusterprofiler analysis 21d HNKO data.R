library(openxlsx)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(clipr)
library(limma)
library(here)
library(tidyverse)
library(UpSetR)

expressions <- read_excel("R/Re-run RNA seq set/21d HNKO effect RNAseq.xlsx")
expressions_significant <- expressions %>%
  filter(FDR < 0.05)

expressions_enztrez= bitr(expressions_significant$ENSEMBL, 
                    fromType = "ENSEMBL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Mm.eg.db",
                    drop = T)


expressions_enztrez_bg = bitr(expressions$ENSEMBL, 
                        fromType = "ENSEMBL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Mm.eg.db",
                        drop = T)

goResults_expressions <- enrichGO(gene = expressions_enztrez$ENTREZID,
                            universe = expressions_enztrez_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP")
dotplot(goResults_expressions)

goResults <- setReadable(goResults_expressions, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults)

#####3d data####
expressions_3d <- read_excel("R/Re-run RNA seq set/3d HNKO effect RNAseq.xlsx")
expressions_significant_3d <- expressions_3d %>%
  filter(FDR < 0.05)

expressions_3d_enztrez= bitr(expressions_significant_3d$ENSEMBL, 
                          fromType = "ENSEMBL", 
                          toType = "ENTREZID", 
                          OrgDb = "org.Mm.eg.db",
                          drop = T)


expressions_3d_enztrez_bg = bitr(expressions_3d$ENSEMBL, 
                              fromType = "ENSEMBL", 
                              toType = "ENTREZID", 
                              OrgDb = "org.Mm.eg.db",
                              drop = T)

goResults_expressions_3d <- enrichGO(gene = expressions_3d_enztrez$ENTREZID,
                                  universe = expressions_3d_enztrez_bg$ENTREZID,
                                  OrgDb = org.Mm.eg.db,
                                  ont = "BP")
dotplot(goResults_expressions_3d)

####6d data ####


goResults_3d <- setReadable(goResults_expressions_3d, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults_3d)

#####Upset plot of edgeR data#####
#import edgeR data from Lars

edgeR_data <- list("HNKO 3d" = NA,
                   "HNKO 6d" = NA,
                   "HNKO 12d" = NA,
                   "HNKO 21d" = NA)

for (i in 1:4){
  edgeR_data[[i]]<- openxlsx::read.xlsx(here::here("Sequencing_time_course/022_JonasTreebak_Anna_edgeR_results.xlsx"),i)
}

edgeR_sig <- edgeR_data  
for (i in 1:4){
  edgeR_sig[[i]]<-edgeR_sig[[i]] %>% 
    dplyr::filter(FDR < 0.05) 
    edgeR_sig[[i]]<-edgeR_sig[[i]]$SYMBOL
} 

order_upset <- c("HNKO 21d", "HNKO 12d", "HNKO 6d","HNKO 3d")
UpSetR::upset(fromList(edgeR_sig),
              sets = order_upset,
              order.by = "freq", 
              keep.order = T,
              text.scale = 3.5
)
grid::grid.text("Genes with effect of genotype", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 48))

#extract significant overlap genes

overlap_genes <- as.data.frame(edgeR_sig[[1]]) %>%
  dplyr::filter(edgeR_sig$KO_at_3d %in% edgeR_sig$KO_at_6d &
                  edgeR_sig$KO_at_3d %in% edgeR_sig$KO_at_12d &
                  edgeR_sig$KO_at_3d %in% edgeR_sig$KO_at_21d
                  )
overlap_genes <- overlap_genes %>% 
  dplyr::filter(!is.na(overlap_genes))
colnames(overlap_genes) <- "Symbol"

overlap_entrez <- bitr(overlap_genes$Symbol, 
                             fromType = "SYMBOL", 
                             toType = "ENTREZID", 
                             OrgDb = "org.Mm.eg.db",
                             drop = T)


overlap_bg = bitr(edgeR_data[[1]]$SYMBOL, 
                                 fromType = "SYMBOL", 
                                 toType = "ENTREZID", 
                                 OrgDb = "org.Mm.eg.db",
                                 drop = T)

goResults_overlap <- enrichGO(gene = overlap_entrez$ENTREZID,
                                     universe = overlap_bg$ENTREZID,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "BP")
enrichplot::dotplot(goResults_overlap)+ggtitle("Overlap Genes GO-terms")
cnet
cpm_matrix <- openxlsx::read.xlsx(here("Sequencing_time_course/cpmData sorted by genotype updated 17.03.20.xlsx"))

#extract 3d only genes
overlap_genes_3d <- as.data.frame(edgeR_sig[[1]]) %>%
  dplyr::filter(!edgeR_sig$KO_at_3d %in% edgeR_sig$KO_at_6d &
                  !edgeR_sig$KO_at_3d %in% edgeR_sig$KO_at_12d &
                  !edgeR_sig$KO_at_3d %in% edgeR_sig$KO_at_21d
  )
overlap_genes_3d <- overlap_genes_3d %>% 
  dplyr::filter(!is.na(overlap_genes_3d))
colnames(overlap_genes_3d) <- "Symbol"

overlap_entrez <- bitr(overlap_genes_3d$Symbol, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)


goResults_overlap_3d <- enrichGO(gene = overlap_entrez$ENTREZID,
                              universe = overlap_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "BP")
#enrichplot::dotplot(goResults_overlap_3d)+ggtitle("Overlap Genes GO-terms")
#no significant GO-terms

#repeat for 21d
overlap_genes_21d <- as.data.frame(edgeR_sig[[4]]) %>%
  dplyr::filter(!edgeR_sig$KO_at_21d %in% edgeR_sig$KO_at_3d &
                  !edgeR_sig$KO_at_21d %in% edgeR_sig$KO_at_12d &
                  !edgeR_sig$KO_at_21d %in% edgeR_sig$KO_at_6d
  )
overlap_genes_21d <- overlap_genes_21d %>% 
  dplyr::filter(!is.na(overlap_genes_21d))
colnames(overlap_genes_21d) <- "Symbol"

overlap_entrez <- bitr(overlap_genes_21d$Symbol, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)


goResults_overlap_21d <- enrichGO(gene = overlap_entrez$ENTREZID,
                                 universe = overlap_bg$ENTREZID,
                                 OrgDb = org.Mm.eg.db,
                                 ont = "BP")
enrichplot::dotplot(goResults_overlap_21d)+ggtitle("Genotype effect - Day 21")

#Day 6
