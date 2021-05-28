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
library(clusterProfiler)

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
upsetRNA <- UpSetR::upset(UpSetR::fromList(edgeR_sig),
              sets = order_upset,
              order.by = "freq", 
              keep.order = T,
              text.scale = 3.5
)

grid::grid.text("Genes with effect of genotype", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 48))
tiff("UpsetRNA.tif", unit = "cm", height = 25, width = 50, res = 600)
upsetRNA

grid::grid.text("Genes with effect of genotype", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))
dev.off()

#extract significant overlap genes


overlap_genes <- as.data.frame(edgeR_sig[[1]]) %>%
  dplyr::filter(edgeR_sig[[1]] %in% edgeR_sig[[2]] &
                  edgeR_sig[[1]] %in% edgeR_sig[[3]] &
                  edgeR_sig[[1]] %in% edgeR_sig[[4]]
                  )
overlap_genes <- overlap_genes %>% 
  dplyr::filter(!is.na(overlap_genes))
colnames(overlap_genes) <- "Symbol"

overlap_entrez <- clusterProfiler::bitr(overlap_genes$Symbol, 
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
rnaGO <- enrichplot::dotplot(goResults_overlap)+ggtitle("Overlap Genes GO-terms")

tiff("GORNA.tif", unit = "cm", height = 10, width = 25, res = 300)
rnaGO
dev.off()


cnet
cpm_matrix <- openxlsx::read.xlsx(here("Sequencing_time_course/cpmData sorted by genotype updated 17.03.20.xlsx"))

#extract 3d only genes
overlap_genes_3d <- as.data.frame(edgeR_sig[[1]]) %>%
  dplyr::filter(!edgeR_sig[[1]] %in% edgeR_sig[[2]] &
                  !edgeR_sig[[1]] %in% edgeR_sig[[3]] &
                  !edgeR_sig[[1]] %in% edgeR_sig[[4]]
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
  dplyr::filter(!edgeR_sig[[4]] %in% edgeR_sig[[2]] &
                  !edgeR_sig[[4]] %in% edgeR_sig[[3]] &
                  !edgeR_sig[[4]] %in% edgeR_sig[[1]]
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

#Day 3 significance test
genes_3d <- as.data.frame(edgeR_sig[[1]])
genes_3d <- genes_3d %>% 
  dplyr::filter(!is.na(genes_3d))
colnames(genes_3d) <- "Symbol"

day3_entrez <- bitr(genes_3d$Symbol, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)


goResults_3d<- enrichGO(gene = day3_entrez$ENTREZID,
                                  universe = overlap_bg$ENTREZID,
                                  OrgDb = org.Mm.eg.db,
                                  ont = "BP")
enrichplot::dotplot(goResults_3d)+ggtitle("Genotype effect - Day 3")
goResults_3d <- setReadable(goResults_3d, keyType = "ENTREZID", OrgDb = org.Mm.eg.db)


goResults_3d_CC<- enrichGO(gene = day3_entrez$ENTREZID,
                        universe = overlap_bg$ENTREZID,
                        OrgDb = org.Mm.eg.db,
                        ont = "CC")
enrichplot::dotplot(goResults_3d_CC)+ggtitle("Genotype effect - Day 3")
#day 21 significance test
genes_21d <- as.data.frame(edgeR_sig[[4]])
genes_21d <- genes_21d %>% 
  dplyr::filter(!is.na(genes_21d))
colnames(genes_21d) <- "Symbol"

day21_entrez <- bitr(genes_21d$Symbol, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Mm.eg.db",
                    drop = T)


goResults_21d<- enrichGO(gene = day21_entrez$ENTREZID,
                        universe = overlap_bg$ENTREZID,
                        OrgDb = org.Mm.eg.db,
                        ont = "BP")
enrichplot::dotplot(goResults_21d)+ggtitle("Genotype effect - Day 21")
goResults_21d <- setReadable(goResults_21d, keyType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrichplot::cnetplot(goResults_21d)

goResults_21d_CC<- enrichGO(gene = day21_entrez$ENTREZID,
                           universe = overlap_bg$ENTREZID,
                           OrgDb = org.Mm.eg.db,
                           ont = "CC")
enrichplot::dotplot(goResults_21d_CC)+ggtitle("Genotype effect - Day 21")

goResults_21d_MF<- enrichGO(gene = day21_entrez$ENTREZID,
                            universe = overlap_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "MF")
enrichplot::dotplot(goResults_21d_MF)+ggtitle("Genotype effect - Day 21")
goResults_21d_MF <- setReadable(goResults_21d_MF, keyType = "ENTREZID", OrgDb = org.Mm.eg.db)
enrichplot::cnetplot(goResults_21d_MF)
#heatmap of calcium ion binding proteins
calcium_binding_proteins <- goResults_21d_MF@result[2,]$geneID
Proteinlist_table <- read.table(text = calcium_binding_proteins, sep = "/") 
Proteinlist_table <- t(Proteinlist_table)
cpm_matrix <- openxlsx::read.xlsx(here::here("Sequencing_time_course/cpm_matrix.xlsx"))
setup <- openxlsx::read.xlsx(here::here("Sequencing_time_course/setup.xlsx"))

cpm_calcium <- cpm_matrix %>% 
  dplyr::filter(SYMBOL %in% Proteinlist_table)

setup <- setup %>% 
  dplyr::mutate(Group = paste(Genotype, Time, sep = "_"))
setup <- setup %>% 
  dplyr::arrange(Time, desc(Genotype))
setup$ID <- as.character(setup$ID)
colnames(cpm_calcium)<-as.character(colnames(cpm_calcium))
setup <- setup %>% 
  dplyr::filter(ID %in% colnames(cpm_calcium))

cpm_calcium_key <- cpm_calcium %>% 
  dplyr::select(ENSEMBL, SYMBOL)
cpm_calcium <- cpm_calcium %>% 
  dplyr::select(-ENSEMBL, -SYMBOL)

cpm_calcium <- cpm_calcium %>% 
  dplyr::select(setup$ID)
cpm_calcium <- as.matrix(cpm_calcium)
rownames(cpm_calcium)<-cpm_calcium_key$SYMBOL
key <- as.data.frame(setup$Group)
colnames(key)<-"Group"
rownames(key) <- setup$ID
key$Group <- factor(key$Group, c("WT_3","KO_3","WT_6","KO_6","WT_12","KO_12","WT_21","KO_21"))

calcium_hm<- pheatmap::pheatmap(cpm_calcium,
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
                            cellwidth = 6,
                            cellheight = 4,
                            annotation_col = key,
                            show_colnames = F,
                            show_rownames = F,
                            main = "Oxidation-reduction process"
)
