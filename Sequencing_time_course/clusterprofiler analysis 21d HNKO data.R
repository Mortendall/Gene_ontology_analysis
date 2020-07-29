library(openxlsx)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(clipr)
library(limma)

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

#3d data
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

goResults_3d <- setReadable(goResults_expressions_3d, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults_3d)


