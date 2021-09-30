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
                                     ont = "MF")
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
                        ont = "MF")
enrichplot::dotplot(goResults_3d_CC)+ggtitle("Genotype effect - Day 3 MF")
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
                            show_rownames = T,
                            main = "Calcium Binding Proteins"
)

#check volcano plot of calcium ion genes for day 3
calcium_edgeR_data <- edgeR_data[[1]] %>% 
  dplyr::filter(SYMBOL %in% Proteinlist_table) 

ggplot2::ggplot(calcium_edgeR_data, aes(x = logFC, y = -log(FDR), label = SYMBOL))+
  geom_point()+
  xlim(-5, 5)+
  geom_text(aes(label = SYMBOL))
m
#heatmap of fatty acid metabolic processes
FA_binding_proteins <- goResults_3d@result[5,]$geneID
Proteinlist_table <- read.table(text = FA_binding_proteins, sep = "/") 
Proteinlist_table <- t(Proteinlist_table)
cpm_matrix <- openxlsx::read.xlsx(here::here("Sequencing_time_course/cpm_matrix.xlsx"))
setup <- openxlsx::read.xlsx(here::here("Sequencing_time_course/setup.xlsx"))

cpm_FA <- cpm_matrix %>% 
  dplyr::filter(SYMBOL %in% Proteinlist_table)

setup <- setup %>% 
  dplyr::mutate(Group = paste(Genotype, Time, sep = "_"))
setup <- setup %>% 
  dplyr::arrange(Time, desc(Genotype))
setup$ID <- as.character(setup$ID)
colnames(cpm_FA)<-as.character(colnames(cpm_FA))
setup <- setup %>% 
  dplyr::filter(ID %in% colnames(cpm_FA))

cpm_FA_key <- cpm_FA %>% 
  dplyr::select(ENSEMBL, SYMBOL)
cpm_FA <- cpm_FA %>% 
  dplyr::select(-ENSEMBL, -SYMBOL)


cpm_FA <- cpm_FA %>% 
  dplyr::select(setup$ID)
cpm_FA <- as.matrix(cpm_FA)
rownames(cpm_FA)<-cpm_FA_key$SYMBOL
key <- as.data.frame(setup$Group)
colnames(key)<-"Group"
rownames(key) <- setup$ID
key$Group <- factor(key$Group, c("WT_3","KO_3","WT_6","KO_6","WT_12","KO_12","WT_21","KO_21"))

FA_hm<- pheatmap::pheatmap(cpm_FA,
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
                                cellwidth = 8,
                                cellheight = 6,
                                annotation_col = key,
                                show_colnames = F,
                                show_rownames = T,
                                main = "Fatty acid metabolic processes genes"
)

 tiff("FA_hm.tif", unit = "cm", height = 12, width = 24, res = 600)
 FA_hm
 dev.off()

#####Showing CPM for individual genes####

cpm_FA_long <- as.data.frame(cpm_FA)
cpm_FA_long <- cpm_FA_long %>% 
  dplyr::mutate(Symbol = rownames(cpm_FA_long))
  cpm_FA_long <- pivot_longer(cpm_FA_long, cols = -Symbol, names_to = "ID", values_to = "CPM")

cpm_joined <- left_join(cpm_FA_long,setup, by = "ID")
cpm_joined$Time <- as.factor(cpm_joined$Time)

cpm_summary <- cpm_joined %>%  
  group_by(Symbol, Group) %>% 
  rstatix::get_summary_stats(CPM, type = "mean_se")

cpm_summary$Group <- factor(cpm_summary$Group, levels = c("WT_3","KO_3","WT_6","KO_6","WT_12","KO_12","WT_21","KO_21"))

ggplot(subset(cpm_summary, Symbol == "Sirt1"), aes(Group, mean))+
  geom_col()

#####Filtering and visualizing ck7, EpcAm and NCAM####
cpm_matrix <- openxlsx::read.xlsx(here::here("Sequencing_time_course/cpm_data/cpmData.xlsx"))
setup <- openxlsx::read.xlsx(here::here("Sequencing_time_course/Setup_RNAseq.xlsx"))
sample_ID <- openxlsx::read.xlsx(here::here("Sequencing_time_course/cpm_data/sampleID.xlsx"))
sample_ID$'#'<-paste("022",sample_ID$'#',sep = "_")
setup <- left_join(setup,sample_ID, by = c("ID"="Sample.ID"))

cpm_candidates <- cpm_matrix %>% 
  dplyr::filter(SYMBOL %in% c("Ncam1", "Epcam", "Krt7"))
cpm_cand_long <-cpm_candidates %>% 
  dplyr::select(-ENSEMBL, -GENENAME) %>% 
  pivot_longer(cols = -SYMBOL, names_to = "ID", values_to = "CPM")
cpm_cand_long$ID <- as.character(cpm_cand_long$ID)
setup <- setup %>% 
  dplyr::mutate(Genotype = case_when(
    Genotype == "KO"~"HNKO",
    Genotype == "WT"~"WT",
    TRUE ~ as.character(Genotype)
  )) %>% 
  dplyr::mutate(Group = paste(Genotype, Time, sep = " "))
cpm_joined <- dplyr::left_join(cpm_cand_long, setup, by = c("ID"="#"))
cpm_joined$Group <- factor(cpm_joined$Group, levels = c("WT 3", "WT 6", "WT 12", "WT 21", "HNKO 3", "HNKO 6", "HNKO 12", "HNKO 21"))

sum_stats <- cpm_joined%>% 
  dplyr::select(SYMBOL, CPM, Group, Genotype) %>% 
  group_by(SYMBOL, Group) %>% 
  rstatix::get_summary_stats(type = "mean_sd")
Ncam_data <-  subset(sum_stats, SYMBOL == "Ncam1")
Ncam <- ggplot2::ggplot()+
  geom_bar(data = subset(sum_stats, SYMBOL == "Ncam1"), aes(x = Group, y = mean), stat = "identity")+
  geom_point(data = subset(cpm_joined, SYMBOL == "Ncam1"), aes(x = Group, y = CPM, color = Genotype))+
  geom_errorbar(data = subset(sum_stats, SYMBOL == "Ncam1"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
  ggtitle("Ncam")

Epcam <- ggplot2::ggplot()+
  geom_bar(data = subset(sum_stats, SYMBOL == "Epcam"), aes(x = Group, y = mean), stat = "identity")+
  geom_point(data = subset(cpm_joined, SYMBOL == "Epcam"), aes(x = Group, y = CPM, color = Genotype))+
  geom_errorbar(data = subset(sum_stats, SYMBOL == "Epcam"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
  ggtitle("Epcam")

Krt7 <- ggplot2::ggplot()+
  geom_bar(data = subset(sum_stats, SYMBOL == "Krt7"), aes(x = Group, y = mean), stat = "identity")+
  geom_point(data = subset(cpm_joined, SYMBOL == "Krt7"), aes(x = Group, y = CPM, color = Genotype))+
  geom_errorbar(data = subset(sum_stats, SYMBOL == "Krt7"), aes(x = Group, y = mean, ymin = mean-sd, ymax = mean+sd), stat = "identity")+
  ggtitle("Krt7")

(Ncam+Epcam)/(Krt7+patchwork::plot_spacer())

