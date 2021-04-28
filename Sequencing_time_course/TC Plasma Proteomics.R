library(openxlsx)
library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)
library(clipr)
library(limma)
library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(Rmisc)
library(gridExtra)
library(UpSetR)
library(pheatmap)
library(here)
?arrange
Plasma_proteomics_dataset <- openxlsx::read.xlsx("Sequencing_time_course/Plasma Proteomics dataset sorted.xlsx")
Plasma_proteomics_dataset <- arrange(Plasma_proteomics_dataset, )
setup_proteomics <- openxlsx::read.xlsx("Sequencing_time_course/Setup_plasma.xlsx")
setup_proteomics <- setup_proteomics %>%
  unite(group, c("Genotype", "Time"), remove = F)
setDT(setup_proteomics)
setup_proteomics <- setup_proteomics %>%
  mutate(ID = as.character(ID))
setkey(setup_proteomics, ID)
View(setup_proteomics)
setup_proteomics <- setup_proteomics %>%
  mutate(sample = ID)%>%
  dplyr::select(-ID)

normalized_proteomics <- normalizeBetweenArrays(log(as.matrix(Plasma_proteomics_dataset[,-c(1:2)])), method = "quantile")
rownames(normalized_proteomics) <- Plasma_proteomics_dataset$Protein_ID

missingSamples_proteomics <- data.table(is.na(normalized_proteomics), keep.rownames = TRUE) %>% 
  melt(measure.vars = colnames(normalized_proteomics), variable.name = "sample")

missingSamples_proteomics <- merge(setup_proteomics, missingSamples_proteomics, by = "sample")
setnames(missingSamples_proteomics, "rn", "Accession")
missingSamples_proteomics <- missingSamples_proteomics %>%
  dplyr::group_by(Accession, group) %>%
  dplyr::mutate(nMissing = sum(value))%>%
  reshape2::dcast(Accession ~ group, value.var = "nMissing", fun.aggregate = mean)

setup_proteomics <- setup_proteomics %>%
  group_by(group)

cutoff <- 2

tooManyMissing_proteomics <- missingSamples_proteomics %>%
  filter(KO_12 > cutoff|
           KO_21 > cutoff | 
           KO_3 > cutoff | 
           KO_6 > cutoff | 
           WT_3 > cutoff| 
           WT_6 > cutoff | 
           WT_12 > cutoff | 
           WT_21 > cutoff)

normalized_proteomics_res <- normalized_proteomics[!(rownames(normalized_proteomics) %in% tooManyMissing_proteomics$Accession), ]
all(colnames(normalized_proteomics_res) == setup_proteomics$sample)
setup_proteomics$sample %in% colnames(normalized_proteomics_res)
colnames(normalized_proteomics_res) %in% setup_proteomics$sample 


View(normalized_proteomics_res)
View(setup_proteomics)



mdsData_proteomics <- plotMDS(normalized_proteomics_res, ndim = 3, plot = FALSE)$cmdscale.out
mdsData_proteomics <- as.data.frame(mdsData_proteomics)
mdsData_proteomics_test <- cbind.data.frame(setup_proteomics, mdsData_proteomics)


ggplot(mdsData_proteomics_test, aes(x = V1, y = V2, colour = group)) +
  geom_point() +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

ggplot(mdsData_proteomics_test, aes(x = V1, y = V3, colour = group)) +
  geom_point() +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

ggplot(mdsData_proteomics_test, aes(x = V2, y = V3, colour = group)) +
  geom_point() +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

design_proteomics <- model.matrix(~ 0 + group, setup_proteomics)

colnames(design_proteomics) <- str_remove_all(colnames(design_proteomics), "group")
fit_proteomics <- lmFit(normalized_proteomics_res, design = design_proteomics, method = "robust")
cont.matrix_proteomics <- makeContrasts(
  Genotype_3 = KO_3 - WT_3,
  Genotype_6 = KO_6 - WT_6,
  Genotype_12 = KO_12 - WT_12,
  Genotype_21 = KO_21 - WT_21,
  HNKO_age_21_3 = KO_21 - KO_3,
  WT_age_21_3 = WT_21 - WT_3,
  Interaction = (KO_21 - WT_21) - (KO_3 - WT_3),
  levels = design_proteomics
)
View(cont.matrix_proteomics)
fit2_proteomics <- contrasts.fit(fit_proteomics, cont.matrix_proteomics)
View(fit2_proteomics)
fit2_proteomics <- eBayes(fit2_proteomics, trend = TRUE, robust = TRUE) 
View(fit2_proteomics) 
resultTables_proteomics <- list(
  Genotype_3 = topTable(fit2_proteomics, coef = "Genotype_3", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  Genotype_6 = topTable(fit2_proteomics, coef = "Genotype_6", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  Genotype_12 = topTable(fit2_proteomics, coef = "Genotype_12", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  Genotype_21 = topTable(fit2_proteomics, coef = "Genotype_21", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  HNKO_age_21_3 = topTable(fit2_proteomics, coef = "HNKO_age_21_3", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  WT_age_21_3 = topTable(fit2_proteomics, coef = "WT_age_21_3", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  Interaction = topTable(fit2_proteomics, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE)
)
View(resultTables_proteomics)
lapply(resultTables_proteomics, setnames, "rn", "Accession")
conv_prot <- Plasma_proteomics_dataset[, 1:2]
colnames(conv_prot)[2] <- "Accession"
conv_prot <- as.data.table(conv_prot)
setkey(conv_prot, Accession)
View(conv_prot)
for (i in resultTables_proteomics){
  i[, Gene:=conv_prot[Accession, Gene_names]]
}
#write.xlsx(resultTables_proteomics, file = "limma_results_proteomics_time_course.xlsx")

#####Load proteomics data #####
resultsTables_proteomics <- list(
  "Genotype_3" = NA,
  "Genotype_6" = NA,
  "Genotype_12" = NA,
  "Genotype_21" = NA
)
for (i in 1:4){
resultsTables_proteomics[[i]] <- openxlsx::read.xlsx(here::here("limma_results_proteomics_time_course.xlsx"), i)
}
#GO analysis
HNKO_3d <- resultTables_proteomics$Genotype_3
View(HNKO_3d)
HNKO_3d_sig <- HNKO_3d %>%
  filter(adj.P.Val<0.05)
View(HNKO_3d_sig)

HNKO_3d_enztrez= bitr(HNKO_3d_sig$Gene, 
                      fromType ="SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db",
                      drop = T)
View(HNKO_3d_enztrez)

HNKO_3d_bg= bitr(HNKO_3d$Gene, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Mm.eg.db",
                 drop = T)
dim(HNKO_3d_bg)
dim(HNKO_3d_enztrez)

goResults_HNKO_3d<- enrichGO(gene = HNKO_3d_enztrez$ENTREZID,
                             universe = HNKO_3d_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")
HNKO_6d <- resultTables_proteomics$Genotype_6
HNKO_6d_sig <- HNKO_6d %>%
  filter(adj.P.Val<0.05)

HNKO_6d_enztrez= bitr(HNKO_6d_sig$Gene, 
                      fromType ="SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db",
                      drop = T)

goResults_HNKO_6d<- enrichGO(gene = HNKO_6d_enztrez$ENTREZID,
                             universe = HNKO_3d_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")
HNKO_12d <- resultTables_proteomics$Genotype_12
HNKO_12d_sig <- HNKO_12d %>%
  filter(adj.P.Val<0.05)

HNKO_12d_entrez = bitr(HNKO_12d_sig$Gene, 
                       fromType ="SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)

goResults_HNKO_12d<- enrichGO(gene = HNKO_12d_entrez$ENTREZID,
                             universe = HNKO_3d_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")
HNKO_21d <- resultTables_proteomics$Genotype_21
HNKO_21d_sig <- HNKO_21d %>%
  filter(adj.P.Val<0.05)

HNKO_21d_entrez = bitr(HNKO_21d_sig$Gene, 
                       fromType ="SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)

goResults_HNKO_21d<- enrichGO(gene = HNKO_21d_entrez$ENTREZID,
                              universe = HNKO_3d_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "BP")

clusterProfiler::dotplot(goResults_HNKO_3d)
clusterProfiler::dotplot(goResults_HNKO_6d)
clusterProfiler::dotplot(goResults_HNKO_12d)
clusterProfiler::dotplot(goResults_HNKO_21d)
View(goResults_HNKO_21d@result)
View(goResults_HNKO_6d@result)

View(goResults_HNKO_6d@result)




goResults_HNKO_6d <- setReadable(goResults_HNKO_6d, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
clusterProfiler::cnetplot(goResults_HNKO_6d)
View(gogoResults_HNKO_6d@result)
View(HNKO_6d_sig)
#No significant GO-terms from 3d

######Upsetplot#####

significant_proteins_upset <- list("HNKO 3d" = resultsTables_proteomics[[1]],
                                   "HNKO 6d" = resultsTables_proteomics[[2]],
                                   "HNKO 12d" = resultsTables_proteomics[[3]],
                                   "HNKO 21d" = resultsTables_proteomics[[4]])
for (i in 1:4){
  significant_proteins_upset [[i]]  <- significant_proteins_upset [[i]] %>% 
    dplyr::filter(adj.P.Val < 0.05) 
  significant_proteins_upset [[i]]  <- significant_proteins_upset [[i]]$Gene
}

order_upset <- c("HNKO 21d", "HNKO 12d", "HNKO 6d","HNKO 3d")
View(significant_proteins_upset)
UpSetR::upset(UpSetR::fromList(significant_proteins_upset),
      sets = order_upset,
      order.by = "freq", 
      keep.order = T,
      text.scale = 2
      )

  grid::grid.text("Plasma Proteomics", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 30))

View(significant_proteins_upset)
library(grid)

#####Extracting the Proteasome subunits #####
goResults_HNKO_6d@result$geneID[1]
Proteinlist <- goResults_HNKO_6d@result$gene[1]
Proteinlist_table <- data.frame(Proteinlist)
View(Proteinlist_table)
Proteinlist_table <- read.table(text = Proteinlist, sep = "/") 
Proteinlist_table <- t(Proteinlist_table)



Candidate_proteins <- Plasma_proteomics_dataset %>%
  filter(Gene_names %in% Proteinlist_table)
View(Candidate_proteins)

Candidate_proteins_hm <- Candidate_proteins %>%
  dplyr::select(-c(Gene_names, Protein_ID))

row.names(Candidate_proteins_hm)<- Candidate_proteins$Gene_names



#####Heatmap#####
pheatmap(Candidate_proteins_hm,
         treeheight_col = 0,
         treeheight_row = 0,
         scale = "row",
         legend = T,
         na_col = "white",
         Colv = NA,
         na.rm = T,
         cluster_cols = F,
         fontsize_row = 8,
         fontsize_col = 11,
         cellwidth = 12,
         cellheight = 7
)

######Extract unique candidates and overlap candidates#####
overlap_proteins <- as.data.frame(significant_proteins_upset[[1]]) %>%
  dplyr::filter(significant_proteins_upset[[1]] %in% significant_proteins_upset[[2]] &
                  significant_proteins_upset[[1]] %in% significant_proteins_upset[[3]] &
                  significant_proteins_upset[[1]] %in% significant_proteins_upset[[4]]
  )
colnames(overlap_proteins) <- "Genes"

Plasma_proteomics_dataset <- openxlsx::read.xlsx("Sequencing_time_course/Plasma Proteomics dataset sorted.xlsx")

setup_proteomics <- openxlsx::read.xlsx("Sequencing_time_course/Setup_plasma.xlsx")
setup_proteomics <- setup_proteomics %>%
  unite(group, c("Genotype", "Time"), remove = F)

setup_heatmap <- setup_proteomics %>% 
  dplyr::arrange(Time, desc(Genotype))

Plasma_proteomics_overlap <- Plasma_proteomics_dataset %>% 
  dplyr::filter(Gene_names %in% overlap_proteins$Genes) %>% 
  dplyr::distinct(Gene_names, .keep_all = T)
rownames(Plasma_proteomics_overlap)<- Plasma_proteomics_overlap$Gene_names
Plasma_proteomics_overlap <- Plasma_proteomics_overlap %>% 
  dplyr::select(-Gene_names, -Protein_ID)
setup_heatmap$ID <- as.character(setup_heatmap$ID)
Plasma_proteomics_overlap <- Plasma_proteomics_overlap %>% 
  dplyr::select(setup_heatmap$ID)

setup_heatmap <- setup_heatmap %>% 
  dplyr::mutate(group = case_when(
    group == "WT_3"~"WT 3",
    group == "KO_3"~"HNKO 3",
    group == "WT_6"~"WT 6",
    group == "KO_6"~ "HNKO 6",
    group == "WT_12"~"WT 12",
    group == "KO_12"~"HNKO 12",
    group == "WT_21"~"WT 21",
    group == "KO_21"~ "HNKO 21"
  ))

heatmap_key <- setup_heatmap %>% 
  dplyr::select(ID, group)
rownames(heatmap_key) <- heatmap_key$ID
heatmap_key <- heatmap_key %>% 
  dplyr::select(-ID)
heatmap_key$group <- factor(heatmap_key$group, levels = c("WT 3", "HNKO 3", "WT 6", "HNKO 6", "WT 12", "HNKO 12", "WT 21", "HNKO 21"))

pheatmap::pheatmap(Plasma_proteomics_overlap,
         treeheight_col = 0,
         treeheight_row = 0,
         scale = "row",
         legend = T,
         na_col = "white",
         Colv = NA,
         na.rm = T,
         cluster_cols = F,
         fontsize_row = 14,
         fontsize_col = 12,
         cellwidth = 10,
         cellheight = 14,
         annotation_col = heatmap_key,
         labels_col = "",
         main = "Proteins w. main effect of genotype at all time points"
)

#####Extract early event proteins####
unique_day3 <- overlap_proteins <- HNKO_3d_sig
