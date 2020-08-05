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
?arrange
Plasma_proteomics_dataset <- read_excel("Sequencing_time_course/Plasma Proteomics dataset sorted.xlsx")
Plasma_proteomics_dataset <- arrange(Plasma_proteomics_dataset, )
setup_proteomics <- read_excel("Sequencing_time_course/Setup_plasma.xlsx")
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
?dcast
cutoff <- 2
View(tooManyMissing_proteomics)
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
View(design_proteomics)
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
write.xlsx(resultTables_proteomics, file = "limma_results_proteomics_time_course.xlsx")


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

HNKO_6d_enztrez= bitr(HNKO_6d_sig$Gene, 
                      fromType ="SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db",
                      drop = T)

goResults_HNKO_6d<- enrichGO(gene = HNKO_6d_enztrez$ENTREZID,
                             universe = HNKO_3d_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")

HNKO_12d_entrez = bitr(HNKO_12d_sig$Gene, 
                       fromType ="SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)

goResults_HNKO_12d<- enrichGO(gene = HNKO_12d_entrez$ENTREZID,
                             universe = HNKO_3d_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")

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
goResults_western <- setReadable(goResults_western, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(goResults_HNKO_3d)


HNKO_6d_sig <- resultTables_proteomics$Genotype_6 %>%
  filter(adj.P.Val<0.05)
HNKO_12d_sig <- resultTables_proteomics$Genotype_12 %>%
  filter(adj.P.Val<0.05)

HNKO_21d_sig <- resultTables_proteomics$Genotype_21 %>%
  filter(adj.P.Val<0.05)

gogoResults_HNKO_6d <- setReadable(goResults_HNKO_6d, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(gogoResults_HNKO_6d)
View(gogoResults_HNKO_6d@result)
View(HNKO_6d_sig)
#No significant GO-terms from 3d

######Upsetplot#####
significant_proteins_upset <- list(HNKO_3d = HNKO_3d_sig$Gene,
                                   HNKO_6d = HNKO_6d_sig$Gene,
                                   HNKO_12d = HNKO_12d_sig$Gene,
                                   HNKO_21d = HNKO_21d_sig$Gene)
order_upset <- c("HNKO_21d", "HNKO_12d", "HNKO_6d","HNKO_3d")
View(significant_proteins_upset)
upset(fromList(significant_proteins_upset),
      sets = order_upset,
      order.by = "freq", 
      keep.order = T,
      )+
  grid.text("Plasma Proteomics", x=0.65, y = 0.95, gp=gpar(fontsize = 30))

View(significant_proteins_upset)
library(grid)

#####Extracting the Proteasome subunits #####
gogoResults_HNKO_6d@result$geneID[1]
Proteinlist <- gogoResults_HNKO_6d@result$gene[1]
Proteinlist_table <- data.frame(Proteinlist)
View(Proteinlist_table)
Proteinlist_table <- read.table(text = Proteinlist, sep = "/") 
Proteinlist_table <- t(Proteinlist_table)

View(Plasma_proteomics_dataset)

Candidate_proteins <- Plasma_proteomics_dataset %>%
  filter(Gene_names %in% Proteinlist_table)
View(Candidate_proteins)
