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
library(readxl)
#####limma analysis#####

Proteomics_dataset <- openxlsx::read.xlsx(here("Sequencing_time_course/Proteomics dataset.xlsx"))
setup_proteomics <- read_excel("R/Re-run RNA seq set/Setup.xlsx")
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

normalized_proteomics <- normalizeBetweenArrays(log(as.matrix(Proteomics_dataset[,-c(1:2)])), method = "quantile")
rownames(normalized_proteomics) <- Proteomics_dataset$Protein_ID

missingSamples_proteomics <- data.table(is.na(normalized_proteomics), keep.rownames = TRUE) %>% 
  melt(measure.vars = colnames(normalized_proteomics), variable.name = "sample")


missingSamples_proteomics <- merge(setup_proteomics, missingSamples_proteomics, by = "sample")
setnames(missingSamples_proteomics, "rn", "Accession")
missingSamples_proteomics <- missingSamples_proteomics %>%
  group_by(Accession, group) %>%
  mutate(nMissing = sum(value))%>%
  reshape2::dcast(Accession ~ group, value.var = "nMissing", fun.aggregate = mean)
View(missingSamples_proteomics)
setup_proteomics <- setup_proteomics %>%
  group_by(group)
?dcast
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

View(tooManyMissing_proteomics)                
normalized_proteomics_res <- normalized_proteomics[!(rownames(normalized_proteomics) %in% tooManyMissing_proteomics$Accession), ]
View(normalized_proteomics_res)

all(colnames(normalized_proteomics_res) == setup_proteomics$sample)

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
conv_prot <- Proteomics_dataset[, 1:2]
colnames(conv_prot)[2] <- "Accession"
conv_prot <- as.data.table(conv_prot)
setkey(conv_prot, Accession)
View(conv_prot)
for (i in resultTables_proteomics){
  i[, Gene:=conv_prot[Accession, Gene_names]]
}
#write.xlsx(resultTables_proteomics, file = "limma_results_proteomics_time_course_liver.xlsx")
getwd()

#####GO analysis#####
resultTables_proteomics <- list("Genotype_3" = NA,
                                "Genotype_6" = NA,
                                "Genotype_12" = NA,
                                "Genotype_21" = NA
                                )
for (i in 1:4){
  resultTables_proteomics[[i]] <- openxlsx::read.xlsx(here("Sequencing_time_course/limma_results_proteomics_time_course_liver.xlsx"),i)
}

proteomics_sig <- resultTables_proteomics 
for (i in 1:4){
  proteomics_sig[[i]]<-proteomics_sig[[i]] %>% 
    dplyr::filter(adj.P.Val < 0.05) 
  proteomics_sig[[i]]<-proteomics_sig[[i]]$Gene
} 

order_upset <- c("Genotype_21", "Genotype_12", "Genotype_6","Genotype_3")
UpSetR::upset(fromList(proteomics_sig),
              sets = order_upset,
              order.by = "freq", 
              keep.order = T,
              text.scale = 2,
)
grid::grid.text("Proteins with effect of genotype", x=0.65, y = 0.95, gp=grid::gpar(fontsize = 24))

#GO for overlap genes
overlap_genes <- as.data.frame(proteomics_sig[[1]]) %>%
  dplyr::filter(proteomics_sig$Genotype_3 %in% proteomics_sig$Genotype_6 &
                  proteomics_sig$Genotype_3 %in% proteomics_sig$Genotype_12 &
                  proteomics_sig$Genotype_3 %in% proteomics_sig$Genotype_21
  )
overlap_genes <- overlap_genes %>% 
  dplyr::filter(!is.na(overlap_genes))
colnames(overlap_genes) <- "Symbol"

overlap_entrez <- bitr(overlap_genes$Symbol, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)


overlap_bg = bitr(resultTables_proteomics[[1]]$Gene, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = "org.Mm.eg.db",
                  drop = T)

goResults_overlap <- enrichGO(gene = overlap_entrez$ENTREZID,
                              universe = overlap_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "BP")
enrichplot::dotplot(goResults_overlap)+ggtitle("Overlap Protein GO-terms")

#old GO analysis
HNKO_3d <- resultTables_proteomics$Genotype_3
View(HNKO_3d)
HNKO_3d_sig <- HNKO_3d %>%
  filter(adj.P.Val<0.05)
HNKO_3d_sig_neg <- HNKO_3d_sig %>%
  filter(logFC<0)
HNKO_3d_enztrez= bitr(HNKO_3d_sig$Gene, 
                          fromType ="SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = "org.Mm.eg.db",
                          drop = T)

HNKO_3d_neg_enztrez= bitr(HNKO_3d_sig_neg$Gene, 
                      fromType ="SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db",
                      drop = T)

HNKO_3d_bg= bitr(HNKO_3d$Gene, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = "org.Mm.eg.db",
                         drop = T)
dim(HNKO_3d_bg)
dim(HNKO_3d_enztrez)
goResults_HNKO_3d_neg<- enrichGO(gene = HNKO_3d_neg_enztrez$ENTREZID,
                              universe = HNKO_3d_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "BP")
goResults_HNKO_3d<- enrichGO(gene = HNKO_3d_enztrez$ENTREZID,
                                 universe = HNKO_3d_bg$ENTREZID,
                                 OrgDb = org.Mm.eg.db,
                                 ont = "BP")
dotplot(goResults_HNKO_3d_neg, showCategory = 10)
goResults_HNKO_3d_neg <- setReadable(goResults_HNKO_3d_neg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(goResults_HNKO_3d_neg)
?dotplot
goResults_western <- setReadable(goResults_western, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(goResults_HNKO_3d)
?cnetplot
View(goResults_HNKO_3d)

#HNKO day 21
day_21_proteins <- proteomics_sig$Genotype_21
day_21_entrez <-  bitr(day_21_proteins, 
                      fromType ="SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db",
                      drop = T)
goResults_HNKO_21d<- enrichGO(gene = day_21_entrez$ENTREZID,
                             universe = overlap_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")
enrichplot::dotplot(goResults_HNKO_21d)
#HNKO day 3
day_3_proteins <- proteomics_sig$Genotype_3
day_3_entrez <-  bitr(day_3_proteins, 
                       fromType ="SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)
goResults_HNKO_3d<- enrichGO(gene = day_3_entrez$ENTREZID,
                              universe = overlap_bg$ENTREZID,
                              OrgDb = org.Mm.eg.db,
                              ont = "BP")
enrichplot::dotplot(goResults_HNKO_3d)
#HNKO unique day 3
unique_3_genes <- as.data.frame(proteomics_sig[[1]]) %>%
  dplyr::filter(!proteomics_sig$Genotype_3 %in% proteomics_sig$Genotype_6 &
                  !proteomics_sig$Genotype_3 %in% proteomics_sig$Genotype_12 &
                  !proteomics_sig$Genotype_3 %in% proteomics_sig$Genotype_21
  )
unique_3_genes <- unique_3_genes %>% 
  dplyr::filter(!is.na(unique_3_genes))
colnames(unique_3_genes) <- "Symbol"

unique_3_entrez <- bitr(unique_3_genes$Symbol, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = "org.Mm.eg.db",
                       drop = T)
goResults_HNKO_3d_unique<- enrichGO(gene = unique_3_entrez$ENTREZID,
                             universe = overlap_bg$ENTREZID,
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")
enrichplot::dotplot(goResults_HNKO_3d_unique)

#unique day 21
unique_21_genes <- as.data.frame(proteomics_sig[[4]]) %>%
  dplyr::filter(!proteomics_sig$Genotype_21 %in% proteomics_sig$Genotype_3 &
                  !proteomics_sig$Genotype_21 %in% proteomics_sig$Genotype_6 &
                  !proteomics_sig$Genotype_21 %in% proteomics_sig$Genotype_12
  )
unique_21_genes <- unique_21_genes %>% 
  dplyr::filter(!is.na(unique_21_genes))
colnames(unique_21_genes) <- "Symbol"

unique_21_entrez <- bitr(unique_21_genes$Symbol, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Mm.eg.db",
                        drop = T)
goResults_HNKO_21d_unique<- enrichGO(gene = unique_21_entrez$ENTREZID,
                                    universe = overlap_bg$ENTREZID,
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP")
enrichplot::dotplot(goResults_HNKO_21d_unique)


#Analysis from Lili
ANOVA_results <- openxlsx::read.xlsx(here("Sequencing_time_course/Two-way_ANOVA_liver_fdr0.05.xlsx"))
bg_proteins <- ANOVA_results$Gene.name
sig_proteins <- ANOVA_results %>% 
  dplyr::filter(qvalue < 0.05) %>% 
  dplyr::select(Gene.name)
main_proteins<- bitr(sig_proteins$Gene.name, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = "org.Mm.eg.db",
                         drop = T)
goResults_main_proteins<- enrichGO(gene = main_proteins$ENTREZID,
                                     universe = overlap_bg$ENTREZID,
                                     OrgDb = org.Mm.eg.db,
                                     ont = "BP")
enrichplot::dotplot(goResults_main_proteins)+ggtitle("Proteins with main effect of genotype")
