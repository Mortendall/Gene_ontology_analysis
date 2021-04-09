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

expressions <- fread("R/Lars Proteomic Data/expressions.csv", header = TRUE)
setup <- fread("R/Lars Proteomic Data/setup.csv")
setnames(setup, c("sample", "Genotype", "Treatment"))
setup[, group:=paste(Genotype, Treatment, sep = "_")]
View(setup)
setup <- setup[sample != "330"]
#330 is the steatotic mouse
expressions <- expressions %>%
  dplyr::select(!"330")

res <- normalizeBetweenArrays(log(as.matrix(expressions[,-c(1:2)])), method = "quantile")
rownames(res) <- expressions$Accession

View(res)
View(expressions)
# Sometimes all samples in a group is missing
missingSamples <- data.table(is.na(res), keep.rownames = TRUE) %>% 
  melt(measure.vars = colnames(res), variable.name = "sample")
View(missingSamples)
setnames(missingSamples, "rn", "Accession")
setup[, sample:=as.character(sample)]
setkey(setup, sample)
missingSamples <- merge(setup, missingSamples, by = "sample")
missingSamples <- missingSamples[, .(nMissing = sum(value)), by = c("Accession", "group")] %>%
  dcast(Accession ~ group, value.var = "nMissing")

setup[, .N, by = "group"]
cutoff <- 2
tooManyMissing <- missingSamples[KO_Control > cutoff | 
                                   KO_NR > cutoff | 
                                   WT_Control > cutoff |
                                   WT_NR > cutoff,
                                 Accession
                                 ]
View(tooManyMissing)
res <- res[!(rownames(res) %in% tooManyMissing), ]

# Now impute
#resImputed <- impute::impute.knn(res)$data

selectedData <- res
#ind <- rowSums(is.na(selectedData)) == 0
#ind <- TRUE

all(colnames(selectedData) == setup$sample)

mdsData <- plotMDS(selectedData, ndim = 3, plot = FALSE)$cmdscale.out
View(mdsData)
View(setup)
all(rownames(mdsData) == setup$sample)
mdsData <- cbind(setup, mdsData)
ggplot(mdsData, aes(x = V1, y = V2, colour = group)) +
  geom_point() +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

ggplot(mdsData, aes(x = V1, y = V3, colour = group)) +
  geom_point() +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")

ggplot(mdsData, aes(x = V2, y = V3, colour = group)) +
  geom_point() +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set1")


design <- model.matrix(~ 0 + group, setup)
View(design)
colnames(design) <- str_remove_all(colnames(design), "group")
fit <- lmFit(selectedData, design = design, method = "robust")

cont.matrix <- makeContrasts(
  Treatment_in_WT = WT_NR - WT_Control,
  Treatment_in_KO = KO_NR - KO_Control,
  KO_in_Control = KO_Control - WT_Control,
  KO_in_NR = KO_NR - WT_NR,
  Interaction = (KO_NR - KO_Control) - (WT_NR - WT_Control),
  levels = design
)
View(cont.matrix)
fit2 <- contrasts.fit(fit, cont.matrix)
View(fit2)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

resultTables <- list(
  Treatment_in_WT = topTable(fit2, coef = "Treatment_in_WT", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  Treatment_in_KO = topTable(fit2, coef = "Treatment_in_KO", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  KO_in_Control = topTable(fit2, coef = "KO_in_Control", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  KO_in_NR = topTable(fit2, coef = "KO_in_NR", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
  Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE)
)
View(resultTables)
lapply(resultTables, setnames, "rn", "Accession")
conv <- expressions[, 1:2]
setkey(conv, Accession)

for (i in resultTables){
  i[, Gene:=conv[Accession, Gene]]
}
#write.xlsx(resultTables, file = "limma_results.xlsx")
getwd()

pD <- data.table(selectedData, keep.rownames = TRUE) %>% 
  melt(measure.vars = colnames(selectedData), variable.name = "sample")
setnames(pD, "rn", "Accession")
setup[, sample:=as.character(sample)]
setkey(setup, sample)
pD <- merge(setup, pD, by = "sample")
pD <- pD[, .(intensity = mean(value)), by = c("Accession", "group")] %>%
  dcast(Accession ~ group, value.var = "intensity")

ggplot(data = pD, aes(x = WT_Control, y = WT_NR, colour = Accession %in% resultTables$Treatment_in_WT[adj.P.Val < 0.05, Accession])) +
  geom_point() +
  geom_point(data = pD[Accession %in% resultTables$Treatment_in_WT[adj.P.Val < 0.05, Accession]]) +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set2")

ggplot(data = pD, aes(x = KO_Control, y = KO_NR, colour = Accession %in% resultTables$Treatment_in_KO[adj.P.Val < 0.05, Accession])) +
  geom_point() +
  geom_point(data = pD[Accession %in% resultTables$Treatment_in_KO[adj.P.Val < 0.05, Accession]]) +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set2")

ggplot(data = pD, aes(x = WT_Control, y = KO_Control, colour = Accession %in% resultTables$KO_in_Control[adj.P.Val < 0.05, Accession])) +
  geom_point() +
  geom_point(data = pD[Accession %in% resultTables$KO_in_Control[adj.P.Val < 0.05, Accession]]) +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set2")

ggplot(data = pD, aes(x = WT_NR, y = KO_NR, colour = Accession %in% resultTables$KO_in_NR[adj.P.Val < 0.05, Accession])) +
  geom_point() +
  geom_point(data = pD[Accession %in% resultTables$KO_in_NR[adj.P.Val < 0.05, Accession]]) +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set2")

ggplot(data = pD, aes(x = WT_NR - WT_Control, y = KO_NR - KO_Control, colour = Accession %in% resultTables$Interaction[adj.P.Val < 0.05, Accession])) +
  geom_point() +
  geom_point(data = pD[Accession %in% resultTables$Interaction[adj.P.Val < 0.05, Accession]]) +
  scale_color_brewer(name = "Significant", type = "qual", palette = "Set2")




pD[(WT_NR - WT_Control) < -1.5 & (KO_NR - KO_Control) > 1 ,]
resultTables$Interaction

tmp <- copy(setup)
tmp[, val:=selectedData["P34927", sample]]
ggplot(tmp, aes(x = group, y = val)) + geom_point(position = position_jitter(width = 0.1))


#this final line does not work
plot(match(resultTables$Interaction$Accession, rownames(selectedData)[indBad]))

#my code for GO
limma_results <- vector(mode = "list", length = 4)
temporary_result <- data.frame(NA, NA, NA, NA)

for (i in 1:4){
  temporary_result <- read.xlsx("C:/Users/tvb217/Documents/R/Lars Proteomic Data/limma_results.xlsx", i)
  limma_results[[i]]<- temporary_result
}

names(limma_results)[[1]]<- "Treatment_in_WT"
names(limma_results)[[2]]<- "Treatment_in_KO"
names(limma_results)[[3]]<- "KO_in_Control"
names(limma_results)[[4]]<- "KO_in_NR"
resultTables <- limma_results

KO_control_list <- resultTables$KO_in_Control
KO_control_list_significant <- KO_control_list %>%
  dplyr::filter(adj.P.Val<0.05) 

KO_list <- KO_control_list_significant[,c(8,2)]
View(KO_control_list)
KO_background <- KO_control_list[,c(8,2)]

eg= bitr(KO_list$Gene, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = "org.Mm.eg.db",
         drop = T)
View(eg)

bg = bitr(KO_background$Gene, 
          fromType = "SYMBOL", 
          toType = "ENTREZID", 
          OrgDb = "org.Mm.eg.db",
          drop = T)
goResults <- enrichGO(gene = eg$ENTREZID,
                      universe = bg$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
enrichplot::dotplot(goResults)
goResults <- setReadable(goResults, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults)

#effect of NR
  NR_control_list <- resultTables$Treatment_in_WT
  NR_control_list_significant <- NR_control_list %>%
    filter(adj.P.Val<0.05) 

NR_WT_enztrez= bitr(NR_control_list_significant$Gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Mm.eg.db",
                    drop = T)


NR_WT_enztrez_bg = bitr(NR_control_list$Gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Mm.eg.db",
                        drop = T)
goResults_NR_WT <- enrichGO(gene = NR_WT_enztrez$ENTREZID,
                            universe = NR_WT_enztrez_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP")
#Treatment in HNKO
NR_HNKO_list <- resultTables$Treatment_in_KO
NR_HNKO_list_significant <- NR_HNKO_list %>%
  filter(adj.P.Val<0.05)

NR_KO_enztrez= bitr(NR_HNKO_list_significant$Gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Mm.eg.db",
                    drop = T)

View(NR_KO_enztrez)
NR_KO_enztrez_bg = bitr(NR_HNKO_list$Gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Mm.eg.db",
                        drop = T)
goResults_NR_KO <- enrichGO(gene = NR_KO_enztrez$ENTREZID,
                            universe = NR_KO_enztrez_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP")
enrichplot::dotplot(goResults_NR_KO)

HNKO_NR_list <- resultTables$KO_in_NR
HNKO_NR_list_significant <- HNKO_NR_list %>%
  filter(adj.P.Val<0.05)
View(HNKO_NR_list_significant)
goResults_NR_KO <- setReadable(goResults_NR_KO, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
cnetplot(goResults_NR_KO)

#Create an upset plot
significant_proteins_upset <- list(Treatment_in_WT = NR_control_list_significant$Gene,
                                   Treatment_in_KO = NR_HNKO_list_significant$Gene,
                                   KO_in_NR = HNKO_NR_list_significant$Gene,
                                   KO_in_control = KO_list$Gene)
View(significant_proteins_upset)
upset(fromList(significant_proteins_upset), order.by = "freq")

#interaction proteins
View(Interaction_list)
Interaction_list <- resultTables$Interaction
Interaction_list_significant <- Interaction_list %>%
  filter(adj.P.Val<0.05)

Interaction_enztrez= bitr(Interaction_list_significant$Gene, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Mm.eg.db",
                    drop = T)

View(Interaction_enztrez)
Interaction_enztrez_bg = bitr(Interaction_list$Gene, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Mm.eg.db",
                        drop = T)
goResults_interaction <- enrichGO(gene = Interaction_enztrez$ENTREZID,
                            universe = Interaction_enztrez_bg$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP")
#No enriched pathways for interaction 
#####Identifying NR response genes#####
NR_overlap <- NR_HNKO_list_significant$Gene %in% NR_control_list_significant$Gene
class(NR_overlap)
NR_overlap <- as.data.frame(NR_overlap)
NR_overlap <- NR_overlap %>%
  mutate(Gene = NR_HNKO_list_significant$Gene)
View(NR_overlap)
NR_overlap_candidates <- NR_overlap %>%
  filter(NR_overlap == T)
View(NR_overlap_candidates)

norm_data <- normalizeBetweenArrays(log(as.matrix(expressions[,-c(1:2)])), method = "quantile")
rownames(norm_data) <- expressions$Gene
candidate_genes <- norm_data
candidate_genes <- as.data.frame(candidate_genes)
candidate_genes <- candidate_genes %>%
  mutate(Name = rownames(norm_data))
View(candidate_genes)

candidate_genes <- candidate_genes %>%
  mutate(Candidate = Name %in% NR_overlap_candidates$Gene)
View(candidate_genes)

candidates_expression <- candidate_genes %>%
  filter(Candidate == T)

class(candidates_expression)
rownames(candidates_expression) <- candidates_expression$Name
View(candidates_expression)
candidates_tidy <- candidates_expression %>%
  pivot_longer(names_to = "sample",
               values_to = "abundance",
               cols= -(Name:Candidate))
candidates_tidy <- candidates_tidy %>%
  dplyr::select(-Candidate)
colnames(candidates_tidy)[1] <- "sample"
View(candidates_tidy)
candidates_tidy$sample <- as.character(candidates_tidy$sample)
setup$sample <- as.character(setup$sample)

joined_candidates <- inner_join(candidates_tidy, setup)
View(joined_candidates)

joined_candidates <- joined_candidates %>%
  dplyr::select(-group)

#generate summary data
summary_candidates <- summarySE(joined_candidates,
                                measurevar = "abundance",
                                groupvars = c("Name", "Genotype", "Treatment"),
                                na.rm = T) %>%
  unite(Group, Genotype, Treatment, sep = "_", remove = F)
summary_candidates <- summary_candidates %>%
  mutate(raw_data = 2^abundance)
View(summary_candidates)
View(expressions)
View(res)
ggplot(summary_candidates,
       aes(x = Group, 
           y = abundance, 
           fill = Group))+
  geom_col()+
  facet_wrap(~Name)

#repeated part of the analysis for raw expression

raw_data <- expressions[,-(2)]
raw_data_candidates <- raw_data %>%
  filter(Gene %in% NR_overlap_candidates$Gene) 
rownames(raw_data_candidates) <- raw_data_candidates$Gene
View(raw_data_candidates)
raw_data_candidates$Gene
raw_data_tidy <- pivot_longer(raw_data_candidates, names_to = "sample",
              values_to = "abundance",
             cols= -Gene)
View(raw_data_tidy)

raw_data_plot <- inner_join(raw_data_tidy, setup)
View(raw_data_plot)
raw_data_plot <- as.data.frame(raw_data_plot)
raw_data_plot <- raw_data_plot %>%
  dplyr::select(-group)

summary_candidates_raw <- summarySE(raw_data_plot,
                                measurevar = "abundance",
                                groupvars = c("Gene", "Genotype", "Treatment"),
                                na.rm = T) %>%
  unite(Group, Genotype, Treatment, sep = "_", remove = F)
View(summary_candidates_raw)

ggplot(summary_candidates_raw,
       aes(x = Group, 
           y = abundance, 
           fill = Group))+
  geom_col()+
  facet_wrap(~Gene)

summary_candidates_raw$Group  <- factor(summary_candidates_raw$Group, levels = c("WT_Control", "WT_NR", "KO_Control", "KO_NR"))

Gcc2 <- 
  ggplot(subset(summary_candidates_raw, Gene == "Gcc2"), 
       aes(x = Group, 
           y = abundance))+
  geom_col()+
  ggtitle("GCC2")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 

Hadh <- ggplot(subset(summary_candidates_raw, Gene == "Hadh"), 
         aes(x = Group, 
             y = abundance))+
  geom_col()+
  ggtitle("Hadh")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 

Nadk <- ggplot(subset(summary_candidates_raw, Gene == "Nadk"), 
       aes(x = Group, 
           y = abundance))+
  geom_col()+
  ggtitle("Nadk")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 

Pah <- ggplot(subset(summary_candidates_raw, Gene == "Pah"), 
              aes(x = Group, 
                  y = abundance))+
  geom_col()+
  ggtitle("Pah")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 

Pebp1 <- ggplot(subset(summary_candidates_raw, Gene == "Pebp1"), 
                aes(x = Group, 
                    y = abundance))+
  geom_col()+
  ggtitle("Pebp1")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 

Safb <- ggplot(subset(summary_candidates_raw, Gene == "Safb"), 
               aes(x = Group, 
                   y = abundance))+
  geom_col()+
  ggtitle("Sabf")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 
Slc47a1 <- ggplot(subset(summary_candidates_raw, Gene == "Slc47a1"), 
       aes(x = Group, 
           y = abundance))+
  geom_col()+
  ggtitle("Slc47a1")+
  geom_errorbar(aes(ymin=abundance-se, 
                    ymax = abundance+se), 
                width = .1) 

grid.arrange(Gcc2, Hadh, Nadk, Pah, Pebp1, Safb, Slc47a1, ncol = 3)
?grid.arrange

