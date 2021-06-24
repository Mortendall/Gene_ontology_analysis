library(here)
library(openxlsx)
library(tidyverse)
library(limma)
library(edgeR)
library(GO.db)
library(org.Mm.eg.db)

cpm_matrix <- openxlsx::read.xlsx(here::here("Sequencing_time_course/cpm_matrix.xlsx"))
setup <- openxlsx::read.xlsx(here::here("Sequencing_time_course/setup.xlsx"))

setup <- setup %>% 
  dplyr::mutate(Group = paste(Genotype, Time, sep = "_"))
setup <- setup %>% 
  dplyr::arrange(Time, desc(Genotype))
setup$ID <- as.character(setup$ID)

colnames(cpm_matrix)<-as.factor(colnames(cpm_matrix))
setup <- setup %>% 
  dplyr::filter(ID %in% colnames(cpm_matrix))

cpm_matrix_key <- cpm_matrix %>% 
  dplyr::select(ENSEMBL, SYMBOL)
cpm_matrix <- cpm_matrix %>% 
  dplyr::select(-ENSEMBL, -SYMBOL)

cpm_matrix <- cpm_matrix %>% 
  dplyr::select(setup$ID)
cpm_matrix <- as.matrix(cpm_matrix)
rownames(cpm_matrix)<-cpm_matrix_key$SYMBOL
key <- as.data.frame(setup$Group)
colnames(key)<-"Group"
rownames(key) <- setup$ID



#####Data loaded, Run analysis for various targets#####
AnnotationDbi::Rkeys(org.Mm.egGO2ALLEGS)<-Go_term
Go_term <- "GO:1903209"

DeepFryer <- function(Go_term, cpm_matrix, setup, org.database){
  #Make design and contrast matrix
  design <- stats::model.matrix(~ 0 + Group, setup)
  colnames(design) <- stringr::str_remove_all(colnames(design), "Group")
  contrasts <- limma::makeContrasts(KO_3 - WT_3, levels = design)
  #extract terms
  # term <- AnnotationDbi::select(GO.db, keys = Go_term, columns = "TERM")
  org.data <- org.Mm.egGO2ALLEGS
  AnnotationDbi::Rkeys(org.data)<-Go_term
  go.genes <- as.list(org.data)
  
  colname_key <- clusterProfiler::bitr(rownames(cpm_matrix), 
                                       fromType = "SYMBOL", 
                                       toType = "ENTREZID", 
                                       OrgDb = "org.Mm.eg.db",
                                       drop = T)
  cpm_entrez <- as.data.frame(cpm_matrix) 
  cpm_entrez <- subset(cpm_entrez, rownames(cpm_entrez)%in%colname_key$SYMBOL)
  cpm_entrez$SYMBOL <- rownames(cpm_entrez)
  cpm_entrez <- left_join(cpm_entrez, colname_key, by ="SYMBOL")
  rownames(cpm_entrez)<-cpm_entrez$ENTREZID
  cpm_entrez <- cpm_entrez %>% 
    dplyr::select(-SYMBOL, -ENTREZID)
  cpm_entrez <- as.matrix(cpm_entrez)
  
  if(all(setup$ID == colnames(cpm_entrez))){
    res <- limma::fry(cpm_entrez, index = go.genes, design = design, contrast = contrasts)
    return(res)
  }
  else{
    return("Your colnames and setup IDs do not match")
  }
  }

#necrosis
DeepFryer("GO:0070265", cpm_matrix = cpm_matrix, setup = setup)

#positive regulation of oxidative stress-induced cell death
DeepFryer("GO:1903209", cpm_matrix = cpm_matrix, setup = setup)

#PPAR signaling
DeepFryer("GO:0035357", cpm_matrix = cpm_matrix, setup = setup)


#autophagy
DeepFryer("GO:0006914", cpm_matrix = cpm_matrix, setup = setup)

#calcium ion binding
DeepFryer("GO:0005509", cpm_matrix = cpm_matrix, setup = setup)

#ECM organization
DeepFryer("GO:0030198", cpm_matrix = cpm_matrix, setup = setup)

#FA metabolism

DeepFryer("GO:0006631", cpm_matrix = cpm_matrix, setup = setup)
