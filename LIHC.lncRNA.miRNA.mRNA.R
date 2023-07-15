setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()
mid.Gene.Type <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol) %>% unique()

Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv")

GeneType.lncRNA <- GeneType %>% filter(Gene.type == "lncRNA")

#### lncipedia LNCipedia database - lncRNA id #####
lncipedia <- read.table("lncipedia_5_2_hc_hg19.gtf.txt",header = F,sep = "\t")
lncipedia2 <- lncipedia %>% dplyr::select(V9)

lncipedia3 <- data.frame()
for (i in 1:dim(lncipedia2)[1]) {
  mid <- (lncipedia2$V9[i] %>% str_split(" ; "))[[1]]
  mid2 <- c()
  for (j in c(1,3,4,9)) {
    id <- (mid[j] %>% str_split(" "))[[1]][2]
    mid2 <- c(mid2, id)
  }
  lncipedia3 <- rbind.data.frame(lncipedia3,mid2)
}
lncipedia4 <- lncipedia3 %>% rename(Gene.ID=1,ENSEMBLE.ID=2,GENE.NAME=3,OTHER=4) %>%
  filter(str_detect(ENSEMBLE.ID,"ENSG000"))
save(lncipedia4,file = "lncipedia.database.Rdata")

lncipedia5 <- lncipedia4 %>% mutate(ENSEMBLE.ID = str_replace_all(ENSEMBLE.ID,"\\..*",""))

#Inter.lncRNA.Up %in% lncipedia5$ENSEMBLE.ID
#Inter.lncRNA.Down %in% lncipedia5$ENSEMBLE.ID
  
######################### TCGA - LIHC - lncRNA + mRNA ####
Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC")

# check whether sample numbers are consistent
length(Phenodata$Sample.ID) # 424
length(colnames(Exp)) # 424

mid.phenodata <- Phenodata %>% filter(Sample.ID %in% colnames(Exp)) %>% 
  filter(Sample.Type == "Solid Tissue Normal" | Sample.Type == "Primary Tumor") %>%
  mutate(Group = case_when(Sample.Type == "Solid Tissue Normal" ~ "Ctrl",
                           Sample.Type == "Primary Tumor" ~ "Tumor")) %>%
  mutate(Group = factor(.$Group,levels = c("Ctrl","Tumor")))

Exp.LIHC <- Exp[,mid.phenodata$Sample.ID]
colData <- mid.phenodata %>% select(Sample.ID,Group) %>% mutate(Sample.ID = factor(.$Sample.ID,levels = colnames(Exp.LIHC))) %>%
  arrange(Sample.ID) %>% unique() %>% remove_rownames() %>% column_to_rownames("Sample.ID")

suppressWarnings(library(DESeq2))
suppressWarnings(library(ggthemes))
dds <- DESeqDataSetFromMatrix(countData = Exp.LIHC, colData = colData, design = ~Group) 
dds <- DESeq(dds)
res <- results(dds,contrast=c("Group","Tumor","Ctrl"))
res <- res[order(res$pvalue),] %>% data.frame() %>% 
  mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange >=1, "Up",ifelse(padj <= 0.01 & log2FoldChange <=-1,"Down","No")))

save(res,file = paste("TCGA-LIHC","_DEG.Rdata",sep = ""))
write.csv(res,file = paste("TCGA-LIHC","_DEG.csv",sep = ""))
res <- read.csv(paste("TCGA-LIHC","_DEG.csv",sep = ""),row.names = 1)

#### lncRNA ####

LIHC.lncRNA <- res[rownames(res) %in% GeneType.lncRNA$Gene.stable.ID,] %>% data.frame()
#LIHC.lncRNA["ENSG00000215424",]
dim(LIHC.lncRNA %>% filter(Threshold != "No"))


LIHC.lncRNA.Up <- LIHC.lncRNA %>% filter(Threshold == "Up")
LIHC.lncRNA.Down <- LIHC.lncRNA %>% filter(Threshold == "Down")

#### mRNA - PCG ####
GeneType.mRNA.PCG <- GeneType %>% filter(Gene.type == "protein_coding")
LIHC.mRNA.PCG <- res[rownames(res) %in% GeneType.mRNA.PCG$Gene.stable.ID,] %>% data.frame()

LIHC.mRNA.PCG.Up <- LIHC.mRNA.PCG%>% filter(Threshold == "Up")
LIHC.mRNA.PCG.Down <- LIHC.mRNA.PCG %>% filter(Threshold == "Down")

LIHC.RNA <- list(LIHC.lncRNA.Up=LIHC.lncRNA.Up,LIHC.lncRNA.Down=LIHC.lncRNA.Down,LIHC.mRNA.PCG.Up=LIHC.mRNA.PCG.Up,LIHC.mRNA.PCG.Down=LIHC.mRNA.PCG.Down)
saveRDS(LIHC.RNA,"LIHC.RNA.rds")


######################### TCGA - LIHC - primary miRNA ######################
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)

Exp <- read.csv("../../miRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.miRNA.Exp.csv",check.names = F,row.names = 1)
Phenodata <- read.csv("../../miRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC")

length(Phenodata$Sample.ID) # 425
length(colnames(Exp)) # 425

mid.phenodata <- Phenodata %>% filter(Sample.ID %in% colnames(Exp)) %>% 
  filter(Sample.Type == "Solid Tissue Normal" | Sample.Type == "Primary Tumor") %>%
  mutate(Group = case_when(Sample.Type == "Solid Tissue Normal" ~ "Ctrl",
                           Sample.Type == "Primary Tumor" ~ "Tumor")) %>%
  mutate(Group = factor(.$Group,levels = c("Ctrl","Tumor")))

Exp.LIHC <- Exp[,mid.phenodata$Sample.ID]
colData <- mid.phenodata %>% select(Sample.ID,Group) %>% mutate(Sample.ID = factor(.$Sample.ID,levels = colnames(Exp.LIHC))) %>%
  arrange(Sample.ID) %>% unique() %>% remove_rownames() %>% column_to_rownames("Sample.ID")

suppressWarnings(library(DESeq2))
suppressWarnings(library(ggthemes))
dds <- DESeqDataSetFromMatrix(countData = Exp.LIHC, colData = colData, design = ~Group) 
dds <- DESeq(dds)
res <- results(dds,contrast=c("Group","Tumor","Ctrl"))
res <- res[order(res$pvalue),] %>% data.frame() %>% 
  mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange >=1, "Up",ifelse(padj <= 0.01 & log2FoldChange <=-1,"Down","No")))

save(res,file = paste("TCGA-LIHC","_Primary.miRNA_DE.Rdata",sep = ""))
write.csv(res,file = paste("TCGA-LIHC","_Primary.miRNA_DE.csv",sep = ""))

######################### TCGA - LIHC - mature miRNA ######################
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)

Exp <- read.csv("../../miRNA.Mature.PanCancer.Exp/PanCancer.TCGA-LIHC.miRNA.Exp.csv",check.names = F)
Phenodata <- read.csv("../../miRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC")

#BiocManager::install("miRBaseVersions.db")
library(miRBaseVersions.db)
#select(miRBaseVersions.db, keys = "MIMAT0000062", keytype = "MIMAT", columns = "*")
#items = select(miRBaseVersions.db, keys = c("MIMAT0000062", "MIMAT0000063","MIMAT0000064","MIMAT0000065","MIMAT0000066"), keytype = "MIMAT", columns = "*")
#res = items[items$VERSION == 22.0, "NAME"]

#Exp$miRNA_region %>% str_remove_all(".*\\,")
items = select(miRBaseVersions.db, keys = unique(Exp$miRNA_region %>% str_remove_all(".*\\,")), keytype = "MIMAT", columns = "*")
mature.miRNA = items[items$VERSION == 22.0,] %>% data.frame()
save(mature.miRNA,file = "mature.miRNA.Rdata")

Exp.mid <- Exp %>% mutate(mature.miRNA = (miRNA_region %>% str_remove_all(".*\\,"))) %>% 
  dplyr::select(-miRNA_region) %>% merge(.,mature.miRNA %>% dplyr::select(ACCESSION,NAME),by.x="mature.miRNA",by.y="ACCESSION") %>%
  mutate(rowname = str_c(miRNA_ID,NAME,sep = ".")) %>%
  dplyr::select(-miRNA_ID,-mature.miRNA,-NAME) %>% remove_rownames() %>% column_to_rownames("rowname")


mid.phenodata <- Phenodata %>% filter(Sample.ID %in% colnames(Exp)) %>% 
  filter(Sample.Type == "Solid Tissue Normal" | Sample.Type == "Primary Tumor") %>%
  mutate(Group = case_when(Sample.Type == "Solid Tissue Normal" ~ "Ctrl",
                           Sample.Type == "Primary Tumor" ~ "Tumor")) %>%
  mutate(Group = factor(.$Group,levels = c("Ctrl","Tumor")))

# 422 372 Primary Tumor  3 Recurrent Tumor 50 Solid Tissue Normal 

Exp.LIHC <- Exp.mid[,mid.phenodata$Sample.ID]
colData <- mid.phenodata %>% dplyr::select(Sample.ID,Group) %>% mutate(Sample.ID = factor(.$Sample.ID,levels = colnames(Exp.LIHC))) %>%
  dplyr::arrange(Sample.ID) %>% unique() %>% remove_rownames() %>% column_to_rownames("Sample.ID")

Exp.LIHC <- as.matrix(Exp.LIHC)
Exp.LIHC[is.na(Exp.LIHC)] = 0

suppressWarnings(library(DESeq2))
suppressWarnings(library(ggthemes))
dds <- DESeqDataSetFromMatrix(countData = Exp.LIHC, colData = colData, design = ~Group) 
dds <- DESeq(dds)
res <- results(dds,contrast=c("Group","Tumor","Ctrl"))
res <- res[order(res$pvalue),] %>% data.frame() %>% 
  mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange >=1, "Up",ifelse(padj <= 0.01 & log2FoldChange <=-1,"Down","No")))

save(res,file = paste("TCGA-LIHC","_Mature.miRNA_DE.Rdata",sep = ""))
write.csv(res,file = paste("TCGA-LIHC","_Mature.miRNA_DE.csv",sep = ""))

######################### TCGA - GSE140845 (PRJNA591153) Total RNA  - > Intersect ###############
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)
### GSE140845
GSE140845.RNA <- read.csv("../../../Metastasis/TissueMicrobiome/PRJNA591152-Bact+Virus/mRNA.lncRNA.csv")
View(GSE140845.RNA)

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.lncRNA <- GeneType %>% filter(Gene.type == "lncRNA")
GeneType.mRNA.PCG <- GeneType %>% filter(Gene.type == "protein_coding")

GSE140845.RNA.lncRNA <- GSE140845.RNA %>% filter(X %in% GeneType.lncRNA$Gene.stable.ID) %>% merge(.,GeneType.lncRNA,by.x="X",by.y = "Gene.stable.ID") %>% data.frame()
GSE140845.RNA.PCG <- GSE140845.RNA %>% filter(X %in% GeneType.mRNA.PCG$Gene.stable.ID) %>% merge(.,GeneType.lncRNA,by.x="X",by.y = "Gene.stable.ID") %>% data.frame()

save(GSE140845.RNA.lncRNA,file = "GSE140845.RNA.lncRNA.Rdata")
save(GSE140845.RNA.PCG,file = "GSE140845.RNA.PCG.Rdata")

### LIHC
load(file = "GSE140845.RNA.lncRNA.Rdata") #GSE140845.RNA.lncRNA
GSE140845.RNA.lncRNA.Up <- GSE140845.RNA.lncRNA %>% filter(Sig == "Up")
GSE140845.RNA.lncRNA.Down <- GSE140845.RNA.lncRNA %>% filter(Sig == "Down")

load(file = "TCGA-LIHC_DEG.Rdata") #res

LIHC.lncRNA <- res[rownames(res) %in% GeneType.lncRNA$Gene.stable.ID,] %>% data.frame()
save(LIHC.lncRNA,file = "LIHC.lncRNA.Rdata") #LIHC.lncRNA

LIHC.lncRNA.Up <- LIHC.lncRNA %>% filter(Threshold == "Up")
LIHC.lncRNA.Down <- LIHC.lncRNA %>% filter(Threshold == "Down")

Inter.lncRNA.Up <- intersect(unique(GSE140845.RNA.lncRNA.Up$X),rownames(LIHC.lncRNA.Up)) #28
Inter.lncRNA.Down <- intersect(unique(GSE140845.RNA.lncRNA.Down$X),rownames(LIHC.lncRNA.Down)) #24

save(GSE140845.RNA.lncRNA.Up,GSE140845.RNA.lncRNA.Down,LIHC.lncRNA.Up,LIHC.lncRNA.Down,Inter.lncRNA.Up,Inter.lncRNA.Down,file = "LIHC.GSE140845.lncRNA.Rdata")

load("LIHC.GSE140845.lncRNA.Rdata")
lncRNA <- GSE140845.RNA.lncRNA.Up[,c(1:9,11)] %>% unique() %>% data.frame() %>% mutate(Cohort="GSE140845") %>% 
  rbind(GSE140845.RNA.lncRNA.Down[,c(1:9,11)] %>% unique() %>% data.frame() %>% mutate(Cohort="GSE140845")) %>%
  dplyr::select(X,logFC,Sig,Gene.type,Gene.name,Cohort) %>%
  magrittr::set_colnames(c("X","log2FoldChange","Threshold","Gene.type","Gene.name","Cohort")) %>%
  rbind(LIHC.lncRNA.Up %>% rownames_to_column("X") %>%
          merge(mid.Gene.Type,by.x="X",by.y="Gene.stable.ID") %>%
          dplyr::select(-HGNC.symbol,-baseMean,-lfcSE,-stat,-pvalue,-padj) %>%
          unique() %>% data.frame() %>% mutate(Cohort="TCGA-LIHC",Gene.type="lncRNA")) %>%
  rbind(LIHC.lncRNA.Down %>% rownames_to_column("X") %>%
          merge(mid.Gene.Type,by.x="X",by.y="Gene.stable.ID") %>%
          dplyr::select(-HGNC.symbol,-baseMean,-lfcSE,-stat,-pvalue,-padj) %>%
          unique() %>% data.frame() %>% mutate(Cohort="TCGA-LIHC",Gene.type="lncRNA"))
colnames(lncRNA)[1]="ID"

load("LIHC.GSE140845.miRNA.Rdata")
miRNA <- GSE140845.miRNA.Up %>% rownames_to_column("X") %>% unique() %>% data.frame() %>% mutate(Gene.type="miRNA",Gene.name=X,Cohort="GSE140845") %>% 
  rbind(GSE140845.miRNA.Down %>% rownames_to_column("X") %>% unique() %>% data.frame() %>% mutate(Gene.type="miRNA",Gene.name=X,Cohort="GSE140845")) %>%
  dplyr::select(X,log2FoldChange,Sig,Gene.type,Gene.name,Cohort) %>%
  magrittr::set_colnames(c("Mature.miRNA","log2FoldChange","Threshold","Gene.type","Gene.name","Cohort")) %>%
  rbind(LIHC.miRNA.Up %>%
          dplyr::select(-baseMean,-lfcSE,-stat,-pvalue,-padj,-Precursor) %>%
          unique() %>% data.frame() %>% mutate(Cohort="TCGA-LIHC",Gene.type="miRNA",Gene.name=Mature.miRNA)) %>%
  rbind(LIHC.miRNA.Down %>%
          dplyr::select(-baseMean,-lfcSE,-stat,-pvalue,-padj,-Precursor) %>%
          unique() %>% data.frame() %>% mutate(Cohort="TCGA-LIHC",Gene.type="miRNA",Gene.name=Mature.miRNA))
colnames(miRNA)[1]="ID"

load("LIHC.GSE140845.mRNA.Rdata")
mRNA <- GSE140845.mRNA.Up %>% 
  merge(mid.Gene.Type,by.x="X",by.y="Gene.stable.ID") %>% 
  dplyr::select(-HGNC.symbol) %>%
  unique() %>% data.frame() %>% mutate(Cohort="GSE140845",Gene.type="mRNA") %>% 
  rbind(GSE140845.mRNA.Down %>% merge(mid.Gene.Type,by.x="X",by.y="Gene.stable.ID") %>% dplyr::select(-HGNC.symbol) %>%
          unique() %>% data.frame() %>% mutate(Cohort="GSE140845",Gene.type="mRNA")) %>%
  dplyr::select(X,logFC,Sig,Gene.type,Gene.name,Cohort) %>%
  magrittr::set_colnames(c("X","log2FoldChange","Threshold","Gene.type","Gene.name","Cohort")) %>%
  rbind(LIHC.mRNA.PCG.Up %>% rownames_to_column("X") %>%
          merge(mid.Gene.Type,by.x="X",by.y="Gene.stable.ID") %>%
          dplyr::select(-HGNC.symbol,-baseMean,-lfcSE,-stat,-pvalue,-padj) %>%
          unique() %>% data.frame() %>% mutate(Cohort="TCGA-LIHC",Gene.type="mRNA")) %>%
  rbind(LIHC.mRNA.PCG.Down %>% rownames_to_column("X") %>%
          merge(mid.Gene.Type,by.x="X",by.y="Gene.stable.ID") %>%
          dplyr::select(-HGNC.symbol,-baseMean,-lfcSE,-stat,-pvalue,-padj) %>%
          unique() %>% data.frame() %>% mutate(Cohort="TCGA-LIHC",Gene.type="mRNA"))
colnames(mRNA)[1]="ID"

data.F <- rbind(lncRNA,miRNA) %>%
  rbind(mRNA)
write.csv(data.F,file = "Supplementary table.csv",row.names = F)

#unique((GSE140845.RNA.lncRNA.Up %>% filter(X %in% Inter.lncRNA.Up))$Gene.name)
#unique((GSE140845.RNA.lncRNA.Down %>% filter(X %in% Inter.lncRNA.Down))$Gene.name)
#GSE140845.RNA.lncRNA.Sig %>% filter(X == "ENSG00000215424")
#LIHC.lncRNA.Up["ENSG00000215424",]
load("LIHC.GSE140845.lncRNA.Rdata")
library(ggvenn)
library("RColorBrewer")
X <- list(GSE140845.Up = GSE140845.RNA.lncRNA.Up$X,
          LIHC.Up = rownames(LIHC.lncRNA.Up),
          LIHC.Down=rownames(LIHC.lncRNA.Down),
          GSE140845.Down = GSE140845.RNA.lncRNA.Down$X)
#saveRDS(X,"LIHC.GSE140845.lncRNA.rds")

p3 <- ggvenn(X,set_name_size = 3,text_size =2,fill_color = brewer.pal(8,"Set1"))
ggsave(p3,filename = "LIHC.GSE140845.lncRNA.Venn.pdf",height = 3,width = 5)

######################### KM Cox Analyais for intersected lncRNA ###################
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)
LIHC.GSE140845.lncRNA <- readRDS("LIHC.GSE140845.lncRNA.rds")

Inter.lncRNA.Up <- intersect(unique(LIHC.GSE140845.lncRNA$GSE140845.Up$X),rownames(LIHC.GSE140845.lncRNA$LIHC.Up)) #28
Inter.lncRNA.Down <- intersect(unique(LIHC.GSE140845.lncRNA$GSE140845.Down$X),rownames(LIHC.GSE140845.lncRNA$LIHC.Down)) #24
#filter primary tumor samples => 371
Pancancer.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC") %>%
  filter(Sample.Type == "Primary Tumor")

#XENA data phenotype and survival datasets
PhenoData <- read.table("TCGA.LIHC.sampleMap_LIHC_clinicalMatrix.txt",sep ="\t",row.names = 1,header = T)
TPMdata <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
TPMdata <- TPMdata[,Pancancer.Phenodata$Sample.ID]

SurvivalData <- read.table("survival_LIHC_survival.txt",header = T,row.names = 1,sep = '\t') 
View(SurvivalData)

#substr(colnames(TPMdata), 0, 15) %in% PhenoData$X_INTEGRATION
######### P:Overall survival -> single factor KM - Up.lncRNA #######
mid.Survival <- SurvivalData[substr(colnames(TPMdata), 0, 15),]
lncRNA.Up.Survival.Data <- TPMdata[Inter.lncRNA.Up,] %>% t() %>% data.frame() %>% rownames_to_column() %>%
  mutate(rowname = substr(.$rowname,0,15)) %>% mutate(rowname = factor(.$rowname,levels = rownames(mid.Survival))) %>%
  arrange(rowname) %>% remove_rownames() %>%
  column_to_rownames("rowname") %>% cbind(mid.Survival)

saveRDS(lncRNA.Up.Survival.Data,"lncRNA.Up.Survival.Data.rds")

lncRNA.Up.Data <- lncRNA.Up.Survival.Data %>% select(starts_with("ENS"))
save(lncRNA.Up.Data,mid.Survival,file = "lncRNA.Up.Survival.Data.Rdata")

library(survival)
library(survminer)

log_rank_p <- apply(lncRNA.Up.Data , 2 , function(gene){
  # gene=exprSet[1,]
  mid.Survival$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(Surv(OS.time, OS)~group,data=mid.Survival)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

if (F) {
  log_rank_p2 <- apply(lncRNA.Up.Data[,1:2] , 2 , function(gene){
    #gene <- lncRNA.Up.Data[,1]
    mid.Survival$group=ifelse(gene>median(gene),'high','low')  
    data.survdiff=survdiff(Surv(OS.time, OS)~group,data=mid.Survival)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    return(p.val)
  })
  
}


lncRNA.Up.OS.KM.logrank.P <- data.frame(log_rank_p) %>% arrange(log_rank_p) 

lncRNA.Up.mid.gene.type1 <- GeneType %>% filter(Gene.stable.ID %in% rownames(lncRNA.Up.OS.KM.logrank.P)) %>% 
  dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol) %>% unique() %>%
  merge(.,lncRNA.Up.OS.KM.logrank.P %>% rownames_to_column(),by.x = "Gene.stable.ID",by.y="rowname")

write.csv(lncRNA.Up.mid.gene.type1,"lncRNA.Up.mid.gene.type,csv",row.names = F)

lncRNA.Up.OS.KM.logrank.P.Sig <- lncRNA.Up.OS.KM.logrank.P  %>% filter(log_rank_p <= 0.05) %>% rownames_to_column()

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
lncRNA.Up.mid.gene.type <- GeneType %>% filter(Gene.stable.ID %in% lncRNA.Up.OS.KM.logrank.P.Sig$rowname) %>% 
  dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol) %>% unique() %>%
  merge(.,lncRNA.Up.OS.KM.logrank.P.Sig,by.x = "Gene.stable.ID",by.y="rowname")

save(lncRNA.Up.OS.KM.logrank.P.Sig,lncRNA.Up.mid.gene.type,file = "lncRNA.Up.OS.KM.Up.logrank.P.Sig.Rdata")

meta <- mid.Survival
library(export)
library(ggthemes)

dir.create("lncRNA.Up.OS.KM.Cox")
for (i in 1:dim(lncRNA.Up.mid.gene.type)[1]) {
  #gene.name <- GeneType %>% filter(Gene.stable.ID == gene) %>% dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol)
  if (lncRNA.Up.mid.gene.type$Gene.name[i] != "") {
    meta[[lncRNA.Up.mid.gene.type$Gene.stable.ID[i]]] = ifelse(lncRNA.Up.Data[,lncRNA.Up.mid.gene.type$Gene.stable.ID[i]]> median(lncRNA.Up.Data[,lncRNA.Up.mid.gene.type$Gene.stable.ID[i]]),'high','low')
    sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~",lncRNA.Up.mid.gene.type$Gene.stable.ID[i],sep = '')), data=meta)
    #pdf(paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 6,width = 6)
    #ggsave(filename = paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 4,width = 4)
    p<-ggsurvplot(sfit1,pval =TRUE, data = meta, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event"ÂàôÂèØÁªòÂà∂cumulative eventsÂõæÔºàfunÈô§‰∫ÜevenÂ§ñÔºåËøòÊúâlogÔºölog transformationÂíåcumhazÔºöcumulative hazard???
                  legend.title = lncRNA.Up.mid.gene.type$Gene.name[i],
                  legend.labs = c("High", "Low"),
                  conf.int.style = "step",
                  xlab = "Time in days",
                  risk.table = "abs_pct",
                  risk.table.y.text.col = T,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  ggtheme = theme_few())
    print(p)
    #dev.off()
    graph2pdf(file=paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 6,width = 5)
    
    p<-ggsurvplot(sfit1,pval =TRUE, data = meta, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event"ÂàôÂèØÁªòÂà∂cumulative eventsÂõæÔºàfunÈô§‰∫ÜevenÂ§ñÔºåËøòÊúâlogÔºölog transformationÂíåcumhazÔºöcumulative hazard???
                  legend.title = lncRNA.Up.mid.gene.type$Gene.name[i],
                  legend.labs = c("High", "Low"),
                  conf.int.style = "step",
                  xlab = "Time in days",
                  risk.table = "abs_pct",
                  risk.table.y.text.col = T,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  fun = "cumhaz",
                  ggtheme = theme_few())
    print(p)
    #dev.off()
    graph2pdf(file=paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.cumhaz.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 6,width = 5)
    
    #ggsave(p,filename = paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 6,width = 4)
    #print(unique(gene.name$Gene.name))
  }else{
    meta[[lncRNA.Up.mid.gene.type$Gene.stable.ID[i]]] = ifelse(lncRNA.Up.Data[,lncRNA.Up.mid.gene.type$Gene.stable.ID[i]]> median(lncRNA.Up.Data[,lncRNA.Up.mid.gene.type$Gene.stable.ID[i]]),'high','low')
    sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~",lncRNA.Up.mid.gene.type$Gene.stable.ID[i],sep = '')), data=meta)
    #pdf(paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.stable.ID[i],".pdf",sep = ''),height = 6,width = 6)
    #ggsave(filename = paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.stable.ID[i],".pdf",sep = ''),height = 4,width = 4)
    p<-ggsurvplot(sfit1,pval =TRUE, data = meta, 
                  surv.median.line = "hv",
                  legend.title = lncRNA.Up.mid.gene.type$Gene.stable.ID[i],
                  conf.int.style = "step",
                  xlab = "Time in days",
                  #break.time.by = 500,
                  risk.table = "abs_pct",
                  risk.table.y.text.col = T,
                  risk.table.y.text = FALSE,
                  legend.labs = c("High", "Low"),
                  #pval = TRUE,
                  conf.int = TRUE,
                  palette = "Set1",
                  ggtheme = theme_bw())
    print(p)
    #dev.off()
    graph2pdf(file=paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.stable.ID[i],".pdf",sep = ''),height = 6,width = 5)
    
    p<-ggsurvplot(sfit1,pval =TRUE, data = meta, 
                  surv.median.line = "hv",
                  legend.title = lncRNA.Up.mid.gene.type$Gene.stable.ID[i],
                  conf.int.style = "step",
                  xlab = "Time in days",
                  #break.time.by = 500,
                  risk.table = "abs_pct",
                  risk.table.y.text.col = T,
                  risk.table.y.text = FALSE,
                  legend.labs = c("High", "Low"),
                  #pval = TRUE,
                  conf.int = TRUE,
                  palette = "Set1",
                  fun = "cumhaz",
                  ggtheme = theme_bw())
    print(p)
    #dev.off()
    graph2pdf(file=paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.cumhaz.",lncRNA.Up.mid.gene.type$Gene.stable.ID[i],".pdf",sep = ''),height = 6,width = 5)
  }
}

######### P:Overall survival -> Cox - Up.lncRNA #########
#Cox ÂõûÂΩíÁöÑÈáçË¶ÅÁªüËÆ°ÊåáÊ†áÔºö È£éÈô©ÊØîÔºà hazard ratio???
#??? ??? HR>1 Êó∂ÔºåËØ¥ÊòéÁ†îÁ©∂ÂØπË±°ÊòØ‰∏Ä‰∏™Âç±Èô©Âõ†Á¥????
#??? ??? HR<1 Êó∂ÔºåËØ¥ÊòéÁ†îÁ©∂ÂØπË±°ÊòØ‰∏Ä‰∏™‰øùÊä§Âõ†Á¥????
#??? ??? HR=1 Êó∂ÔºåËØ¥ÊòéÁ†îÁ©∂ÂØπË±°ÂØπÁîüÂ≠òÊó∂Èó¥‰∏çËµ∑‰ΩúÁî????

lnc.RNA.Up.cox_results <-apply(lncRNA.Up.Data , 2 , function(gene){
  #gene= lncRNA.Up.Data[,1]
  meta$gene = gene
  #ÂèØÁõ¥Êé•‰ΩøÁî®ËøûÁª≠ÂûãÂèòÈáè
  m = coxph(Surv(OS.time, OS) ~ gene, data =  meta)
  #‰πüÂèØ‰ΩøÁî®‰∫åÂàÜÁ±ªÂèò???
  #meta$group=ifelse(gene>median(gene),'high','low') 
  #m=coxph(Surv(time, event) ~ group, data =  meta)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, 
                     se = se, z = beta/se, 
                     p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, 
                     HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)), 6)
  
  return(tmp['gene',]) 
  #return(tmp['grouplow',])#‰∫åÂàÜÁ±ªÂèò???
})

lnc.RNA.Up.cox_results <- lnc.RNA.Up.cox_results %>% t() %>% data.frame()
lnc.RNA.Up.cox_results %>% rownames_to_column() %>% 
  merge.data.frame(.,Inter.lncRNA.GeneCard,by.x="rowname",by.y="lncRNA.ENSEMBLE.ID") %>% write.csv(.,file="lnc.RNA.Up.cox_results.csv",row.names=F)

save(lnc.RNA.Up.cox_results,file = "lnc.RNA.Up.cox_results.Rdata")
saveRDS(lnc.RNA.Up.cox_results,"lnc.RNA.Up.cox_results.rds")

lnc.RNA.Up.cox_results.Sig <- lnc.RNA.Up.cox_results %>% 
  filter(rownames(.) %in% lncRNA.Up.mid.gene.type$Gene.stable.ID) %>% rownames_to_column() %>%
  merge(.,lncRNA.Up.mid.gene.type,by.x = "rowname",by.y = "Gene.stable.ID")
  
saveRDS(lnc.RNA.Up.cox_results.Sig,"lnc.RNA.Up.cox_results.Sig.rds")


######### P:Overall survival -> KM - Down.lncRNA ########

mid.Survival <- SurvivalData[substr(colnames(TPMdata), 0, 15),]
lncRNA.Down.Survival.Data <- TPMdata[Inter.lncRNA.Down,] %>% t() %>% data.frame() %>% rownames_to_column() %>%
  mutate(rowname = substr(.$rowname,0,15)) %>% mutate(rowname = factor(.$rowname,levels = rownames(mid.Survival))) %>%
  arrange(rowname) %>% remove_rownames() %>%
  column_to_rownames("rowname") %>% cbind(mid.Survival)

saveRDS(lncRNA.Down.Survival.Data,"lncRNA.Down.Survival.Data.rds")
lncRNA.Down.Survival.Data <- readRDS("lncRNA.Down.Survival.Data.rds")

lncRNA.Down.Data <- lncRNA.Down.Survival.Data %>% select(starts_with("ENS"))
save(lncRNA.Down.Data,mid.Survival,file = "lncRNA.Down.Survival.Data.Rdata")
load("lncRNA.Down.Survival.Data.Rdata")

library(survival)
library(survminer)

log_rank_p <- apply(lncRNA.Down.Data , 2 , function(gene){
  # gene=exprSet[1,]
  mid.Survival$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(Surv(OS.time, OS)~group,data=mid.Survival)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

lncRNA.Down.OS.KM.logrank.P <- data.frame(log_rank_p) %>% arrange(log_rank_p) 
lncRNA.Down.OS.KM.logrank.P %>% rownames_to_column() %>% merge.data.frame(.,Inter.lncRNA.GeneCard,by.x="rowname",by.y="lncRNA.ENSEMBLE.ID") %>%
  write.csv(.,file = "lncRNA.Down.OS.KM.logrank.P.csv",row.names = F)

lncRNA.Down.OS.KM.logrank.P.Sig <- lncRNA.Down.OS.KM.logrank.P %>% filter(log_rank_p <= 0.05) %>% rownames_to_column()

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
lncRNA.Down.mid.gene.type <- GeneType %>% filter(Gene.stable.ID %in% lncRNA.Down.OS.KM.logrank.P.Sig$rowname) %>% 
  dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol) %>% unique() %>%
  merge(.,lncRNA.Down.OS.KM.logrank.P.Sig,by.x = "Gene.stable.ID",by.y="rowname")

save(lncRNA.Down.OS.KM.logrank.P.Sig,lncRNA.Down.mid.gene.type,file = "lncRNA.Down.OS.KM.Down.logrank.P.Sig.Rdata")

meta <- mid.Survival
library(export)
library(ggthemes)

dir.create("lncRNA.Down.OS.KM.Cox")
for (i in 1:dim(lncRNA.Down.mid.gene.type)[1]) {
  #gene.name <- GeneType %>% filter(Gene.stable.ID == gene) %>% dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol)
  if (lncRNA.Down.mid.gene.type$Gene.name[i] != "") {
    meta[[lncRNA.Down.mid.gene.type$Gene.stable.ID[i]]] = ifelse(lncRNA.Down.Data[,lncRNA.Down.mid.gene.type$Gene.stable.ID[i]]> median(lncRNA.Down.Data[,lncRNA.Down.mid.gene.type$Gene.stable.ID[i]]),'high','low')
    sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~",lncRNA.Down.mid.gene.type$Gene.stable.ID[i],sep = '')), data=meta)
    P<-ggsurvplot(sfit1,pval =TRUE, data = meta, #risk.table = TRUE,
                  surv.median.line = "hv",
                  conf.int.style = "step",
                  xlab = "Time in days",
                  risk.table = "abs_pct",
                  risk.table.y.text.col = T,
                  risk.table.y.text = FALSE,
               legend.title = lncRNA.Down.mid.gene.type$Gene.name[i],
               legend.labs = c("High", "Low"),
               conf.int = TRUE,
               palette = "Set1",
               ggtheme = theme_few())
    
    print(p)
    graph2pdf(file=paste("lncRNA.Down.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Down.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 6,width = 5)
    #ggsave(p,filename = paste("lncRNA.Down.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Down.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 4,width = 4)
    print(unique(gene.name$Gene.name))
  }else{
    meta[[lncRNA.Down.mid.gene.type$Gene.stable.ID[i]]] = ifelse(lncRNA.Down.Data[,lncRNA.Down.mid.gene.type$Gene.stable.ID[i]]> median(lncRNA.Down.Data[,lncRNA.Down.mid.gene.type$Gene.stable.ID[i]]),'high','low')
    sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~",lncRNA.Down.mid.gene.type$Gene.stable.ID[i],sep = '')), data=meta)
    p<-ggsurvplot(sfit1,pval =TRUE, data = meta, #risk.table = TRUE,
                  surv.median.line = "hv",
               legend.title = lncRNA.Down.mid.gene.type$Gene.stable.ID[i],
               legend.labs = c("High", "Low"),
               conf.int.style = "step",
               xlab = "Time in days",
               risk.table = "abs_pct",
               risk.table.y.text.col = T,
               risk.table.y.text = FALSE,
               #pval = TRUE,
               conf.int = TRUE,
               palette = "Set1",
               ggtheme = theme_bw())
    print(p)
    graph2pdf(file=paste("lncRNA.Down.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Down.mid.gene.type$Gene.stable.ID[i],".pdf",sep = ''),height = 6,width = 5)
  }
}




######### P:Overall survival -> single factor Cox - Down.lncRNA #########
lnc.RNA.Down.cox_results <-apply(lncRNA.Down.Data , 2 , function(gene){
  #gene= lncRNA.Down.Data[,1]
  meta$gene = gene
  #ÂèØÁõ¥Êé•‰ΩøÁî®ËøûÁª≠ÂûãÂèòÈáè
  m = coxph(Surv(OS.time, OS) ~ gene, data =  meta)
  #‰πüÂèØ‰ΩøÁî®‰∫åÂàÜÁ±ªÂèò???
  #meta$group=ifelse(gene>median(gene),'high','low') 
  #m=coxph(Surv(time, event) ~ group, data =  meta)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, 
                     se = se, z = beta/se, 
                     p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, 
                     HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)), 6)
  
  return(tmp['gene',]) 
  #return(tmp['grouplow',])#‰∫åÂàÜÁ±ªÂèò???
})

lnc.RNA.Down.cox_results <- lnc.RNA.Down.cox_results %>% t() %>% data.frame()
lnc.RNA.Down.cox_results %>% rownames_to_column() %>% merge.data.frame(.,Inter.lncRNA.GeneCard,by.x="rowname",by.y="lncRNA.ENSEMBLE.ID") %>%
  write.csv(.,file = "lncRNA.Down.OS.cox_results.csv",row.names = F)

lnc.RNA.Down.cox_results %>% rownames_to_column() %>% 
  merge.data.frame(.,Inter.lncRNA.GeneCard,by.x="rowname",by.y="lncRNA.ENSEMBLE.ID") %>% write.csv(.,file="lnc.RNA.Down.cox_results.csv",row.names=F)


save(lnc.RNA.Down.cox_results,file = "lnc.RNA.Down.cox_results.Rdata")
saveRDS(lnc.RNA.Down.cox_results,"lnc.RNA.Down.cox_results.rds")

lnc.RNA.Down.cox_results.Sig <- lnc.RNA.Down.cox_results %>% filter(p <= 0.05) %>% 
  filter(rownames(.) %in% lncRNA.Down.mid.gene.type$Gene.stable.ID) %>% rownames_to_column() %>%
  merge(.,lncRNA.Down.mid.gene.type,by.x = "rowname",by.y = "Gene.stable.ID")

saveRDS(lnc.RNA.Down.cox_results,"lnc.RNA.Down.cox_results.Sig.rds")

######### P:Overall survival -> multi factor Cox - DE COX.lncRNA #########
### down -> no; up -> lnc.RNA.Up.cox_results.Sig
mid.Survival <- SurvivalData[substr(colnames(TPMdata), 0, 15),] #TPMdata only contain the primary samples
Exp.data <- TPMdata[lnc.RNA.Up.cox_results.Sig$rowname,] %>% t() %>% data.frame() %>% 
  rownames_to_column() %>% mutate(rowname = substr(rowname, 0, 15)) %>%
  mutate(rowname = factor(.$rowname,levels = c(rownames(mid.Survival)))) %>%
  arrange(rowname)

Exp.data.Surv <- cbind(Exp.data,mid.Survival) %>% data.frame()
saveRDS(Exp.data.Surv,"lncRNA.Multi.Cox.Exp.data.Surv.rds")

library(survival)
library(survminer)
library(My.stepwise)

Input <- as.formula(paste("Surv(OS.time, OS == 1) ~ ",paste(lnc.RNA.Up.cox_results.Sig$rowname,collapse = "+"),sep = ""))

res.cox <- coxph(Input, data =  Exp.data.Surv)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
write.table(file="lncRNA.multivariate_cox_result.txt",multi_res,quote=F,sep="\t")

ggforest(res.cox,  #coxphÂæóÂà∞ÁöÑCoxÂõûÂΩíÁªìÊûú
         data = Exp.data.Surv,  #Êï∞ÊçÆ???
         main = 'Hazard ratio of lncRNA',  #Ê†áÈ¢ò
         #cpositions = c(0.05, 0.35, 0.35),  #Ââç‰∏âÂàóË∑ù???
         fontsize = 1, #Â≠ó‰ΩìÂ§ßÂ∞è
         refLabel = 'reference', #Áõ∏ÂØπÂèòÈáèÁöÑÊï∞ÂÄºÊ†áÁ≠æÔºå‰πüÂèØÊîπ‰∏∫1
         noDigits = 3 #‰øùÁïôHRÂÄº‰ª•???95%CIÁöÑÂ∞èÊï∞‰Ωç???
)

######################### TCGA - GSE140845 (PRJNA591153) miRNA  ############
### miRNA Expression data -> RPM
setwd("K:/TCGA/Anlysis/LIHC")
LIHC.miRNA.exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)

### GSE140845 => log2FoldChange 1 => padj 0.05
GSE140845.miRNA <- read.csv("../../../Metastasis/TissueMicrobiome/PRJNA591152-Bact+Virus/GSE140845.miRNA.DEG.csv",row.names = 1)
GSE140845.miRNA.Up <- na.omit(GSE140845.miRNA) %>% filter(Sig == "Up")
GSE140845.miRNA.Down <- na.omit(GSE140845.miRNA) %>% filter(Sig == "Down")

### LIHC => log2FoldChange 2 => padj 0.01
LIHC.miRNA <- read.csv("../Pancancer/TCGA-LIHC_DE_Mature.miRNA.csv")
LIHC.miRNA <- LIHC.miRNA %>% separate(X, c("Precursor", "Mature.miRNA"),sep = "\\.")
LIHC.miRNA.Up <- LIHC.miRNA %>% filter(Threshold == "Up")
LIHC.miRNA.Down <- LIHC.miRNA %>% filter(Threshold == "Down")


load(file = "LIHC.GSE140845.miRNA.Rdata")
library(ggvenn)
library("RColorBrewer")
X <- list(GSE140845.Up = rownames(GSE140845.miRNA.Up),
          LIHC.Up = LIHC.miRNA.Up$Mature.miRNA,
          LIHC.Down= LIHC.miRNA.Down$Mature.miRNA,
          GSE140845.Down = rownames(GSE140845.miRNA.Down))

p3 <- ggvenn(X,set_name_size = 3,text_size =2,fill_color = brewer.pal(8,"Set1"))
ggsave(p3,filename = "LIHC.GSE140845.miRNA.Venn.pdf",height = 3,width = 5)

Inter.miRNA.Up <- intersect(rownames(GSE140845.miRNA.Up),LIHC.miRNA.Up$Mature.miRNA) #6
Inter.miRNA.Down <- intersect(unique(LIHC.miRNA.Down$Mature.miRNA),rownames(GSE140845.miRNA.Down)) #2

Inter.miRNA <- c(Inter.miRNA.Up,Inter.miRNA.Down)
save(Inter.miRNA.Up,Inter.miRNA.Down,file = "LIHC.GSE140845.miRNA.Intersect.Rdata")
save(GSE140845.miRNA.Up,GSE140845.miRNA.Down,LIHC.miRNA.Up,LIHC.miRNA.Down,Inter.miRNA.Up,Inter.miRNA.Down,Inter.miRNA, file = "LIHC.GSE140845.miRNA.Rdata")

######### P:Overall survival -> KM - miRNA ##########
SurvivalData <- read.table("survival_LIHC_survival.txt",header = T,row.names = 1,sep = '\t') 
# confirm sampletype to filter primary tumor
PhenoData <- read.csv("../../miRNA.Mature.Phenotype.csv") %>% dplyr::filter(Project.ID == "TCGA-LIHC") %>%
  dplyr::select(Sample.ID,Sample.Type) %>% 
  filter(Sample.Type == "Primary Tumor")

#substr(PhenoData$Sample.ID, 0, 15) %in% rownames(SurvivalData)
mid.Survival <- SurvivalData[substr(PhenoData$Sample.ID, 0, 15),]

RPMdata <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
RPMdata <- RPMdata %>% rownames_to_column() %>% separate(rowname, c("Precursor", "Mature.miRNA"),sep = "\\.") %>%
  dplyr::filter(Mature.miRNA %in% Inter.miRNA) %>% 
  dplyr::select(-Precursor) %>% remove_rownames() %>%
  column_to_rownames("Mature.miRNA") %>%
  dplyr::select(PhenoData$Sample.ID) %>% t() %>%
  as.data.frame(.,make.names=F)

mid.Survival <- mid.Survival %>% t() %>% as.data.frame(.,make.names=F) %>% dplyr::select(substr(PhenoData$Sample.ID, 0, 15)) %>%
  t() %>% as.data.frame(.,make.names=F)
mid.Survival$OS <- as.numeric(mid.Survival$OS)
mid.Survival$OS.time <- as.numeric(mid.Survival$OS.time)

save(RPMdata, mid.Survival, file = "LIHC.GSE140845.miRNA.KM.Data.Rdata")

library(survival)
library(survminer)

log_rank_p <- apply(RPMdata , 2 , function(gene){
  # gene=exprSet[1,]
  mid.Survival$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(Surv(OS.time, OS)~group,data=mid.Survival)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
# no significant p value for miRNA KM analysis

library(export)
library(ggthemes)

meta <- mid.Survival

dir.create("miRNA.OS.KM.Cox")
RPMdata2 <- RPMdata %>% data.frame()
for (i in 1:dim(RPMdata)[2]) {
  meta[[colnames(RPMdata2)[i]]] = ifelse(as.numeric(as.character(RPMdata2[,i]))> median(as.numeric(as.character(RPMdata2[,i]))),'high','low')
  sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","",colnames(RPMdata2)[i],"",sep = '')), data=meta)
  #pdf(paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 6,width = 6)
  #ggsave(filename = paste("lncRNA.Up.OS.KM.Cox/lncRNA.OS.KM.",lncRNA.Up.mid.gene.type$Gene.name[i],".pdf",sep = ''),height = 4,width = 4)
  p<-ggsurvplot(sfit1,pval =TRUE, data = meta, #risk.table = TRUE,
                surv.median.line = "hv", #fun = "event"ÂàôÂèØÁªòÂà∂cumulative eventsÂõæÔºàfunÈô§‰∫ÜevenÂ§ñÔºåËøòÊúâlogÔºölog transformationÂíåcumhazÔºöcumulative hazard???
                legend.title = colnames(RPMdata)[i],
                legend.labs = c("High", "Low"),
                conf.int.style = "step",
                xlab = "Time in days",
                risk.table = "abs_pct",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE,
                conf.int = TRUE,
                palette = "Set1",
                ggtheme = theme_few())
  print(p)
  #dev.off()
  graph2pdf(file=paste("miRNA.OS.KM.Cox/miRNA.OS.KM.",colnames(RPMdata)[i],".pdf",sep = ''),height = 6,width = 5)
  
  p<-ggsurvplot(sfit1,pval =TRUE, data = meta, #risk.table = TRUE,
                #surv.median.line = "hv", #fun = "event"ÂàôÂèØÁªòÂà∂cumulative eventsÂõæÔºàfunÈô§‰∫ÜevenÂ§ñÔºåËøòÊúâlogÔºölog transformationÂíåcumhazÔºöcumulative hazard???
                legend.title = colnames(RPMdata)[i],
                legend.labs = c("High", "Low"),
                conf.int.style = "step",
                xlab = "Time in days",
                risk.table = "abs_pct",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE,
                conf.int = TRUE,
                palette = "Set1",
                fun = "cumhaz",
                ggtheme = theme_few())
  print(p)
  #dev.off()
  graph2pdf(file=paste("miRNA.OS.KM.Cox/miRNA.OS.KM.cumhaz.",colnames(RPMdata)[i],".pdf",sep = ''),height = 6,width = 5)
  
}

######### P:Overall survival -> Cox - miRNA ##########
meta <- mid.Survival
miRNA.cox_results <-apply(RPMdata , 2 , function(gene){
  #gene= lncRNA.Up.Data[,1]
  meta$gene = gene
  #ÂèØÁõ¥Êé•‰ΩøÁî®ËøûÁª≠ÂûãÂèòÈáè
  m = coxph(Surv(OS.time, OS) ~ gene, data =  meta)
  #‰πüÂèØ‰ΩøÁî®‰∫åÂàÜÁ±ªÂèò???
  #meta$group=ifelse(gene>median(gene),'high','low') 
  #m=coxph(Surv(time, event) ~ group, data =  meta)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, 
                     se = se, z = beta/se, 
                     p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, 
                     HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)), 6)
  
  return(tmp['gene',]) 
  #return(tmp['grouplow',])#‰∫åÂàÜÁ±ªÂèò???
})

save(miRNA.cox_results,file = "miRNA.cox_results.Rdata")
saveRDS(miRNA.cox_results,"miRNA.cox_results.rds")

######################### TCGA - GSE140845 (PRJNA591153) mRNA  ############
#log2FC 2 -> padj 0.01
LIHC.lncRNA.mRNA <- readRDS("LIHC.RNA.rds")

GSE140845.mRNA <- read.csv("../../../Metastasis/TissueMicrobiome/PRJNA591152-Bact+Virus/mRNA.lncRNA.csv")

GSE140845.mRNA.Up <- na.omit(GSE140845.mRNA) %>% filter(Sig == "Up")
GSE140845.mRNA.Down <- na.omit(GSE140845.mRNA) %>% filter(Sig == "Down")

LIHC.mRNA.PCG.Up <- LIHC.lncRNA.mRNA$LIHC.mRNA.PCG.Up
LIHC.mRNA.PCG.Down <- LIHC.lncRNA.mRNA$LIHC.mRNA.PCG.Down

load("LIHC.GSE140845.mRNA.Rdata")
library(ggvenn)
library("RColorBrewer")
X <- list(GSE140845.Up = GSE140845.mRNA.Up$X,
          LIHC.Up = rownames(LIHC.mRNA.PCG.Up),
          LIHC.Down= rownames(LIHC.mRNA.PCG.Down),
          GSE140845.Down = GSE140845.mRNA.Down$X)

p3 <- ggvenn(X,set_name_size = 3,text_size =2,fill_color = brewer.pal(8,"Set1"))
ggsave(p3,filename = "LIHC.GSE140845.mRNA.Venn.pdf",height = 3,width = 5)

Inter.mRNA.Up <- intersect(rownames(LIHC.mRNA.PCG.Up),GSE140845.mRNA.Up$X) #150
Inter.mRNA.Down <- intersect(rownames(LIHC.mRNA.PCG.Down),GSE140845.mRNA.Down$X) #130

Inter.mRNA <- c(Inter.mRNA.Up,Inter.mRNA.Down)

save(GSE140845.mRNA.Up,GSE140845.mRNA.Down,LIHC.mRNA.PCG.Up,LIHC.mRNA.PCG.Down,Inter.mRNA.Up, Inter.mRNA.Down,Inter.mRNA,file = "LIHC.GSE140845.mRNA.Rdata")


######### KM logRank + UNICOX => DElncRNA DEmiRNA DEmRNA ###########
#### Keep consistent with LASSO model
## Datasets => Primary tumor
load(file = "LIHC.GSE140845.miRNA.Intersect.Rdata")
load(file = "LIHC.GSE140845.lncRNA.Rdata")
load(file = "LIHC.GSE140845.mRNA.Rdata")

library(export)
library(ggthemes)
library(survminer)
library(survival)

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
#mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

miRNA.Exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
#miRNA.Exp.log2.TPM <- log2(miRNA.Exp+1)

mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC") %>%
  filter(Sample.ID %in% colnames(mRNA.Exp)) %>% filter(Sample.Type == "Primary Tumor") #371

miRNA.Phenodata <- read.csv("../../miRNA.Mature.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC") %>%
  filter(Sample.ID %in% colnames(miRNA.Exp)) %>% filter(Sample.Type == "Primary Tumor") #371

#miRNA.Phenodata$Sample.ID %in% mRNA.Phenodata$Sample.ID
Intet.Sample.mRNA.miRNA <- intersect(miRNA.Phenodata$Sample.ID,mRNA.Phenodata$Sample.ID) #367

SurvivalData <- read.table("survival_LIHC_survival.txt",header = T,row.names = 1,sep = '\t') 
mid.surv <- SurvivalData[substr(Intet.Sample.mRNA.miRNA,1,15),] %>% dplyr::select(OS,OS.time) %>% na.omit() %>%
  filter(OS.time != 0)#361
Inter.Sample.361 <- Intet.Sample.mRNA.miRNA[substr(Intet.Sample.mRNA.miRNA,1,15) %in% rownames(mid.surv)]

save(Inter.Sample.361,mid.surv,file = "Inter.miRNA.mRNA.Samples.361.Rdata")

####  DElncRNA  ####
load("Inter.miRNA.mRNA.Samples.361.Rdata")
dir.create("361.lncRNA.361.KM.Cox")
lncRNA.361 <- mRNA.Exp[c(Inter.lncRNA.Down,Inter.lncRNA.Up),Inter.Sample.361]

rt <- cbind.data.frame(mid.surv,t.data.frame(lncRNA.361))

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #coxÂàÜÊûê
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  #KMÂàÜÊûê
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  if(length(table(group))==1) return(NULL) #ÂéªÊéâË°®ËææÈáèÊó†ÂèòÂåñÁöÑÂü∫???
  diff=survdiff(Surv(OS.time, OS) ~group,data = rt)
  
  rt2 <- rt %>% dplyr::select(OS.time, OS) %>% mutate(group=group)
  
  sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","group",sep = '')), data=rt2)
  
  pValue=1-pchisq(diff$chisq,df=1)
  
  p<-ggsurvplot(sfit1,pval =TRUE, data = rt2, 
                surv.median.line = "hv",
                legend.title = i,
                conf.int.style = "step",
                xlab = "Time in days",
                #break.time.by = 500,
                risk.table = "abs_pct",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE,
                legend.labs = c("High", "Low"),
                #pval = TRUE,
                conf.int = TRUE,
                palette = "Set1",
                ggtheme = theme_bw())
  print(p)
  #dev.off()
  graph2pdf(file=paste("361.lncRNA.361.KM.Cox/lncRNA.OS.KM.",i,".pdf",sep = ''),height = 6,width = 5)
  
  outTab=rbind(outTab,
               cbind(id=i,
                     KM=pValue,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv")

outTab <- merge.data.frame(outTab,Inter.lncRNA.GeneCard,by.x = "id",by.y="lncRNA.ENSEMBLE.ID")

write.csv(outTab,file = "361.lncRNA.KM.COX.csv",row.names = F)

####  DEmRNA  ####
load("Inter.miRNA.mRNA.Samples.361.Rdata")
dir.create("361.mRNA.361.KM.Cox")
mRNA.361 <- mRNA.Exp[c(Inter.mRNA.Down,Inter.mRNA.Up),Inter.Sample.361]

rt <- cbind.data.frame(mid.surv,t.data.frame(mRNA.361))

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #coxÂàÜÊûê
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  #KMÂàÜÊûê
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  if(length(table(group))==1) return(NULL) #ÂéªÊéâË°®ËææÈáèÊó†ÂèòÂåñÁöÑÂü∫???
  diff=survdiff(Surv(OS.time, OS) ~group,data = rt)
  
  rt2 <- rt %>% dplyr::select(OS.time, OS) %>% mutate(group=group)
  
  sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","group",sep = '')), data=rt2)
  
  pValue=1-pchisq(diff$chisq,df=1)
  if (F) {
    p<-ggsurvplot(sfit1,pval =TRUE, data = rt2, 
                  surv.median.line = "hv",
                  legend.title = i,
                  conf.int.style = "step",
                  xlab = "Time in days",
                  #break.time.by = 500,
                  risk.table = "abs_pct",
                  risk.table.y.text.col = T,
                  risk.table.y.text = FALSE,
                  legend.labs = c("High", "Low"),
                  #pval = TRUE,
                  conf.int = TRUE,
                  palette = "Set1",
                  ggtheme = theme_bw())
    print(p)
    #dev.off()
    #graph2pdf(file=paste("361.mRNA.361.KM.Cox/mRNA.OS.KM.",i,".pdf",sep = ''),height = 6,width = 5)
    
  }

  outTab=rbind(outTab,
               cbind(id=i,
                     KM=pValue,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

mid.Gene.Type <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol) %>% unique()
outTab2 <- merge.data.frame(outTab,mid.Gene.Type,by.x = "id",by.y="Gene.stable.ID") %>% unique()

write.csv(outTab2,file = "361.mRNA.KM.COX.csv",row.names = F)


####  DEmiRNA  ####
load("Inter.miRNA.mRNA.Samples.361.Rdata")
dir.create("361.miRNA.361.KM.Cox")
miRNA.Index <- str_remove_all(rownames(miRNA.Exp),".*\\.") %in% c(Inter.miRNA.Up,Inter.miRNA.Down)
miRNA.361 <- miRNA.Exp[miRNA.Index,Inter.Sample.361]

rt <- cbind.data.frame(mid.surv,t.data.frame(miRNA.361))

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #coxÂàÜÊûê
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  #KMÂàÜÊûê
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  if(length(table(group))==1) return(NULL) #ÂéªÊéâË°®ËææÈáèÊó†ÂèòÂåñÁöÑÂü∫???
  diff=survdiff(Surv(OS.time, OS) ~group,data = rt)
  
  rt2 <- rt %>% dplyr::select(OS.time, OS) %>% mutate(group=group)
  
  sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","group",sep = '')), data=rt2)
  
  pValue=1-pchisq(diff$chisq,df=1)
  name = str_remove_all(i,".*\\.")
  p<-ggsurvplot(sfit1,pval =TRUE, data = rt2, 
                surv.median.line = "hv",
                legend.title = name,
                conf.int.style = "step",
                xlab = "Time in days",
                #break.time.by = 500,
                risk.table = "abs_pct",
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE,
                legend.labs = c("High", "Low"),
                #pval = TRUE,
                conf.int = TRUE,
                palette = "Set1",
                ggtheme = theme_bw())
  print(p)
  #dev.off()
  graph2pdf(file=paste("361.miRNA.361.KM.Cox/miRNA.OS.KM.",name,".pdf",sep = ''),height = 6,width = 5)
  
  outTab=rbind(outTab,
               cbind(id=name,
                     KM=pValue,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

write.csv(outTab,file = "361.miRNA.KM.COX.csv",row.names = F)




######################### ceRNA - Database #######################
Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv") 
load("LIHC.GSE140845.lncRNA.Rdata")
load("LIHC.GSE140845.mRNA.Rdata")
load("LIHC.GSE140845.miRNA.Rdata")
mRNA.ENSEMBLE.NAME <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name) %>% filter(Gene.stable.ID %in% Inter.mRNA) %>% unique()

#### lncRNA - miRNA -> lnc2base database ####
lnc2miRNA <- read.csv("../../../Metastasis/LncBase2.MicroRNA.txt",header = T,sep = "\t")

lnc.miRNA <- data.frame()
for (i in 1:dim(lnc2miRNA)[1]) {
  miRNA <- lnc2miRNA$Set[i]
  lncRNA <- (lnc2miRNA$LncRNA[i] %>% str_split(";"))[[1]]
  for (j in 1:length(lncRNA)) {
    lnc.miRNA <- rbind.data.frame(lnc.miRNA,c(miRNA,lncRNA[j]))
  }
}

lnc.miRNA2 <- lnc.miRNA %>% rename(miRNA=1,lncRNA=2) %>% arrange(miRNA,lncRNA) %>% unique()
saveRDS(lnc.miRNA2,"ceRNA.lnc2Base.lncRNA_miRNA.rds")

Inter.lncRNA.GeneCard$lncRNA.ID[Inter.lncRNA.GeneCard$lncRNA.ID %in% lnc.miRNA2$lncRNA]
#"MCM3AP-AS1"  "SNHG20"      "TMEM147-AS1" "CDKN2B-AS1"  "B4GALT1-AS1"

lnc2Base.DElncRNA.DEmiRNA <- lnc.miRNA2 %>% filter(lncRNA %in% Inter.lncRNA.GeneCard$lncRNA.ID) %>% 
  filter(miRNA %in% Inter.miRNA) %>% merge.data.frame(.,Inter.lncRNA.GeneCard,by.x= "lncRNA",by.y="lncRNA.ID") %>%
  dplyr::select(miRNA,lncRNA.ENSEMBLE.ID,lncRNA,Sig) %>%
  rename(miRNA.ID=1,lncRNA.ENSEMBLE.ID=2,lncRNA.NAME=3,lncRNA.Sig=4) %>%
  mutate(miRNA.Sig = if_else(miRNA.ID %in% Inter.miRNA.Up,"Up","Down")) %>%
  dplyr::select(miRNA.ID,lncRNA.ENSEMBLE.ID,lncRNA.NAME,miRNA.Sig,lncRNA.Sig)

save(lnc2Base.DElncRNA.DEmiRNA,file = "lnc2Base.DElncRNA.DEmiRNA.rds")

#### lncRNA - miRNA + mRNA - miRNA -> ENCORI database ####
setwd("K:/TCGA/Anlysis/LIHC/ceRNA-ENCORI-Database")
### download lncRNA_miRNA_interaction 
family=read.table("refData/hg19_all_fimaly.txt",sep="\t")
miRNA=unlist(strsplit(as.character(family$V4),",")) %>% sort() %>% unique()

dir.create("lncRNA-miRNA")
for(mir in miRNA){
  file=paste("lncRNA-miRNA/",mir,".txt",sep="")
  link=paste("http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=lncRNA&miRNA=",mir,"&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all",sep="")
  download.file(link,file)
  Sys.sleep(2)
}

lncRNA_files <- list.files(path="lncRNA-miRNA", full.names=TRUE)

library(plyr)
lncRNA.list <- llply(lncRNA_files, function(x)read.table(x,header=T,sep="\t",comment.char ="#",stringsAsFactors=F))
combind_lncRNA=do.call(rbind,lncRNA.list)
write.table(file="ENCORI_lncRNA_miRNA_interaction.txt",combind_lncRNA,quote=F,sep="\t",row.names = F)


### download mRNA_miRNA_interaction
dir.create("mRNA-miRNA")
for(mir in miRNA){
  file=paste("mRNA-miRNA/",mir,".txt",sep="")
  link=paste("http://starbase.sysu.edu.cn/api/miRNATarget/?assembly=hg19&geneType=mRNA&miRNA=",mir,"&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=1&program=None&target=all",sep="")
  download.file(link,file)
  Sys.sleep(2)
}

lncRNA_files <- list.files(path="mRNA-miRNA", full.names=TRUE)

library(plyr)
lncRNA.list <- llply(lncRNA_files, function(x)read.table(x,header=T,sep="\t",stringsAsFactors=F,comment.char = "#",row.names =NULL))
combind_mRNA=do.call(rbind,lncRNA.list)
combind_mRNA2 <- combind_mRNA[,-5]
colnames(combind_mRNA2)[1:4] <- colnames(combind_mRNA)[2:5]
write.table(file="ENCORI_mRNA_miRNA_interaction.txt",combind_mRNA2,quote=F,sep="\t",row.names = F)

ENCORI.lncRNA.miRNA <- read.table("ceRNA-ENCORI-Database/ENCORI_lncRNA_miRNA_interaction.txt",sep = "\t",header = T)
ENCORI.mRNA.miRNA <- read.table("ceRNA-ENCORI-Database/ENCORI_mRNA_miRNA_interaction.txt",sep = "\t",header = T)

Inter.lncRNA.GeneCard$lncRNA.ENSEMBLE.ID[Inter.lncRNA.GeneCard$lncRNA.ENSEMBLE.ID %in% ENCORI.lncRNA.miRNA$geneID]
Inter.lncRNA.GeneCard$lncRNA.ID[Inter.lncRNA.GeneCard$lncRNA.ID %in% ENCORI.lncRNA.miRNA$geneName]
#"PKP4-AS1"    "MCM3AP-AS1"  "SNHG20"      "TMEM147-AS1" "CDKN2B-AS1"  "LINC01831"   "LINC02027"

ENCORI.DElncRNA.DEmiRNA <- ENCORI.lncRNA.miRNA %>% dplyr::select(miRNAname,geneID,geneName,geneType) %>%
  filter(geneID %in% (Inter.lncRNA.GeneCard$lncRNA.ENSEMBLE.ID[Inter.lncRNA.GeneCard$lncRNA.ENSEMBLE.ID %in% ENCORI.lncRNA.miRNA$geneID])) %>%
  filter(miRNAname %in% Inter.miRNA) %>% dplyr::select(-geneType) %>% 
  mutate(miRNA.Sig = if_else(miRNAname %in% Inter.miRNA.Up,"Up","Down")) %>%
  mutate(lncRNA.Sig = if_else(geneID %in% Inter.lncRNA.Up,"Up","Down")) %>% unique() %>%
  rename(miRNA.ID=1,lncRNA.ENSEMBLE.ID=2,lncRNA.NAME=3,miRNA.Sig=4,lncRNA.Sig=5) %>%
  dplyr::select(miRNA.ID,lncRNA.ENSEMBLE.ID,lncRNA.NAME,miRNA.Sig,lncRNA.Sig)

saveRDS(ENCORI.DElncRNA.DEmiRNA,file = "ENCORI.DElncRNA.DEmiRNA.rds")

# "hsa-miR-10b-5p"  "hsa-miR-183-5p"  "hsa-miR-216b-5p" "hsa-miR-424-5p"  "hsa-miR-552-3p" 
# "ENSG00000236144" "ENSG00000215068" "ENSG00000273329" "ENSG00000240498" "ENSG00000242553" "ENSG00000215424" "ENSG00000272068"

#### lncRNA - miRNA + mRNA - miRNA -> miRWalk2.0 database ####
load("ceRNA-miRWalk-Database/hsaLncRNA-EnsemblGID.rdata") #hsa

miRWalk2.lncRNA.miRNA <- data.frame()
for (i in Inter.miRNA) {
  for (j in 1:length(hsa[[i]])) {
    miRWalk2.lncRNA.miRNA <- rbind.data.frame(miRWalk2.lncRNA.miRNA,c(i,(hsa[[i]])[j]))
  }
}

miRWalk2.DElncRNA.DEmiRNA <- miRWalk2.lncRNA.miRNA %>% rename(miRNA.ID=1, lncRNA.ENSEMBLE.ID=2) %>%
  filter(lncRNA.ENSEMBLE.ID %in% Inter.lncRNA.GeneCard$lncRNA.ENSEMBLE.ID) %>% 
  merge(.,Inter.lncRNA.GeneCard,by = "lncRNA.ENSEMBLE.ID") %>% 
  mutate(miRNA.Sig = if_else(miRNA.ID %in% Inter.miRNA.Up,"Up","Down")) %>%
  rename(lncRNA.ENSEMBLE.ID=1,miRNA.ID=2,lncRNA.NAME=3,lncRNA.Sig=4,miRNA.Sig=5) %>%
  dplyr::select(miRNA.ID,lncRNA.ENSEMBLE.ID,lncRNA.NAME,miRNA.Sig,lncRNA.Sig)

saveRDS(miRWalk2.DElncRNA.DEmiRNA,file = "miRWalk2.DElncRNA.DEmiRNA.rds")

###### lncRNA - miRNA => lnc2base v3.0 + ENCORI + miRWalk2.0 ####
Integrate.DElncRNA.DEmiRNA <- rbind(rbind(lnc2Base.DElncRNA.DEmiRNA %>% mutate(Database = "lnc2Base"),
                                          ENCORI.DElncRNA.DEmiRNA %>% mutate(Database = "ENCORI")),
                                    miRWalk2.DElncRNA.DEmiRNA %>% mutate(Database = "miRWalk2")) %>% data.frame()
Integrate.UpDElncRNA.DownDEmiRNA <- Integrate.DElncRNA.DEmiRNA %>% filter(lncRNA.Sig == "Up" & miRNA.Sig == "Down") %>%
  arrange(miRNA.ID,lncRNA.ENSEMBLE.ID) %>% unique()

write.csv(Integrate.UpDElncRNA.DownDEmiRNA,file = "Integrate.UpDElncRNA.DownDEmiRNA.csv",row.names = F)

Integrate.UpDElncRNA.DownDEmiRNA <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.csv")
Integrate.UpDElncRNA.DownDEmiRNA.Unique <- Integrate.UpDElncRNA.DownDEmiRNA %>% dplyr::select(-Database) %>% arrange(miRNA.ID,lncRNA.ENSEMBLE.ID) %>% unique()

## lncRNA down -> miRNA Up => mRNA up
Integrate.DownDElncRNA.UpDEmiRNA <- Integrate.DElncRNA.DEmiRNA %>% filter(lncRNA.Sig == "Down" & miRNA.Sig == "Up") %>%
  arrange(miRNA.ID,lncRNA.ENSEMBLE.ID) %>%
  unique()

save(Integrate.DElncRNA.DEmiRNA,
     Integrate.UpDElncRNA.DownDEmiRNA,
     Integrate.UpDElncRNA.DownDEmiRNA.Unique,
     Integrate.DownDElncRNA.UpDEmiRNA,
     file = "Integrate.DElncRNA.DEmiRNA.Rdata")

load("Integrate.DElncRNA.DEmiRNA.Rdata")

Integrate.UpDElncRNA.DownDEmiRNA.Count <- Integrate.UpDElncRNA.DownDEmiRNA %>% 
  dplyr::group_by(miRNA.ID,lncRNA.ENSEMBLE.ID,lncRNA.NAME,miRNA.Sig,lncRNA.Sig) %>%
  summarise(lncRNA.miRNA.Count = n())

Integrate.DownDElncRNA.UpDEmiRNA.Count <- Integrate.DownDElncRNA.UpDEmiRNA %>% 
  dplyr::group_by(miRNA.ID,lncRNA.ENSEMBLE.ID,lncRNA.NAME,miRNA.Sig,lncRNA.Sig) %>%
  summarise(lncRNA.miRNA.Count = n())

save(Integrate.UpDElncRNA.DownDEmiRNA.Count,Integrate.DownDElncRNA.UpDEmiRNA.Count,file = "Integrate.DElncRNA.DEmiRNA.Count.Rdata")

## lncRNA Up + miRNA Down + mRNA Up <====> lncRNA Down + miRNA Up + mRNA Down


#### mRNA - miRNA -> ENCORI database ####
ENCORI.mRNA.miRNA <- read.csv("ceRNA-ENCORI-Database/ENCORI_mRNA_miRNA_interaction.txt",header = T,sep = "\t")
ENCORI.DEmRNA.DEmiRNA <- ENCORI.mRNA.miRNA %>% dplyr::select(miRNAname,geneID,geneName) %>%
  rename(miRNA.ID=1,mRNA.ENSEMBLE.ID=2,mRNA.NAME=3) %>% filter(miRNA.ID %in% Inter.miRNA) %>%
  filter(mRNA.ENSEMBLE.ID %in% Inter.mRNA) %>% arrange(miRNA.ID,mRNA.ENSEMBLE.ID) %>% unique() %>%
  mutate(miRNA.Sig = if_else(miRNA.ID %in% Inter.miRNA.Up,"Up","Down")) %>%
  mutate(mRNA.Sig = if_else(mRNA.ENSEMBLE.ID %in% Inter.mRNA.Up,"Up","Down"))

saveRDS(ENCORI.DEmRNA.DEmiRNA,file = "ENCORI.DEmRNA.DEmiRNA.rds")

ENCORI.DEmRNA.DEmiRNA <- readRDS("ENCORI.DEmRNA.DEmiRNA.rds")

#### mRNA - miRNA -> miRWalk2.0 Validated database ####
load("ceRNA-miRWalk-Database/hsa-vtm-gene.rdata.Rdata")#id
mRNA.ENSEMBLE.NAME <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name) %>% filter(Gene.stable.ID %in% Inter.mRNA) %>% unique()

miRWalk2.mRNA.miRNA <- data.frame()
for (i in Inter.miRNA) {
  for (j in 1:length(id[[i]])) {
    miRWalk2.mRNA.miRNA <- rbind.data.frame(miRWalk2.mRNA.miRNA,c(i,(id[[i]])[j]))
  }
}

miRWalk2.DEmRNA.DEmiRNA <- miRWalk2.mRNA.miRNA %>% rename(miRNA.ID=1,mRNA.NAME=2) %>% filter(miRNA.ID!=mRNA.NAME) %>%
  filter(miRNA.ID %in% Inter.miRNA) %>% filter(mRNA.NAME %in% mRNA.ENSEMBLE.NAME$Gene.name) %>%
  merge.data.frame(.,mRNA.ENSEMBLE.NAME,by.x = "mRNA.NAME",by.y = "Gene.name") %>% 
  dplyr::select(miRNA.ID,Gene.stable.ID,mRNA.NAME) %>%
  rename(miRNA.ID=1,mRNA.ENSEMBLE.ID=2,mRNA.NAME=3) %>% 
  mutate(miRNA.Sig = if_else(miRNA.ID %in% Inter.miRNA.Up,"Up","Down")) %>%
  mutate(mRNA.Sig = if_else(mRNA.ENSEMBLE.ID %in% Inter.mRNA.Up,"Up","Down"))
  
saveRDS(miRWalk2.DEmRNA.DEmiRNA,file = "miRWalk2.DEmRNA.DEmiRNA.rds")

miRWalk2.DEmRNA.DEmiRNA <- readRDS("miRWalk2.DEmRNA.DEmiRNA.rds")

if (F) {
  files <- list.files(path="ceRNA-miRWalk-Database/miRNA.mRNA/", pattern = "[txt]",full.names = T)
  miRWalk2.mRNA.miRNA.Predict <- data.frame()
  library(tidyverse)
  for (x in files) {
    testdata <- read.table(x,sep = "\t",header = T,row.names = 1,check.names = F)
    testdata.mid <- testdata %>% dplyr::select(-EntrezID,-NmiRs) %>% rownames_to_column() %>%
      gather(miRNA.ID,Interaction,-rowname) %>% mutate(miRNA.ID = str_remove_all(miRNA.ID,"\\_.*")) %>% 
      filter(Interaction == 1) %>% rename(mRNA.NAME=1,miRNA.ID=2) %>% arrange(miRNA.ID,mRNA.NAME)  %>% 
      group_by(miRNA.ID,mRNA.NAME) %>% summarise(Count = n()) %>%
      filter(miRNA.ID %in% Inter.miRNA) %>% filter(mRNA.NAME %in% mRNA.ENSEMBLE.NAME$Gene.name) %>% unique()
    
    miRWalk2.mRNA.miRNA.Predict <- rbind.data.frame(miRWalk2.mRNA.miRNA.Predict,testdata.mid)
  }
  
  miRWalk2.mRNA.miRNA.Predict <- read.csv("ceRNA-miRWalk-Database/miRWalk2.mRNA.miRNA.Predict.csv") %>% dplyr::select(-Interaction) %>%
    dplyr::select(miRNA.ID,mRNA.NAME)
  
  miRWalk2.mRNA.miRNA <- rbind.data.frame(miRWalk2.mRNA.miRNA,miRWalk2.mRNA.miRNA.Predict) %>%
    merge.data.frame(.,mRNA.ENSEMBLE.NAME,by.x = "mRNA.NAME",by.y = "Gene.name") %>% dplyr::select(miRNA.ID,Gene.stable.ID,mRNA.NAME)
}


#### mRNA - miRNA -> miRWalk3.0 database -> query DEmiRNA ####
files <- list.files(path="ceRNA-miRWalk-Database/miRWalk3/", pattern = "[csv]",full.names = T)

miRWalk3.mRNA.miRNA <- data.frame()
for (x in files) {
  testdata <- read.csv(x) %>% rename(miRNA.ID=1,mRNA.NAME=3) %>%
    dplyr::select(miRNA.ID,mRNA.NAME) %>% filter(miRNA.ID %in% Inter.miRNA) %>%
    filter(mRNA.NAME %in% mRNA.ENSEMBLE.NAME$Gene.name) %>% unique() %>%
    merge.data.frame(.,mRNA.ENSEMBLE.NAME,by.x = "mRNA.NAME",by.y = "Gene.name") %>% 
    dplyr::select(miRNA.ID,Gene.stable.ID,mRNA.NAME) %>%
    rename(miRNA.ID=1,mRNA.ENSEMBLE.ID=2,mRNA.NAME=3) %>% 
    mutate(miRNA.Sig = if_else(miRNA.ID %in% Inter.miRNA.Up,"Up","Down")) %>%
    mutate(mRNA.Sig = if_else(mRNA.ENSEMBLE.ID %in% Inter.mRNA.Up,"Up","Down"))
  
  miRWalk3.mRNA.miRNA <- rbind.data.frame(miRWalk3.mRNA.miRNA,testdata)
}
miRWalk3.DEmRNA.DEmiRNA <- miRWalk3.mRNA.miRNA

saveRDS(miRWalk3.DEmRNA.DEmiRNA,"miRWalk3.DEmRNA.DEmiRNA.rds")

miRWalk3.DEmRNA.DEmiRNA  <- readRDS("miRWalk3.DEmRNA.DEmiRNA.rds")

#### mRNA - miRNA -> multiMIR R package -> Cannor access ####
#3‰∏™Â∑≤È™åËØÅÁöÑmiRNA-Èù∂Âü∫Âõ†Áõ∏‰∫íÂÖ≥Á≥ªÁöÑÊï∞ÊçÆÂ∫ìÔºàmirecords", "mirtarbase", "tarbase"ÔºâÔºõ
#8‰∏™È¢ÑÊµãÊï∞ÊçÆÂ∫ì???"diana_microt", "elmmo", "microcosm", "miranda", "mirdb", "pictar", "pita", "targetscan"ÔºâÔºõ
#3‰∏™miRNA‰∏éÁñæÁóÖÂíåËçØÁâ©Áõ∏ÂÖ≥Êï∞ÊçÆÂ∫ìÔºà"mir2disease", "pharmaco_mir", "phenomir"???
library(multiMiR)
#BiocManager::install("multiMiR",ask = F,update = F)
#multimir_switchDBVersion(db_version = "2.0.0")
multimir_switchDBVersion(db_version = curr_vers)


#### mRNA - miRNA -> miRNAtap R package ####
#BiocManager::install("miRNAtap",ask = F,update = F)
library(miRNAtap)
library(org.Hs.eg.db)
library(clusterProfiler)

miRNAtap.DEmiRNA.mRNA <- data.frame()
for (miRNA in Inter.miRNA) {
  miRNAtap.DEmiRNA = getPredictedTargets(miRNA, species = 'hsa',method = 'geom', min_src = 2)
  map_dt <- bitr(rownames(miRNAtap.DEmiRNA), fromType = "ENTREZID",toType = c("ENSEMBL","SYMBOL"),OrgDb = org.Hs.eg.db)
  map_dt <- map_dt %>% data.frame() %>% mutate(miRNA.ID = miRNA)
  miRNAtap.DEmiRNA.mRNA <- rbind.data.frame(miRNAtap.DEmiRNA.mRNA,map_dt)
}

miRNAtap.DEmiRNA.DEmRNA <- miRNAtap.DEmiRNA.mRNA %>% dplyr::select(miRNA.ID,ENSEMBL,SYMBOL) %>% 
  rename(mRNA.ENSEMBLE.ID=2,mRNA.NAME=3) %>% mutate(miRNA.Sig = if_else(miRNA.ID %in% Inter.miRNA.Up,"Up","Down")) %>%
  mutate(mRNA.Sig = if_else(mRNA.ENSEMBLE.ID %in% Inter.mRNA.Up,"Up","Down"))

saveRDS(miRNAtap.DEmiRNA.DEmRNA,"miRNAtap.DEmiRNA.DEmRNA.rds")

miRNAtap.DEmiRNA.DEmRNA <- readRDS("miRNAtap.DEmiRNA.DEmRNA.rds")

###### mRNA - miRNA -> ENCORI + miRWalk2.0 + miRWalk3.0 + miRNAtap #######################
Integrate.DEmiRNA.DEmRNA <- rbind(rbind(rbind(ENCORI.DEmRNA.DEmiRNA %>% mutate(Database = "ENCORI"),
                                              miRWalk2.DEmRNA.DEmiRNA %>% mutate(Database = "miRWalk2")),
                                        miRWalk3.DEmRNA.DEmiRNA %>% mutate(Database = "miRWalk3")),
                                  miRNAtap.DEmiRNA.DEmRNA %>% mutate(Database = "miRNAtap"))

Integrate.UpDEmiRNA.DownDEmRNA <- Integrate.DEmiRNA.DEmRNA %>% data.frame() %>% filter(miRNA.Sig == "Up" & mRNA.Sig == "Down")
Integrate.DownDEmiRNA.UpDEmRNA <- Integrate.DEmiRNA.DEmRNA %>% data.frame() %>% filter(miRNA.Sig == "Down" & mRNA.Sig == "Up")  

save(Integrate.UpDEmiRNA.DownDEmRNA,Integrate.DownDEmiRNA.UpDEmRNA,file = "Integrate.DEmiRNA.DEmRNA.True.Rdata")

Integrate.DEmiRNA.DEmRNA.True <- rbind(Integrate.UpDEmiRNA.DownDEmRNA,Integrate.DownDEmiRNA.UpDEmRNA)
write.csv(Integrate.DEmiRNA.DEmRNA.True,"Integrate.DEmiRNA.DEmRNA.True.csv",row.names = F)


Integrate.UpDEmiRNA.DownDEmRNA.Count <- Integrate.UpDEmiRNA.DownDEmRNA %>% 
  group_by(miRNA.ID,mRNA.ENSEMBLE.ID,mRNA.NAME,miRNA.Sig,mRNA.Sig) %>%
  summarise(Count=n()) %>% filter(miRNA.mRNA.Count > 1)
  
Integrate.DownDEmiRNA.UpDEmRNA.Count <- Integrate.DownDEmiRNA.UpDEmRNA %>% 
  group_by(miRNA.ID,mRNA.ENSEMBLE.ID,mRNA.NAME,miRNA.Sig,mRNA.Sig) %>%
  summarise(Count=n()) %>% filter(miRNA.mRNA.Count > 1)
  
save(Integrate.UpDEmiRNA.DownDEmRNA.Count,Integrate.DownDEmiRNA.UpDEmRNA.Count,file = "Integrate.DEmiRNA.DEmRNA.True.Count.Rdata")


###### CeRNA constructing -1 => Database Integrate => lncRNA:miRNA -> 1 Database + miRNA:mRNA -> 2 Database ######
load("Integrate.DElncRNA.DEmiRNA.Count.Rdata")
load("Integrate.DEmiRNA.DEmRNA.True.Count.Rdata")

Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA <- merge.data.frame(Integrate.DownDElncRNA.UpDEmiRNA.Count,
                                                                Integrate.UpDEmiRNA.DownDEmRNA.Count,by="miRNA.ID")

write.csv(Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA,"ceRNA.Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.csv",row.names = F)


Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA <- merge.data.frame(Integrate.UpDElncRNA.DownDEmiRNA.Count,
                                                                Integrate.DownDEmiRNA.UpDEmRNA.Count,by="miRNA.ID")

write.csv(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA,"ceRNA.Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.csv",row.names = F)

save(Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA,Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA,file = "ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata")

###### Cytoscape => 16UplncRNA 28DEmRNA 2DownDEmiRNA + 21DownDElncRNA 33DownDEmRNA 6 UpDEmiRNA #####
load("ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata") #Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA,Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA
Data.set <- matrix(nrow = 0,ncol = 6)
colnames(Data.set) = c("miRNA","Gene","Biotype","miRNA.Sig","Gene.Sig","Interaction.Count")
mid1 <- data.frame(miRNA = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$miRNA.ID,
                   Gene = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$lncRNA.NAME,
                   Biotype = "lncRNA",
                   miRNA.Sig = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$miRNA.Sig.x,
                   Gene.Sig = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$lncRNA.Sig,
                   Interaction.Count = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$lncRNA.miRNA.Count)

mid2 <- data.frame(miRNA = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$miRNA.ID,
                   Gene = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$mRNA.NAME,
                   Biotype = "mRNA",
                   miRNA.Sig = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$miRNA.Sig.x,
                   Gene.Sig = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$mRNA.Sig,
                   Interaction.Count = Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$Count)

mid3 <- data.frame(miRNA = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$miRNA.ID,
                   Gene = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$lncRNA.NAME,
                   Biotype = "lncRNA",
                   miRNA.Sig = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$miRNA.Sig.x,
                   Gene.Sig = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$lncRNA.Sig,
                   Interaction.Count = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$lncRNA.miRNA.Count)

mid4 <- data.frame(miRNA = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$miRNA.ID,
                   Gene = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$mRNA.NAME,
                   Biotype = "mRNA",
                   miRNA.Sig = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$miRNA.Sig.x,
                   Gene.Sig = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$mRNA.Sig,
                   Interaction.Count = Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$Count)


Data.set <- rbind(rbind(rbind(rbind(Data.set,mid1),mid2),mid3),mid4) %>% data.frame() %>% arrange(miRNA,Gene) %>% unique()
write.csv(Data.set,"ceRNA.DElncRNA.DEmiRNA.DEmRNA.Cytoscape.csv",row.names = F)


######################## P:DElncRNA DEmiRNA DEmRNA => Pearson Correlation -> Including Normal Samples ################
setwd("K:/TCGA/Anlysis/LIHC")
load("ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata") #Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA; Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA

library(ggcor)
library(dplyr)
library(ggplot2)
library(psych)
library(RColorBrewer)
library(export)

Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
Exp.log2.TPM <- log2(Exp+1)
Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC")
mid.Phenodata <- Phenodata %>% filter(Sample.ID %in% colnames(Exp)) %>% 
  filter(Sample.Type == "Solid Tissue Normal" | Sample.Type == "Primary Tumor") #421
mid.Exp.log2.TPM <- Exp.log2.TPM %>% dplyr::select(mid.Phenodata$Sample.ID)

for (miRNA in unique(Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$miRNA.ID)) {
  LINE <- Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA %>% filter(miRNA.ID == miRNA)
  mRNA <- LINE %>% dplyr::select(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% arrange(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% unique()
  lncRNA <- LINE %>% dplyr::select(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% arrange(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% unique()
  
  lncRNA.Exp <- mid.Exp.log2.TPM[lncRNA$lncRNA.ENSEMBLE.ID,] %>% t()
  mRNA.Exp <- mid.Exp.log2.TPM[mRNA$mRNA.ENSEMBLE.ID,] %>% t()
  
  lnc.mRNA.ccor <- corr.test(lncRNA.Exp,mRNA.Exp, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=TRUE)
  
  rho <- lnc.mRNA.ccor$r
  colnames(rho) = mRNA$mRNA.NAME
  rownames(rho) = lncRNA$lncRNA.NAME
  
  padj <- lnc.mRNA.ccor$p.adj
  colnames(padj) = mRNA$mRNA.NAME
  rownames(padj) = lncRNA$lncRNA.NAME
  
  rho2 <- rho %>% data.frame(check.names = F) %>% rownames_to_column() %>% gather(mRNA.Symbol,Pearson.Correlation,-rowname) %>%
    rename(lncRNA.Symbol=1)
  
  padj2 <- padj %>% data.frame(check.names = F) %>% rownames_to_column() %>% gather(mRNA.Symbol,Pearson.FDR,-rowname) %>%
    rename(lncRNA.Symbol=1)
  
  rho.padj <- cbind.data.frame(rho2,data.frame(Pearson.FDR=padj2$Pearson.FDR)) %>%
    mutate(r = cut(abs(Pearson.Correlation), breaks = c(-Inf, 0.3, 0.6, Inf),
                   labels = c("0~0.3", "0.3~0.6", ">=0.6"),
                   right = FALSE),
           p.adj = cut(Pearson.FDR, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
                         right = FALSE)) %>%
    mutate(r=factor(r,levels = c("0~0.3", "0.3~0.6", ">=0.6"))) %>%
    mutate(Cor.Direction = if_else(Pearson.Correlation >= 0,"Positive","Negative")) %>%
    mutate(Cor.Direction = factor(Cor.Direction,levels = c("Positive","Negative")))
  
  mRNA.ccor <- corr.test(mRNA.Exp, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=TRUE)
  
  Rho = mRNA.ccor$r
  colnames(Rho) = rownames(Rho) = mRNA$mRNA.NAME
  
  colnames(mRNA.Exp) = mRNA$mRNA.NAME
  
  p1<-quickcor(mRNA.Exp, type = "upper") + geom_square(inherit.aes = T) +
    anno_link(rho.padj, mapping = aes(colour = p.adj, size = r,linetype = Cor.Direction)) +
    scale_size_manual(values = c(0.1, 0.5, 1)) 
  p2<-p1+scale_fill_gradient2(midpoint = 0, low = brewer.pal(8,"Set1")[2], mid = "white",
                              high = brewer.pal(8,"Set1")[1], space = "Lab" )
  p3 <- p2+scale_color_manual(values=c(brewer.pal(8,"Dark2")[4:2],"grey60"))
  p4 <- p3 + scale_linetype_manual(values = c(1,2))
  p5 <- p4+guides(size=guide_legend(title="Pearson's CC",override.aes=list(colour="grey35"),order=1),
            colour=guide_legend(title="Pearson's FDR",override.aes = list(size=3),order=2),
            linetype = guide_legend(title="Correlativity",override.aes = list(size=0.6),order=3),
            fill=guide_colorbar(title="Pearson's CC",order=4))
  
  print(p5)
  graph2pdf(file=paste("lncRNA.mRNA.Correlation.DownDElncRNA.UpDEmiRNA.DownDEmRNA.",miRNA,".pdf",sep = ""),height = 6,width = 12)
  
  save(lnc.mRNA.ccor,mRNA.ccor,rho.padj,Rho,file = paste("lncRNA.mRNA.Correlation.DownDElncRNA.UpDEmiRNA.DownDEmRNA.",miRNA,".Rdata",sep = ""))
}

for (miRNA in unique(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$miRNA.ID)) {
  LINE <- Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA %>% filter(miRNA.ID == miRNA)
  mRNA <- LINE %>% dplyr::select(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% arrange(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% unique()
  lncRNA <- LINE %>% dplyr::select(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% arrange(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% unique()
  
  lncRNA.Exp <- mid.Exp.log2.TPM[lncRNA$lncRNA.ENSEMBLE.ID,] %>% t()
  mRNA.Exp <- mid.Exp.log2.TPM[mRNA$mRNA.ENSEMBLE.ID,] %>% t()
  
  lnc.mRNA.ccor <- corr.test(lncRNA.Exp,mRNA.Exp, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=TRUE)
  
  rho <- lnc.mRNA.ccor$r
  colnames(rho) = mRNA$mRNA.NAME
  rownames(rho) = lncRNA$lncRNA.NAME
  
  padj <- lnc.mRNA.ccor$p.adj
  colnames(padj) = mRNA$mRNA.NAME
  rownames(padj) = lncRNA$lncRNA.NAME
  
  rho2 <- rho %>% data.frame(check.names = F) %>% rownames_to_column() %>% gather(mRNA.Symbol,Pearson.Correlation,-rowname) %>%
    rename(lncRNA.Symbol=1)
  
  padj2 <- padj %>% data.frame(check.names = F) %>% rownames_to_column() %>% gather(mRNA.Symbol,Pearson.FDR,-rowname) %>%
    rename(lncRNA.Symbol=1)
  
  rho.padj <- cbind.data.frame(rho2,data.frame(Pearson.FDR=padj2$Pearson.FDR)) %>%
    mutate(r = cut(abs(Pearson.Correlation), breaks = c(-Inf, 0.3, 0.6, Inf),
                   labels = c("0~0.3", "0.3~0.6", ">=0.6"),
                   right = FALSE),
           p.adj = cut(Pearson.FDR, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE)) %>%
    mutate(r=factor(r,levels = c("0~0.3", "0.3~0.6", ">=0.6"))) %>%
    mutate(Cor.Direction = if_else(Pearson.Correlation >= 0,"Positive","Negative")) %>%
    mutate(Cor.Direction = factor(Cor.Direction,levels = c("Positive","Negative")))
  
  mRNA.ccor <- corr.test(mRNA.Exp, use = "pairwise",method="pearson",adjust="fdr", alpha=.05,ci=TRUE)
  
  Rho = mRNA.ccor$r
  colnames(Rho) = rownames(Rho) = mRNA$mRNA.NAME
  
  colnames(mRNA.Exp) = mRNA$mRNA.NAME
  
  p1<-quickcor(mRNA.Exp, type = "upper") + geom_square(inherit.aes = T) +
    anno_link(rho.padj, mapping = aes(colour = p.adj, size = r,linetype = Cor.Direction)) +
    scale_size_manual(values = c(0.1, 0.5, 1)) 
  p2<-p1+scale_fill_gradient2(midpoint = 0, low = brewer.pal(8,"Set1")[2], mid = "white",
                              high = brewer.pal(8,"Set1")[1], space = "Lab" )
  p3 <- p2+scale_color_manual(values=c(brewer.pal(8,"Dark2")[4:2],"grey60"))
  p4 <- p3 + scale_linetype_manual(values = c(1,2))
  p5 <- p4+guides(size=guide_legend(title="Pearson's CC",override.aes=list(colour="grey35"),order=1),
                  colour=guide_legend(title="Pearson's FDR",override.aes = list(size=3),order=2),
                  linetype = guide_legend(title="Correlativity",override.aes = list(size=0.6),order=3),
                  fill=guide_colorbar(title="Pearson's CC",order=4))
  
  print(p5)
  graph2pdf(file=paste("lncRNA.mRNA.Correlation.UpDElncRNA.DownDEmiRNA.UpDEmRNA.",miRNA,".pdf",sep = ""),height = 6,width = 12)
  
  save(lnc.mRNA.ccor,mRNA.ccor,rho.padj,Rho,file = paste("lncRNA.mRNA.Correlation.UpDElncRNA.DownDEmiRNA.UpDEmRNA.",miRNA,".Rdata",sep = ""))
}

load(file = paste("lncRNA.mRNA.Correlation.UpDElncRNA.DownDEmiRNA.UpDEmRNA.","hsa-miR-424-5p",".Rdata",sep = ""))


######################## U:DElncRNA DEmiRNA DEmRNA => Pearson Correlation -> Only Tumor Samples ################
setwd("K:/TCGA/Anlysis/LIHC")
#load("ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata")
library(ggcor)
library(dplyr)
library(ggplot2)
library(psych)
library(RColorBrewer)
library(export)

load(file = "LIHC.GSE140845.miRNA.Intersect.Rdata")
load(file = "LIHC.GSE140845.lncRNA.Rdata")
load(file = "LIHC.GSE140845.mRNA.Rdata")
load("ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata") #Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA,Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA
load("Inter.miRNA.mRNA.Samples.361.Rdata")

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1) %>% dplyr::select(Inter.Sample.361)


mid.data <- data.frame(SYMBOL=c("MCM2","CBX2","CEP55","MCM3AP-AS1","DUXAP8","lnc-RGS5-1","CDKN2B-AS1"),
                       ENSEMBL = c("ENSG00000073111","ENSG00000173894","ENSG00000138180",
                                   "ENSG00000215424","ENSG00000206195","ENSG00000232995",
                                   "ENSG00000240498"))

Target.Symbol <- mRNA.Exp.log2.TPM[mid.data$ENSEMBL,]
rownames(Target.Symbol)=mid.data$SYMBOL
write.csv(t(Target.Symbol),file = "Source.Data.1H.csv",row.names = T)
library(corrplot) 
library(ggcor)
library(RColorBrewer)

p <- quickcor(t(Target.Symbol), cor.test = TRUE)+  ## cor.test 
  geom_square(data = get_data(type = "lower", show.diag = FALSE))+ 
  geom_mark(data = get_data(type = "upper", show.diag = FALSE), size =2.5)+ 
  geom_abline(slope = -1, intercept = 8)+
  scale_fill_gradientn(colors = brewer.pal(9,"Reds"))
ggsave(p,filename = "Hub.CeRNA.Cor.pdf",height = 4,width = 4)

miRNA.Exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
miRNA.Exp.log2.TPM <- log2(miRNA.Exp+1) %>% dplyr::select(Inter.Sample.361)

### miRNA mRNA => 361 Tumur samples 
GeneType[GeneType$HGNC.symbol =="MCM3AP-AS1",] %>% unique() #ENSG00000215424
#GeneType[GeneType$HGNC.symbol =="MAP2",] %>% unique() #ENSG00000078018
GeneType[GeneType$HGNC.symbol =="MCM2",] %>% unique() #ENSG00000073111
GeneType[GeneType$HGNC.symbol =="LRRC1",] %>% unique() #ENSG00000137269

GeneType[GeneType$HGNC.symbol =="DUXAP8",] %>% unique() #ENSG00000206195

GeneType[GeneType$HGNC.symbol =="lnc-RGS5-1",] %>% unique() #ENSG00000232995

GeneType[GeneType$HGNC.symbol =="CKAP2L",] %>% unique() #ENSG00000169607
GeneType[GeneType$HGNC.symbol =="RMI2",] %>% unique() #ENSG00000175643


GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180

GeneType[GeneType$HGNC.symbol =="ESR1",] %>% unique() #ENSG00000091831
GeneType[GeneType$HGNC.symbol =="CPEB3",] %>% unique() #ENSG00000107864

# CDKN2B-AS1 ENSG00000240498


######### CeRNA constructing -2 => Positive Rho 0.3 of lncRNA and mRNA #############
###### ggalluvial ######
setwd("K:/TCGA/Anlysis/LIHC")
load("ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata") #Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA,Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA
#### Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA => 11 lncRNA 6 miRNA 17 mRNA ####
CorData <- read.csv("Correlation.Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.csv")
CorData2 <- CorData %>% filter(Rho >= 0.3)
DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final <- merge(Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA,CorData2,by=c("miRNA.ID","lncRNA.ENSEMBLE.ID","mRNA.ENSEMBLE.ID"))
write.csv(DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final,"Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv",row.names = F)

#Figure.DownDElncRNA.UpDEmiRNA.DownDEmRNA <- rbind(DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final[,c(1,2,4,5,6,7)] %>% as.matrix(),DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final[,c(1,3,8,5,10,11)] %>% as.matrix())
library(ggalluvial)
plot.data <- DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final %>% dplyr::select(lncRNA.NAME,miRNA.ID,mRNA.NAME) %>% mutate(Freq=1)
df_lodes <- to_lodes_form(plot.data,key ="x", value = "stratum", id = "alluvium",axes =1:3)
mycol3=colorRampPalette(c("#00abef","#64b036","#ffe743","#64b036","#00abef"))(36)
mycol <-rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),3)
library(RColorBrewer)
mycol <- c(brewer.pal(8,"Set1"),brewer.pal(11,"Paired"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"))

p<-ggplot(df_lodes,aes(x = x, stratum =stratum, alluvium = alluvium,
                    fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0.2, knot.pos = 0.1) +
  geom_stratum(alpha = .9,color="grey20",width = 1/5) +
  geom_text(stat = "stratum", size =1.5,color="black") +
  scale_fill_manual(values = mycol) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())+
  guides(fill = FALSE)

ggsave(p,filename = "Figure.DownDElncRNA.UpDEmiRNA.DownDEmRNA.ggalluvial.2.pdf",height = 5,width = 4)

#### Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA => 14 lncRNA 2 miRNA 25 mRNA => 246 ####
setwd("K:/TCGA/Anlysis/LIHC")
CorData <- read.csv("Correlation.Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.csv")
CorData2 <- CorData %>% filter(Rho >= 0.3)
Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final <- merge(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA,CorData2,by=c("miRNA.ID","lncRNA.ENSEMBLE.ID","mRNA.ENSEMBLE.ID"))
write.csv(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final,"Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv",row.names = F)

#Figure.Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA <- rbind(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final[,c(1,2,4,5,6,7)] %>% as.matrix(),DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final[,c(1,3,8,5,10,11)] %>% as.matrix())
library(ggalluvial)
plot.data <- Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final %>% dplyr::select(lncRNA.NAME,miRNA.ID,mRNA.NAME) %>% mutate(Freq=1)
df_lodes <- to_lodes_form(plot.data,key ="x", value = "stratum", id = "alluvium",axes =1:3)
mycol3=colorRampPalette(c("#00abef","#64b036","#ffe743","#64b036","#00abef"))(36)
#mycol <-rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),3)
library(RColorBrewer)
mycol <- c(brewer.pal(8,"Set1"),brewer.pal(11,"Paired"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(7,"Accent"))

p<-ggplot(df_lodes,aes(x = x, stratum =stratum, alluvium = alluvium,
                       fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0.2, knot.pos = 0.1) +
  geom_stratum(alpha = .9,color="grey20",width = 1/5) +
  geom_text(stat = "stratum", size =1.5,color="black") +
  scale_fill_manual(values = mycol) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())+
  guides(fill = FALSE)

ggsave(p,filename = "Figure.Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.ggalluvial.pdf",height = 5,width = 7)






#### Cytoscape ####
Dataset.mid1 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")
Dataset.mid1.miRNA.lncRNA <- Dataset.mid1 %>% dplyr::select(miRNA.ID,lncRNA.NAME,lncRNA.Sig,lncRNA.miRNA.Count) %>%
  dplyr::rename(miRNA=1,Symbol=2,Sig=3,Interaction=4)  %>% mutate(Pos="P",Type="lncRNA") %>% arrange(miRNA,Symbol) %>% unique()

Dataset.mid1.miRNA.mRNA <- Dataset.mid1 %>% dplyr::select(miRNA.ID,mRNA.NAME,mRNA.Sig,Count) %>%
  dplyr::rename(miRNA=1,Symbol=2,Sig=3,Interaction=4) %>% mutate(Pos="P",Type="mRNA") %>% arrange(miRNA,Symbol) %>% unique()

rbind.data.frame(Dataset.mid1.miRNA.lncRNA,Dataset.mid1.miRNA.mRNA) %>%
  write.table(file = "Cytoscape.FilterRho.ceRNA.UpDownUp.txt",sep = "\t",row.names = F,quote = F)

Dataset.mid1 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")
Dataset.mid1.miRNA.lncRNA <- Dataset.mid1 %>% dplyr::select(miRNA.ID,lncRNA.NAME,lncRNA.Sig,lncRNA.miRNA.Count) %>%
  dplyr::rename(miRNA=1,Symbol=2,Sig=3,Interaction=4)  %>% mutate(Pos="P",Type="lncRNA") %>% arrange(miRNA,Symbol) %>% unique()

Dataset.mid1.miRNA.mRNA <- Dataset.mid1 %>% dplyr::select(miRNA.ID,mRNA.NAME,mRNA.Sig,Count) %>%
  dplyr::rename(miRNA=1,Symbol=2,Sig=3,Interaction=4) %>% mutate(Pos="P",Type="mRNA") %>% arrange(miRNA,Symbol) %>% unique()
rbind.data.frame(Dataset.mid1.miRNA.lncRNA,Dataset.mid1.miRNA.mRNA)  %>%
  write.table(file = "Cytoscape.FilterRho.ceRNA.DownUpDown.txt",sep = "\t",row.names = F,quote = F)



##### Boxplot + Violin plot for DElncRNA DEmRNA DEmiRNA ######
setwd("K:/TCGA/Anlysis/LIHC")
dir.create("DE.Boxplot")
Cor.ceRNA1 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")
Cor.ceRNA2 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")

mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor" | Sample.Type == "Solid Tissue Normal")
miRNA.Phenodata <- read.csv("../../miRNA.Mature.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor" | Sample.Type == "Solid Tissue Normal")

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

miRNA.Exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
miRNA.Exp.log2.TPM <- log2(miRNA.Exp+1)

#Inter.Samples <- intersect(intersect(intersect(mRNA.Phenodata$Sample.ID,miRNA.Phenodata$Sample.ID),colnames(mRNA.Exp)),colnames(miRNA.Exp))
#415 Samples including tumor and normal samples

# miRNA 
for (miRNA in unique(c(Cor.ceRNA1$miRNA.ID,Cor.ceRNA2$miRNA.ID))) {
  Index <- which((rownames(miRNA.Exp.log2.TPM)%>% str_remove_all(".*\\.")) == miRNA)
  middata <- data.frame(miRNA.Exp.log2.TPM[Index,]) %>% t.data.frame() %>% data.frame(check.names = F) %>%
    mutate(Sample.ID = colnames(miRNA.Exp.log2.TPM)) %>% 
    merge(.,miRNA.Phenodata %>% dplyr::select(Sample.ID,Sample.Type),by="Sample.ID") %>%
    mutate(Group = if_else(Sample.Type == "Primary Tumor","HCC","Ctrl"))
  
  colnames(middata)[2] = miRNA
  
  p<-ggviolin(middata, "Group", miRNA, fill = "Group",
           palette = brewer.pal(8,"Set1")[1:2],
           add = "boxplot", add.params = list(fill = "white"))+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",miRNA,sep = ""))
  
  ggsave(p,filename = paste("DE.Boxplot/",miRNA,".Violin.Boxplot.pdf",sep = ""),width = 2.5,height = 2.5)
}

# mRNA 
for (mRNA in unique(c(Cor.ceRNA1$mRNA.ENSEMBLE.ID,Cor.ceRNA2$mRNA.ENSEMBLE.ID))) {
  Index <- which((rownames(mRNA.Exp.log2.TPM) == mRNA))
  middata <- data.frame(mRNA.Exp.log2.TPM[Index,]) %>% t.data.frame() %>% data.frame(check.names = F) %>%
    mutate(Sample.ID = colnames(mRNA.Exp.log2.TPM)) %>% 
    merge(.,mRNA.Phenodata %>% dplyr::select(Sample.ID,Sample.Type),by="Sample.ID") %>%
    mutate(Group = if_else(Sample.Type == "Primary Tumor","HCC","Ctrl"))
  
  p<-ggviolin(middata, "Group", mRNA, fill = "Group",
              palette = brewer.pal(8,"Set1")[1:2],
              add = "boxplot", add.params = list(fill = "white"))+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",mRNA,sep = ""))
  
  ggsave(p,filename = paste("DE.Boxplot/mRNA.",mRNA,".Violin.Boxplot.pdf",sep = ""),width = 2.5,height = 2.5)
}

# lncRNA 
for (lncRNA in unique(c(Cor.ceRNA1$lncRNA.ENSEMBLE.ID,Cor.ceRNA2$lncRNA.ENSEMBLE.ID))) {
  Index <- which((rownames(mRNA.Exp.log2.TPM) == lncRNA))
  middata <- data.frame(mRNA.Exp.log2.TPM[Index,]) %>% t.data.frame() %>% data.frame(check.names = F) %>%
    mutate(Sample.ID = colnames(mRNA.Exp.log2.TPM)) %>% 
    merge(.,mRNA.Phenodata %>% dplyr::select(Sample.ID,Sample.Type),by="Sample.ID") %>%
    mutate(Group = if_else(Sample.Type == "Primary Tumor","HCC","Ctrl"))
  
  p<-ggviolin(middata, "Group", lncRNA, fill = "Group",
              palette = brewer.pal(8,"Set1")[1:2],
              add = "boxplot", add.params = list(fill = "white"))+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA,sep = ""))
  
  ggsave(p,filename = paste("DE.Boxplot/lncRNA.",lncRNA,".Violin.Boxplot.pdf",sep = ""),width = 2.5,height = 2.5)
}

##### Pheatmap for top 20 Inter DE genes or miRNA #####
load(file = "LIHC.GSE140845.miRNA.Intersect.Rdata")
load(file = "LIHC.GSE140845.lncRNA.Rdata")
load(file = "LIHC.GSE140845.mRNA.Rdata")
#### lncRNA => 421 = 371+50 ####

mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor" | Sample.Type == "Solid Tissue Normal")
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1) #%>% as.matrix() %>% t() %>% data.frame(check.names = F)

Inter.mRNA.Samples <- intersect(mRNA.Phenodata$Sample.ID,colnames(mRNA.Exp.log2.TPM)) #421
mRNA.Phenodata.mid <- mRNA.Phenodata %>% filter(Sample.ID %in% Inter.mRNA.Samples) %>% arrange(Sample.Type)

load(file = "TCGA-LIHC_DEG.Rdata") #res
Data.set1 <- res[rownames(res) %in% c(Inter.lncRNA.Up,Inter.lncRNA.Down),] %>% data.frame() %>% rownames_to_column("ENSEMBLE.ID") %>%
  arrange(padj) %>% top_n(-21,padj) %>% filter(ENSEMBLE.ID != "ENSG00000255980") %>% 
  arrange(Threshold) 

Data.set2 <- mRNA.Exp.log2.TPM[Data.set1$ENSEMBLE.ID,mRNA.Phenodata.mid$Sample.ID]

Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv")
mid.Gene.Type3 <- Inter.lncRNA.GeneCard  %>% filter(lncRNA.ENSEMBLE.ID %in% rownames(Data.set2)) %>% 
  mutate(lncRNA.ENSEMBLE.ID = factor(lncRNA.ENSEMBLE.ID,levels=rownames(Data.set2))) %>%
  arrange(lncRNA.ENSEMBLE.ID)

rownames(Data.set2) = mid.Gene.Type3$lncRNA.ID

annotation_col = data.frame(Sample.Type = factor(rep(c("Primary Tumor","Solid Tissue Normal"),c(371,50))))
rownames(annotation_col) = colnames(Data.set2)
annotation_row = data.frame(Sig = factor(rep(c("Down","Up"),c(6,14))))
rownames(annotation_row) = rownames(Data.set2)

anno_color=list(Sample.Type = c(`Primary Tumor`=brewer.pal(8,"Set1")[1],`Solid Tissue Normal`=brewer.pal(8,"Set1")[2]),
                Sig = c(Up = brewer.pal(8,"Set1")[5],Down = brewer.pal(8,"Set1")[3]))

pheatmap(Data.set2,
         annotation_col = annotation_col,
         #annotation_row = annotation_row, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),#colorRampPalette(brewer.pal(11, "RdYlBu"))(50), #
         border_color = "black",
         annotation_colors = anno_color,
         #cellwidth = 7,cellheight = 5,
         #fontsize_row = 5,fontsize_col = 6,
         #clustering_distance_cols = "euclidean",
         #clustering_method = "complete",cluster_cols = T,
         annotation_row = annotation_row,
         cluster_rows = T,cluster_cols = T,show_colnames = F)
graph2pdf(file="Pheatmap.DE.lncRNA.pdf",height = 4,width = 7)

#### mRNA => 421 = 371+50 ####
load(file = "TCGA-LIHC_DEG.Rdata")

Data.set1 <- res[rownames(res) %in% c(Inter.mRNA.Up,Inter.mRNA.Down),] %>% data.frame() %>% rownames_to_column("ENSEMBLE.ID") %>%
  arrange(padj) %>% top_n(-20,padj) %>% filter(ENSEMBLE.ID != "ENSG00000255980") %>% 
  arrange(Threshold) 

Data.set2 <- mRNA.Exp.log2.TPM[Data.set1$ENSEMBLE.ID,mRNA.Phenodata.mid$Sample.ID]

mid.Gene.Type3 <- GeneType  %>% filter(Gene.stable.ID %in% rownames(Data.set2)) %>% 
  dplyr::select(Gene.stable.ID,Gene.name) %>% unique() %>%
  mutate(Gene.stable.ID = factor(Gene.stable.ID,levels=rownames(Data.set2))) %>%
  arrange(Gene.stable.ID)

rownames(Data.set2) = mid.Gene.Type3$Gene.name

annotation_col = data.frame(Sample.Type = factor(rep(c("Primary Tumor","Solid Tissue Normal"),c(371,50))))
rownames(annotation_col) = colnames(Data.set2)
annotation_row = data.frame(Sig = factor(rep(c("Down","Up"),c(1,19))))
rownames(annotation_row) = rownames(Data.set2)

anno_color=list(Sample.Type = c(`Primary Tumor`=brewer.pal(8,"Set1")[1],`Solid Tissue Normal`=brewer.pal(8,"Set1")[2]),
                Sig = c(Up = brewer.pal(8,"Set1")[5],Down = brewer.pal(8,"Set1")[3]))

pheatmap(Data.set2,
         annotation_col = annotation_col,
         #annotation_row = annotation_row, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),#colorRampPalette(brewer.pal(11, "RdYlBu"))(50), #
         border_color = "black",
         annotation_colors = anno_color,
         #cellwidth = 7,cellheight = 5,
         #fontsize_row = 5,fontsize_col = 6,
         #clustering_distance_cols = "euclidean",
         #clustering_method = "complete",cluster_cols = T,
         annotation_row = annotation_row,
         cluster_rows = T,cluster_cols = T,show_colnames = F)
graph2pdf(file="Pheatmap.DE.mRNA.pdf",height = 4,width = 7)


#### miRNA => 422 = 372+50 ####
miRNA.Phenodata <- read.csv("../../miRNA.Mature.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor" | Sample.Type == "Solid Tissue Normal")
miRNA.Exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
miRNA.Exp.log2.TPM <- log2(miRNA.Exp+1)#%>% dplyr::select(colnames(.) %in% miRNA.Phenodata$Sample.ID)

Inter.miRNA.Samples <- intersect(miRNA.Phenodata$Sample.ID,colnames(miRNA.Exp.log2.TPM)) #422
miRNA.Phenodata.mid <- miRNA.Phenodata %>% filter(Sample.ID %in% Inter.miRNA.Samples) %>% arrange(Sample.Type)

Index <- rownames(miRNA.Exp.log2.TPM) %>% str_remove_all(".*\\.") %in% c(Inter.miRNA.Down,Inter.miRNA.Up)
Data.set <- miRNA.Exp.log2.TPM[Index,] %>% dplyr::select(miRNA.Phenodata.mid$Sample.ID)
rownames(Data.set) <- rownames(Data.set) %>% str_remove_all(".*\\.")
Data.set2 <- Data.set %>% rownames_to_column() %>% mutate(rowname = factor(.$rowname,levels=c(Inter.miRNA.Down,Inter.miRNA.Up))) %>%
  arrange(rowname) %>% column_to_rownames("rowname")

library(pheatmap)
library(export)
annotation_col = data.frame(Sample.Type = factor(rep(c("Primary Tumor","Solid Tissue Normal"),c(372,50))))
rownames(annotation_col) = colnames(Data.set2)
annotation_row = data.frame(Sig = factor(rep(c("Down","Up"),c(2,6))))
rownames(annotation_row) = rownames(Data.set2)

anno_color=list(Sample.Type = c(`Primary Tumor`=brewer.pal(8,"Set1")[1],`Solid Tissue Normal`=brewer.pal(8,"Set1")[2]),
                Sig = c(Up = brewer.pal(8,"Set1")[5],Down = brewer.pal(8,"Set1")[3]))

pheatmap(Data.set2,
         annotation_col = annotation_col,
         #annotation_row = annotation_row, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = "black",
         annotation_colors = anno_color,
         annotation_row = annotation_row,
         #cellwidth = 7,cellheight = 5,
         #fontsize_row = 5,fontsize_col = 6,
         #clustering_distance_cols = "euclidean",
         #clustering_method = "complete",cluster_cols = T,
         cluster_rows = T,cluster_cols = T,show_colnames = F)
graph2pdf(file="Pheatmap.DE.miRNA.pdf",height = 4,width = 7)



###### Unique KM COX ######
lncRNA.361 <- read.csv("361.lncRNA.KM.COX.csv")
miRNA.361 <- read.csv("361.miRNA.KM.COX.csv")
mRNA.361 <- read.csv("361.mRNA.KM.COX.csv")
Cor.ceRNA1 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")
Cor.ceRNA2 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")

#### lncRNA => 24 ceRNA => 6 Up ####

lncRNA.361.km.cox.sig <- lncRNA.361 %>% filter(KM <= 0.05 & pvalue <= 0.05) %>% 
  filter((id %in% Cor.ceRNA1$lncRNA.ENSEMBLE.ID) | (id %in% Cor.ceRNA2$lncRNA.ENSEMBLE.ID))
write.csv(lncRNA.361.km.cox.sig,file = "lncRNA.361.km.cox.sig.csv",row.names = F)

mRNA.361.km.cox.sig <- mRNA.361 %>% filter(KM <= 0.05 & pvalue <= 0.05) %>% 
  filter((id %in% Cor.ceRNA1$mRNA.ENSEMBLE.ID) | (id %in% Cor.ceRNA2$mRNA.ENSEMBLE.ID))
write.csv(mRNA.361.km.cox.sig,file = "mRNA.361.km.cox.sig.csv",row.names = F)

#### mRNA => 42 ceRNA => 17 Up + 3 Down ####

mid.data <- merge(mRNA.361.km.cox.sig,res%>%rownames_to_column("id"),by="id")
write.csv(mid.data,file = "mRNA.361.km.cox.sig.csv",row.names = F)

####



###### mRNA function #######
#### GO ####
ceRNA1 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")
ceRNA2 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")

library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
#Dm_msigdbr <- msigdbr(species="Homo sapiens")
#head(Dm_msigdbr, 2) %>% as.data.frame
#DmGO <- msigdbr(species="Drosophila melanogaster",category="C5") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)
#em <- enricher(ceRNA1$mRNA.NAME,TERM2GENE=DmGO[,c(1,3)])
#data.frame(em) %>% View()

ego <- enrichGO(ceRNA1$mRNA.ENSEMBLE.ID, OrgDb=org.Hs.eg.db, ont='ALL',
                  pAdjustMethod='BH', pvalueCutoff=0.05, 
                  qvalueCutoff=0.2, keyType='ENSEMBL')

write.csv(ego,file = "Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.EnrichGO.csv")

ego <- enrichGO(ceRNA2$mRNA.ENSEMBLE.ID, OrgDb=org.Hs.eg.db, ont='ALL',
                pAdjustMethod='BH', pvalueCutoff=0.05, 
                qvalueCutoff=0.2, keyType='ENSEMBL')
write.csv(ego,file = "Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.EnrichGO.csv")

## Figure barplot
data <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.EnrichGO.csv")
data.mid <- data[1:12,]
library(ggthemes)
p<-ggplot(data.mid,aes(reorder(Description,-log10(qvalue)),y=-log10(qvalue),fill = Description))+
  geom_col()+theme_few()+
  coord_flip()+
  theme(legend.position = "none")+
  labs(x="")+
  scale_y_continuous(expand = c(0,0))+
  geom_hline(yintercept = c(-log10(0.05)))

ggsave(p,filename = "Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.EnrichGO.pdf",width = 4,height = 3)

## 
data <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.EnrichGO.csv")
#data.mid <- data[1:12,]
library(ggthemes)
p<-ggplot(data,aes(reorder(Description,-log10(qvalue)),y=-log10(qvalue),fill = Description))+
  geom_col()+theme_few()+
  coord_flip()+
  theme(legend.position = "none")+
  labs(x="")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  geom_hline(yintercept = c(-log10(0.05)))

ggsave(p,filename = "Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.EnrichGO.pdf",width = 6,height = 2.5)

#### GSEA FOR H	hallmark gene sets ####
library(msigdbr)
Dm_msigdbr <- msigdbr(species="Homo sapiens")
head(Dm_msigdbr, 2) %>% as.data.frame
DmGO <- msigdbr(species="Homo sapiens",category="H") %>% dplyr::select(gs_name, entrez_gene, gene_symbol)

### ceRNA1
#em <- enricher(ceRNA1$mRNA.NAME,TERM2GENE=DmGO[,c(1,3)])
#data.frame(em) %>% View()
Part.Hallmark <- DmGO %>% dplyr::filter(gene_symbol %in% ceRNA1$mRNA.NAME)
Part.Hallmark <- Part.Hallmark[1:17,]
Part.Hallmark.Plot.Count <- Part.Hallmark %>% data.frame() %>% group_by(gs_name) %>% summarise(Count = n())
#Part.Hallmark.Plot.Label <- Part.Hallmark %>% data.frame() %>% group_by(gs_name) %>% mutate(Label = paste())
p<-ggplot(Part.Hallmark.Plot.Count,aes(reorder(gs_name,Count),Count,fill=gs_name))+
  geom_col()+theme_few()+
  coord_flip()+
  theme(legend.position = "none")+
  labs(x="",y="Gene count") +
  scale_y_continuous(expand = c(0,0),limits = c(0,3.5))
ggsave(p,filename = "Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.msigdb.H.pdf",width = 6.5,height = 3)

### ceRNA2
#em <- enricher(ceRNA2$mRNA.NAME,TERM2GENE=DmGO[,c(1,3)])
#data.frame(em) %>% View()
Part.Hallmark <- DmGO %>% dplyr::filter(gene_symbol %in% ceRNA2$mRNA.NAME)
#Part.Hallmark <- Part.Hallmark[1:17,]
Part.Hallmark.Plot.Count <- Part.Hallmark %>% data.frame() %>% group_by(gs_name) %>% summarise(Count = n())
p<-ggplot(Part.Hallmark.Plot.Count,aes(reorder(gs_name,Count),Count,fill=gs_name))+
  geom_col()+theme_few()+
  coord_flip()+
  theme(legend.position = "none")+
  labs(x="",y="Gene count") +
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5))
ggsave(p,filename = "Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.msigdb.H.pdf",width = 5.5,height = 2.5)

#### PPI ####
#string database

######################## Prognostic model -> Candidate ##################
###### Prognostic model => filter out gene ########
#### Datasets Process #####
load("ceRNA.DElncRNA.DEmiRNA.DEmRNA.Rdata") #Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA; Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(Exp+1)

miRNA.Exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
miRNA.Exp.log2.TPM <- log2(miRNA.Exp+1)

mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC") %>%
  filter(Sample.ID %in% colnames(mRNA.Exp)) %>% filter(Sample.Type == "Primary Tumor") #371

miRNA.Phenodata <- read.csv("../../miRNA.Mature.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC") %>%
  filter(Sample.ID %in% colnames(miRNA.Exp)) %>% filter(Sample.Type == "Primary Tumor") #371

#miRNA.Phenodata$Sample.ID %in% mRNA.Phenodata$Sample.ID

Intet.Sample.mRNA.miRNA <- intersect(miRNA.Phenodata$Sample.ID,mRNA.Phenodata$Sample.ID) #367


SurvivalData <- read.table("survival_LIHC_survival.txt",header = T,row.names = 1,sep = '\t') 
mid.surv <- SurvivalData[substr(Intet.Sample.mRNA.miRNA,1,15),] %>% dplyr::select(OS,OS.time) %>% na.omit() %>%
  filter(OS.time != 0)#361

#Target gene
DEmiRNA <- unique(c(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA$miRNA.ID,Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA$miRNA.ID))
DElncRNA<- rbind.data.frame(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA %>% dplyr::select(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% arrange(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% unique(),
                            Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA %>% dplyr::select(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% arrange(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% unique()) %>%
  arrange(lncRNA.ENSEMBLE.ID,lncRNA.NAME) %>% unique() #38

DEmRNA <- rbind.data.frame(Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA %>% dplyr::select(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% arrange(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% unique(),
                           Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA %>% dplyr::select(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% arrange(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% unique())%>%
  arrange(mRNA.ENSEMBLE.ID,mRNA.NAME) %>% unique() #61

DEmRNA.DElncRNA.Exp <- mRNA.Exp.log2.TPM[c(DElncRNA$lncRNA.ENSEMBLE.ID,DEmRNA$mRNA.ENSEMBLE.ID),] %>%
  dplyr::select(Intet.Sample.mRNA.miRNA[substr(Intet.Sample.mRNA.miRNA,1,15) %in% rownames(mid.surv)]) %>% t() %>% data.frame(check.names = F)

MIRNAINDEX <- (rownames(miRNA.Exp.log2.TPM) %>%str_remove_all(".*\\.")) %in% DEmiRNA
DEmiRNA.Exp <- miRNA.Exp.log2.TPM[MIRNAINDEX,]%>%
  dplyr::select(Intet.Sample.mRNA.miRNA[substr(Intet.Sample.mRNA.miRNA,1,15) %in% rownames(mid.surv)]) %>% t() %>% data.frame(check.names = F)

DElncRNA.DEmiRNA.DEmRNA.Exp <- cbind.data.frame(DEmRNA.DElncRNA.Exp,DEmiRNA.Exp)

save(DElncRNA.DEmiRNA.DEmRNA.Exp,mid.surv,file = "Prognostic.model.Rdata")

#### glmnet => LASSO model ######
library(glmnet)

load(file = "Prognostic.model.Rdata")

lncRNA.361.km.cox.sig <- read.csv("lncRNA.361.km.cox.sig.csv")
mRNA.361.km.cox.sig <- read.csv("mRNA.361.km.cox.sig.csv")

DElncRNA.DEmiRNA.DEmRNA.Exp <- DElncRNA.DEmiRNA.DEmRNA.Exp %>% dplyr::select(c(lncRNA.361.km.cox.sig$id,mRNA.361.km.cox.sig$id,contains("miR")))

x <- as.matrix(DElncRNA.DEmiRNA.DEmRNA.Exp) #gene => column
y <- Surv(as.double(mid.surv$OS.time),as.double(mid.surv$OS))

Select.Variables <- c()
for (i in 1:1000) {
  #fit <- glmnet(x, y, family="cox")
  #plot(fit,xvar = "lambda")
  #using 10-fold CV to select lambda:
  cv.fit =cv.glmnet(x, y, family="cox", nfolds=10)
  #cv.fit$lambda.min
  fit_regL=cv.fit
  fit_regL_CVdev=cv.fit
  
  # final Lasso model:
  model_lasso_min <- glmnet(x=x, y=y,family="cox",lambda=cv.fit$lambda.min)
  fit_regL_coef <- coef(fit_regL, s=fit_regL_CVdev$lambda.min)
  ### Plot Cross-Validation LASSO model
  
  #par(mfrow=c(1,2))
  #plot(fit_regL_CVdev,las=1, main="Lasso fit CV")
  #abline(v=log(fit_regL_CVdev$lambda.min), col="orange", lwd=2)
  #text(log(fit_regL_CVdev$lambda.min), 1.4, paste("lambda.min=",round(fit_regL_CVdev$lambda.min,4),"\n",length(fit_regL_coef@i), " genes" ,sep=""), col="orange", cex=0.75, pos=4)
  
  ### Plot lambda fit
  #plot(fit_regL_CVdev, xvar="lambda", main="Lasso coefficient")
  #abline(v=log(fit_regL_CVdev$lambda.min), col="orange", lwd=2)
  ###############
  #summary(cv.fit)
  #print(cv.fit)
  #cv.fit$lambda.min
  #cv.fit$lambda.1se
  #coef(cv.fit, s = "lambda.min") #----------candidates
  #coef.min= coef(cv.fit, s = "lambda.min")
  select.varialbes = rownames(as.data.frame(which(coef(cv.fit, s = "lambda.min")[,1]!=0)))
  #select.varialbes = select.varialbes[-1] #remove intercept
  #select.varialbes
  Select.Variables<-c(Select.Variables,select.varialbes)
}

SV.LASSO.COX <- table(Select.Variables) %>% data.frame()

DEmiRNA.LASSO.COX <- SV.LASSO.COX[str_detect(SV.LASSO.COX$Select.Variables,"miR"),] %>%
  mutate(Gene.type = "miRNA",Gene.name = str_remove_all(Select.Variables,".*\\.")) %>% 
  filter(Freq >= 10)

DEmRNA.DElncRNA.DEmiRNA.LASSO.COX <- SV.LASSO.COX %>% merge.data.frame(.,GeneType,by.x="Select.Variables",by.y="Gene.stable.ID") %>% 
  dplyr::select(Select.Variables,Freq,Gene.type,Gene.name) %>%
  unique() %>% filter(Freq >= 10) %>% rbind(DEmiRNA.LASSO.COX)

save(DEmRNA.DElncRNA.DEmiRNA.LASSO.COX,file = "DEmRNA.DElncRNA.DEmiRNA.LASSO.COX.ValidateCeRNA.Rdata")
write.csv(DEmRNA.DElncRNA.DEmiRNA.LASSO.COX,"DEmRNA.DElncRNA.DEmiRNA.LASSO.COX.ValidateCeRNA.csv")

#### Random forest => Importance for DElncRNA,DEmiRNA,DEmRNA ########
setwd("K:/TCGA/Anlysis/LIHC")
load("Prognostic.model.Rdata")
library(survival)
library(randomForestSRC)

UseData <- cbind.data.frame(DElncRNA.DEmiRNA.DEmRNA.Exp,mid.surv)

##
## Minimal depth variable selection
## survival analysis
## use larger node size which is better for minimal depth
## 

pbc.obj <- rfsrc(Surv(OS.time,OS) ~ ., UseData, nodesize = 20, importance = TRUE)
vs.pbc <- var.select(object = pbc.obj, conservative = "high")
topvars <- vs.pbc$topvars

Vectors <- max.subtree(pbc.obj)$topvars %>% data.frame() %>% dplyr::rename(ID=1)

# the above is equivalent to
max.subtree(pbc.obj)$topvars
# different levels of conservativeness
#var.select(object = pbc.obj, conservative = "low")
#vs.pbc<-var.select(object = pbc.obj, conservative = "medium")
#vs.pbc$topvars
var.select(object = pbc.obj, conservative = "high")

DEmiRNA.RF <- Vectors[str_detect(Vectors$ID,"miR"),] %>% data.frame() %>% dplyr::rename(ID=1) %>%
  mutate(Gene.type = "miRNA",Gene.name = str_remove_all(ID,".*\\."))

DEmRNA.DElncRNA.DEmiRNA.RF <- Vectors %>% merge.data.frame(.,GeneType,by.x="ID",by.y="Gene.stable.ID") %>% 
  dplyr::select(ID,Gene.type,Gene.name) %>%
  unique() %>% rbind(DEmiRNA.RF)

save(DEmRNA.DElncRNA.DEmiRNA.RF,file = "DEmRNA.DElncRNA.DEmiRNA.RF.ValidateCeRNA.Rdata")
write.csv(DEmRNA.DElncRNA.DEmiRNA.RF,"DEmRNA.DElncRNA.DEmiRNA.RF.ValidateCeRNA.csv",row.names = F)

if (F) {
  ## 
  ## Variable hunting high-dimensional example
  ## using weights from univarate cox p-values
  ## nrep is small for illustration; typical values are nrep = 100
  ##
  ###
  if (library("survival", logical.return = TRUE)){
    cox.weights <- function(rfsrc.f, rfsrc.data) {
      event.names <- all.vars(rfsrc.f)[1:2]
      p <- ncol(rfsrc.data) - 2
      event.pt <- match(event.names, names(rfsrc.data))
      xvar.pt <- setdiff(1:ncol(rfsrc.data), event.pt)
      sapply(1:p, function(j) {
        cox.out <- coxph(rfsrc.f, rfsrc.data[, c(event.pt, xvar.pt[j])])
        pvalue <- summary(cox.out)$coef[5]
        if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
      })
    }
    data_rf=UseData  #ËæìÂÖ•Êï∞ÊçÆ
    rfsrc.f <- as.formula(Surv(OS.time, OS) ~ .)  #ËæìÂÖ•yÔºåËÆ∞ÂæóÊîπPFSÂíåPFS_TÁöÑÂàó???
    cox.wts <- cox.weights(rfsrc.f, data_rf)
    b.obj <- rfsrc(rfsrc.f, data_rf , nsplit = 10, xvar.wt = cox.wts,importance = "random",na.action ="na.impute",ntree = 1000)
    vh.breast.cox <- var.select(rfsrc.f, data_rf, method = "vh", nstep = 5,nrep = 10, xvar.wt = cox.wts,conservative = "high")
  }
  
  random.plot=plot(b.obj)
  vimp=as.data.frame(b.obj$importance)
  vimp$symbol=rownames(vimp)
  vimp=vimp[order(vimp$`b.obj$importance`,decreasing = T),]
  vimp <- vimp %>% filter(`b.obj$importance` > 0)
  write.csv(vimp,"DEmRNA.DElncRNA.DEmiRNA.RF.Cox.Coef.csv",row.names = T)
}


###### Prognostic model => Confirm Candidate gene #########
Cor.ceRNA1 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")
Cor.ceRNA2 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")

LASSO <- read.csv("DEmRNA.DElncRNA.DEmiRNA.LASSO.COX.ValidateCeRNA.csv",row.names = 1)
RF <- read.csv("DEmRNA.DElncRNA.DEmiRNA.RF.ValidateCeRNA.csv")
RF <- RF %>% merge(.,mid.Gene.Type,by.x = "ID",by.y="Gene.stable.ID",all.x=T)

library(ggvenn)
library("RColorBrewer")
X <- list(LASSO.COX = LASSO$Select.Variables,
          RSF = RF$ID)

p3 <- ggvenn(X,set_name_size = 3,text_size =3,fill_color = brewer.pal(8,"Set2")[c(1,2)])

ggsave(p3,filename = "DEmRNA.DElncRNA.DEmiRNA.RSF.LASSOCox.Venn.ValidateCeRNA.pdf",height = 3,width = 5)

Intersect.V <- LASSO[LASSO$Select.Variables %in% intersect(LASSO$Select.Variables,RF$ID),]
write.csv(Intersect.V,file = "DEmRNA.DElncRNA.DEmiRNA.RSF.LASSOCox.Intersect.ValidateCeRNA.csv",row.names = F)

####  Multi Cox => HR for DElncRNA,DEmiRNA,DEmRNA => Interseted 7 DE ####

#### Clinical Relevance 
meta <- read.csv("clinical.cart.2021-08-19/clinical.tsv",sep = "\t",header = T)
#meta <- column_to_rownames(meta,var = "case_submitter_id")
meta=meta[,colnames(meta) %in% c("case_submitter_id",
                                 "gender",
                                 "ajcc_pathologic_t",
                                 "ajcc_pathologic_n",
                                 "ajcc_pathologic_m",
                                 "age","Stage")] %>% unique()

## 
load(file = "Prognostic.model.Rdata") #DElncRNA.DEmiRNA.DEmRNA.Exp,mid.surv

mid.surv2 <-  meta[meta$case_submitter_id %in% substr(rownames(mid.surv),1,12),] %>% data.frame() %>%
  mutate(case_submitter_id = factor(case_submitter_id,levels=substr(rownames(mid.surv),1,12))) %>%
  arrange(case_submitter_id) %>% cbind(mid.surv)

Intersect.V <- read.csv("DEmRNA.DElncRNA.DEmiRNA.RSF.LASSOCox.Intersect.ValidateCeRNA.csv")
#
ID = unique(c(LASSO$Select.Variables,RF$ID))
Gene.Name <- unique(c(LASSO$Gene.name,RF$Gene.name))
Intersect.V.Total <- data.frame(Select.Variables = ID,Gene.name = Gene.Name)
rownames(Intersect.V.Total) = ID

Export.Data <<- data.frame()
for (i in 2:15) {
  print(i)
  mid <- as.data.frame(t(combn(Intersect.V.Total$Select.Variables,i)))
  for (j in 1:dim(mid)[1]) {
    print(j)
    DElncRNA.DEmiRNA.DEmRNA.Exp.Target <- DElncRNA.DEmiRNA.DEmRNA.Exp[,as.character(mid[j,])]
    
    mid.name <- Intersect.V.Total[as.character(mid[j,]),]
    
    colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target) = mid.name$Gene.name
    colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target) = str_replace_all(colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target),"-",".")
    
    
    mid.surv3 <- mid.surv2[substr(rownames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target),1,15),]
    
    ## multi cox dataset
    Usedata <- cbind(mid.surv3,DElncRNA.DEmiRNA.DEmRNA.Exp.Target)
    
    ## formaula
    formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target),collapse = "+")))
    multi_variate_cox <- coxph(formula,data=Usedata)
    
    ## check PH hypothesis
    ph_hypo_multi <- cox.zph(multi_variate_cox)
    ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
    
    if (sum(ph_hypo_table[,3] <= 0.05) == i) {
      Export.Data <- rbind.data.frame(Export.Data,c(paste(as.character(mid[j,]),collapse = ","),NA,NA,NA,NA))
    }else{
      formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(rownames(ph_hypo_table)[ph_hypo_table[,3] > 0.05],collapse = "+")))
      multi_variate_cox_2 <- coxph(formula,data=Usedata)
      
      ## consider the correlation
      library(rms)
      #vif <- rms::vif(multi_variate_cox_2)
      #sqrt(vif) < 2
      
      #candidate_genes_for_cox2 <- names(vif[sqrt(vif) < 2])
      
      #formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(candidate_genes_for_cox2,collapse = "+")))
      #multi_variate_cox_2 <- coxph(formula,data=Usedata)
      
      
      # C-index=0.5ÂàôÂÆåÂÖ®ÈöèÊú∫Ôºå=1‰∏∫È¢ÑÊµã‰∏éÂÆûÈôÖÂÆåÂÖ®‰∏Ä???
      #ggforest(model = multi_variate_cox_2,data = Usedata,main = "Hazard ratios of candidate genes",fontsize = 1)
      
      # 
      riskscore <- function(survival_cancer_df,candidate_genes_for_cox,cox_report){
        library(dplyr)
        risk_score_table <- survival_cancer_df %>% dplyr::select(candidate_genes_for_cox)
        for (each_sig_gene in colnames(risk_score_table)) {
          risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
        }
        risk_score_table <- cbind(risk_score_table, 'total_risk_score'=exp(rowSums(risk_score_table))) %>%
          cbind(survival_cancer_df[,c('OS.time','OS')])
        risk_score_table <- risk_score_table[,c('OS.time','OS',candidate_genes_for_cox,'total_risk_score')]
        risk_score_table
      }
      
      candidate_genes_for_cox2 <- c(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05])
      risk_score_table_multi_cox2 <- riskscore(Usedata,candidate_genes_for_cox2,multi_variate_cox_2)
      
      ## ROC 1 3 5
      library('survivalROC')
      multi_ROC <- function(time_vector, risk_score_table){
        
        single_ROC <- function(single_time){
          for_ROC <- survivalROC(Stime = risk_score_table$OS.time,
                                 status = risk_score_table$OS,
                                 marker = risk_score_table$total_risk_score,
                                 predict.time = single_time,method = 'KM')
          data.frame('True_positive'=for_ROC$TP,'False_positive'=for_ROC$FP,'Cut_values'=for_ROC$cut.values,
                     'Time_point'=rep(single_time, length(for_ROC$TP)),'AUC'=rep(for_ROC$AUC,length(for_ROC$TP)))
        }
        multi_ROC_list <- lapply(time_vector, single_ROC)
        do.call(rbind, multi_ROC_list)
      }
      # calculate 1,3,5
      for_multi_ROC <- multi_ROC(time_vector = c(365*c(1,3,5)), risk_score_table = risk_score_table_multi_cox2)
      AUC.1 <- (for_multi_ROC[for_multi_ROC$Time_point==365,])$AUC %>% unique()
      AUC.3 <- (for_multi_ROC[for_multi_ROC$Time_point==1095,])$AUC %>% unique()
      AUC.5 <- (for_multi_ROC[for_multi_ROC$Time_point==1825,])$AUC %>% unique()
      
      Export.Data <- rbind.data.frame(Export.Data,c(paste(as.character(mid[j,]),collapse = ","),multi_variate_cox_2$concordance["concordance"],AUC.1,AUC.3,AUC.5))
    }
  }
}

#### confirm the best combinations ####
#Export.Data <- read.csv("Export.Data.csv") %>% dplyr::rename(Combinations=1, C.index=2, AUC.1=3, AUC.3=4, AUC.5=5)
#Export.Data[is.na(Export.Data)] = 0

#Export.Data.Target <- Export.Data[which(Export.Data$AUC.1 == max(Export.Data$AUC.1)),]

DElncRNA.DEmiRNA.DEmRNA.Exp.Target <- DElncRNA.DEmiRNA.DEmRNA.Exp[,Intersect.V$Select.Variables]#[,unlist(as.character(Export.Data.Target[2,1]) %>% str_split(","))]

#mid.name <- Intersect.V.Total[unlist(as.character(Export.Data.Target[2,1]) %>% str_split(",")),]

colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target) = Intersect.V$Gene.name#mid.name$Gene.name
colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target) = str_replace_all(colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target),"-",".")

mid.surv3 <- mid.surv2[substr(rownames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target),1,15),]

## multi cox dataset
Usedata <- cbind(mid.surv3,DElncRNA.DEmiRNA.DEmRNA.Exp.Target)

## formaula
formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(colnames(DElncRNA.DEmiRNA.DEmRNA.Exp.Target),collapse = "+")))
multi_variate_cox <- coxph(formula,data=Usedata)

## check PH hypothesis
ph_hypo_multi <- cox.zph(multi_variate_cox)
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]

formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(rownames(ph_hypo_table)[ph_hypo_table[,3] > 0.05],collapse = "+")))
multi_variate_cox_2 <- coxph(formula,data=Usedata)
summary(multi_variate_cox_2)
## consider the correlation
library(rms)
#vif <- rms::vif(multi_variate_cox_2)
#sqrt(vif) < 2
#candidate_genes_for_cox2 <- names(vif[sqrt(vif) < 2])
#formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(candidate_genes_for_cox2,collapse = "+")))
#multi_variate_cox_2 <- coxph(formula,data=Usedata)

library(survival)
library(survminer)
# C-index=0.5ÂàôÂÆåÂÖ®ÈöèÊú∫Ôºå=1‰∏∫È¢ÑÊµã‰∏éÂÆûÈôÖÂÆåÂÖ®‰∏Ä???
ggforest(model = multi_variate_cox_2,data = Usedata,main = "Hazard ratios of candidate genes",fontsize = 1)

# 
riskscore <- function(survival_cancer_df,candidate_genes_for_cox,cox_report){
  library(dplyr)
  risk_score_table <- survival_cancer_df %>% dplyr::select(candidate_genes_for_cox)
  for (each_sig_gene in colnames(risk_score_table)) {
    risk_score_table[,each_sig_gene] <- risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table, 'total_risk_score'=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c('OS.time','OS')])
  risk_score_table <- risk_score_table[,c('OS.time','OS','total_risk_score')]
  risk_score_table
}

candidate_genes_for_cox2 <- c(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(Usedata,candidate_genes_for_cox2,multi_variate_cox_2)

risk_score_table_multi_cox2$RiskScore <- if_else(risk_score_table_multi_cox2$total_risk_score > median(risk_score_table_multi_cox2$total_risk_score),"High","Low")

save(risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table,file = "Survival.Signature.RiskScore.Rdata")
load("Survival.Signature.RiskScore.Rdata")

tdmultiCoxSum <- summary(multi_variate_cox_2)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CIH=tdmultiCoxSum$conf.int[,"upper .95"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"])
outResult=cbind(id=row.names(outResult),outResult)
write.csv(outResult,file = "SouceData.Figure2B.csv")
###

fit<-survfit(Surv(OS.time,OS)~RiskScore,data=risk_score_table_multi_cox2)
#p
write.csv(risk_score_table_multi_cox2,file = "SourceData.Risk.csv")
data.survdiff=survdiff(Surv(OS.time,OS)~RiskScore,data=risk_score_table_multi_cox2)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

library(survminer)
library(survival)
ggsurvplot(fit,
           pval=TRUE,#Ê∑ªÂä†P???
           conf.int = TRUE,#Â¢ûÂä†ÁΩÆ‰ø°Âå∫Èó¥
           risk.tabel.col="strata",  #change risk tabel color by groups
           linetype = "strata",      #change line type by groups
           surv.median.line = "hv",   #specify median survival Â¢ûÂä†‰∏≠‰ΩçÁîüÂ≠òÊó∂Èó¥
           ggtheme = ggthemes::theme_few(),      #change ggplots theme
           palette = "Set1",
           #palette = c("#E7B800","#2E9FDF"),#Ëá™ÂÆö‰πâË∞ÉËâ≤Êùø
           risk.table = TRUE,#add risk tabelÁªòÂà∂Á¥ØËÆ°È£éÈô©Êõ≤Á∫ø
           xlab = "Follow up time(d)", # ÊåáÂÆöxËΩ¥Ê†á???
           #legend = c(0.8,0.75), # ÊåáÂÆöÂõæ‰æã‰ΩçÁΩÆ
           legend.title = "", # ËÆæÁΩÆÂõæ‰æãÊ†áÈ¢ò
           legend.labs = c("Riskscore high", "Riskscore low"), # ÊåáÂÆöÂõæ‰æãÂàÜÁªÑÊ†áÁ≠æ
           break.x.by = 1000
           #fun="event"
)
library(export)
graph2pdf(file = "Survival.Signiture.KM.COX.pdf",height = 5,width=4)

res.cox<-coxph(Surv(OS.time,OS)~total_risk_score,data=risk_score_table_multi_cox2)
res.cox.summary <- summary(res.cox)


#### time dependent survival ROC 1 3 5 year ####
library('survivalROC')
multi_ROC <- function(time_vector, risk_score_table){
  
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$OS.time,
                           status = risk_score_table$OS,
                           marker = risk_score_table$total_risk_score,
                           predict.time = single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP,'False_positive'=for_ROC$FP,'Cut_values'=for_ROC$cut.values,
               'Time_point'=rep(single_time, length(for_ROC$TP)),'AUC'=rep(for_ROC$AUC,length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
# calculate 1,3,5
for_multi_ROC <- multi_ROC(time_vector = c(365*c(1,3,5)), risk_score_table = risk_score_table_multi_cox2)
write.csv(for_multi_ROC,file = "SourceData.F2D.csv")
require(plotROC)
library(ggthemes)
for_multi_ROC$Time_point <- factor(for_multi_ROC$Time_point)
p <- ggplot(for_multi_ROC,aes(False_positive,True_positive,label = Cut_values,color = Time_point))+
  geom_roc(labels = F,stat = 'identity', n.cuts = 0)+
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_few()
ggsave(p,filename = "Survival.Signiture.ROC.pdf",width = 3.5,height = 2)

###### Nomogram + calibration + C-index #########
#options(unzip ='internal')
#devtools::install_github('yikeshu0611/ggDCA')
#devtools::install_github('yikeshu0611/do')
#library(ggDCA) > 4.1 R version
save(multi_variate_cox_2,Usedata,file = "DCA.Rdata")
#dca1 <- dca(multi_variate_cox_2,new.data = NULL,times=c(12,36,60))
library(survival)
library(rms)
library(caret)
library(PredictABEL)
library(survivalROC)
library(RColorBrewer)
formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(rownames(ph_hypo_table)[ph_hypo_table[,3] > 0.05],collapse = "+")))

res.cox1 <- cph(formula,data=Usedata,surv=T,x=TRUE, y=TRUE,time.inc=365) #0.69
cal1 <- calibrate(res.cox1, cmethod='KM', method="boot",u=365,m=60,B=1000)

res.cox3 <- cph(formula,data=Usedata,surv=T,x=TRUE, y=TRUE,time.inc=1095) #0.71
cal3 <- calibrate(res.cox3, cmethod='KM', method="boot",u=1095,m=60,B=1000)

res.cox5 <- cph(formula,data=Usedata,surv=T,x=TRUE, y=TRUE,time.inc=1825) #0.71
cal5 <- calibrate(res.cox5, cmethod='KM', method="boot",u=1825,m=60,B=1000)

#### TCGA - calibration ####
pdf("Survival.Signiture.calibration.pdf",height = 5,width = 5)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", 
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[1], pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[2], pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[3], pch = 16)

abline(0,1, lwd = 2, lty = 3, col = brewer.pal(8,"Set1")[4])
legend("topleft", #Âõæ‰æãÁöÑ‰Ωç???
       legend = c("1-year","3-year","5-year"), #Âõæ‰æãÊñáÂ≠ó
       col = brewer.pal(8,"Set1")[1:3], #Âõæ‰æãÁ∫øÁöÑÈ¢úËâ≤Ôºå‰∏éÊñáÂ≠óÂØπÂ∫î
       lwd = 2,#Âõæ‰æã‰∏≠Á∫øÁöÑÁ≤ó???
       cex = 1.2,#Âõæ‰æãÂ≠ó‰ΩìÂ§ßÂ∞è
       bty = "n")#‰∏çÊòæÁ§∫Âõæ‰æãËæπ???
dev.off()

#### TCGA - external validation (2 folds) ####
formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(rownames(ph_hypo_table)[ph_hypo_table[,3] > 0.05],collapse = "+")))

set.seed(12345)
folds <- createFolds(Usedata$OS,k=2)
for (i in 1:2) {
  train_data <- Usedata[-folds[[i]],]
  test_data <- Usedata[folds[[i]],]
  f <- cph(formula,data=train_data,surv=T,x=TRUE, y=TRUE,time.inc=365)
  fev1 <- cph(Surv(OS.time, OS) ~ predict(f, newdata=test_data), x=T, y=T, surv=T, data=test_data, time.inc=365)
  cal1 <- calibrate(fev1, cmethod="KM", method="boot", u=365, m=30, B=1000)
  
  fev3 <- cph(Surv(OS.time, OS) ~ predict(f, newdata=test_data), x=T, y=T, surv=T, data=test_data, time.inc=365*3)
  cal3 <- calibrate(fev3, cmethod="KM", method="boot", u=365*3, m=30, B=1000)
  
  fev5 <- cph(Surv(OS.time, OS) ~ predict(f, newdata=test_data), x=T, y=T, surv=T, data=test_data, time.inc=365*5)
  cal5 <- calibrate(fev5, cmethod="KM", method="boot", u=365*5, m=30, B=1000)
  
  
  pdf(paste("Survival.Signiture.calibration.TCGA.Folds.",i,".pdf",sep = ""),height = 5,width = 5)
  plot(cal1,lwd = 2,lty = 0,errbar.col = brewer.pal(8,"Set1")[1],
       bty = "l", 
       xlim = c(0,1),ylim= c(0,1),
       xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
       col = brewer.pal(8,"Set1")[1],
       cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
  lines(cal1[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[1], pch = 16)
  mtext("")
  
  plot(cal3,lwd = 2,lty = 0,errbar.col = brewer.pal(8,"Set1")[2],
       xlim = c(0,1),ylim= c(0,1),col = brewer.pal(8,"Set1")[2],add = T)
  lines(cal3[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[2], pch = 16)
  
  plot(cal5,lwd = 2,lty = 0,errbar.col = brewer.pal(8,"Set1")[3],
       xlim = c(0,1),ylim= c(0,1),col = brewer.pal(8,"Set1")[3],add = T)
  lines(cal5[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[3], pch = 16)
  
  abline(0,1, lwd = 2, lty = 3, col = brewer.pal(8,"Set1")[4])
  legend("topleft", #Âõæ‰æãÁöÑ‰Ωç???
         legend = c("1-year","3-year","5-year"), #Âõæ‰æãÊñáÂ≠ó
         col = brewer.pal(8,"Set1")[1:3], #Âõæ‰æãÁ∫øÁöÑÈ¢úËâ≤Ôºå‰∏éÊñáÂ≠óÂØπÂ∫î
         lwd = 2,#Âõæ‰æã‰∏≠Á∫øÁöÑÁ≤ó???
         cex = 1.2,#Âõæ‰æãÂ≠ó‰ΩìÂ§ßÂ∞è
         bty = "n")#‰∏çÊòæÁ§∫Âõæ‰æãËæπ???
  dev.off()
  
}


#### TCGA - Intra - nomogram ####
dd <- datadist(Usedata)
options(datadist = "dd")
coxpbc<-cph(formula = formula,data=Usedata,x=T,y=T,surv = T,na.action=na.delete)
surv<-Survival(coxpbc) 
med <- Quantile(f2)
surv2<-function(x) surv(365,x)
surv3<-function(x) surv(1095,x)
surv4<-function(x) surv(1825,x)
#surv5 <- function(x) median(lp=x)
x<-nomogram(coxpbc,fun = list(surv2,surv3,surv4,surv5),lp=T,
            funlabel = c('1-year survival probability','3-year survival probability','5-year survival probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("Survival.Signiture.nomogram.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

#### TCGA - dependent time ROC -> external validation (2 folds) ####
#risk_score_table_multi_cox2
dir.create("TimeROC")
TimeAUC <- data.frame()
TimeROC <- data.frame()
for (j in 1:50) {
  folds <- createFolds(Usedata$OS,k=2)
  for (i in 1:2) {
    train_data <- Usedata[-folds[[i]],]
    test_data <- Usedata[folds[[i]],]
    
    res.cox1 <- coxph(formula,data=train_data)
    #validate(res.cox1, method="boot", B=1000, dxy=T)
    fvad <- coxph(Surv(OS.time, OS) ~ predict(res.cox1, newdata=test_data), data=test_data)
    risk_score <- predict(fvad, type="risk",newdata=test_data)
    test_data$risk_score <- risk_score
    #for_multi_ROC <- multi_ROC(time_vector = c(365*c(1,3,5)), risk_score_table = risk_score_table_multi_cox2)
    roc1=survivalROC(Stime=test_data$OS.time, status=test_data$OS, marker = test_data$risk_score, predict.time =365*1, method="KM") #0.73
    TimeAUC <- rbind.data.frame(TimeAUC,c(roc1$AUC,1,j,i))
    TimeROC <- rbind.data.frame(TimeROC,data.frame(TP=roc1$TP,FP=roc1$FP,Cut_values = roc1$cut.values) %>% mutate(Times=j,Folds=i,Year=1))
    roc3=survivalROC(Stime=test_data$OS.time, status=test_data$OS, marker = test_data$risk_score, predict.time =365*3, method="KM") #0.72
    TimeAUC <- rbind.data.frame(TimeAUC,c(roc3$AUC,3,j,i))
    TimeROC <- rbind.data.frame(TimeROC,data.frame(TP=roc3$TP,FP=roc3$FP,Cut_values = roc3$cut.values) %>% mutate(Times=j,Folds=i,Year=3))
    roc5=survivalROC(Stime=test_data$OS.time, status=test_data$OS, marker = test_data$risk_score, predict.time =365*5, method="KM") #
    TimeAUC <- rbind.data.frame(TimeAUC,c(roc5$AUC,5,j,i))
    TimeROC <- rbind.data.frame(TimeROC,data.frame(TP=roc5$TP,FP=roc5$FP,Cut_values = roc5$cut.values) %>% mutate(Times=j,Folds=i,Year=5))
  }
}
#ggplot(longtest, aes(d = D, m = M, color = name)) + geom_roc() + style_roc()
TimeAUC <- TimeAUC %>% dplyr::rename(AUC=1,Year=2,Times=3,Folds=4)
write.csv(TimeAUC,file = "TimeROC/TimeAUC.TCGA.Folds.csv",row.names = F)
write.csv(TimeROC,file = "TimeROC/TimeROC.TCGA.Folds.csv",row.names = F)

## figure time AUC mean sd ## 
#Year1 0.7370186 0.03908213
#Year3 0.6908085 0.04146177
#Year5 0.6662195 0.05107059
library(ggthemes)
TimeAUC$Year <- factor(TimeAUC$Year,levels = c(1,3,5))
p<-ggplot(TimeAUC,aes(x=Year,y=AUC,fill=Year))+
  geom_boxplot()+
  geom_jitter(height = 0.01,width = 0.1)+
  theme_few()+
  theme(legend.position = "none")
ggsave(p,filename = "TimeROC/TimeAUC.TCGA.Folds.pdf",height = 2.5,width = 2.5)

Target.data <- TimeROC %>% dplyr::filter(Times==2 & Folds==1)
Target.data$Year <- factor(Target.data$Year,levels = c(1,3,5))
library(plotROC)
p <-ggplot(Target.data,aes(FP,TP,label = Cut_values,color = Year))+
  geom_roc(labels = F,stat = 'identity', n.cuts = 0)+
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_few()+
  labs(x="False positive",y="True positive")

ggsave(p,filename = "TimeROC/Representative.ROC.TCGA.Folds.pdf",width = 3.5,height = 2.8)

###### Validation -> ICGC + PCAWG Validation ######
#### LICA-FR -> Samples -> No OK ####
#exp_seq.*.tsc.gz ÔºöÁ±ª‰ººÁ≥ªÊï∞Áü©ÈòµÂΩ¢ÂºèÔºåÁªô‰∫ÜÂéüÂßãcount/ÂΩí‰∏ÄÂåñcount
#sample.*.tsv ??? ÂåÖÂê´‰∫ÜÊ†∑???/ÊçêËµ†‰∫∫‰ø°???
#specimen.PRAD-CA.tsvÔºöÊ†áËÆ∞Âá∫Ê†∑Êú¨ÊòØÂê¶‰∏∫ÁôåÁóáËøòÊòØÁôå???.
#donor : survival data

### build the cox model based TCGA count dataset => 361 samples
LIHC.TCGA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",row.names = 1,check.names = F)
LIHC.TCGA.Exp.TPM = log2(LIHC.TCGA.Exp+1)

#Inter.Sample.361
#five genes
Signature.Five.Gene <- Intersect.V[-c(1,7),]
write.csv(Signature.Five.Gene,file = "Signature.Five.Gene.csv",row.names = F)

load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
load("Inter.miRNA.mRNA.Samples.361.Rdata")
LIHC.TCGA.Exp.Target <- LIHC.TCGA.Exp.TPM[Signature.Five.Gene$Select.Variables,Inter.Sample.361] %>% data.frame(check.names = F) %>%
  t() %>% data.frame() %>% rownames_to_column("ID") %>% mutate(Sample.ID=substr(ID,1,15)) %>%
  mutate(Sample.ID = factor(Sample.ID,levels = rownames(risk_score_table_multi_cox2))) %>% 
  dplyr::arrange(Sample.ID) %>%
  cbind.data.frame(risk_score_table_multi_cox2[,1:2])

#
library(rms)
library(tidyverse)

### process data 
setwd("K:/TCGA/Anlysis/LIHC")
LICA.Exp <- read.csv("../../ICGC/LICA-FR/LICA-FR.csv")
LICA.Exp$gene_id <- LICA.Exp$gene_id %>% str_remove_all("\\..*")
LICA.Exp <- LICA.Exp %>% remove_rownames() %>% column_to_rownames("gene_id")

LICA.specimen <- read.csv("../../ICGC/LICA-FR/specimen.LICA-FR.tsv",header =T,row.names = 1,sep = "\t") %>% dplyr::filter(specimen_type == "Primary tumour - solid tissue")
#LICA.sample <- read.csv("../../ICGC/LICA-FR/sample.LICA-FR.tsv",header =T,row.names = 1,sep = "\t")
LICA.donor <- read.csv("../../ICGC/LICA-FR/donor.LICA-FR.tsv",header =T,row.names = 1,sep = "\t")
LICA.donor2 <- LICA.donor %>% dplyr::select(donor_vital_status,donor_survival_time) %>% na.omit() %>%
  mutate(OS=ifelse(donor_vital_status == "alive",0,1)) %>% mutate(OS.time = donor_survival_time)
  
#Êï¥ÁêÜÂêéË°®ËææÁü©ÈòµÔºöË°åÂêç‰∏∫Âü∫Âõ†ÂêçÔºõÂàóÂêçÂåÖÂê´‰∫ÜÊ†∑Êú¨-ÊçêËµ†???-ÊòØÂê¶‰∏∫ÁôåÊóÅÔºàCancer-C ;Normal-N???
## ÂèëÁé∞ÁîüÂ≠òÊï∞ÊçÆ‰∏éË°®ËææÊï∞ÊçÆÂá†‰πéÂØπÂ∫î‰∏ç???
Index1 <- (colnames(LICA.Exp) %>% str_remove_all("\\..*")) %in% rownames(LICA.specimen)
Col.mid <- unlist(colnames(LICA.Exp) %>% str_split("\\."))
Index2 <- Col.mid[Col.mid %>% str_detect("DO")] %in% rownames(LICA.donor2)

### PCAWG - LINC.JP C-index 0.832 => 1 YEAR + 3 YEAR
pcawg.specimen <- read.csv("../../PCAWG/pcawg.sp_specimen_type.txt",header = T,sep = "\t") %>% filter(dcc_specimen_type == "Primary tumour - solid tissue" | 
                                                                                                  dcc_specimen_type == "Primary tumour - other" | 
                                                                                                  dcc_specimen_type =="Primary tumour" ) #SP
pcawg.project <- read.csv("../../PCAWG/pcawg.project_code_sp.txt",header = T,sep = "\t") %>% dplyr::filter(dcc_project_code == "LICA-FR") %>%
  dplyr::filter(icgc_specimen_id %in% pcawg.specimen$icgc_specimen_id)

pcawg.survival <- read.csv("../../PCAWG/pcawg.survival_sp.txt",header = T,sep = "\t") %>% dplyr::select(xena_sample,donor_vital_status,donor_survival_time) %>%
  na.omit() %>% mutate(OS=ifelse(donor_vital_status == "alive",0,1)) %>% mutate(OS.time = donor_survival_time)

pcawg.RNA <- read.csv("../../PCAWG/pcawg.tophat_star_fpkm_uq.v2_aliquot_gl.sp.log.txt",header = T,row.names = 1,sep = "\t")
#log2(fpkm-uq+0.001) => log2(TPM+1)
pcawg.RNA.2 <- 2^pcawg.RNA - 0.001

FPKMUQtoTPM <- function(x) {
  return(exp(log(x) - log(sum(x)) + log(1e6)))
}
pcawg.RNA.TPM <- apply(pcawg.RNA.2, 2,FPKMUQtoTPM)
write.csv(pcawg.RNA.TPM,file = "PCAWG/pcawg.RNA.TPM.csv")
pcawg.RNA.TPM <- read.csv("../../PCAWG/pcawg.RNA.TPM.csv",row.names = 1)
pcawg.RNA.TPM.log <- log2(pcawg.RNA.TPM+1)


rownames(pcawg.RNA) <- rownames(pcawg.RNA) %>% str_remove_all("\\..*")

#### LIHC-US -> 54  LICA-FR -> 6  LINC-JP -> 31 => 0  LIRI-JP => 269 => 67 ####
#### LIRI-JP -> Only 67  ####
LIRI.JP <- read.csv("../../PCAWG/pcawg.project_code_sp.txt",header = T,sep = "\t") %>% dplyr::filter(dcc_project_code == "LIRI-JP") %>%
  dplyr::filter(icgc_specimen_id %in% pcawg.specimen$icgc_specimen_id) %>%
  dplyr::filter(icgc_specimen_id %in% pcawg.survival$xena_sample)

Inter.Samples <- intersect(LIRI.JP$icgc_specimen_id, colnames(pcawg.RNA))

LIRI.JP.Target <- pcawg.RNA.TPM.log[Signature.Five.Gene$Select.Variables,Inter.Samples] %>% t() %>% data.frame(check.names = F) %>%
  rownames_to_column("ID") %>% merge.data.frame(pcawg.survival,by.x="ID",by.y="xena_sample")

### 
formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(Signature.Five.Gene$Select.Variables,collapse = "+")))

f <- cph(formula,data=LIHC.TCGA.Exp.Target,surv=T,x=TRUE, y=TRUE,time.inc=365)
fev1 <- cph(Surv(OS.time, OS) ~ predict(f, newdata=LIRI.JP.Target), x=T, y=T, surv=T, data=LIRI.JP.Target, time.inc=365) #0.664/2+0.5 => 0.832
cal1 <- calibrate(fev1, cmethod="KM", method="boot", u=365, m=20, B=1000)

fev3 <- cph(Surv(OS.time, OS) ~ predict(f, newdata=LIRI.JP.Target), x=T, y=T, surv=T, data=LIRI.JP.Target, time.inc=365*3) ##0.832
cal3 <- calibrate(fev3, cmethod="KM", method="boot", u=365*3, m=20, B=1000)

#fev5 <- cph(Surv(OS.time, OS) ~ predict(f, newdata=LIRI.JP.Target), x=T, y=T, surv=T, data=LIRI.JP.Target, time.inc=365*5)
#cal5 <- calibrate(fev5, cmethod="KM", method="boot", u=365*5, m=10, B=1000)

pdf(paste("Survival.Signiture.calibration.TCGA.To.LIRI-JP.pdf",sep = ""),height = 5,width = 5)
plot(cal1,lwd = 2,lty = 0,errbar.col = brewer.pal(8,"Set1")[1],
     bty = "l", 
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = brewer.pal(8,"Set1")[1],
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[1], pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = brewer.pal(8,"Set1")[2],
     xlim = c(0,1),ylim= c(0,1),col = brewer.pal(8,"Set1")[2],add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[2], pch = 16)

if (F) {
  plot(cal5,lwd = 2,lty = 0,errbar.col = brewer.pal(8,"Set1")[3],
       xlim = c(0,1),ylim= c(0,1),col = brewer.pal(8,"Set1")[3],add = T)
  lines(cal5[,c('mean.predicted',"KM")],
        type = 'b', lwd = 1, col = brewer.pal(8,"Set1")[3], pch = 16)
}

abline(0,1, lwd = 2, lty = 3, col = brewer.pal(8,"Set1")[4])
legend("topleft", #Âõæ‰æãÁöÑ‰Ωç???
       legend = c("1-year","3-year"), #Âõæ‰æãÊñáÂ≠ó
       col = brewer.pal(8,"Set1")[1:2], #Âõæ‰æãÁ∫øÁöÑÈ¢úËâ≤Ôºå‰∏éÊñáÂ≠óÂØπÂ∫î
       lwd = 2,#Âõæ‰æã‰∏≠Á∫øÁöÑÁ≤ó???
       cex = 1.2,#Âõæ‰æãÂ≠ó‰ΩìÂ§ßÂ∞è
       bty = "n")#‰∏çÊòæÁ§∫Âõæ‰æãËæπ???
dev.off()

### timeROC
res.cox1 <- coxph(formula,data=LIHC.TCGA.Exp.Target)
#validate(res.cox1, method="boot", B=1000, dxy=T)
fvad <- coxph(Surv(OS.time, OS) ~ predict(res.cox1, newdata=LIRI.JP.Target), data=LIRI.JP.Target) #0.832
risk_score <- predict(fvad, type="risk",newdata=LIRI.JP.Target)
LIRI.JP.Target$risk_score <- risk_score

TimeROC <- data.frame()
#for_multi_ROC <- multi_ROC(time_vector = c(365*c(1,3,5)), risk_score_table = risk_score_table_multi_cox2)
roc1=survivalROC(Stime=LIRI.JP.Target$OS.time, status=LIRI.JP.Target$OS, marker = LIRI.JP.Target$risk_score, predict.time =365*1, method="KM") #0.9035492 + survival 0.8938336
#TimeAUC <- rbind.data.frame(TimeAUC,c(roc1$AUC,1,j,i))
TimeROC <- rbind.data.frame(TimeROC,data.frame(TP=roc1$TP,FP=roc1$FP,Cut_values = roc1$cut.values) %>% mutate(Year=1))
roc3=survivalROC(Stime=LIRI.JP.Target$OS.time, status=LIRI.JP.Target$OS, marker = LIRI.JP.Target$risk_score, predict.time =365*3, method="KM") #0.8757653 + Survival 0.7885328
#TimeAUC <- rbind.data.frame(TimeAUC,c(roc3$AUC,3,j,i))
TimeROC <- rbind.data.frame(TimeROC,data.frame(TP=roc3$TP,FP=roc3$FP,Cut_values = roc3$cut.values) %>% mutate(Year=3))
roc5=survivalROC(Stime=LIRI.JP.Target$OS.time, status=LIRI.JP.Target$OS, marker = LIRI.JP.Target$risk_score, predict.time =365*5, method="KM") #0.7853117 + Survival 0.5456575
#TimeAUC <- rbind.data.frame(TimeAUC,c(roc5$AUC,5,j,i))
TimeROC <- rbind.data.frame(TimeROC,data.frame(TP=roc5$TP,FP=roc5$FP,Cut_values = roc5$cut.values) %>% mutate(Year=5))
TimeROC$Year <- factor(TimeROC$Year,levels = c(1,3,5))
write.csv(TimeROC,file = "SourceData.F2G.csv")

p <-ggplot(TimeROC,aes(FP,TP,label = Cut_values,color = Year))+
  geom_roc(labels = F,stat = 'identity', n.cuts = 0)+
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_few()+
  labs(x="False positive",y="True positive")

ggsave(p,filename = "Survival.Signature.TimeROC.TCGA.To.LIRI-JP.pdf",width = 3.5,height = 2.8)

#### LIRI-JP -> HCCCDB  ####
library(tidyverse)
Exp.Samples <- read.csv("../../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18.sample.txt",sep = "\t",check.names = F,row.names = 1) %>%
  t() %>% data.frame(check.names = F) %>% 
  filter(TYPE=="HCC") %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID,TYPE,PATIENT_ID)
Exp.Clinical <- read.csv("../../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18.patient.txt",sep = "\t",check.names = F,row.names = 1) %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("PATIENT_ID") %>%
  mutate(OS=if_else(STATUS=="Dead",1,0)) %>%
  mutate(OS.time=SUR) %>%
  dplyr::select(PATIENT_ID,PATIENT,OS,OS.time) %>%
  merge(Exp.Samples,by="PATIENT_ID")

Exp <- read.csv("../../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18_mRNA_level3.txt",sep = "\t",check.names = F)
Exp.Data <- Exp[,-1] %>% data.frame(check.names = F) %>%
  remove_rownames() %>%
  column_to_rownames("Symbol") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("SampleID") %>%
  filter(SampleID %in% Exp.Clinical$SampleID) %>%
  column_to_rownames("SampleID") %>%
  t()

## CBX2 CEP55
middata <- Exp.Data[c("CBX2","CEP55"),] %>% t() %>%
  data.frame() %>% 
  rownames_to_column("SampleID") %>%
  merge(Exp.Clinical,by="SampleID")
library(survival)
library(survminer)
middata$OS <- as.numeric(as.character(middata$OS))
middata$OS.time <- as.numeric(as.character(middata$OS.time))
middata$Group <- if_else(middata$CEP55 > median(middata$CEP55),"High","Low")
sfit1=survfit(Surv(OS.time,OS)~Group, data=middata)
p<-ggsurvplot(sfit1,pval =TRUE, data = middata, #risk.table = TRUE,
              #surv.median.line = "hv", #fun = "event"‘Úø…ªÊ÷∆cumulative eventsÕº£®fun≥˝¡ÀevenÕ‚£¨ªπ”–log£∫log transformation∫Õcumhaz£∫cumulative hazard???
              legend.title = "CEP55",
              legend.labs = c("High", "Low"),
              conf.int.style = "step",
              xlab = "Time in months",
              #risk.table = "abs_pct",
              risk.table.y.text.col = F,
              risk.table.y.text = FALSE,
              conf.int = TRUE,
              palette = "Set1",
              ggtheme = ggthemes::theme_few())
print(p)
#dev.off()
#graph2pdf(file=paste("miRNA.OS.KM.Cox/miRNA.OS.KM.",colnames(RPMdata)[i],".pdf",sep = ''),height = 6,width = 5)


###### Risk Score with TMN => 361 ######
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table

meta <- read.csv("clinical.cart.2021-08-19/clinical.tsv",sep = "\t",header = T)
#meta <- column_to_rownames(meta,var = "case_submitter_id")
meta=meta[,colnames(meta) %in% c("case_submitter_id",
                                 "race",
                                 "gender",
                                 "ajcc_pathologic_t",
                                 "ajcc_pathologic_n",
                                 "ajcc_pathologic_m",
                                 "age","Stage")] %>% unique()

library(tidyverse)
risk.score.meta <- risk_score_table_multi_cox2[substr(rownames(risk_score_table_multi_cox2),1,12) %in% meta$case_submitter_id,] %>% 
  rownames_to_column("sampleid") %>% mutate(case_submitter_id = substr(sampleid,1,12)) %>%
  merge.data.frame(.,meta,by="case_submitter_id")

save(risk.score.meta,file = "risk.score.meta.Rdata") #risk.score.meta

risk.score.meta.mid <- risk.score.meta %>% filter(Stage != "")

p<-ggboxplot(risk.score.meta.mid, "Stage", "total_risk_score", fill = "Stage",
             palette = "Set1",
             add = "jitter")+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Risk score")+
  stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                     ref.group = ".all")
ggsave(p,filename = "Survival.Signature.RiskScore.Stage.pdf",height = 2.5,width = 2.5)


risk.score.meta.mid <- risk.score.meta %>% filter(T != "")
risk.score.meta.mid$T <- factor(risk.score.meta.mid$T,levels = paste("T",c("X",1,2,3,4),sep = ""))

p<-ggboxplot(risk.score.meta.mid, "T", "total_risk_score", fill = "T",
             palette = "Set1",
             add = "jitter")+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Risk score")+
  stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                     ref.group = ".all")
ggsave(p,filename = "Survival.Signature.RiskScore.ajcc_pathologic_t.pdf",height = 2.5,width = 2.5)


risk.score.meta.mid <- risk.score.meta %>% filter(ajcc_pathologic_m != "")

p<-ggboxplot(risk.score.meta.mid, "ajcc_pathologic_m", "total_risk_score", fill = "ajcc_pathologic_m",
             palette = "Set1",
             add = "jitter")+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Risk score")+
  stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                     ref.group = ".all")
ggsave(p,filename = "Survival.Signature.RiskScore.ajcc_pathologic_m.pdf",height = 2.5,width = 2.5)

risk.score.meta.mid <- risk.score.meta %>% filter(ajcc_pathologic_n != "")

p<-ggboxplot(risk.score.meta.mid, "ajcc_pathologic_n", "total_risk_score", fill = "ajcc_pathologic_n",
             palette = "Set1",
             add = "jitter")+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="",y="Risk score")+
  stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                     ref.group = ".all")
ggsave(p,filename = "Survival.Signature.RiskScore.ajcc_pathologic_n.pdf",height = 2.5,width = 2.5)






###### Risk Score dependent or independent signature #######
setwd("K:/TCGA/Anlysis/LIHC")
#load("risk.score.meta.Rdata")
#write.csv(risk.score.meta,file = "risk.score.meta.csv",row.names = F)
risk.score.meta <- read.csv("risk.score.meta.csv")
risk.score.meta$Stage <- risk.score.meta$Stage %>% str_remove_all(".* ")
risk.score.meta$Stage <- factor(risk.score.meta$Stage,levels = c("I","II","III","IV"))
risk.score.meta$gender <- factor(risk.score.meta$gender,levels = c("male","female"))
colnames(risk.score.meta)[10:12] = c("M","N","T")
risk.score.meta$T <- factor(risk.score.meta$T,levels =paste("T",c(1,2,3,4,"X"),sep = ""))
risk.score.meta$M <- factor(risk.score.meta$M,levels = paste("M", c("0","1","X"),sep = ""))
risk.score.meta$N <- factor(risk.score.meta$N,levels = paste("N", c("0","1","X"),sep = ""))

formula <- as.formula(paste0("Surv(OS.time,OS)~",paste(c("age","gender","T","M","N","total_risk_score"),collapse = "+")))
multi_variate_cox <- coxph(formula,data=risk.score.meta)

cox_sum <- summary(multi_variate_cox)
HR<- round(cox_sum$coefficients[,2],5) 
PValue<- round(cox_sum$coefficients[,5],15) 
lower<-round(cox_sum$conf.int[,3],5)
upper<-round(cox_sum$conf.int[,4],5)
#4.∂‡“ÚÀÿΩ·π˚”≈ªØ≤¢≥…±Ì£∫mul_cox1
HRa<-paste(HR," (", lower,'-',upper,")",sep = "")
mul_results<- data.frame("HRa"=HRa,"P"=PValue,"HR"=HR,"lower"=lower,"upper"=upper)
Variants=data.frame(rownames(mul_results))
rownames(mul_results)=NULL
mul_results<-cbind(Variants,mul_results)
colnames(mul_results)=c("Variants","Hazard Ratio (95%CI)","P-value","","","")
write.csv(mul_results,file = "SourceData.F2F.csv")
## check PH hypothesis
#ph_hypo_multi <- cox.zph(multi_variate_cox)
#ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
ggforest(multi_variate_cox,  #coxphÂæóÂà∞ÁöÑCoxÂõûÂΩíÁªìÊûú
         data = risk.score.meta,  #Êï∞ÊçÆ???
         main = 'Hazard ratio',  #Ê†áÈ¢ò
         #cpositions = c(0.05, 0.35, 0.35),  #Ââç‰∏âÂàóË∑ù???
         fontsize = 1, #Â≠ó‰ΩìÂ§ßÂ∞è
         refLabel = 'reference', #Áõ∏ÂØπÂèòÈáèÁöÑÊï∞ÂÄºÊ†áÁ≠æÔºå‰πüÂèØÊîπ‰∏∫1
)

library(export)
graph2pdf(file="Survival.Signature.Forest.Independent.pdf",height=6,width=8)

###### Expression level with TMN #####
dir.create("Stage.Expression")
meta <- read.csv("clinical.cart.2021-08-19/clinical.tsv",sep = "\t",header = T)
#meta <- column_to_rownames(meta,var = "case_submitter_id")
meta=meta[,colnames(meta) %in% c("case_submitter_id",
                                 "race",
                                 "gender",
                                 "ajcc_pathologic_t",
                                 "ajcc_pathologic_n",
                                 "ajcc_pathologic_m",
                                 "age","Stage")] %>% unique()


mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor")

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

Inter.PrimaryTumor.Samples <- intersect(colnames(mRNA.Exp.log2.TPM),mRNA.Phenodata$Sample.ID) #371

Cor.ceRNA1 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")
Cor.ceRNA2 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")

#### lncRNA 371 Samples ####
Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv")

mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[unique(c(Cor.ceRNA1$lncRNA.ENSEMBLE.ID,Cor.ceRNA2$lncRNA.ENSEMBLE.ID)),Inter.PrimaryTumor.Samples] %>%
  t() %>% data.frame(check.names = F) %>% mutate(case_submitter_id = substr(rownames(.),1,12)) %>%
  merge(.,meta,by="case_submitter_id")

write.csv(mRNA.Exp.log2.TPM.Target,"TMN.lncRNA.Exp.log2.TPM.Target.csv")

library(ggpubr)
library(ggthemes)
for (lncRNA in unique(c(Cor.ceRNA1$lncRNA.ENSEMBLE.ID,Cor.ceRNA2$lncRNA.ENSEMBLE.ID))) {
  ## stage
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(Stage != "")
  mRNA.Exp.log2.TPM.Target.mid$Stage <- factor(mRNA.Exp.log2.TPM.Target.mid$Stage,levels = paste("Stage", c("I","II","III","IV")))
  
  lncRNA.name <- (Inter.lncRNA.GeneCard %>% filter(lncRNA.ENSEMBLE.ID==lncRNA))$lncRNA.ID
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "Stage", lncRNA, fill = "Stage",
           palette = "Set1",
           add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                       ref.group = ".all") #comparisons = combn(levels(mRNA.Exp.log2.TPM.Target.mid$Stage), 2, simplify =FALSE))
  
                       
  ggsave(p,filename = paste("Stage.Expression/lncRNA.",lncRNA,".",lncRNA.name,".Stage.pdf",sep = ""),width = 3,height = 3)
  
  ## t
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(ajcc_pathologic_t != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t,levels = paste("T",c(1,2,"2a","2b",3,"3a","3b",4,"X"),sep = ""))
  
  lncRNA.name <- (Inter.lncRNA.GeneCard %>% filter(lncRNA.ENSEMBLE.ID==lncRNA))$lncRNA.ID
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_t", lncRNA, fill = "ajcc_pathologic_t",
              palette = "Set1",
              add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
                       #comparisons = combn(levels(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t), 2, simplify =FALSE))
  
  ggsave(p,filename = paste("Stage.Expression/lncRNA.",lncRNA,".",lncRNA.name,".ajcc_pathologic_t.pdf",sep = ""),width = 3,height = 3)
  
  ## m
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(ajcc_pathologic_m != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_m <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_m,levels = paste("M", c("X","0","1"),sep = ""))
  
  lncRNA.name <- (Inter.lncRNA.GeneCard %>% filter(lncRNA.ENSEMBLE.ID==lncRNA))$lncRNA.ID
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_m", lncRNA, fill = "ajcc_pathologic_m",
              palette = "Set1",
              add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  
  ggsave(p,filename = paste("Stage.Expression/lncRNA.",lncRNA,".",lncRNA.name,".ajcc_pathologic_m.pdf",sep = ""),width = 3,height = 3)
  
  ## n
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(ajcc_pathologic_n != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_n <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_n,levels = paste("N", c("X","0","1"),sep = ""))
  
  lncRNA.name <- (Inter.lncRNA.GeneCard %>% filter(lncRNA.ENSEMBLE.ID==lncRNA))$lncRNA.ID
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_n", lncRNA, fill = "ajcc_pathologic_n",
              palette = "Set1",
              add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  
  ggsave(p,filename = paste("Stage.Expression/lncRNA.",lncRNA,".",lncRNA.name,".ajcc_pathologic_n.pdf",sep = ""),width = 3,height = 3)
}


#### mRNA 371 Samples ####
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[unique(c(Cor.ceRNA1$mRNA.ENSEMBLE.ID,Cor.ceRNA2$mRNA.ENSEMBLE.ID)),Inter.PrimaryTumor.Samples] %>%
  t() %>% data.frame(check.names = F) %>% mutate(case_submitter_id = substr(rownames(.),1,12)) %>%
  merge(.,meta,by="case_submitter_id") 

write.csv(mRNA.Exp.log2.TPM.Target,"TMN.mRNA.Exp.log2.TPM.Target.csv")

library(ggpubr)
library(ggthemes)
for (lncRNA in unique(c(Cor.ceRNA1$mRNA.ENSEMBLE.ID,Cor.ceRNA2$mRNA.ENSEMBLE.ID))) {
  ## stage
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(Stage != "")
  mRNA.Exp.log2.TPM.Target.mid$Stage <- factor(mRNA.Exp.log2.TPM.Target.mid$Stage,levels = paste("Stage", c("I","II","III","IV")))
  
  lncRNA.name <- GeneType[GeneType$Gene.stable.ID == lncRNA,]$Gene.name %>% unique()
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "Stage", lncRNA, fill = "Stage",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                       ref.group = ".all") #comparisons = combn(levels(mRNA.Exp.log2.TPM.Target.mid$Stage), 2, simplify =FALSE))
  
  
  ggsave(p,filename = paste("Stage.Expression/mRNA.",lncRNA,".",lncRNA.name,".Stage.pdf",sep = ""),width = 3,height = 3)
  
  ## t
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(ajcc_pathologic_t != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t,levels = paste("T",c(1,2,"2a","2b",3,"3a","3b",4,"X"),sep = ""))
  
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_t", lncRNA, fill = "ajcc_pathologic_t",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  #comparisons = combn(levels(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t), 2, simplify =FALSE))
  
  ggsave(p,filename = paste("Stage.Expression/mRNA.",lncRNA,".",lncRNA.name,".ajcc_pathologic_t.pdf",sep = ""),width = 3,height = 3)
  
  ## m
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(ajcc_pathologic_m != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_m <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_m,levels = paste("M", c("X","0","1"),sep = ""))

  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_m", lncRNA, fill = "ajcc_pathologic_m",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  
  ggsave(p,filename = paste("Stage.Expression/mRNA.",lncRNA,".",lncRNA.name,".ajcc_pathologic_m.pdf",sep = ""),width = 3,height = 3)
  
  ## n
  mRNA.Exp.log2.TPM.Target.mid <- mRNA.Exp.log2.TPM.Target %>% filter(ajcc_pathologic_n != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_n <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_n,levels = paste("N", c("X","0","1"),sep = ""))
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_n", lncRNA, fill = "ajcc_pathologic_n",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA.name,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  
  ggsave(p,filename = paste("Stage.Expression/mRNA.",lncRNA,".",lncRNA.name,".ajcc_pathologic_n.pdf",sep = ""),width = 3,height = 3)
}

#### miRNA 367 Samples ####
miRNA.Exp <- read.csv("../Pancancer/TCGA-LIHC_Mature_miRNA.RPM.csv",check.names = F,row.names = 1)
miRNA.Exp.log2.TPM <- log2(miRNA.Exp+1)
miRNA.Phenodata <- read.csv("../../miRNA.Mature.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor")

Inter.PrimaryTumor.Samples <- intersect(colnames(mRNA.Exp.log2.TPM),miRNA.Phenodata$Sample.ID)

Index <- (rownames(miRNA.Exp.log2.TPM) %>% str_remove_all(".*\\.")) %in% unique(c(Cor.ceRNA1$miRNA.ID,Cor.ceRNA2$miRNA.ID))
miRNA.Exp.log2.TPM.Mid <- miRNA.Exp.log2.TPM[Index,Inter.PrimaryTumor.Samples] %>%
  t() %>% data.frame(check.names = F) %>% mutate(case_submitter_id = substr(rownames(.),1,12)) %>%
  merge(.,meta,by="case_submitter_id")

write.csv(miRNA.Exp.log2.TPM.Mid,file = "TMN.miRNA.Exp.log2.TPM.Mid.csv")


library(ggpubr)
library(ggthemes)
for (lncRNA in unique(c(Cor.ceRNA1$miRNA.ID,Cor.ceRNA2$miRNA.ID))) {
  miRNA.indexname <- str_detect(colnames(miRNA.Exp.log2.TPM.Mid),lncRNA)
  miRNA.name <- colnames(miRNA.Exp.log2.TPM.Mid)[miRNA.indexname]
  
  ## stage
  mRNA.Exp.log2.TPM.Target.mid <- miRNA.Exp.log2.TPM.Mid %>% filter(Stage != "")
  mRNA.Exp.log2.TPM.Target.mid$Stage <- factor(mRNA.Exp.log2.TPM.Target.mid$Stage,levels = paste("Stage", c("I","II","III","IV")))
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "Stage", miRNA.name, fill = "Stage",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,
                       ref.group = ".all") #comparisons = combn(levels(mRNA.Exp.log2.TPM.Target.mid$Stage), 2, simplify =FALSE))
  
  
  ggsave(p,filename = paste("Stage.Expression/miRNA.",lncRNA,".Stage.pdf",sep = ""),width = 3,height = 3)
  
  ## t
  mRNA.Exp.log2.TPM.Target.mid <- miRNA.Exp.log2.TPM.Mid %>% filter(ajcc_pathologic_t != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t,levels = paste("T",c(1,2,"2a","2b",3,"3a","3b",4,"X"),sep = ""))
  
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_t",miRNA.name, fill = "ajcc_pathologic_t",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  #comparisons = combn(levels(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_t), 2, simplify =FALSE))
  
  ggsave(p,filename = paste("Stage.Expression/miRNA.",lncRNA,".ajcc_pathologic_t.pdf",sep = ""),width = 3,height = 3)
  
  ## m
  mRNA.Exp.log2.TPM.Target.mid <- miRNA.Exp.log2.TPM.Mid %>% filter(ajcc_pathologic_m != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_m <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_m,levels = paste("M", c("X","0","1"),sep = ""))
  
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_m", miRNA.name, fill = "ajcc_pathologic_m",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  
  ggsave(p,filename = paste("Stage.Expression/miRNA.",lncRNA,".ajcc_pathologic_m.pdf",sep = ""),width = 3,height = 3)
  
  ## n
  mRNA.Exp.log2.TPM.Target.mid <- miRNA.Exp.log2.TPM.Mid %>% filter(ajcc_pathologic_n != "")
  mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_n <- factor(mRNA.Exp.log2.TPM.Target.mid$ajcc_pathologic_n,levels = paste("N", c("X","0","1"),sep = ""))
  
  p<-ggboxplot(mRNA.Exp.log2.TPM.Target.mid, "ajcc_pathologic_n", miRNA.name, fill = "ajcc_pathologic_n",
               palette = "Set1",
               add = "jitter")+
    theme_few()+
    theme(legend.position = "none")+
    labs(x="",y=paste("log2(TPM+1) ",lncRNA,sep = ""))+
    stat_compare_means(method = "t.test",label = "p.signif",hide.ns = TRUE,ref.group = ".all")
  
  ggsave(p,filename = paste("Stage.Expression/miRNA.",lncRNA,".ajcc_pathologic_n.pdf",sep = ""),width = 3,height = 3)
}



######### mRNA High/Low stratified ########
#### LRRC1  MAP2  MCM2 RMI2  CKAP2L CEP55 CBX2 DUXAP8 #Down => ESR1  CPEB3   ANXA10
### UALCAN
# DNMT1 DNMT3A DNMT3B
# CBX2 => 1.962340E-03
# LRRC1 => 1.223570E-02
ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  scale_y_log10()+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL"),
                                                c("HCC","HL")),method = "t.test")
### DiseaseMeth version 2.0
#CBX2
setwd("C:/Users/51331/Downloads")
data <- read.csv("DiseaseMeth.CBX2_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "CBX2.NM_005189.chr17.77749976.77752476", 
          fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("CBX methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.CBX2.Boxplot.pdf",height = 2.5,width = 2.5)

data <- read.csv("DiseaseMeth.ESR1_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "ESR1.NM_000125.chr6.152126813.152129313", 
             fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("CBX methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.ESR1.Boxplot.pdf",height = 2.5,width = 2.5)

data <- read.csv("DiseaseMeth.CEP55_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "CEP55.NM_018131.chr10.95254388.95256888", 
             fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("CBX methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.CEP55.Boxplot.pdf",height = 2.5,width = 2.5)

data <- read.csv("DiseaseMeth.CKAP2L_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "CKAP2L.NM_001304361.chr2.113521754.113524254", 
             fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("CBX methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.CKAP2L.Boxplot.pdf",height = 2.5,width = 2.5)

data <- read.csv("DiseaseMeth.LRRC1_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "LRRC1.NM_018214.chr6.53657777.53660277", 
             fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("LRRC1 methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.LRRC1.Boxplot.pdf",height = 2.5,width = 2.5)

data <- read.csv("DiseaseMeth.MCM2_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "MCM2.NR_073375.chr3.127315199.127317699", 
             fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("MCM2 methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.MCM2.Boxplot.pdf",height = 2.5,width = 2.5)

###
data <- read.csv("DiseaseMeth.RMI_C2239176_450k_20210901193434-205259.txt",header = T,row.names = 1,sep = "")
data$isControl <- factor(data$isControl,levels = c(0,1),labels = c("HCC","Ctrl"))

p<-ggboxplot(data, "isControl", "RMI2.NM_152308.chr16.11437294.11439794", 
             fill = "isControl",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("RMI2 methylation value",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = "DiseaseMeth.RMI2.Boxplot.pdf",height = 2.5,width = 2.5)


#### methylation ####
#ualcan
#MEXPRESS

#### CBX2 Low/High ####
#DNMT1 ENSG00000130816
#DNMT3A ENSG00000119772
#DNMT3B ENSG00000088305
#CBX2 ENSG00000173894
#CEP55 ENSG00000138180
#RMI2 ENSG00000175643
#MCM2 ENSG00000073111
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
load("Inter.miRNA.mRNA.Samples.361.Rdata")
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000073111",
                                                "ENSG00000130816",
                                                "ENSG00000119772",
                                                "ENSG00000088305",
                                                "ENSG00000173894",
                                                "ENSG00000138180",
                                                "ENSG00000175643",
                                                "ENSG00000107864",
                                                "ENSG00000091831",
                                                "ENSG00000168036"),Inter.Sample.361] %>%
  t() %>% data.frame(check.names = F)

mRNA.Exp.log2.TPM.Target$CBX2Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000173894 > median(mRNA.Exp.log2.TPM.Target$ENSG00000173894),"CBX2High","CBX2Low")
mRNA.Exp.log2.TPM.Target$CBX2Group <- factor(mRNA.Exp.log2.TPM.Target$CBX2Group,levels = c("CBX2Low","CBX2High"))
write.csv(mRNA.Exp.log2.TPM.Target,"SourceData.F5DE.csv")
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CBX2Group", "ENSG00000130816", 
          fill = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

p12 <- ggplot(mRNA.Exp.log2.TPM.Target,aes(CBX2Group,ENSG00000168036,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="log2(TPM+1) CTNNB1",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "t.test",angle=-90,label.x=2)
ggsave(p12,filename = "CBX2.Group.CTNNB1.pdf",height = 1.5,width = 4)

ENSG00000138180

p13 <- ggplot(mRNA.Exp.log2.TPM.Target,aes(CBX2Group,ENSG00000138180,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="log2(TPM+1) CEP55",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "t.test",angle=-90,label.x=2)
ggsave(p13,filename = "CBX2.Group.CEP55.pdf",height = 1.5,width = 4)


ggsave(p,filename = "stratified/CBX2.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CBX2Group", "ENSG00000088305", 
             fill = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CBX2.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CBX2Group", "ENSG00000119772", 
             fill = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CBX2.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)

#### CEP55 Low/High ####
mRNA.Exp.log2.TPM.Target$CEP55Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000138180 > median(mRNA.Exp.log2.TPM.Target$ENSG00000138180),"CEP55High","CEP55Low")
mRNA.Exp.log2.TPM.Target$CEP55Group <- factor(mRNA.Exp.log2.TPM.Target$CEP55Group,levels = c("CEP55High","CEP55Low"))
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CEP55Group", "ENSG00000130816", 
             fill = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CEP55.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CEP55Group", "ENSG00000088305", 
             fill = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CEP55.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CEP55Group", "ENSG00000119772", 
             fill = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CEP55.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)

#### RMI2 Low/High ####
mRNA.Exp.log2.TPM.Target$RMI2Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000175643 > median(mRNA.Exp.log2.TPM.Target$ENSG00000175643),"RMI2High","RMI2Low")
mRNA.Exp.log2.TPM.Target$RMI2Group <- factor(mRNA.Exp.log2.TPM.Target$RMI2Group,levels = c("RMI2High","RMI2Low"))
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "RMI2Group", "ENSG00000130816", 
             fill = "RMI2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/RMI2.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "RMI2Group", "ENSG00000088305", 
             fill = "RMI2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/RMI2.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "RMI2Group", "ENSG00000119772", 
             fill = "RMI2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/RMI2.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)


#### ESR1 Low/High #### ENSG00000091831
mRNA.Exp.log2.TPM.Target$ESR1Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000091831 > median(mRNA.Exp.log2.TPM.Target$ENSG00000091831),"ESR1High","ESR1Low")
mRNA.Exp.log2.TPM.Target$ESR1Group <- factor(mRNA.Exp.log2.TPM.Target$ESR1Group,levels = c("ESR1High","ESR1Low"))
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "ESR1Group", "ENSG00000130816", 
             fill = "ESR1Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/ESR1.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "ESR1Group", "ENSG00000088305", 
             fill = "ESR1Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/ESR1.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "ESR1Group", "ENSG00000119772", 
             fill = "ESR1Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/ESR1.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)
#### SER1 Low/High ####
mRNA.Exp.log2.TPM.Target$ESR1Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000091831 > median(mRNA.Exp.log2.TPM.Target$ENSG00000091831),"ESR1High","ESR1Low")
mRNA.Exp.log2.TPM.Target$ESR1Group <- factor(mRNA.Exp.log2.TPM.Target$ESR1Group,levels = c("ESR1High","ESR1Low"))
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "ESR1Group", "ENSG00000130816", 
             fill = "ESR1Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/ESR1.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "ESR1Group", "ENSG00000088305", 
             fill = "ESR1Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/ESR1.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "ESR1Group", "ENSG00000119772", 
             fill = "ESR1Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/ESR1.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)


#### CEPB3 ####
mRNA.Exp.log2.TPM.Target$CPEB3Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000107864 > median(mRNA.Exp.log2.TPM.Target$ENSG00000107864),"CPEB3High","CPEB3Low")
mRNA.Exp.log2.TPM.Target$CPEB3Group <- factor(mRNA.Exp.log2.TPM.Target$CPEB3Group,levels = c("CPEB3High","CPEB3Low"))
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CPEB3Group", "ENSG00000130816", 
             fill = "CPEB3Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CPEB3.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CPEB3Group", "ENSG00000088305", 
             fill = "CPEB3Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CPEB3.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CPEB3Group", "ENSG00000119772", 
             fill = "CPEB3Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CPEB3.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)


#### MCM2 ENSG00000073111 ####
mRNA.Exp.log2.TPM.Target$MCM2Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000073111 > median(mRNA.Exp.log2.TPM.Target$ENSG00000073111),"MCM2High","MCM2Low")
mRNA.Exp.log2.TPM.Target$MCM2Group <- factor(mRNA.Exp.log2.TPM.Target$MCM2Group,levels = c("MCM2High","MCM2Low"))
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "MCM2Group", "ENSG00000130816", 
             fill = "MCM2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT1",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/MCM2.stratified.DNMT1.Expression.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "MCM2Group", "ENSG00000088305", 
             fill = "MCM2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3B",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/MCM2.stratified.DNMT3B.Expression.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "MCM2Group", "ENSG00000119772", 
             fill = "MCM2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("log2(TPM+1) DNMT3A",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/MCM2.stratified.DNMT3A.Expression.pdf",height = 2.5,width = 2.5)
######### Samples cluster -> According 5 gene => inconsistent ###########
#BiocManager::install("ConsensusClusterPlus")

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[Inter.Five.Gene$Select.Variables,Inter.Sample.361]

load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))

mRNA.Exp.log2.TPM.Target2 <- mRNA.Exp.log2.TPM.Target %>% t() %>% data.frame(check.names = F) %>% mutate(Sample.ID = substr(rownames(.),1,15)) %>%
  merge(.,risk_score_table_multi_cox2,by="Sample.ID")

mRNA.Exp.log2.TPM.Target3 <- mRNA.Exp.log2.TPM.Target2 %>% dplyr::select(Sample.ID:ENSG00000206195) %>%
  remove_rownames() %>% column_to_rownames("Sample.ID") %>% t() %>% data.frame(check.names = F)
rownames(mRNA.Exp.log2.TPM.Target3) = Inter.Five.Gene$Gene.name

annotation_col = data.frame(RiskScore = factor(mRNA.Exp.log2.TPM.Target2$RiskScore))
rownames(annotation_col) = colnames(mRNA.Exp.log2.TPM.Target3)
annotation_row = data.frame(Type = factor(rep(c("mRNA","lncRNA"),c(4,1))))
rownames(annotation_row) = rownames(mRNA.Exp.log2.TPM.Target3)

anno_color=list(RiskScore = c(High=brewer.pal(8,"Set1")[1],Low=brewer.pal(8,"Set1")[2]),
                Type= c(mRNA = brewer.pal(8,"Set1")[5],lncRNA = brewer.pal(8,"Set1")[3]))

library(pheatmap)
pheatmap(mRNA.Exp.log2.TPM.Target3,
         annotation_col = annotation_col,
         #annotation_row = annotation_row, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),#colorRampPalette(brewer.pal(11, "RdYlBu"))(50), #
         border_color = "black",
         annotation_colors = anno_color,
         #cellwidth = 7,cellheight = 5,
         #fontsize_row = 5,fontsize_col = 6,
         #clustering_distance_cols = "euclidean",
         #clustering_method = "complete",cluster_cols = T,
         annotation_row = annotation_row,
         cluster_rows = F,cluster_cols = T,show_colnames = F)
library(export)
graph2pdf(file="Pheatmap.Five.Gene.Riskgroup.pdf",height = 3,width = 5)

library(ConsensusClusterPlus)
dir.create("ConsensusClusterPlus")
title="./ConsensusClusterPlus/"

dataset <- mRNA.Exp.log2.TPM.Target3 <- mRNA.Exp.log2.TPM.Target2 %>% dplyr::select(Sample.ID:ENSG00000206195) %>%
  remove_rownames() %>% column_to_rownames("Sample.ID") %>% as.matrix() %>% t()
results <- ConsensusClusterPlus(dataset, maxK = 3,
                                reps = 50, pItem = 0.8,
                                pFeature = 1,  
                                clusterAlg = "hc", 
                                distance = "pearson",
                                title = title,
                                plot = "png") 

cluster <- results[[2]]$consensusClass %>% data.frame(check.names = F) %>% dplyr::rename(cluster=1) %>% rownames_to_column("Sample.ID")
write.csv(cluster,"ConsensusClusterPlus/Cluster.csv",row.names = F)

cluster.risk <- merge(cluster,mRNA.Exp.log2.TPM.Target2,by="Sample.ID")
tableR <- table(cluster.risk$cluster,cluster.risk$RiskScore) %>% as.matrix()
chisq.test(tableR)
fisher.test(tableR)



######################## TIDE + IPS ########
#The patient‚Äôs response to immunotherapy was inferred by the tumor immune dysfunction and exclusion (TIDE) score and immunophenoscore (IPS) TCIA #
#a patient‚Äôs IPS can be derived in an unbiased manner using machine learning by considering the four major categories of genes that determine immunogenicity (effector cells, immunosuppressive cells, MHC molecules, and immunomodulators)
#MHC molecular(MHC); Effector cells(EC); Immune Checkpoints(CP); Immunsuppressive Cells(SC)
#### TIDE ####
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))

TIDE.Data <- read.csv("../../TIDE/TCGA-LIHC.TIDE.csv") %>% mutate(Sample.ID = substr(Patient,1,15)) %>%
  merge(.,risk_score_table_multi_cox2,by="Sample.ID")

p<-ggboxplot(TIDE.Data, "RiskScore", "Dysfunction", fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("TIDE dysfunction",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = paste("TIDE.dysfunction.Boxplot.pdf",sep = ""),width = 3,height = 3)

p<-ggboxplot(TIDE.Data, "RiskScore", "Exclusion", fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("TIDE Exclusion",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = paste("TIDE.Exclusion.Boxplot.pdf",sep = ""),width = 3,height = 3)

#### IPS ####
#The IPS analysis suggested that the low-risk patients possessed a higher IPS score and a higher immunoreactivity phenotype, which were correlated with better immunotherapy response
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = substr(rownames(.),1,12))

IPS <- read.csv("../../TIDE/TCIA-ClinicalData.tsv",header = T,sep = "\t") %>% dplyr::select(barcode,ips_ctla4_neg_pd1_neg) %>%
  merge(.,risk_score_table_multi_cox2,by.x="barcode",by.y="Sample.ID")
#ips_ctla4_neg_pd1_neg
p<-ggboxplot(IPS, "RiskScore", "ips_ctla4_neg_pd1_neg", fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("TCIA IPS score",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = paste("TIDE.IPS.Boxplot.pdf",sep = ""),width = 3,height = 3)

setwd("K:/TCGA/Anlysis/LIHC")
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = substr(rownames(.),1,12)) %>%
  rownames_to_column("Sample")

ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

ImmunoPhenoScore = data.frame()
for (project in c("TCGA-LIHC")) {
  print(project)
  filename = paste("../mRNA.PanCancer.Exp/",project,".mRNA.TPM.csv",sep = "")
  dat <- read.csv(filename,check.names = F,row.names = 1)
  
  library('biomaRt')
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),values=rownames(dat),mart= mart)
  write.csv(G_list,"biomaRt.GeneID.tansfer.csv",row.names = F)
  
  dat.2 <- dat %>% data.frame(check.names = F) %>%
    rownames_to_column("ensembl_gene_id") %>% data.frame(check.names = F) %>%
    merge(.,G_list[,1:2],by="ensembl_gene_id") %>% 
    data.frame(check.names = F) %>%
    na.omit() %>% dplyr::select(-ensembl_gene_id)
  
  expr_mean=aggregate(.~hgnc_symbol,mean,data=dat.2)
  dat.3 <- expr_mean %>% remove_rownames() %>% column_to_rownames("hgnc_symbol")
  dat.4 <- log2(dat.3+1)
  
  sample_names<-names(dat.4)
  ### Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  IPSG<-read.table("../TIDE/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
  
  unique_ips_genes<-as.vector(unique(IPSG$NAME))
  
  IPS<-NULL
  MHC<-NULL
  CP<-NULL
  EC<-NULL
  SC<-NULL
  AZ<-NULL
  
  # Gene names in expression file
  GVEC<-row.names(dat.4)
  # Genes names in IPS genes file
  VEC<-as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind<-which(is.na(match(VEC,GVEC)))
  # List genes missing or differently named
  MISSING_GENES<-VEC[ind]
  dat<-IPSG[ind,]
  if (length(MISSING_GENES)>0) {
    cat("differently named or missing genes: ",MISSING_GENES,"\n")
  }
  for (x in 1:length(ind)) {
    print(IPSG[ind,])
  }
  
  for (i in 1:length(sample_names)) {	
    GE<-dat.4[[i]]
    mGE<-mean(GE)
    sGE<-sd(GE)
    Z1<-(dat.4[as.vector(IPSG$GENE),i]-mGE)/sGE
    W1<-IPSG$WEIGHT
    WEIGHT<-NULL
    MIG<-NULL
    k<-1
    for (gen in unique_ips_genes) {
      MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
      k<-k+1
    }
    WG<-MIG*WEIGHT
    MHC[i]<-mean(WG[1:10])
    CP[i]<-mean(WG[11:20])
    EC[i]<-mean(WG[21:24])
    SC[i]<-mean(WG[25:26])
    AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
    IPS[i]<-ipsmap(AZ[i])
  }
  ImmunoPhenoScore <- data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS,ProjectID = project) %>%
    rbind.data.frame(ImmunoPhenoScore)
}

write.csv(ImmunoPhenoScore,"IPS.LIHC.csv",row.names = F)
# MHC molecules (MHC)
# Immunomodulators (CP)
# Effector cells (EC)
# Suppressor cells (SC)
ImmunoPhenoScore <- read.csv("IPS.LIHC.csv")
ImmunoPhenoScore2 <- ImmunoPhenoScore %>% mutate(Sample = substr(SAMPLE,1,15)) %>%
  merge(.,risk_score_table_multi_cox2,by.x="Sample")

#ImmunoPhenoScore2$MHC.Zscore <- scale(ImmunoPhenoScore2$MHC)
library(ggpubr)
p <- ggboxplot(ImmunoPhenoScore2, "RiskScore", "MHC", fill = "RiskScore",
          palette = "Set1",add = "jitter")+
  ggthemes::theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none",axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        axis.title = element_text(size=20))+
  labs(x="",y=paste("Antigen presentation",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T,size=6)+
  scale_x_discrete(limits=c("High","Low"),labels = c("High risk","Low risk"))

ggsave(p,filename = "IPS.MHC.pdf",height = 4,width = 4)
#
#ImmunoPhenoScore2$EC.Zscore <- scale(ImmunoPhenoScore2$EC)
library(ggpubr)
p2 <- ggboxplot(ImmunoPhenoScore2, "RiskScore", "EC", fill = "RiskScore",
          palette = "Set1",add = "jitter")+
  ggthemes::theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none",axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        axis.title = element_text(size=20))+
  labs(x="",y=paste("Effector cells",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T,size=6)+
  scale_x_discrete(limits=c("High","Low"),labels = c("High risk","Low risk"))
ggsave(p2,filename = "IPS.Effector.pdf",height = 4,width = 4)
##
#ImmunoPhenoScore2$SC.Zscore <- scale(ImmunoPhenoScore2$SC)
library(ggpubr)
p3 <- ggboxplot(ImmunoPhenoScore2, "RiskScore", "SC", fill = "RiskScore",
          palette = "Set1",add = "jitter")+
  ggthemes::theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none",axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        axis.title = element_text(size=20))+
  labs(x="",y=paste("Suppressor cells",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T,size=6)+
  scale_x_discrete(limits=c("High","Low"),labels = c("High risk","Low risk"))#+
  #scale_y_continuous(limits = c(-5,3))
ggsave(p3,filename = "IPS.Suppressor.pdf",height = 4,width = 4)
#
#ImmunoPhenoScore2$CP.Zscore <- scale(ImmunoPhenoScore2$CP)
library(ggpubr)
p4<-ggboxplot(ImmunoPhenoScore2, "RiskScore", "CP", fill = "RiskScore",
          palette = "Set1",add = "jitter")+
  ggthemes::theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none",axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        axis.title = element_text(size=20))+
  labs(x="",y=paste("Checkpionts",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T,size=6)+
  scale_x_discrete(limits=c("High","Low"),labels = c("High risk","Low risk"))#+
  #scale_y_continuous(limits = c(-5,2))
ggsave(p4,filename = "IPS.Checkpionts.pdf",height = 4,width = 4)
#
#ImmunoPhenoScore2$IPS.Zscore <- scale(ImmunoPhenoScore2$IPS)
library(ggpubr)
p5<-ggboxplot(ImmunoPhenoScore2, "RiskScore", "IPS", fill = "RiskScore",
          palette = "Set1",add = "jitter")+
  ggthemes::theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none",axis.text.x = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        axis.title = element_text(size=20))+
  labs(x="",y=paste("IPS Z-score",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T,size=6)+
  scale_x_discrete(limits=c("High","Low"),labels = c("High risk","Low risk"))
  #scale_y_continuous(limits = c(-5,2))
ggsave(p5,filename = "IPS.IPS.pdf",height = 4,width = 4)

write.csv(ImmunoPhenoScore2,"IPS.LIHC2.csv",row.names = F)
######################## ESTIMATE => TME + TIMER2 (All Tools) => TIIC #######
#Furthermore, tumor microenvironment (TME) was evaluated using ESTIMATE algorithm. Tumor-infiltrating immune cells (TIICs) were inferred using TIMER2.
#ËÇøÁò§ÂæÆÁéØÂ¢É‰∏≠ÔºåÂÖçÁñ´ÁªÜËÉûÂíåÂü∫Ë¥®ÁªÜËÉûÊòØ‰∏§Áßç‰∏ªË¶ÅËÇ∫ËÇøÁò§ÁªÜËÉûÁ±ªÂûãÔºõËÄåESTIMATEÔºåÂà©Áî®Ë°®ËææË∞±Êï∞ÊçÆÊù•È¢ÑÊµãÂü∫Ë¥®ÁªÜËÉûÂíåÂÖçÁñ´ÁªÜËÉûËØÑÂàÜÔºåËøõËÄåÈ¢ÑÊµãËøô‰∏§ÁßçÁªÜËÉûÁöÑÂê´???

##BiocManager::install("estimate")
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
#mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor")
mRNA.Exp.Target <- mRNA.Exp[,Inter.Sample.361]

dat=log2(edgeR::cpm(mRNA.Exp.Target)+1)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),values=rownames(dat),mart= mart)
write.csv(G_list,"biomaRt.GeneID.tansfer.csv",row.names = F)

dat.2 <- dat %>% data.frame(check.names = F) %>%
  rownames_to_column("ensembl_gene_id") %>% data.frame(check.names = F) %>%
  merge(.,G_list[,1:2],by="ensembl_gene_id") %>% 
  data.frame(check.names = F) %>%
  na.omit() %>% dplyr::select(-ensembl_gene_id)

expr_mean=aggregate(.~hgnc_symbol,mean,data=dat.2)
dat.3 <- expr_mean %>% remove_rownames() %>% column_to_rownames("hgnc_symbol")

library(estimate)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina") ## Ê≥®ÊÑèËøô‰∏™platformÂèÇÊï∞
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}

pro='TCGA.LIHC'
scores=estimate(dat.3,pro)

scores=read.table("TCGA.LIHC_estimate_score.gct",skip = 2,header = T,check.names = F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
head(scores)
rownames(scores) <- rownames(scores) %>% str_replace_all("\\.","-")
write.csv(scores,"TCGA-LIHC.361Samples.ESTIMATE.csv",row.names = T)

###### ESTIMATE => calculate the correlation => Risk, Gene, Immune score ######
dir.create("ESTIMATE.Correlation")

mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor") %>%
  filter(Sample.ID %in% Inter.Sample.361)

## risk data
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))

## meta data
meta <- read.csv("clinical.cart.2021-08-19/clinical.tsv",sep = "\t",header = T)
#meta <- column_to_rownames(meta,var = "case_submitter_id")
meta=meta[,colnames(meta) %in% c("case_submitter_id",
                                 "race",
                                 "gender",
                                 "ajcc_pathologic_t",
                                 "ajcc_pathologic_n",
                                 "ajcc_pathologic_m",
                                 "age","Stage")] %>% 
  unique() %>% filter(case_submitter_id %in% mRNA.Phenodata$Case.ID)

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

Cor.ceRNA1 <- read.csv("Integrate.DownDElncRNA.UpDEmiRNA.DownDEmRNA.Final.csv")
Cor.ceRNA2 <- read.csv("Integrate.UpDElncRNA.DownDEmiRNA.UpDEmRNA.Final.csv")

scores <- read.csv("TCGA-LIHC.361Samples.ESTIMATE.csv")
#### lncRNA ####
Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv")

mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[unique(c(Cor.ceRNA1$lncRNA.ENSEMBLE.ID,Cor.ceRNA2$lncRNA.ENSEMBLE.ID)),Inter.Sample.361] %>%
  t() %>% data.frame(check.names = F) %>% mutate(SampleID = rownames(.)) %>%
  mutate(case_submitter_id = substr(SampleID,1,12)) %>%
  merge(.,meta,by="case_submitter_id") %>%
  merge(.,scores,by.x = "SampleID",by.y="X") %>%
  mutate(Sample.ID = substr(SampleID,1,15)) %>%
  merge(.,risk_score_table_multi_cox2,by="Sample.ID")
  

write.csv(mRNA.Exp.log2.TPM.Target,"ESTIMATE.lncRNA.Exp.log2.TPM.Target.csv")

library(ggpubr)
library(ggthemes)
Cortab <- data.frame()
for (lncRNA in unique(c(Cor.ceRNA1$lncRNA.ENSEMBLE.ID,Cor.ceRNA2$lncRNA.ENSEMBLE.ID))) {
  lncRNA.name <- (Inter.lncRNA.GeneCard %>% filter(lncRNA.ENSEMBLE.ID==lncRNA))$lncRNA.ID
  
  
  p<- ggscatter(mRNA.Exp.log2.TPM.Target, x= lncRNA, y = "ImmuneScore", color="RiskScore",
                palette = brewer.pal(8,"Set1")[1:2],#"#426671", size =2,#color = "group", palette = "Set1",
                add = "reg.line",
                add.params = list(color = "#764C29", fill = "#E7E1D7"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson",label.x = 1, label.sep = "\n"),
                ggtheme = theme_few()) +
    #stat_cor(method = "pearson",label.x = 1)+
    labs(x = paste("log2(TPM+1) ",lncRNA.name,sep = ""), 
         y = paste("Immune score",sep = ""))
  
  ggsave(p,filename = paste("ESTIMATE.Correlation/lncRNA.",lncRNA,".",lncRNA.name,".ImmuneScore.pdf",sep = ""),width = 4.3,height = 3)
  
  Cor <- cor.test(mRNA.Exp.log2.TPM.Target[[lncRNA]],mRNA.Exp.log2.TPM.Target$ImmuneScore)
  Cortab <- rbind(Cortab,c(lncRNA,lncRNA.name,"ImmuneScore",Cor$estimate,Cor$p.value,"lncRNA"))
}

#### mRNA ####
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[unique(c(Cor.ceRNA1$mRNA.ENSEMBLE.ID,Cor.ceRNA2$mRNA.ENSEMBLE.ID)),Inter.Sample.361] %>%
  t() %>% data.frame(check.names = F) %>% mutate(SampleID = rownames(.)) %>%
  mutate(case_submitter_id = substr(SampleID,1,12)) %>%
  merge(.,meta,by="case_submitter_id") %>%
  merge(.,scores,by.x = "SampleID",by.y="X") %>%
  mutate(Sample.ID = substr(SampleID,1,15)) %>%
  merge(.,risk_score_table_multi_cox2,by="Sample.ID")
write.csv(mRNA.Exp.log2.TPM.Target,"ESTIMATE.mRNA.Exp.log2.TPM.Target.csv")
for (lncRNA in unique(c(Cor.ceRNA1$mRNA.ENSEMBLE.ID,Cor.ceRNA2$mRNA.ENSEMBLE.ID))) {
  lncRNA.name <- (GeneType %>% filter(Gene.stable.ID == lncRNA))$HGNC.symbol %>% unique()
  
  p<- ggscatter(mRNA.Exp.log2.TPM.Target, x= lncRNA, y = "ImmuneScore", color="RiskScore",
                palette = brewer.pal(8,"Set1")[1:2],#"#426671", size =2,#color = "group", palette = "Set1",
                add = "reg.line",
                add.params = list(color = "#764C29", fill = "#E7E1D7"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson",label.x = 1, label.sep = "\n"),
                ggtheme = theme_few()) +
    #stat_cor(method = "pearson",label.x = 1)+
    labs(x = paste("log2(TPM+1) ",lncRNA.name,sep = ""), 
         y = paste("Immune score",sep = ""))
  
  ggsave(p,filename = paste("ESTIMATE.Correlation/mRNA.",lncRNA,".",lncRNA.name,".ImmuneScore.pdf",sep = ""),width = 4.3,height = 3)
  
  Cor <- cor.test(mRNA.Exp.log2.TPM.Target[[lncRNA]],mRNA.Exp.log2.TPM.Target$ImmuneScore)
  Cortab <- rbind(Cortab,c(lncRNA,lncRNA.name,"ImmuneScore",Cor$estimate,Cor$p.value,"mRNA"))
}

colnames(Cortab) <- c("ENSENBL","SYMBOL","IMMNUESCORE","RHO","PVALUE")
write.csv(Cortab,"ESTIMATE.Correlation/Immunescore.lncRNA.mRNA.Pearson.csv",row.names = F)
Cortab <- read.csv("ESTIMATE.Correlation/Immunescore.lncRNA.mRNA.Pearson.csv")
Cortab %>% filter(RHO >= 0.1) # CEP55

#ENSENBL  SYMBOL IMMNUESCORE       RHO       PVALUE  NA.
#1 ENSG00000087258   GNAO1 ImmuneScore 0.1144507 2.969233e-02 mRNA
#2 ENSG00000164619   BMPER ImmuneScore 0.1892280 2.996744e-04 mRNA
#3 ENSG00000129993 CBFA2T3 ImmuneScore 0.2932940 1.357831e-08 mRNA
#4 ENSG00000154928   EPHB1 ImmuneScore 0.1438933 6.167199e-03 mRNA
#5 ENSG00000113594    LIFR ImmuneScore 0.1515110 3.909162e-03 mRNA
#6 ENSG00000137269   LRRC1 ImmuneScore 0.1009010 5.544643e-02 mRNA
#7 ENSG00000138180   CEP55 ImmuneScore 0.1812904 5.375830e-04 mRNA

#### Risk + Immune score ####
p<- ggscatter(mRNA.Exp.log2.TPM.Target, x= "total_risk_score", y = "ImmuneScore", color="RiskScore",
              palette = brewer.pal(8,"Set1")[1:2],#"#426671", size =2,#color = "group", palette = "Set1",
              add = "reg.line",
              add.params = list(color = "#764C29", fill = "#E7E1D7"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
              cor.coeff.args = list(method = "pearson",label.x = 1, label.sep = "\n"),
              ggtheme = theme_few()) +
  #stat_cor(method = "pearson",label.x = 1)+
  labs(x = paste("total_risk_score",sep = ""), 
       y = paste("Immune score",sep = ""))

ggsave(p,filename = paste("ESTIMATE.Correlation/","risk_score",".ImmuneScore.pdf",sep = ""),width = 4.3,height = 3)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "RiskScore", "ImmuneScore", fill = "RiskScore",
          palette = "Set1",
          add = "jitter")+
  theme_few()+
  theme(legend.position = "none")+
  labs(x="RiskScore",y=paste("ImmuneScore",sep = ""))+
  stat_compare_means(method = "t.test")
ggsave(p,filename = paste("ESTIMATE.Correlation/","risk_score",".ImmuneScore.Boxplot.pdf",sep = ""),width = 3,height = 3)


###### TIMER2 =>calculate the correlation => Risk, Gene, Immune score #####
# immunedeconv  => quantiseq, timer, cibersort, cibersort_abs, mcp_counter, xcell, epic
# TIMER2 
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))

Timer2 <- read.csv("../../TIMER2/infiltration_estimation_for_tcga.csv",check.names = F)
#TIMER CIBERSORT CIBERSORT-ABS QUANTISEQ MCPCOUNTER XCELL EPIC
risk_score.Immune <- merge.data.frame(risk_score_table_multi_cox2,Timer2,by.x="Sample.ID",by.y="cell_type")
dim(risk_score.Immune)
write.csv(risk_score.Immune,"Survival.Signature.Riskscore.Immune.TIMER2.csv",row.names = F)
risk_score.Immune <- risk_score.Immune %>% mutate(RiskScore = factor(RiskScore,levels = c("High","Low"))) %>%
  dplyr::arrange(RiskScore) %>% mutate(Sample.ID = factor(Sample.ID,levels=Sample.ID))

dir.create("TIMER2")
#### TIMER ####
risk_score.Immune.TIMER <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore",ends_with("TIMER"))) %>%
  gather(Immune.Cell,Percent,`B cell_TIMER`:`Myeloid dendritic cell_TIMER`)
risk_score.Immune.TIMER2 <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore",ends_with("TIMER")))

p<-ggplot(risk_score.Immune.TIMER, aes( x = Sample.ID,y=Percent,fill = Immune.Cell))+
  #geom_colÂíågeom_barËøô‰∏§Êù°ÂëΩ‰ª§ÈÉΩÂèØ‰ª•ÁªòÂà∂Â†ÜÂè†Êü±ÂΩ¢???
  geom_col(position = 'fill', width = 0.6)+
  #geom_bar(position = "fill", stat = "identity", width = 0.6) 
  theme_few()+   #ËÆæÁΩÆÂõ∫ÂÆö‰∏ªÈ¢ò‰∏∫‰º†ÁªüÁöÑÁôΩËâ≤ËÉåÊôØÂíåÊ∑±ÁÅ∞Ëâ≤ÁöÑÁΩëÊ†ºÁ∫ø
  scale_fill_manual(values=brewer.pal(8,"Set1"))+  #Ëá™ÂÆö‰πâÁöÑ‰øÆÊîπÂ°´ÂÖÖÈ¢úËâ≤
  scale_y_continuous(expand = c(0,0))+# Ë∞ÉÊï¥yËΩ¥Â±ûÊÄßÔºå‰ΩøÊü±Â≠ê‰∏éXËΩ¥ÂùêÊ†áÊé•???
  labs(x="",y="Immune cell(%)")+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,
                             vjust = 0.5),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12)
  )+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_vline(xintercept = c(180))+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) #‰øÆÊîπÂõæ‰æãÁöÑÊ°ÜÂ§ßÂ∞è

ggsave(p,filename = "TIMER2.RiskScore.Barplot.TIMER.pdf",height = 3,width = 7)

for (immune in unique(risk_score.Immune.TIMER$Immune.Cell)) {
  p<-ggboxplot(risk_score.Immune.TIMER2, "RiskScore", immune, fill = "RiskScore",palette = "Set1",add = "jitter")+
    theme_few()+
    #facet_wrap(~Immune.Cell)+
    theme(legend.position = "none")+
    labs(x="",y=paste("TIMER predicted score ",sep = ""))+
    stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
  ggsave(p,filename = paste("TIMER2/TIMER2.",immune,".RiskScore.Boxplot.TIMER.pdf",sep = ""),height = 2.5,width = 2.5)
}

#### CIBERSORT -> More accurate ####
risk_score.Immune.CIBERSORT <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore","OS.time", "OS" ,"total_risk_score",ends_with("CIBERSORT"))) %>%
  gather(Immune.Cell,Percent,`B cell naive_CIBERSORT`:`Neutrophil_CIBERSORT`)
risk_score.Immune.CIBERSORT2 <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore","OS.time", "OS" ,"total_risk_score",ends_with("CIBERSORT")))

colnames(risk_score.Immune.CIBERSORT2) = str_remove_all(colnames(risk_score.Immune.CIBERSORT2),"_CIBERSORT")
risk_score.Immune.CIBERSORT3 <- risk_score.Immune.CIBERSORT2 %>% gather(Celltype,Preidicted,`B cell naive`:`Neutrophil`)


p<-ggboxplot(risk_score.Immune.CIBERSORT3, "Celltype","Preidicted", color = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "top")+
  labs(x="",y=paste("CIBERSORT predicted score ",sep = ""))+
  stat_compare_means(aes(group = RiskScore),method = "wilcox.test",label="p.format",hide.ns = T,label.y = 0.75)+
  coord_flip()+
  scale_y_sqrt(limits = c(0,0.8))+
  theme(axis.text.y = element_text(size = 12,face = "bold"),axis.text = element_text(colour = "black"))

ggsave(p,filename = "Figure/RiskScore.CIBERSORT.pdf",height = 8,width = 5.5)

### Immune Cell COX ###
library("survival")
library("survminer")
Results <- data.frame()
# KM
for (cell in unique(risk_score.Immune.CIBERSORT3$Celltype)) {
  test.middata <- risk_score.Immune.CIBERSORT3 %>% filter(Celltype == cell) %>% dplyr::select("OS.time", "OS","Celltype","Preidicted") %>%
    mutate(Group = ifelse(Preidicted > median(Preidicted),"High","Low"))
  data.survdiff=survdiff(Surv(OS.time, OS)~Group,data=test.middata)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  
  m = coxph(Surv(OS.time, OS) ~ Preidicted, data =  test.middata)
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, 
                     se = se, z = beta/se, 
                     p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, 
                     HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)), 6)
  Results <- rbind(Results,c(cell,p.val,tmp["Preidicted",]))
}

colnames(Results) = c("CellType","KM.P","coef","se","z","COX.P","HR","HRse","HRz","HRp","HRCI.LL","HRCI.UL")
write.csv(Results,file = "Figure/CIBERSORT.ImmuneCell.Surviaval.csv",row.names = F)

library(pheatmap)
pheatmap.data <- Results %>% dplyr::select(CellType,HR) %>% mutate(CellType=factor(CellType,levels = rev(CellType))) %>% arrange(CellType) %>%
  remove_rownames() %>% column_to_rownames("CellType") %>% mutate(HR = as.numeric(HR)) %>%
  mutate(HR = ifelse(HR > 1,"Risky","Protective")) %>% as.matrix()

pheatmap.Sig <- Results %>% dplyr::select(CellType,KM.P,COX.P) %>% mutate(CellType=factor(CellType,levels = rev(CellType))) %>% arrange(CellType) %>%
  remove_rownames() %>% column_to_rownames("CellType") %>% mutate(KM.P = as.numeric(KM.P),COX.P = as.numeric(COX.P)) %>%
  mutate(HR = ifelse(KM.P <= 0.05 | COX.P <= 0.05,"*","")) %>% dplyr::select(HR) %>% as.matrix()

library(RColorBrewer)
pheatmap(pheatmap.data, 
         display_numbers = pheatmap.Sig, cluster_cols = F,cluster_rows = F,
         #annotation_row = annotation_row, 
         color = colorRampPalette(c(brewer.pal(9,"Set1")[2], "white", brewer.pal(9,"Set1")[1])),
         border_color = "black",cellwidth = 8,cellheight = 6,fontsize_row = 5,fontsize_col = 7)

library(ComplexHeatmap)
colors <- structure(brewer.pal(9,"Set1")[1:2],names=c("Risky","Protective"))
Heatmap(pheatmap.data,col = colors,border_gp = gpar(col = "black"),rect_gp = gpar(col = "black"),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(pheatmap.Sig[i,j], x = x, y = y,gp = gpar(fontsize = 20,col="white"))})
graph2pdf(file=paste("Figure/CIBERSORT.ImmunceCell.Survival.pdf",sep = ''),height = 7.8,width = 4)

# Macrophage M0 + T cell CD4+ memory resting
cell = "T cell CD4+ memory resting" 
test.middata <- risk_score.Immune.CIBERSORT3 %>% filter(Celltype == cell) %>% dplyr::select("OS.time", "OS","Celltype","Preidicted") %>%
  mutate(Group = ifelse(Preidicted > median(Preidicted),"High","Low"))
sfit1=survfit(Surv(OS.time, OS)~Group,data=test.middata)
p<-ggsurvplot(sfit1,pval =TRUE, data = test.middata, #risk.table = TRUE,
              surv.median.line = "hv", 
              legend.title = cell,
              legend.labs = c("High", "Low"),
              conf.int.style = "step",
              xlab = "Time in days",
              #risk.table = "abs_pct",
              risk.table.y.text.col = T,
              risk.table.y.text = FALSE,
              conf.int = TRUE,
              palette = "Set1",
              ggtheme = theme_few())

print(p)
#dev.off()
library(export)
graph2pdf(file=paste("Figure/CIBERSORT.TcellCD4+memoryresting.KM.pdf",sep = ''),height = 2.5,width = 2)

### CBX2 Immune cell Correlation ###
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
colnames(mRNA.Exp.log2.TPM) = substr(colnames(mRNA.Exp.log2.TPM),1,15)
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000173894"),unique(risk_score.Immune.CIBERSORT3$Sample.ID)] %>% t() %>%
  data.frame(check.names = F) %>% rownames_to_column("Sample.ID") %>% merge(.,risk_score.Immune.CIBERSORT3,by="Sample.ID")

library(psych)
Results.Cor <- data.frame()
# KM
for (cell in unique(risk_score.Immune.CIBERSORT3$Celltype)) {
  test.middata <- mRNA.Exp.log2.TPM.Target %>% filter(Celltype == cell)
  #partial.r(data=test.middata, x="Preidicted", y="ENSG00000173894")
  test <- cor.test(test.middata$ENSG00000173894,test.middata$Preidicted,method = "spearman")
  Results.Cor <- rbind(Results.Cor,c(cell,test$estimate,test$p.value))
}
colnames(Results.Cor) =c("CellType","Rho","Pvalue")
write.csv(Results.Cor,"Figure/CBX2.ImmuneCellCor.csv",row.names = F)

Results.Cor <- read.csv("Figure/CBX2.ImmuneCellCor.csv")
Results.Cor <- Results.Cor %>% mutate(Label = ifelse(Pvalue > 0.05,"",ifelse(Pvalue >= 0.01,"*",ifelse(Pvalue >= 0.001,"**","***")))) %>% 
  mutate(Direction = ifelse(Rho > 0,"P","N"))
Results.Cor$CellType = factor(Results.Cor$CellType,levels = Results.Cor$CellType)
p<-ggplot(Results.Cor,aes(x=CellType,y=Rho,fill=Direction)) + geom_bar(stat="identity", position=position_dodge())+
  theme_few() + coord_flip()+
  scale_fill_manual(values = brewer.pal(9,"Set1")[2:1])+
  labs(x="",y="Spearman's correaltion efficience")+
  theme(legend.position = "top")+
  geom_text(aes(label=Label),angle=90,size=7)+
  ylim(-0.2,0.2)

ggsave(p,filename = "Figure/CIBERSORT.ImmuneCell.CBX2.Spearmen.pdf",height = 8,width = 3.5)

### CEP55 Immune cell Correlation ###
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
colnames(mRNA.Exp.log2.TPM) = substr(colnames(mRNA.Exp.log2.TPM),1,15)
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000138180"),unique(risk_score.Immune.CIBERSORT3$Sample.ID)] %>% t() %>%
  data.frame(check.names = F) %>% rownames_to_column("Sample.ID") %>% merge(.,risk_score.Immune.CIBERSORT3,by="Sample.ID")

library(psych)
Results.Cor <- data.frame()
# KM
for (cell in unique(risk_score.Immune.CIBERSORT3$Celltype)) {
  test.middata <- mRNA.Exp.log2.TPM.Target %>% filter(Celltype == cell)
  #partial.r(data=test.middata, x="Preidicted", y="ENSG00000138180")
  test <- cor.test(test.middata$ENSG00000138180,test.middata$Preidicted,method = "spearman")
  Results.Cor <- rbind(Results.Cor,c(cell,test$estimate,test$p.value))
}
colnames(Results.Cor) =c("CellType","Rho","Pvalue")
write.csv(Results.Cor,"Figure/CEP55.ImmuneCellCor.csv",row.names = F)

Results.Cor <- read.csv("Figure/CEP55.ImmuneCellCor.csv")
Results.Cor <- Results.Cor %>% mutate(Label = ifelse(Pvalue > 0.05,"",ifelse(Pvalue >= 0.01,"*",ifelse(Pvalue >= 0.001,"**","***")))) %>% 
  mutate(Direction = ifelse(Rho > 0,"P","N"))
Results.Cor$CellType = factor(Results.Cor$CellType,levels = Results.Cor$CellType)
p<-ggplot(Results.Cor,aes(x=CellType,y=Rho,fill=Direction)) + geom_bar(stat="identity", position=position_dodge())+
  theme_few() + coord_flip()+
  scale_fill_manual(values = brewer.pal(9,"Set1")[2:1])+
  labs(x="",y="Spearman's correaltion efficience")+
  theme(legend.position = "top")+
  geom_text(aes(label=Label),angle=90,size=7)#+
  #ylim(-0.2,0.2)

ggsave(p,filename = "Figure/CIBERSORT.ImmuneCell.CEP55.Spearmen.pdf",height = 8,width = 3.5)



####
p<-ggplot(risk_score.Immune.CIBERSORT, aes( x = Sample.ID,y=Percent,fill = Immune.Cell))+
  #geom_colÂíågeom_barËøô‰∏§Êù°ÂëΩ‰ª§ÈÉΩÂèØ‰ª•ÁªòÂà∂Â†ÜÂè†Êü±ÂΩ¢???
  geom_col(position = 'fill', width = 0.6)+
  #geom_bar(position = "fill", stat = "identity", width = 0.6) 
  theme_few()+   #ËÆæÁΩÆÂõ∫ÂÆö‰∏ªÈ¢ò‰∏∫‰º†ÁªüÁöÑÁôΩËâ≤ËÉåÊôØÂíåÊ∑±ÁÅ∞Ëâ≤ÁöÑÁΩëÊ†ºÁ∫ø
  scale_fill_manual(values=c(brewer.pal(8,"Set1"),brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")))+  #Ëá™ÂÆö‰πâÁöÑ‰øÆÊîπÂ°´ÂÖÖÈ¢úËâ≤
  scale_y_continuous(expand = c(0,0))+# Ë∞ÉÊï¥yËΩ¥Â±ûÊÄßÔºå‰ΩøÊü±Â≠ê‰∏éXËΩ¥ÂùêÊ†áÊé•???
  labs(x="",y="Immune cell(%)")+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,
                             vjust = 0.5),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12)
  )+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_vline(xintercept = c(180))+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1)) #‰øÆÊîπÂõæ‰æãÁöÑÊ°ÜÂ§ßÂ∞è

ggsave(p,filename = "TIMER2.RiskScore.Barplot.CIBERSORT.pdf",height = 3,width = 10)

#dir.create("TIMER2")
for (immune in unique(risk_score.Immune.CIBERSORT$Immune.Cell)) {
  p<-ggboxplot(risk_score.Immune.CIBERSORT2, "RiskScore", immune, fill = "RiskScore",palette = "Set1",add = "jitter")+
    theme_few()+
    #facet_wrap(~Immune.Cell)+
    theme(legend.position = "none")+
    labs(x="",y=paste("CIBERSORT predicted score ",sep = ""))+
    stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
  ggsave(p,filename = paste("TIMER2/TIMER2.",immune,".RiskScore.Boxplot.CIBERSORT.pdf",sep = ""),height = 2.5,width = 2.5)
}

#### XCELL ####
risk_score.Immune.XCELL2 <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore",ends_with("XCELL")))
ggboxplot(risk_score.Immune.XCELL2, "RiskScore", "microenvironment score_XCELL", fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("XCELL predicted score ",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

p<-ggboxplot(risk_score.Immune.XCELL2, "RiskScore", "immune score_XCELL", fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("XCELL predicted score ",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = paste("XCELL.risk_score",".ImmuneScore.Boxplot.pdf",sep = ""),width = 3,height = 3)


risk_score.Immune.XCELL <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore",ends_with("XCELL"))) %>%
  dplyr::select(Sample.ID:`T cell regulatory (Tregs)_XCELL`) %>%
  gather(Immune.Cell,Percent,`Myeloid dendritic cell activated_XCELL`:`T cell regulatory (Tregs)_XCELL`)

p<-ggplot(risk_score.Immune.XCELL, aes( x = Sample.ID,y=Percent,fill = Immune.Cell))+
  #geom_colÂíågeom_barËøô‰∏§Êù°ÂëΩ‰ª§ÈÉΩÂèØ‰ª•ÁªòÂà∂Â†ÜÂè†Êü±ÂΩ¢???
  geom_col(position = 'fill', width = 0.6)+
  #geom_bar(position = "fill", stat = "identity", width = 0.6) 
  theme_few()+   #ËÆæÁΩÆÂõ∫ÂÆö‰∏ªÈ¢ò‰∏∫‰º†ÁªüÁöÑÁôΩËâ≤ËÉåÊôØÂíåÊ∑±ÁÅ∞Ëâ≤ÁöÑÁΩëÊ†ºÁ∫ø
  scale_fill_manual(values=c(brewer.pal(8,"Set1"),brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(8,"Set2")))+  #Ëá™ÂÆö‰πâÁöÑ‰øÆÊîπÂ°´ÂÖÖÈ¢úËâ≤
  scale_y_continuous(expand = c(0,0))+# Ë∞ÉÊï¥yËΩ¥Â±ûÊÄßÔºå‰ΩøÊü±Â≠ê‰∏éXËΩ¥ÂùêÊ†áÊé•???
  labs(x="",y="Immune cell(%)")+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,
                             vjust = 0.5),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12)
  )+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_vline(xintercept = c(180))
  #guides(fill=guide_legend(keywidth = 1, keyheight = 1)) #‰øÆÊîπÂõæ‰æãÁöÑÊ°ÜÂ§ßÂ∞è

ggsave(p,filename = "TIMER2.RiskScore.Barplot.XCELL.pdf",height = 8,width = 12)

for (immune in unique(risk_score.Immune.XCELL$Immune.Cell)) {
  p<-ggboxplot(risk_score.Immune.XCELL2, "RiskScore", immune, fill = "RiskScore",palette = "Set1",add = "jitter")+
    theme_few()+
    #facet_wrap(~Immune.Cell)+
    theme(legend.position = "none")+
    labs(x="",y=immune)+
    stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
  ggsave(p,filename = paste("TIMER2/TIMER2.",immune,".RiskScore.Boxplot.XCELL.pdf",sep = ""),height = 2.5,width = 2.5)
}

#### MCPCOUNTER => cytotoxicity score -> More accurate ####
library(ggpubr)
library(ggthemes)
risk_score.Immune.MCPCOUNTER <- risk_score.Immune %>% dplyr::select(c("Sample.ID","RiskScore",ends_with("MCPCOUNTER")))
colnames(risk_score.Immune.MCPCOUNTER) = str_remove_all(colnames(risk_score.Immune.MCPCOUNTER),"_MCPCOUNTER")

p<-ggboxplot(risk_score.Immune.MCPCOUNTER, "RiskScore", "cytotoxicity score_MCPCOUNTER", fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("MCPCOUNTER predicted score ",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)

ggsave(p,filename = paste("MCPCOUNTER.risk_score",".cytotoxicity score_MCPCOUNTER.Boxplot.pdf",sep = ""),width = 3,height = 3)

setwd("K:/TCGA/Anlysis/LIHC")
risk_score.Immune.MCPCOUNTER.New <- risk_score.Immune.MCPCOUNTER %>% gather(Celltype,Preidicted,`T cell`:`Cancer associated fibroblast`)

p<-ggboxplot(risk_score.Immune.MCPCOUNTER.New, "Celltype","Preidicted", color = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "top")+
  labs(x="",y=paste("MCPCOUNTER predicted score ",sep = ""))+
  stat_compare_means(aes(group = RiskScore),method = "wilcox.test",label="p.format",hide.ns = T,label.y = 25)+
  coord_flip()+
  scale_y_sqrt()+
  theme(axis.text.y = element_text(size = 12,face = "bold"),axis.text = element_text(colour = "black"))

ggsave(p,filename = "Figure/RiskScore.MCPCOUNTER.pdf",height = 5,width = 5.5)


####### CBX2 CEP55 correlation with immune cell ####
##GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
##GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180
Timer2 <- read.csv("../../TIMER2/infiltration_estimation_for_tcga.csv",check.names = F) #%>% dplyr::select("cell_type",ends_with("XCELL"))

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
load("Inter.miRNA.mRNA.Samples.361.Rdata")

Timer2 <- Timer2 %>% dplyr::filter(cell_type %in% substr(Inter.Sample.361,1,15)) %>%
  mutate(cell_type = factor(cell_type,levels = substr(Inter.Sample.361,1,15))) %>% dplyr::arrange(cell_type)
colnames(Timer2) <- colnames(Timer2) %>% str_replace_all("/",".") %>%
  str_replace_all(" ",".")
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000173894","ENSG00000138180"),Inter.Sample.361] %>% t() %>%
  data.frame(check.names = F) %>% cbind.data.frame(Timer2)



library(ggthemes)
cor.data <- data.frame()
for (gene in c("ENSG00000173894","ENSG00000138180")) {
  for (immunecell in colnames(Timer2)[2:120]) {
    if (gene == "ENSG00000173894") {
      gene.name = "CBX2"
    }else{gene.name = "CEP55"}
    p<- ggscatter(mRNA.Exp.log2.TPM.Target, x= gene, y = immunecell, #color=brewer.pal(8,"Set1")[1],
                  #palette = brewer.pal(8,"Set1")[1],#"#426671", size =2,#color = "group", palette = "Set1",
                  add = "reg.line",
                  add.params = list(color = brewer.pal(8,"Set1")[1], fill = "#E7E1D7"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson",label.x = 1, label.sep = "\n"),
                  ggtheme = theme_few()) +
      #stat_cor(method = "pearson",label.x = 1)+
      labs(x = paste("log2(TPM+1) ",gene.name,sep = ""), 
           y = paste(immunecell,sep = ""))
    
    test <- cor.test(mRNA.Exp.log2.TPM.Target[[gene]],mRNA.Exp.log2.TPM.Target[[immunecell]],method = "pearson")
    cor.data <- rbind.data.frame(cor.data,c(gene.name,immunecell,test$p.value,test$estimate))
    
    immunecell <- str_replace(immunecell,"/",".")
    ggsave(p,filename = paste("TIMER2.CBX2.CEP55/",gene.name,".",immunecell,".pdf",sep = ""),width = 2.5,height = 2.5)
    
  }
}  

cor.data <- cor.data %>% dplyr::rename(Gene=1,ImmuneCell=2,pvalue=3,Rho=4)
write.csv(cor.data,"TIMER2.CBX2.CEP55/Corrlation.csv",row.names = F)

####@ CIBERSORT #######
cor.data <- read.csv("TIMER2.CBX2.CEP55/Corrlation.csv")

cor.data2 <- cor.data %>% 
  filter(!str_detect(ImmuneCell,"score")) %>%
  filter(!str_detect(ImmuneCell,"uncharacterized")) %>%
  mutate(Method=str_remove(ImmuneCell,".*_")) %>%
  arrange(ImmuneCell) %>%
  mutate(Size=abs(Rho)) %>%
  filter(Method=="CIBERSORT") %>%
  mutate(ImmuneCell=(str_remove_all(ImmuneCell,"_CIBERSORT") %>% str_replace_all("\\."," ")))

cor.data3 <- cor.data2 %>% filter(pvalue <= 0.05)

Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p21<-ggplot()+
  geom_point(data=cor.data2,aes(x=Gene,y=ImmuneCell,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(Gene,ImmuneCell,size=Size,alpha=2,fill=Rho,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p21,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.CIBERSORT.pdf",width = 4,height =5)

####@ XCELL #######
cor.data <- read.csv("TIMER2.CBX2.CEP55/Corrlation.csv")

cor.data2 <- cor.data %>% 
  filter(!str_detect(ImmuneCell,"score")) %>%
  filter(!str_detect(ImmuneCell,"uncharacterized")) %>%
  mutate(Method=str_remove(ImmuneCell,".*_")) %>%
  arrange(ImmuneCell) %>%
  mutate(Size=abs(Rho)) %>%
  filter(Method=="XCELL") %>%
  mutate(ImmuneCell=(str_remove_all(ImmuneCell,"_XCELL") %>% str_replace_all("\\."," ")))

cor.data3 <- cor.data2 %>% filter(pvalue <= 0.05)

Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p22<-ggplot()+
  geom_point(data=cor.data2,aes(x=Gene,y=ImmuneCell,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(Gene,ImmuneCell,size=Size,alpha=2,fill=Rho,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p22,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.XCELL.pdf",width = 4,height =6)

####@ EPIC #######
cor.data <- read.csv("TIMER2.CBX2.CEP55/Corrlation.csv")

cor.data2 <- cor.data %>% 
  filter(!str_detect(ImmuneCell,"score")) %>%
  filter(!str_detect(ImmuneCell,"uncharacterized")) %>%
  mutate(Method=str_remove(ImmuneCell,".*_")) %>%
  arrange(ImmuneCell) %>%
  mutate(Size=abs(Rho)) %>%
  filter(Method=="EPIC") %>%
  mutate(ImmuneCell=(str_remove_all(ImmuneCell,"_EPIC") %>% str_replace_all("\\."," ")))

cor.data3 <- cor.data2 %>% filter(pvalue <= 0.05)

Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p23<-ggplot()+
  geom_point(data=cor.data2,aes(x=Gene,y=ImmuneCell,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(Gene,ImmuneCell,size=Size,alpha=2,fill=Rho,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p23,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.EPIC.pdf",width = 4,height =2)

####@ MCPCOUNTER #######
cor.data <- read.csv("TIMER2.CBX2.CEP55/Corrlation.csv")

cor.data2 <- cor.data %>% 
  filter(!str_detect(ImmuneCell,"score")) %>%
  filter(!str_detect(ImmuneCell,"uncharacterized")) %>%
  mutate(Method=str_remove(ImmuneCell,".*_")) %>%
  arrange(ImmuneCell) %>%
  mutate(Size=abs(Rho)) %>%
  filter(Method=="MCPCOUNTER") %>%
  mutate(ImmuneCell=(str_remove_all(ImmuneCell,"_MCPCOUNTER") %>% str_replace_all("\\."," ")))

cor.data3 <- cor.data2 %>% filter(pvalue <= 0.05)

Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p25<-ggplot()+
  geom_point(data=cor.data2,aes(x=Gene,y=ImmuneCell,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(Gene,ImmuneCell,size=Size,alpha=2,fill=Rho,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p25,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.MCPCOUNTER.pdf",width = 4,height =3)

####@ TIMER #######
cor.data <- read.csv("TIMER2.CBX2.CEP55/Corrlation.csv")

cor.data2 <- cor.data %>% 
  filter(!str_detect(ImmuneCell,"score")) %>%
  filter(!str_detect(ImmuneCell,"uncharacterized")) %>%
  mutate(Method=str_remove(ImmuneCell,".*_")) %>%
  arrange(ImmuneCell) %>%
  mutate(Size=abs(Rho)) %>%
  filter(Method=="TIMER") %>%
  mutate(ImmuneCell=(str_remove_all(ImmuneCell,"_TIMER") %>% str_replace_all("\\."," ")))

cor.data3 <- cor.data2 %>% filter(pvalue <= 0.05)

Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p26<-ggplot()+
  geom_point(data=cor.data2,aes(x=Gene,y=ImmuneCell,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(Gene,ImmuneCell,size=Size,alpha=2,fill=Rho,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p26,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.TIMER.pdf",width = 3.5,height =2.5)

####@ QUANTISEQ #######
cor.data <- read.csv("TIMER2.CBX2.CEP55/Corrlation.csv")

cor.data2 <- cor.data %>% 
  filter(!str_detect(ImmuneCell,"score")) %>%
  filter(!str_detect(ImmuneCell,"uncharacterized")) %>%
  mutate(Method=str_remove(ImmuneCell,".*_")) %>%
  arrange(ImmuneCell) %>%
  mutate(Size=abs(Rho)) %>%
  filter(Method=="QUANTISEQ") %>%
  mutate(ImmuneCell=(str_remove_all(ImmuneCell,"_QUANTISEQ") %>% str_replace_all("\\."," ")))

cor.data3 <- cor.data2 %>% filter(pvalue <= 0.05)

Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p27<-ggplot()+
  geom_point(data=cor.data2,aes(x=Gene,y=ImmuneCell,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(Gene,ImmuneCell,size=Size,alpha=2,fill=Rho,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p27,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.QUANTISEQ.pdf",width = 4,height =4)

####@ ImmuCellAI #######
library(readxl)
cor.data <- read_xlsx("TIMER2.CBX2.CEP55/ImmuneAndExprTable.Cor.ImmuCellAI.xlsx")
cor.data2 <- cor.data %>% mutate(Size=abs(cor)) %>% filter(!str_detect(cell_type,"Score"))
cor.data3 <- cor.data2 %>% filter(p_value <= 0.05)
#Anno.data <- cor.data2 %>% dplyr::select(-Gene) %>% unique() %>% mutate(Count=1)

p30<-ggplot()+
  geom_point(data=cor.data2,aes(x=symbol,y=cell_type,fill=cor,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,4))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=cor.data3,aes(symbol,cell_type,size=Size,alpha=2,fill=cor,stroke=1.5),color="black",shape=21)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(p30,filename = "TIMER2.CBX2.CEP55/CBX2.CEP55.ImmuneCell.ImmuCellAI.pdf",width = 3.5,height =4.5)


######################## lncRNA location ######################
lnclocation <- read_xlsx("lncRNA.lncLocator.xlsx")

library(ggpubr)
library(RColorBrewer)
library(ggthemes)
p<-ggplot(lnclocation,aes(Location,lncLocator,fill = Location))+
  geom_col()+
  scale_fill_manual(values = brewer.pal(8,"Set1")[1:5])+
  facet_wrap(~lncRNA)+
  labs(x="",y="lncLocator score")+
  theme_few()+
  theme(axis.text = element_text(colour = "black"))+
  theme(legend.position = "none")

ggsave(p,filename = "lncLocator.location.pdf",width = 7,height = 4)






######################## DElncRNA DEmiRNA DEmRNA Validated with Cohorts ###################
###### miRNA Validate ######
load("../LIHC.GSE140845.miRNA.Rdata")
dir.create("miRNA.Validate")

setwd("K:/TCGA/Anlysis/LIHC/miRNA.Validate")

library(AnnoProbe)
library(GEOquery)
library(ggplot2) 
library(ggstatsplot) 
library(reshape2)
library(patchwork)
library(limma)
library(tidyverse)
library(ggthemes)

#### miRNA - GSE115016 -> GEL -> Limma ####

dir.create("miRNA.GSE115016")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE115016'
eSet <- getGEO(gse, 
               destdir = 'miRNA.GSE115016')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "miRNA.GSE115016/GSE115016.Rdata")

load("miRNA.GSE115016/GSE115016.Rdata")
group_list=ifelse(str_detect(pd$title,"normal"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "miRNA.GSE115016"))
Anno_data2 <- Anno_data[,1:4]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data2,by.x = "probe_id",by.y="ID") %>%
  filter(str_detect(`Transcript ID(Array Design)`,"hsa"))

write.csv(hsa.deg,file = "miRNA.GSE115016/GSE115016.de.csv",row.names = F)

Target.gene <- hsa.deg %>% filter(`Transcript ID(Array Design)` %in% Inter.miRNA) %>% 
  mutate(`Transcript ID(Array Design)` = factor(.$`Transcript ID(Array Design)`,levels=Inter.miRNA))

p<-ggplot(Target.gene,aes(x=`Transcript ID(Array Design)`,y=logFC,fill=`Transcript ID(Array Design)`))+
  geom_col() + 
  geom_text(aes(x=`Transcript ID(Array Design)`,y=logFC,label = sprintf("%2.3f", adj.P.Val)))+
  theme_few() + 
  coord_flip() + labs(x="")+
  theme(legend.position = "none")

ggsave(p,filename = "miRNA.GSE115016/Target.logFC.pdf",height = 3,width = 6)


#### miRNA - GSE140370 -> Hiseq2500 -> DEseq2 ####
dir.create("miRNA.GSE140370")

Dta <- data.frame(nrow = 2586)
colname <- c(paste("HCC",1:3,sep = ""),paste("Ctrl",1:3,sep = ""))
library(readxl)
for (i in 1:6) {
  Dta <- read_xlsx("miRNA.GSE140370/GSE140370_Small_RNA_Sequences_non-normalized.xlsx",sheet = i) %>% arrange(miRNA) %>%
    dplyr::select(RawReads) %>%
    cbind.data.frame(Dta,.)
}

test <- read_xlsx("miRNA.GSE140370/GSE140370_Small_RNA_Sequences_non-normalized.xlsx",sheet = i) %>% arrange(miRNA)
colnames(Dta) = c("nrow",colname)
Dta2 <- Dta %>% dplyr::select(-nrow) %>% mutate(ID = test$miRNA) %>% remove_rownames() %>% column_to_rownames("ID")

DESEQ2_Function <- function(CountData=CountData,colData=colData,Filename="Test",Type="mRNA",FoldChange = 1, P.adj = 0.05){ #colData => Group
  suppressWarnings(library(DESeq2))
  suppressWarnings(library(ggthemes))
  dds <- DESeqDataSetFromMatrix(countData = CountData, colData = colData, design = ~Group) 
  dds <- DESeq(dds)
  res <- results(dds,contrast=c("Group","Tumor","Ctrl"))
  res <- res[order(res$pvalue),] %>% data.frame() %>% na.omit() %>%
    mutate(Threshold = ifelse(padj <= P.adj & log2FoldChange >=FoldChange, "Up",ifelse(padj <= P.adj & log2FoldChange <=-FoldChange,"Down","No")))
  
  save(res,file = paste(Filename,"_","DE_",Type,".Rdata",sep = ''))
  write.csv(res,file = paste(Filename,"_","DE_",Type,".csv",sep = ''))
  res <- na.omit(res)
  
  this_tile <- paste("Up:",sum(res$Threshold == "Up"),"\t","    ","Down:",sum(res$Threshold == "Down"),sep = '')
  
  p<- ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=Threshold))+  
    geom_point(alpha=0.8, size=2) +
    scale_colour_manual(values = c("Up"= "#FC4E07", "Down"="#00AFBB",  "No"= "grey60")) +
    #xlim(c(-8, 8)) +
    geom_vline(xintercept=c(-FoldChange,FoldChange),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(P.adj),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",y="-log10 (padj)") +
    theme_few()+
    ggtitle(this_tile)+
    theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    #geom_text_repel(data = subset(hs_data, hs_data$Pvalue < 0.000001 & abs(hs_data$log2FC) >= 3),aes(label = ID),size = 3,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
    #geom_text_repel(aes(x=log2FoldChange,y=-log10(padj),color=Sig,label=Label),size=3)
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(p,filename = paste(Filename,"_","DE_",Type,".pdf",sep = ''),width = 4,height = 4)
}

ColData <- data.frame(Sample=colnames(Dta2),Group = rep(c("Tumor","Ctrl"),c(3,3))) %>% remove_rownames() %>% column_to_rownames("Sample")
DESEQ2_Function(CountData = Dta2,colData = ColData,Filename="GSE140370",Type="miRNA")

Dta3 <- read.csv("miRNA.GSE140370/GSE140370_DE_miRNA.csv")
Target.gene <- Dta3 %>% filter(X %in% Inter.miRNA) %>% 
  mutate(X = factor(.$X,levels=Inter.miRNA))

p<-ggplot(Target.gene,aes(x=X,y=log2FoldChange,fill=X))+
  geom_col() + 
  geom_text(aes(x=X,y=log2FoldChange,label = sprintf("%2.2e", padj)))+
  theme_few() + 
  coord_flip() + labs(x="")+
  theme(legend.position = "none")

ggsave(p,filename = "miRNA.GSE140370/GSE140370.Target.logFC.pdf",height = 3,width = 6)


#### miRNA - GSE22058 -> CEL -> GEO2R => didn't detect ####
dir.create("miRNA.GSE22058")

Dat1 <- read.table("miRNA.GSE22058/GSE22058.top.table.tsv",header = T,sep = "\t")

Target.gene <- Dat1 %>% filter(miRNA_ID %in% Inter.miRNA) %>% 
  mutate(miRNA_ID = factor(.$miRNA_ID,levels=Inter.miRNA))

#### miRNA - GSE10694 -> GPR -> GEO2R => didn't detect ####
dir.create("miRNA.GSE10694")

Dat1 <- read.table("miRNA.GSE10694/GSE10694.top.table.tsv",header = T,sep = "\t")

Target.gene <- Dat1 %>% filter(miRNA_ID %in% Inter.miRNA) %>% 
  mutate(miRNA_ID = factor(.$miRNA_ID,levels=Inter.miRNA))


#### miRNA - GSE63046 -> Hiseq2500 (RPM) -> Limma ####
dir.create("miRNA.GSE63046")

Dat1 <- read.table("miRNA.GSE63046/GSE63046_liver_RPM_miRs.txt",header = T,row.names = 1,sep = "\t",check.names = F)

group_list=ifelse(str_detect(colnames(Dat1),"N"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Dat1=normalizeBetweenArrays(log2(Dat1+1))

design=model.matrix(~group_list)
fit=lmFit(Dat1,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

write.csv(deg,file = "miRNA.GSE63046/GSE63046.de.csv",row.names = F)

Target.gene <- deg %>% rownames_to_column() %>% filter(rowname %in% Inter.miRNA) %>% 
  mutate(rowname = factor(.$rowname,levels=Inter.miRNA))

p<-ggplot(Target.gene,aes(x=rowname,y=logFC,fill=rowname))+
  geom_col() + 
  geom_text(aes(x=rowname,y=logFC,label = sprintf("%2.2e", adj.P.Val)))+
  theme_few() + 
  coord_flip() + labs(x="")+
  theme(legend.position = "none")

ggsave(p,filename = "miRNA.GSE63046/GSE63046.Target.logFC.pdf",height = 3,width = 6)



###### lncRNA Validate ######
Inter.lncRNA.GeneCard <- read.csv("Inter.lncRNA.GeneCard.csv")
setwd("K:/TCGA/Anlysis/LIHC/miRNA.Validate")

library(AnnoProbe)
library(GEOquery)
library(ggplot2) 
library(ggstatsplot) 
library(reshape2)
library(patchwork)
library(limma)
library(tidyverse)
library(ggthemes)

#### lncRNA - GSE125469 - FPKM -> limma ########
dir.create("lncRNA.GSE125469")
exp <- read_xlsx("lncRNA.GSE125469/GSE125469_Expression_Gene.xlsx",col_names = TRUE)
exp <- exp %>% dplyr::select(c("Track_id",contains("-"))) %>% remove_rownames() %>% column_to_rownames("Track_id")
exp2 <- sapply(exp, function(x) {as.numeric(x)})
rownames(exp2) = rownames(exp)
exp2=normalizeBetweenArrays(log2(exp2+1))

group_list=ifelse(str_detect(colnames(exp2),"C1"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

design=model.matrix(~group_list)
fit=lmFit(exp2,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

write.csv(deg,file = "lncRNA.GSE125469/GSE125469.de.csv",row.names = T)

hsa.deg <- deg %>% mutate(lncRNA.ENSEMBLE.ID=(rownames(.) %>% str_remove_all("\\..*"))) %>% merge.data.frame(Inter.lncRNA.GeneCard,.,by="lncRNA.ENSEMBLE.ID")
write.csv(hsa.deg,file = "lncRNA.GSE125469/GSE125469.InterlncRNA.de.csv",row.names = T)

#### lncRNA - GSE54238 -> Pair -> limma -> LOC153684 Up ####

dir.create("lncRNA.GSE54238")

#   Differential expression analysis with limma
gse <- 'GSE54238'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE54238')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE54238/GSE54238.Rdata")

load("lncRNA.GSE54238/GSE54238.Rdata")
pd <- pd %>% filter(str_detect(title,"NL") | str_detect(title,"HCC"))

group_list=ifelse(str_detect(pd$title,"NL"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

exp = exp[,match(rownames(pd),colnames(exp))]

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE54238"))
#Anno_data <- Anno_data %>% filter(TYPE == "noncoding")
Anno_data <- Anno_data[,1:6]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE54238/GSE54238.de.csv",row.names = F)

Target.gene <- hsa.deg %>% filter(GENESYMBOL %in% Inter.lncRNA.GeneCard$lncRNA.ID)

#### lncRNA - GSE27462 -> -> limma  ####
dir.create("lncRNA.GSE27462")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE27462'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE27462')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE27462/GSE27462.Rdata")

load("lncRNA.GSE27462/GSE27462.Rdata")
group_list=ifelse(str_detect(pd$source_name_ch1,"non-tumor"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE27462"))
Anno_data[1:8,1:12]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE27462/GSE27462.de.csv",row.names = F)


#### lncRNA - GSE49713 -> -> limma ####
dir.create("lncRNA.GSE49713")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE49713'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE49713')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE49713/GSE49713.Rdata")

load("lncRNA.GSE49713/GSE49713.Rdata")
group_list=ifelse(str_detect(pd$source_name_ch1,"non-tumor"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE49713"))
Anno_data[1:8,1:12]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE49713/GSE49713.de.csv",row.names = F)


#### lncRNA - GSE58043 -> -> limma ####
dir.create("lncRNA.GSE58043")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE58043'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE58043')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE58043/GSE58043.Rdata")

load("lncRNA.GSE58043/GSE58043.Rdata")
group_list=ifelse(str_detect(pd$source_name_ch1,"non-tumor"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE58043"))
Anno_data[1:8,1:8]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE58043/GSE58043.de.csv",row.names = F)
#### lncRNA - GSE64633 -> -> limma -> Useless ####
dir.create("lncRNA.GSE64633")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE64633'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE64633')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE64633/GSE64633.Rdata")

load("lncRNA.GSE64633/GSE64633.Rdata")
group_list=ifelse(str_detect(pd$source_name_ch1,"normal"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE64633"))
Anno_data[1:8,1:8]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE64633/GSE64633.de.csv",row.names = F)
#### lncRNA - GSE55191 -> -> limma -> Useless ####
dir.create("lncRNA.GSE55191")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE55191'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE55191')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE55191/GSE55191.Rdata")

load("lncRNA.GSE55191/GSE55191.Rdata")
group_list=ifelse(str_detect(pd$source_name_ch1,"normal"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE55191"))
Anno_data[1:8,1:8]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE55191/GSE55191.de.csv",row.names = F)
#### lncRNA - GSE70880 -> -> limma -> Useless ####
dir.create("lncRNA.GSE70880")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE70880'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE70880')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE70880/GSE70880.Rdata")

load("lncRNA.GSE70880/GSE70880.Rdata")

pd <- pd %>% filter(str_detect(source_name_ch1,"Liver"))


group_list=ifelse(str_detect(pd$source_name_ch1,"Normal"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE70880"))
Anno_data[1:8,1:4]

exp <- exp[,pd$geo_accession]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE70880/GSE70880.de.csv",row.names = F)

#### lncRNA - GSE67260 -> -> limma -> Useless ####
dir.create("lncRNA.GSE67260")

#   Differential expression analysis with limma
# Ëé∑ÂèñË°®ËææÈáèÁü©???
gse <- 'GSE67260'
eSet <- getGEO(gse, 
               destdir = 'lncRNA.GSE67260')

#(1)ÊèêÂèñË°®ËææÁü©Èòµexp
exp <- exprs(eSet[[1]])
exp=normalizeBetweenArrays(log2(exp+1))

#(2)ÊèêÂèñ‰∏¥Â∫ä‰ø°ÊÅØ
pd <- pData(eSet[[1]])

#(3)Ë∞ÉÊï¥pdÁöÑË°åÂêçÈ°∫Â∫è‰∏éexpÂàóÂêçÂÆåÂÖ®‰∏Ä???
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

#(4)ÊèêÂèñËäØÁâáÂπ≥Âè∞ÁºñÂè∑
gpl <- eSet[[1]]@annotation
save(gse,pd,exp,gpl,file = "lncRNA.GSE67260/GSE67260.Rdata")

load("lncRNA.GSE67260/GSE67260.Rdata")
group_list=ifelse(str_detect(pd$source_name_ch1,"normal"),"Ctrl","HCC")
group_list = factor(group_list,
                    levels = c("Ctrl","HCC"))

Anno_data <- Table(getGEO(gpl,destdir = "lncRNA.GSE67260"))
Anno_data[1:8,1:3]

design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
hsa.deg <- deg %>% mutate(probe_id=rownames(.)) %>% merge.data.frame(.,Anno_data,by.x = "probe_id",by.y="ID")

write.csv(hsa.deg,file = "lncRNA.GSE67260/GSE67260.de.csv",row.names = F)





###### mRNA Validate #######
setwd("K:/TCGA/Anlysis/LIHC/miRNA.Validate")
library(AnnoProbe)
library(GEOquery)
library(ggplot2) 
library(ggstatsplot) 
library(reshape2)
library(patchwork)
library(limma)
library(tidyverse)
library(ggthemes) #http://lifeome.net/database/hccdb/home.html
#### mRNA - GSE144269 -> 70+70 -> DEseq2 => FoldChange = 2, P.adj = 0.01 ####
dir.create("mRNA.GSE144269")

DESEQ2_Function <- function(CountData=CountData,colData=colData,Filename="Test",Type="mRNA",FoldChange = 1, P.adj = 0.05){ #colData => Group
  suppressWarnings(library(DESeq2))
  suppressWarnings(library(ggthemes))
  dds <- DESeqDataSetFromMatrix(countData = CountData, colData = colData, design = ~Group) 
  dds <- DESeq(dds)
  res <- results(dds,contrast=c("Group","Tumor","Ctrl"))
  res <- res[order(res$pvalue),] %>% data.frame() %>% 
    mutate(Threshold = ifelse(padj <= P.adj & log2FoldChange >=FoldChange, "Up",ifelse(padj <= P.adj & log2FoldChange <=-FoldChange,"Down","No")))
  
  save(res,file = paste(Filename,"_","DE_",Type,".pdf",sep = ''))
  write.csv(res,file = paste(Filename,"_","DE_",Type,".csv",sep = ''))
  res <- na.omit(res)
  
  this_tile <- paste("Up:",sum(res$Threshold == "Up"),"\t","    ","Down:",sum(res$Threshold == "Down"),sep = '')
  
  p<- ggplot(res,aes(x=log2FoldChange,y=-log10(padj),color=Threshold))+  
    geom_point(alpha=0.8, size=2) +
    scale_colour_manual(values = c("Up"= "#FC4E07", "Down"="#00AFBB",  "No"= "grey60")) +
    #xlim(c(-8, 8)) +
    geom_vline(xintercept=c(-FoldChange,FoldChange),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(P.adj),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",y="-log10 (padj)") +
    theme_few()+
    ggtitle(this_tile)+
    theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    #geom_text_repel(data = subset(hs_data, hs_data$Pvalue < 0.000001 & abs(hs_data$log2FC) >= 3),aes(label = ID),size = 3,box.padding = unit(0.5, "lines"),point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
    #geom_text_repel(aes(x=log2FoldChange,y=-log10(padj),color=Sig,label=Label),size=3)
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(p,filename = paste(Filename,"_","DE_",Type,".pdf",sep = ''),width = 4,height = 4)
}

Data.set <- read.table("mRNA.GSE144269/GSE144269_RSEM_GeneCounts.txt",header = T,sep = "\t",row.names = 1)
Group = if_else(str_detect(colnames(Data.set),"A"),"Tumor","Ctrl")
colData <- data.frame(Group = Group) %>% mutate(Group = factor(.$Group,levels = c("Ctrl","Tumor")))
rownames(colData) = colnames(Data.set)

DESEQ2_Function(Data.set,colData,Filename="mRNA.GSE144269/GSE144269.DEG",Type="mRNA",FoldChange = 2, P.adj = 0.01)

#Res <- read.csv("mRNA.GSE144269/GSE144269.DEG_DE_mRNA.csv")

#### mRNA - GSE148355 -> 15+54 -> DEseq2 => FoldChange = 2, P.adj = 0.01 ####
dir.create("mRNA.GSE148355")

Data.set <- read.table("mRNA.GSE148355/GSE148355_1.catholic.readCount.RNAseqId.txt",header = T,row.names = 1,sep = "\t",check.names = F)
rownames(Data.set) = paste(rownames(Data.set)%>%str_remove_all("\\..*"),Data.set$Symbol,sep = "_")
Data.set <- Data.set %>% dplyr::select(contains("Cat_N") | contains("Cat_TG"))
Group = if_else(str_detect(colnames(Data.set),"Cat_TG"),"Tumor","Ctrl")
colData <- data.frame(Group = Group) %>% mutate(Group = factor(.$Group,levels = c("Ctrl","Tumor")))
rownames(colData) = colnames(Data.set)
DESEQ2_Function(Data.set,colData,Filename="mRNA.GSE148355/GSE148355.DEG",Type="mRNA",FoldChange = 2, P.adj = 0.01)

#### mRNA - GSE138485 -> 32+32 -> limma + voom => FoldChange = 2, P.adj = 0.01  ####
dir.create("mRNA.GSE138485")

library(R.utils)
library(plyr)

files = list.files("mRNA.GSE138485/",pattern = "[gz]",recursive = F,full.names = T)
for (file in files) {
  gunzip(file, remove = T)
}

library(readxl)
files = list.files("mRNA.GSE138485/",pattern = "[xls$]",recursive = F,full.names = T)
files <- files[str_detect(files,"xls")]
Group <- c()
Data.set <- data.frame()
colnames.V = c()
for (file in files) {
  name = (str_split(file,"/")[[1]])[2] %>% str_remove_all("_.*")
  colnames.V <- c(colnames.V,name)
  if(str_detect(file,"Non_Tumor")){
    Group <- c(Group,"Ctrl")
  }else{
    Group <- c(Group,"Tumor")
  }
  middata <- read_xls(file)
  middata <- middata %>% mutate(rowname = str_remove_all(GeneSymbol,"\\..*")) %>% mutate(rowname = paste(rowname,GeneName,sep = "_")) %>%
    remove_rownames() %>% column_to_rownames("rowname") %>% dplyr::select(TPM) %>% plyr::rename(c("TPM" = name))
  Data.set <- rbind.data.frame(Data.set,t.data.frame(middata))
}

Data.set <- t.data.frame(Data.set)
write.csv(Data.set,file = "mRNA.GSE138485/GSE138485.TPM.csv",row.names = T)

Data.set <- read.csv("mRNA.GSE138485/GSE138485.TPM.csv",row.names = 1)

ColData <- data.frame(Group = Group) %>% mutate(Group = factor(Group,levels = c("Ctrl","Tumor")))
rownames(ColData) = colnames(Data.set)

design <- model.matrix(~factor(Group))
colnames(design)=levels(factor(Group))
rownames(design)=colnames(Data.set)
## voom
Data.set2 <- voom(Data.set,design,normalize="quantile")

fit <- lmFit(Data.set2,design)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom = na.omit(tempOutput)
write.csv(DEG_voom,file = "mRNA.GSE138485/GSE138485.DEG.csv")

data <- read.csv("mRNA.GSE138485/GSE138485.DEG.csv") 
data$EMS = data$X %>% str_remove_all("_.*")
write.csv(data,file = "mRNA.GSE138485/GSE138485.DEG.csv")


#### mRNA - GSE94660  -> 21+21 -> limma + voom => FoldChange = 2, P.adj = 0.01 ####
dir.create("mRNA.GSE94660")

Data.set <- read.csv("mRNA.GSE94660/GSE94660_RPKM_normalized.txt",header = T,sep = "\t")

Group <- if_else(str_detect(colnames(Data.set),"A"),"Tumor","Ctrl")
Group <- factor(Group,levels = c("Ctrl","Tumor"))

design <- model.matrix(~factor(Group))
colnames(design)=levels(factor(Group))
rownames(design)=colnames(Data.set)
cont.matrix <- makeContrasts(contrasts = c('Tumor-Ctrl'), levels = design)
## voom
Data.set2 <- voom(Data.set,design,normalize="quantile")

fit <- lmFit(Data.set2,design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom2 = na.omit(tempOutput)
write.csv(DEG_voom2,file = "mRNA.GSE94660/GSE94660.DEG.csv")


#### mRNA - GSE25599  -> 10+10 -> limma + voom => FoldChange = 2, P.adj = 0.01 ####
dir.create("mRNA.GSE25599")

Data.set <- read.table("mRNA.GSE25599/GSE25599_rpkm_gene_all.txt",header = T,sep = "\t") %>% na.omit() %>% remove_rownames() %>% column_to_rownames("id") 

Group <- if_else(str_detect(colnames(Data.set),"c"),"Tumor","Ctrl")
Group <- factor(Group,levels = c("Ctrl","Tumor"))

design <- model.matrix(~factor(Group))
colnames(design)=levels(factor(Group))
rownames(design)=colnames(Data.set)
cont.matrix <- makeContrasts(contrasts = c('Tumor-Ctrl'), levels = design)
## voom
Data.set2 <- voom(Data.set,design,normalize="quantile")

fit <- lmFit(Data.set2,design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom2 = na.omit(tempOutput) %>% data.frame() %>% rownames_to_column() %>% merge.data.frame(.,mid.Gene.Type,by.x = "rowname",by.y = "Gene.stable.ID")
write.csv(DEG_voom2,file = "mRNA.GSE25599/GSE25599.DEG.csv")


#### mRNA - GSE55758  -> 8+8 -> DEseq2 => FoldChange = 2, P.adj = 0.01 ####
dir.create("mRNA.GSE55758")

files <- list.files("mRNA.GSE55758/",pattern = "[gz]",full.names = T)
for (file in files) {
  gunzip(file,remove = T)
}

files <- list.files("mRNA.GSE55758/",pattern = '[txt]',full.names = T)
files <- files[str_detect(files,"txt")]
Data.set <- data.frame()
for (file in files) {
  middata <- read.table(file,header = T,sep = "\t",check.names = F)
  middata <- middata[,c(1,3,4)] %>% data.frame() %>% arrange(geneID) %>% remove_rownames() %>% column_to_rownames("geneID") %>% t() %>% data.frame(check.names = F )
  #rownames(middata) <- paste(c("Ctrl","Tumor"),i,sep = "")
  Data.set <- plyr::rbind.fill(Data.set,middata)
}

Data.set <- Data.set %>% t() %>% data.frame(check.names = F)
colnames(Data.set) = paste(c("Ctrl","Tumor"),rep(1:8,each = 2),sep = "")
write.csv(Data.set,file = "mRNA.GSE55758/GSE55758.Count.csv")

Data.set <- read.csv("mRNA.GSE55758/GSE55758.Count.csv",row.names = 1)
Data.set[is.na(Data.set)] = 0

colData <- data.frame(Group=if_else(str_detect(colnames(Data.set),"Ctrl"),"Ctrl","Tumor")) %>%
  mutate(Group = factor(.$Group,levels = c("Ctrl","Tumor")))
rownames(colData) <- colnames(Data.set)

DESEQ2_Function(Data.set,colData,Filename="mRNA.GSE55758/GSE55758.DEG",Type="mRNA",FoldChange = 2, P.adj = 0.01)

Data.set <- read.csv("mRNA.GSE55758/GSE55758.DEG_DE_mRNA.csv",check.names = F,row.names = 1)
mid.Gene.Type2 <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name,HGNC.symbol,NCBI.gene..formerly.Entrezgene..ID) %>% unique()

Data.set2 <- Data.set %>% rownames_to_column() %>% merge.data.frame(.,mid.Gene.Type2,by.x = "rowname",by.y = "NCBI.gene..formerly.Entrezgene..ID")
write.csv(Data.set2,file = "mRNA.GSE55758/GSE55758.DEG.Final.csv",row.names = F)


#### mRNA - GSE124535 -> 35+35 -> limma + voom => FoldChange = 2, P.adj = 0.01 ####
dir.create("mRNA.GSE124535")

Data.set <- read.table("mRNA.GSE124535/GSE124535_HCC.RNA-seq.35.samples.fpkm.txt",header = T,row.names = 1,sep = "\t")
rownames(Data.set) <- paste(rownames(Data.set),Data.set$gene_symbol,sep = "_")
Data.set <- Data.set %>% dplyr::select(-gene_symbol,-Chr,-Biotype)

Group <- if_else(str_detect(colnames(Data.set),"P"),"Ctrl","Tumor")

Group <- factor(Group,levels = c("Ctrl","Tumor"))

design <- model.matrix(~factor(Group))
colnames(design)=levels(factor(Group))
rownames(design)=colnames(Data.set)
cont.matrix <- makeContrasts(contrasts = c('Tumor-Ctrl'), levels = design)
## voom
Data.set2 <- voom(Data.set,design,normalize="quantile")

fit <- lmFit(Data.set2,design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
DEG_voom2 = na.omit(tempOutput) %>% mutate(Ensm = str_remove_all(rownames(.),"_.*"))

write.csv(DEG_voom2,file = "mRNA.GSE124535/GSE124535.DEG.csv")

#### mRNA - GSE77314  -> 50+50 -> limma + voom => FoldChange = 2, P.adj = 0.01 => ÁªìÊûú‰∏çÂ•Ω ####
dir.create("mRNA.GSE77314")

library(readxl)
Data.set <- read.csv("mRNA.GSE77314/GSE77314_expression.csv",check.names = F) 
Gene.Name <- Data.set$gene_short_name
Data.set <- Data.set %>% dplyr::select(-gene_short_name)

Group <- if_else(str_detect(colnames(Data.set),"N"),"Ctrl","Tumor")

Group <- factor(Group,levels = c("Ctrl","Tumor"))

design <- model.matrix(~factor(Group))
colnames(design)=levels(factor(Group))
rownames(design)=colnames(Data.set)
cont.matrix <- makeContrasts(contrasts = c('Tumor-Ctrl'), levels = design)
## voom
Data.set2 <- voom(Data.set,design,normalize="quantile")

fit <- lmFit(Data.set2,design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit)
tempOutput = topTable(fit2, coef=2, n=Inf)
tempOutput$Gene.Name <- Gene.Name
DEG_voom2 = na.omit(tempOutput)
write.csv(DEG_voom2,file = "mRNA.GSE77314/GSE77314.DEG.csv")

#### Ingegrated mRNA Datasets ####
setwd("K:/TCGA/Anlysis/LIHC/miRNA.Validate")
Gene <- c("LRRC1","MAP2","MCM2","RMI2","CKAP2L","CPEB3","ESR1","MCM3AP-AS1",Inter.Five.Gene$Gene.name)
Gene.ID <- GeneType %>% dplyr::select(Gene.stable.ID,HGNC.symbol) %>% unique() %>%
  filter(HGNC.symbol %in% Gene)

Gene.ID <- Gene.ID[2:13,] #ENSG00000271672
setwd("K:/TCGA/Anlysis/LIHC/miRNA.Validate")
for (GSE in c("GSE144269","GSE148355","GSE138485","GSE94660","GSE25599","GSE55758","GSE124535","GSE77314")) {}

GSE144269.DEG <- read.csv("mRNA.GSE144269/GSE144269.DEG_DE_mRNA.csv") %>% 
  mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange >=1, "Up",ifelse(padj <= 0.01 & log2FoldChange <=-1,"Down","No"))) %>%
  mutate(ENS = str_remove_all(X,"|.*")) %>% mutate(ENS = str_remove_all(ENS,"\\..*")) %>% 
  dplyr::select(ENS,Threshold) %>%
  filter(ENS %in% Gene.ID$Gene.stable.ID) %>% mutate(Dataset = "GSE144269")

GSE148355.DEG <- read.csv("mRNA.GSE148355/GSE148355.DEG_DE_mRNA.csv") %>% 
  mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange >=1, "Up",ifelse(padj <= 0.01 & log2FoldChange <=-1,"Down","No"))) %>%
  mutate(ENS = str_remove_all(X,"_.*")) %>% mutate(ENS = str_remove_all(ENS,"\\..*")) %>% 
  dplyr::select(ENS,Threshold) %>%
  filter(ENS %in% Gene.ID$Gene.stable.ID) %>% mutate(Dataset = "GSE148355")

GSE138485.DEG <- read.csv("mRNA.GSE138485/GSE138485.DEG.csv")%>% 
  mutate(Threshold = ifelse(adj.P.Val <= 0.01 & logFC >=1, "Up",ifelse(adj.P.Val <= 0.01 & logFC <=-1,"Down","No"))) %>%
  dplyr::select(X,Threshold) %>% dplyr::rename(ENS=1) %>%
  filter(ENS %in% Gene.ID$Gene.stable.ID) %>% mutate(Dataset = "GSE138485")

GSE94660.DEG <- read.csv("mRNA.GSE94660/GSE94660.DEG.csv")%>% 
  mutate(Threshold = ifelse(adj.P.Val <= 0.01 & logFC >=1, "Up",ifelse(adj.P.Val <= 0.01 & logFC <=-1,"Down","No"))) %>%
  dplyr::select(X,Threshold) %>% dplyr::rename(ENS=1) %>%
  filter(ENS %in% Gene.ID$HGNC.symbol) %>% mutate(Dataset = "GSE94660") %>%
  merge(.,Gene.ID,by.x="ENS",by.y="HGNC.symbol") %>% dplyr::select(Gene.stable.ID,Threshold,Dataset) %>%
  dplyr::rename(ENS=1)

GSE25599.DEG <- read.csv("mRNA.GSE25599/GSE25599.DEG.csv") %>% 
  mutate(Threshold = ifelse(adj.P.Val <= 0.01 & logFC >=1, "Up",ifelse(adj.P.Val <= 0.01 & logFC <=-1,"Down","No"))) %>%
  dplyr::select(rowname,Threshold) %>% dplyr::rename(ENS=1) %>%
  filter(ENS %in% Gene.ID$Gene.stable.ID) %>% mutate(Dataset = "GSE25599")

GSE55758.DEG <- read.csv("mRNA.GSE55758/GSE55758.DEG.Final.csv") %>% 
  mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange >=1, "Up",ifelse(padj <= 0.01 & log2FoldChange <=-1,"Down","No"))) %>%
  mutate(ENS = str_remove_all(Gene.stable.ID,"_.*")) %>% mutate(ENS = str_remove_all(ENS,"\\..*")) %>% 
  dplyr::select(ENS,Threshold) %>%
  filter(ENS %in% Gene.ID$Gene.stable.ID) %>% mutate(Dataset = "GSE55758")

GSE124535.DEG <- read.csv("mRNA.GSE124535/GSE124535.DEG.csv") %>% 
  mutate(Threshold = ifelse(adj.P.Val <= 0.01 & logFC >=1, "Up",ifelse(adj.P.Val <= 0.01 & logFC <=-1,"Down","No"))) %>%
  dplyr::select(Ensm,Threshold) %>% dplyr::rename(ENS=1) %>%
  filter(ENS %in% Gene.ID$Gene.stable.ID) %>% mutate(Dataset = "GSE124535")

Cor.data <- rbind.data.frame(GSE144269.DEG,GSE148355.DEG) %>% 
  rbind.data.frame(GSE138485.DEG) %>%
  rbind.data.frame(GSE94660.DEG) %>%
  rbind.data.frame(GSE25599.DEG) %>%
  rbind.data.frame(GSE55758.DEG) %>%
  rbind.data.frame(GSE124535.DEG)

library(ComplexHeatmap)
matrix <- Cor.data %>% data.frame(check.names = F) %>% 
  mutate(Representative=if_else(Threshold == "Up",1,if_else(Threshold == "Down",-1,0))) %>%
  merge(.,Gene.ID,by.x = "ENS",by.y="Gene.stable.ID") %>% dplyr::select(Dataset,HGNC.symbol,Threshold) %>%
  spread(Dataset,Threshold,fill = NA) %>% remove_rownames() %>%
  column_to_rownames("HGNC.symbol")

colors <- structure(brewer.pal(9,"Set1")[c(1,2,9)],names=c('Up','Down','No')) 
Heatmap(matrix, col = colors,border = "black",na_col = "grey60",rect_gp = gpar(col = "white", lwd = 2),heatmap_legend_param = list(title = "Significance"))
graph2pdf(file="mRNA.11genes.Sig.Datasets.pdf",height=4,width=6)


######################################################
########### CBX2 CEP55 corrleated gene ####################
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

Pheno <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID=="TCGA-LIHC") %>% filter(Sample.Type == "Primary Tumor")

Cor.data <- mRNA.Exp.log2.TPM %>% dplyr::select(Pheno$Sample.ID) %>% 
  rownames_to_column("ENSEMBL")
#GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
#GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180
CBX2.CEP55.Exp <- Cor.data %>% 
  filter(ENSEMBL == "ENSG00000173894" | ENSEMBL == "ENSG00000138180") %>%
  column_to_rownames("ENSEMBL")

Other.Exp <- Cor.data %>% 
  filter(ENSEMBL != "ENSG00000173894") %>%
  filter(ENSEMBL != "ENSG00000138180") %>%
  column_to_rownames("ENSEMBL")

cor <- data.frame()
for (i in 1:nrow(CBX2.CEP55.Exp)) {
  for (j in 1:nrow(Other.Exp)) {
    test<-cor.test(as.numeric(as.character(CBX2.CEP55.Exp[i,])),
                   as.numeric(as.character(Other.Exp[j,])),
                   method = "pearson",
                   exact = F,
                   alternative = "two.sided")
    cor <- data.frame(Gene1=rownames(CBX2.CEP55.Exp)[i],
                      Gene2=rownames(Other.Exp)[j],
                      Rho=test$estimate,
                      Pvalue=test$p.value) %>%
      rbind.data.frame(cor)
  }
}
write.table(cor,file="TCGA-LIHC.CBX2CEP55.Pearson.txt",sep="\t",row.names=F,quote=F)

######## GSEA ######
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
###### TCGA LIHC ####
setwd("K:/TCGA/Anlysis/LIHC")
Pearson.data <- read.csv("TCGA-LIHC.CBX2CEP55.Pearson.txt",sep = "\t")
Pearson.data <- na.omit(Pearson.data)

KEGG_df <- read.gmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO_df <- read.gmt("../../msigdb/c5.go.v7.5.1.symbols.gmt")
Hallmark <- read.gmt("../../msigdb/h.all.v7.5.1.symbols.gmt")

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()

#### CBX2 ####
Pearson.data.CBX2 <- Pearson.data %>%
  filter(Gene1=="ENSG00000173894") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  dplyr::select(Rho,Gene.name) %>%
  na.omit()

ge = Pearson.data.CBX2$Rho
names(ge) = Pearson.data.CBX2$Gene.name
ge = sort(ge,decreasing = T)

library(fgsea)
CBX.Pearson.GO <- GSEA(ge, TERM2GENE = GO_df,eps = 1e-100)
save(CBX.Pearson.GO,file = "CBX.Pearson.GO.Rdata")
CBX.Pearson.KEGG <- GSEA(ge, TERM2GENE = KEGG_df,eps = 1e-100)
save(CBX.Pearson.KEGG,file = "CBX.Pearson.KEGG.Rdata")
CBX.Pearson.HALLMARK <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100)
save(CBX.Pearson.HALLMARK,file = "CBX.Pearson.HALLMARK.Rdata")

setwd("K:/TCGA/Anlysis/LIHC")
load(file="CBX.Pearson.GO.Rdata")
load(file="CBX.Pearson.KEGG.Rdata")

TERM1 <- c("KEGG_CELL_CYCLE",
          "KEGG_DNA_REPLICATION",
          "KEGG_RIBOSOME",
          "KEGG_SPLICEOSOME")

# plot
library(GseaVis)
lapply(TERM1, function(x){
  gseaNb(object = CBX.Pearson.KEGG,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         subPlot = 2,
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')

load(file="CBX.Pearson.HALLMARK.Rdata")

TERM2 <- c("HALLMARK_E2F_TARGETS",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
          "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_MITOTIC_SPINDLE",
          "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_MYC_TARGETS_V2",
          "HALLMARK_APOPTOSIS",
          "HALLMARK_ANGIOGENESIS",
          "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

middata <- CBX.Pearson.HALLMARK@result %>%
  filter(ID %in% TERM2) %>%
  dplyr::select(ID,NES,p.adjust) %>%
  mutate(Project="TCGA") %>%
  mutate(SYMBOL="CBX2")



#### CEP55 ####
Pearson.data.CEP55 <- Pearson.data %>%
  filter(Gene1=="ENSG00000138180") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  dplyr::select(Rho,Gene.name) %>%
  na.omit()

ge = Pearson.data.CEP55$Rho
names(ge) = Pearson.data.CEP55$Gene.name
ge = sort(ge,decreasing = T)

library(fgsea)
CEP55.Pearson.GO <- GSEA(ge, TERM2GENE = GO_df,eps = 1e-100)
save(CEP55.Pearson.GO,file = "CEP55.Pearson.GO.Rdata")
CEP55.Pearson.KEGG <- GSEA(ge, TERM2GENE = KEGG_df,eps = 1e-100)
save(CEP55.Pearson.KEGG,file = "CEP55.Pearson.KEGG.Rdata")
CEP55.Pearson.HALLMARK <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100)
save(CEP55.Pearson.HALLMARK,file = "CEP55.Pearson.HALLMARK.Rdata")

load(file = "CEP55.Pearson.KEGG.Rdata")
TERM1 <- c("KEGG_CELL_CYCLE",
           "KEGG_DNA_REPLICATION",
           "KEGG_RIBOSOME",
           "KEGG_SPLICEOSOME")

# plot
library(GseaVis)
lapply(TERM1, function(x){
  gseaNb(object = CEP55.Pearson.KEGG,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         subPlot = 2,
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')


load(file="CEP55.Pearson.HALLMARK.Rdata")


TERM2 <- c("HALLMARK_E2F_TARGETS",
           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_G2M_CHECKPOINT",
           "HALLMARK_MITOTIC_SPINDLE",
           "HALLMARK_MYC_TARGETS_V1",
           "HALLMARK_MYC_TARGETS_V2",
           "HALLMARK_APOPTOSIS",
           "HALLMARK_ANGIOGENESIS",
           "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

middata <- CEP55.Pearson.HALLMARK@result %>%
  filter(ID %in% TERM2) %>%
  dplyr::select(ID,NES,p.adjust) %>%
  mutate(Project="TCGA") %>%
  mutate(SYMBOL="CEP55") %>%
  rbind.data.frame(middata)

###### ICGC LIRI JP ####
setwd("K:/TCGA/Anlysis/LIHC")
Pearson.data <- read.csv("TCGA-LIHC.CBX2CEP55.Pearson.txt",sep = "\t")
Pearson.data <- na.omit(Pearson.data)

KEGG_df <- read.gmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO_df <- read.gmt("../../msigdb/c5.go.v7.5.1.symbols.gmt")
Hallmark <- read.gmt("../../msigdb/h.all.v7.5.1.symbols.gmt")

#### CBX2 ####
Pearson.data.CBX2 <- Pearson.data %>%
  filter(Gene1=="CBX2") #%>%
  #merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  #na.omit()

ge = Pearson.data.CBX2$Rho
names(ge) = Pearson.data.CBX2$Gene2
ge = sort(ge,decreasing = T)

library(fgsea)
CBX2.Pearson.GO <- GSEA(ge, TERM2GENE = GO_df,eps = 1e-100)
save(CBX2.Pearson.GO,file = "CBX2.ICGC-LIRI-JP.Pearson.GO.Rdata")
CBX2.Pearson.KEGG <- GSEA(ge, TERM2GENE = KEGG_df,eps = 1e-100)
save(CBX2.Pearson.KEGG,file = "CBX2.ICGC-LIRI-JP.Pearson.KEGG.Rdata")
CBX2.Pearson.HALLMARK <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100)
save(CBX2.Pearson.HALLMARK,file = "CBX2.ICGC-LIRI-JP.Pearson.HALLMARK.Rdata")

load(file = "CBX2.ICGC-LIRI-JP.Pearson.KEGG.Rdata")
TERM1 <- c("KEGG_CELL_CYCLE",
           "KEGG_DNA_REPLICATION",
           #"KEGG_RIBOSOME",
           "KEGG_SPLICEOSOME")

# plot
library(GseaVis)
lapply(TERM1, function(x){
  gseaNb(object = CBX2.Pearson.KEGG,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         subPlot = 2,
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')

load(file = "CBX2.ICGC-LIRI-JP.Pearson.HALLMARK.Rdata")
middata <- CBX2.Pearson.HALLMARK@result %>%
  filter(ID %in% TERM2) %>%
  dplyr::select(ID,NES,p.adjust) %>%
  mutate(Project="ICGC") %>%
  mutate(SYMBOL="CBX2") %>%
  rbind.data.frame(middata)

#### CEP55 ####
Pearson.data.CEP55 <- Pearson.data %>%
  filter(Gene1=="CEP55") #%>%
  #merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  #na.omit()

ge = Pearson.data.CEP55$Rho
names(ge) = Pearson.data.CEP55$Gene2
ge = sort(ge,decreasing = T)

library(fgsea)
CEP55.Pearson.GO <- GSEA(ge, TERM2GENE = GO_df,eps = 1e-100)
save(CEP55.Pearson.GO,file = "CEP55.ICGC-LIRI-JP.Pearson.GO.Rdata")
CEP55.Pearson.KEGG <- GSEA(ge, TERM2GENE = KEGG_df,eps = 1e-100)
save(CEP55.Pearson.KEGG,file = "CEP55.ICGC-LIRI-JP.Pearson.KEGG.Rdata")
CEP55.Pearson.HALLMARK <- GSEA(ge, TERM2GENE = Hallmark,eps = 1e-100)
save(CEP55.Pearson.HALLMARK,file = "CEP55.ICGC-LIRI-JP.Pearson.HALLMARK.Rdata")

load(file = "CEP55.ICGC-LIRI-JP.Pearson.KEGG.Rdata")
TERM1 <- c("KEGG_CELL_CYCLE",
           "KEGG_DNA_REPLICATION",
           "KEGG_RIBOSOME",
           "KEGG_SPLICEOSOME")

# plot
library(GseaVis)
lapply(TERM1, function(x){
  gseaNb(object = CEP55.Pearson.KEGG,
         geneSetID = x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         subPlot = 2,
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')

load(file = "CEP55.ICGC-LIRI-JP.Pearson.HALLMARK.Rdata")
middata <- CEP55.Pearson.HALLMARK@result %>%
  filter(ID %in% TERM2) %>%
  dplyr::select(ID,NES,p.adjust) %>%
  mutate(Project="ICGC") %>%
  mutate(SYMBOL="CEP55") %>%
  rbind.data.frame(middata)
middata$ID <- str_remove_all(middata$ID,"HALLMARK_")
middata$ID <- factor(middata$ID,levels = rev(sort(unique(middata$ID))))
write.csv(middata,file = "SourceData.F4F.csv")
p <- ggplot(middata,aes(Project,ID,color=NES,size=-log10(p.adjust)))+
  geom_point()+
  facet_wrap(~SYMBOL)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colors = brewer.pal(11,"RdBu")[4:1])+
  theme(strip.text = element_text(size=15))

ggsave(p,filename = "Pearson.CBX2.CEP55.GSEA.HALLMARK.pdf",height = 3.5,width = 6)

######## |Rho| > 0.4 Overlapped ####
###### CBX2 ######
setwd("K:/TCGA/Anlysis/LIHC")
TCGA.Pearson <- read.csv("TCGA-LIHC.CBX2CEP55.Pearson.txt",sep = "\t")
ICGC.Pearson <- read.csv("ICGC-LIRI-JP.CBX2CEP55.Pearson.txt",sep = "\t")

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()

### 729 genes
TCGA.Pearson.CBX2 <- TCGA.Pearson %>% 
  filter(Gene1=="ENSG00000173894") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  filter(abs(Rho) >= 0.5)

TCGA.Pearson.CBX2.2 <- TCGA.Pearson %>% 
  filter(Gene1=="ENSG00000173894") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  na.omit() %>% filter(abs(Rho) >= 0.4)

TCGA.Pearson.CBX2.P <- TCGA.Pearson.CBX2.2 %>% filter(Rho >= 0.5)
TCGA.Pearson.CBX2.N <- TCGA.Pearson.CBX2.2 %>% filter(Rho <= -0.4)

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

Pheno <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID=="TCGA-LIHC") %>% filter(Sample.Type == "Primary Tumor")

Cor.data <- mRNA.Exp.log2.TPM %>% dplyr::select(Pheno$Sample.ID) %>% t() %>%
  data.frame()

library(ggstatsplot)
library(RColorBrewer)
p23 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000083807,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) SLC27A5")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p23,filename = "TCGA.Pearson.CBX2.N.SLC27A5.pdf",height = 3,width = 3)

p24 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000100652,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) SLC10A1")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p24,filename = "TCGA.Pearson.CBX2.N.SLC10A1.pdf",height = 3,width = 3)

p25 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000141485,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) SLC13A5")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p25,filename = "TCGA.Pearson.CBX2.N.SLC13A5.pdf",height = 3,width = 3)


p26 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000166840,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) GLYATL1")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p26,filename = "TCGA.Pearson.CBX2.N.GLYATL1.pdf",height = 3,width = 3)


p27 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000204653,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) ASPDH")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p27,filename = "TCGA.Pearson.CBX2.N.ASPDH.pdf",height = 3,width = 3)


p28 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000129596,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) CDO1")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p28,filename = "TCGA.Pearson.CBX2.N.CDO1.pdf",height = 3,width = 3)


library(clusterProfiler)
TCGA.Pearson.CBX2.P.KEGG <- enrichKEGG(
  na.omit(TCGA.Pearson.CBX2.P$NCBI.gene..formerly.Entrezgene..ID),
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
save(TCGA.Pearson.CBX2.P.KEGG,file = "TCGA.Pearson.CBX2.P.KEGG.Rdata")

TCGA.Pearson.CBX2.P.GO <- enrichGO(
  gene = TCGA.Pearson.CBX2.P$Gene2,
  keyType="ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
save(TCGA.Pearson.CBX2.P.GO,file = "TCGA.Pearson.CBX2.P.GO.Rdata")

### 78 genes
ICGC.Pearson.CBX2 <- ICGC.Pearson %>% na.omit() %>%
  filter(Gene1=="CBX2") %>% filter(abs(Rho) >= 0.5)
### Overlapped 55 genes
Inter.Genes <- merge(TCGA.Pearson.CBX2,ICGC.Pearson.CBX2,
                     by.x="Gene.name",by.y="Gene2")

save(Inter.Genes,file = "TCGA.ICGC.Pearson.CBX2.Inter.Rdata")

#### Corrleation examples ####

  #rownames_to_column("ENSEMBL")

Inter.Genes2 <- Inter.Genes %>% 
  dplyr::select(Gene.name,Rho.x,Rho.y,Pvalue.x,Pvalue.y) %>%
  magrittr::set_colnames(c("Gene","TCGA-LIHC.Rho","ICGC-LIRI-JP.Rho",
                           "TCGA-LIHC.Pvalue","ICGC-LIRI-JP.Pvalue"))
  
Inter.Genes.FigData <- data.frame(Gene=c(Inter.Genes2$Gene,Inter.Genes2$Gene),
                                  Rho=c(Inter.Genes2$`TCGA-LIHC.Rho`,
                                        Inter.Genes2$`ICGC-LIRI-JP.Rho`),
                                  Pvalue=c(Inter.Genes2$`TCGA-LIHC.Pvalue`,
                                           Inter.Genes2$`ICGC-LIRI-JP.Pvalue`)) %>%
  mutate(Dataset2=rep(c("TCGA-LIHC","ICGC-LIRI-JP"),c(55,55)))
  

P22 <- ggplot(Inter.Genes.FigData,aes(x=Gene,
                               y=Dataset2,size=-log10(Pvalue),
                               fill=Rho))+
  #geom_bar(stat="identity",position = position_dodge(0.9)) +
  geom_point(shape=21,color="black")+
  #coord_flip()+
  #ggthemes::theme_few() +
  ggthemes::theme_few()+
  theme(legend.position = "top",
        legend.text = element_text(size=8),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle = 90),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="",x="")+
  scale_fill_gradientn("Pearson'CC",colors = c("white", "#E41A1C"))
  
  #scale_color_gradientn(colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))
ggsave(P22,filename = "TCGA.ICGC.Pearson.CBX2.Inter.Dotplot.pdf",height = 3,width = 10)

###### CEP55 ######
TCGA.Pearson <- read.csv("TCGA-LIHC.CBX2CEP55.Pearson.txt",sep = "\t")
ICGC.Pearson <- read.csv("ICGC-LIRI-JP.CBX2CEP55.Pearson.txt",sep = "\t")

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.type,Gene.name,NCBI.gene..formerly.Entrezgene..ID) %>% unique()

### 3949 genes
TCGA.Pearson.CEP55 <- TCGA.Pearson %>% 
  filter(Gene1=="ENSG00000138180") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  filter(abs(Rho) >= 0.5)

TCGA.Pearson.CEP55.2 <- TCGA.Pearson %>% 
  filter(Gene1=="ENSG00000138180") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  na.omit() %>% filter(abs(Rho) >= 0.4)

TCGA.Pearson.CEP55.P <- TCGA.Pearson.CEP55.2 %>% filter(Rho >= 0.5)
TCGA.Pearson.CEP55.N <- TCGA.Pearson.CEP55.2 %>% filter(Rho <= -0.6)

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

Pheno <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Project.ID=="TCGA-LIHC") %>% filter(Sample.Type == "Primary Tumor")

Cor.data <- mRNA.Exp.log2.TPM %>% dplyr::select(Pheno$Sample.ID) %>% t() %>%
  data.frame()

p30 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000138180,
  y=ENSG00000083807,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) SLC27A5")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p30,filename = "TCGA.Pearson.CEP55.N.SLC27A5.pdf",height = 3,width = 3)


p31 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000138180,
  y=ENSG00000166816,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) LDHD")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p31,filename = "TCGA.Pearson.CEP55.N.LDHD.pdf",height = 3,width = 3)

p32 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000138180,
  y=ENSG00000196616,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) ADH1B")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p32,filename = "TCGA.Pearson.CEP55.N.ADH1B.pdf",height = 3,width = 3)

#### Corrleation examples #####
TCGA.Pearson.CEP55 <- TCGA.Pearson %>% 
  filter(Gene1=="ENSG00000138180") %>%
  merge(GeneType.Unique,by.x="Gene2",by.y="Gene.stable.ID") %>%
  #dplyr::select(Rho,Gene.name) %>%
  filter(abs(Rho) >= 0.5)
ICGC.Pearson.CEP55 <- ICGC.Pearson %>% na.omit() %>%
  filter(Gene1=="CEP55") %>% filter(abs(Rho) >= 0.5)

Inter.Genes <- merge(TCGA.Pearson.CEP55,ICGC.Pearson.CEP55,
                     by.x="Gene.name",by.y="Gene2")

save(Inter.Genes,file = "TCGA.ICGC.Pearson.CEP55.Inter0.5.Rdata")

Inter.Genes2 <- Inter.Genes %>% 
  dplyr::select(Gene.name,Rho.x,Rho.y,Pvalue.x,Pvalue.y) %>%
  magrittr::set_colnames(c("Gene","TCGA-LIHC.Rho","ICGC-LIRI-JP.Rho",
                           "TCGA-LIHC.Pvalue","ICGC-LIRI-JP.Pvalue"))

Inter.Genes.FigData <- data.frame(Gene=c(Inter.Genes2$Gene,Inter.Genes2$Gene),
                                  Rho=c(Inter.Genes2$`TCGA-LIHC.Rho`,
                                        Inter.Genes2$`ICGC-LIRI-JP.Rho`),
                                  Pvalue=c(Inter.Genes2$`TCGA-LIHC.Pvalue`,
                                           Inter.Genes2$`ICGC-LIRI-JP.Pvalue`)) %>%
  mutate(Dataset2=rep(c("TCGA-LIHC","ICGC-LIRI-JP"),c(362,362))) %>%
  arrange(Gene,Dataset2,Rho)
Select.Gene <- Inter.Genes.FigData %>%
  filter(Rho >= 0.85 | Rho <= -0.6)

Inter.Genes.FigData2 <- Inter.Genes.FigData %>% filter(Gene %in% Select.Gene$Gene)


P35 <- ggplot(Inter.Genes.FigData2,aes(x=Gene,
                                      y=Dataset2,size=-log10(Pvalue),
                                      fill=Rho))+
  #geom_bar(stat="identity",position = position_dodge(0.9)) +
  geom_point(shape=21,color="black")+
  #coord_flip()+
  #ggthemes::theme_few() +
  ggthemes::theme_few()+
  theme(legend.position = "top",
        legend.text = element_text(size=8),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="",x="")+
  #scale_fill_gradientn("Pearson'CC",colors = c("white", "#E41A1C"))
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0) #,limits=c(-1,1)

#scale_color_gradientn(colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))
ggsave(P35,filename = "TCGA.ICGC.Pearson.CEP55.Inter.Dotplot.pdf",height = 3,width = 10)



######## CBX2 CEP55 AFP #####
p35 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000138180,
  y=ENSG00000081051,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) AFP")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p35,filename = "TCGA.Pearson.CEP55.AFP.pdf",height = 3,width = 3)


p36 <- ggscatterstats(
  data=Cor.data,
  x=ENSG00000173894,
  y=ENSG00000081051,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) AFP")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p36,filename = "TCGA.Pearson.CBX2.AFP.pdf",height = 3,width = 3)

ICGC <- read.csv("../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18_mRNA_level3.txt",sep = "\t")
Samples <- read.csv("../../../HepG2-CBX2/Validate/ICGC-LIRI-JP/HCCDB18.sample.txt",sep = "\t") %>%
  remove_rownames() %>% column_to_rownames("SAMPLE_ID") %>%
  t() %>% data.frame(check.names = F) %>% filter(TYPE == "HCC")
ICGC2 <- ICGC %>% dplyr::select(c("Symbol",rownames(Samples))) %>% remove_rownames() %>%
  column_to_rownames("Symbol") %>% t() %>% data.frame()

p38 <- ggscatterstats(
  data=ICGC2,
  x=CBX2,
  y=AFP,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) AFP")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p38,filename = "ICGC.Pearson.CBX2.AFP.pdf",height = 3,width = 3)

p39 <- ggscatterstats(
  data=ICGC2,
  x=CEP55,
  y=AFP,
  point.args = list(size = 3, alpha = 1, stroke = 0, na.rm = TRUE),
  smooth.line.args = list(size = 1.5, color = brewer.pal(8,"Set1")[2], method = "lm", formula = y ~ x,
                          na.rm = TRUE),
  type = "r",
)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) AFP")+
  ggpubr::stat_cor(method = "pearson",color=brewer.pal(8,"Set1")[1],label.x = 2,size=3)

ggsave(p39,filename = "ICGC.Pearson.CEP55.AFP.pdf",height = 3,width = 3)

########## Stem index -> mRNAsi + mDNAsi #########################
#install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser)
synLogin(email = "476418002@qq.com", password = "zhouqian123.")
synRNA <- synGet( "syn2701943", downloadLocation = "../../PCBC-StemIndex/" )
# ËØªÂÖ•Ë°®ËææÊï∞ÊçÆ
library(tidyverse)

exp <- read_delim(file = synRNA$path) %>%
  # ÂéªÈô§ Ensembl ID ÁöÑÂêéÁºÄ
  separate(col = "tracking_id", sep = "\\.", into = c("Ensembl_ID", "suffix")) %>%
  dplyr::select(-suffix) %>%
  column_to_rownames("Ensembl_ID") %>%
  as.matrix()

###‰ΩøÁî® synTableQuery ÂáΩÊï∞Êù•Ëé∑ÂèñÂÖÉÊï∞ÊçÆÔºå‰ΩøÁî? SQL ËØ≠Ê≥ïÈÄâÊã© syn3156503 Êï∞ÊçÆ‰∏≠ÁöÑ UID Âí? Diffname_short Âà?
synMeta <- synTableQuery("SELECT UID, Diffname_short FROM syn3156503")
#Â§ÑÁêÜÊï∞ÊçÆ
metaInfo <- synMeta$asDataFrame() %>%
  dplyr::select(UID, Diffname_short) %>%
  column_to_rownames("UID") %>%
  filter(!is.na(Diffname_short))

X <- exp
y <- metaInfo[colnames(X), ]
names(y) <- colnames(X)

## 2. ÊûÑÂª∫Ê®°Âûã
#gelnet(X, y, l1, l2)
#X: Ë°å‰∏∫Ê†∑Êú¨ÔºåÂàó‰∏∫Âü∫Âõ? => transpose(X.sc)
#y: ËøôÈáå‰ΩøÁî®ÂçïÁ±ªÊ®°Âûã => NULL
#l1: L1-Ê≠£ÂàôÂåñÊÉ©ÁΩöÁ≥ªÊï? => 0
#l2: L2-Ê≠£ÂàôÂåñÊÉ©ÁΩöÁ≥ªÊï? => 1

# ÂØπÊï∞ÊçÆËøõË°åÂùáÂÄº‰∏≠ÂøÉÂåñ
X <- exp
m <- apply(X, 1, mean)
X <- X - m
# Â∞ÜÊ†∑Êú¨ÂàÜ‰∏∫Âπ≤ÁªÜËÉûÁªÑÂíåÈùûÂπ≤ÁªÜËÉûÁª?
sc <- which(y == "SC")
X.sc <- X[, sc]
X.or <- X[, -sc]

library(gelnet)
model.RNA <- gelnet(t(X.sc), NULL, 0, 1)
save(X, y, model.RNA, file = "../..//PCBC-StemIndex/model.rda")

##### Preidict mRNAsi for LIHC #####
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
load("Inter.miRNA.mRNA.Samples.361.Rdata")
mRNA.Exp.361 <- mRNA.Exp[,Inter.Sample.361]

common <- intersect(names(model.RNA$w), rownames(mRNA.Exp.361))
X <- mRNA.Exp.361[common, ]
w <- model.RNA$w[common]
#ÂØπ‰∫é RNA Ë°®ËææÊï∞ÊçÆÔºå‰ΩøÁî? spearman ËÆ°ÁÆóÊùÉÈáç‰∏éË°®ËææÂÄº‰πãÈó¥ÁöÑÁõ∏ÂÖ≥ÊÄßÊù•Ë°°ÈáèÊ†∑Êú¨ÁöÑÂπ≤ÊÄßÊåáÊï∞ÔºåÂπ∂ËøõË°åÊ†áÂáÜÂåñ‰ΩøÂÖ∂ËêΩÂú® [0,1] ‰πãÈó¥
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)

StemIndex <- data.frame(Score=score) %>% rownames_to_column("SampleID")
write.csv(StemIndex,"LIHC.361.csv",row.names = F)

StemIndex <- read.csv("LIHC.361.csv")

#GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
#GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000173894","ENSG00000138180"),StemIndex$SampleID] %>%
  t() %>% data.frame(check.names = F) %>% cbind(StemIndex) %>%
  mutate(CBX2Group = ifelse(ENSG00000173894 >= median(ENSG00000173894),"CBX2high","CBX2low")) %>%
  mutate(CEP55Group = ifelse(ENSG00000138180 >= median(ENSG00000138180),"CEP55high","CEP55low"))

mRNA.Exp.log2.TPM.Target$CBX2Group <- factor(mRNA.Exp.log2.TPM.Target$CBX2Group,levels = c("CBX2high","CBX2low"))
mRNA.Exp.log2.TPM.Target$CEP55Group <- factor(mRNA.Exp.log2.TPM.Target$CEP55Group,levels = c("CEP55high","CEP55low"))
write.csv(mRNA.Exp.log2.TPM.Target,file = "SouceData.F5FG.csv")
p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CBX2Group", "Score", 
             fill = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("Stem score",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/Stemindex.CBX2.pdf",height = 2.5,width = 2.5)

p131 <- ggplot(mRNA.Exp.log2.TPM.Target,aes(CBX2Group,Score,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="mRNAsi score",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test")
ggsave(p131,filename = "CBX2.Group.mRNAsi.pdf",height = 2.2,width = 2.5)


p132 <- ggplot(mRNA.Exp.log2.TPM.Target,aes(CEP55Group,Score,fill=CEP55Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CEP55Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CEP55Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CEP55Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="mRNAsi score",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test")
ggsave(p132,filename = "CEP55.Group.mRNAsi.pdf",height = 2.2,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CEP55Group", "Score", 
             fill = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("Stem score",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/Stemindex.CEP55.pdf",height = 2.5,width = 2.5)

p13 <- ggplot(mRNA.Exp.log2.TPM.Target,aes(CBX2Group,ENSG00000138180,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="log2(TPM+1) CEP55",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "t.test",angle=-90,label.x=2)
ggsave(p13,filename = "CBX2.Group.CEP55.pdf",height = 1.5,width = 4)

library(ggpubr)
library(RColorBrewer)
library(ggthemes)
for (gene in c("ENSG00000173894","ENSG00000138180","ENSG00000215424","ENSG00000206195","ENSG00000232995","ENSG00000240498")) {
  p<-ggscatter(mRNA.Exp.log2.TPM.Target, x= gene, y = "Score", #color=brewer.pal(8,"Set1")[1],
            #palette = brewer.pal(8,"Set1")[1],#"#426671", size =2,#color = "group", palette = "Set1",
            add = "reg.line",
            add.params = list(color = brewer.pal(8,"Set1")[1], fill = "#E7E1D7"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson",label.x = 1, label.sep = "\n"),
            ggtheme = theme_few()) +
    #stat_cor(method = "pearson",label.x = 1)+
    labs(x = paste("log2(TPM+1) ",gene,sep = ""), 
         y = paste("Stem score",sep = ""))
  ggsave(p,filename = paste("TIMER2.CBX2.CEP55/StemScore.",gene,".pdf",sep = ""),width = 2.5,height = 2.5)
}

## risk <=> stem
#save(risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table,file = "Survival.Signature.RiskScore.Rdata")
load("Survival.Signature.RiskScore.Rdata")
StemIndex <- StemIndex %>% mutate(Sample = substr(SampleID,1,15))
Mid.dataset <- risk_score_table_multi_cox2[StemIndex$Sample,] %>% data.frame(check.names = F) %>% cbind(StemIndex)
write.csv(Mid.dataset,file = "SourceData.F5H.csv")
Mid.dataset$RiskScore <- factor(Mid.dataset$RiskScore,levels = c("High","Low"))
p<-ggboxplot(Mid.dataset, "RiskScore", "Score", 
             fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("Stem score",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/Stemscore.Riskgroup.pdf",height = 2.5,width = 2.5)


#### stem score => KM ####
library(survival)
library(survminer)
group=ifelse(Mid.dataset$Score > median(Mid.dataset$Score),"high","low")
if(length(table(group))==1) return(NULL) 
diff=survdiff(Surv(OS.time, OS) ~group,data = Mid.dataset)

rt2 <- Mid.dataset %>% dplyr::select(OS.time, OS) %>% mutate(group=group)

sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","group",sep = '')), data=rt2)

pValue=1-pchisq(diff$chisq,df=1)

p<-ggsurvplot(sfit1,pval =TRUE, data = rt2, 
              surv.median.line = "hv",
              legend.title = "Stem",
              conf.int.style = "step",
              xlab = "Time in days",
              #break.time.by = 500,
              risk.table = "abs_pct",
              risk.table.y.text.col = T,
              risk.table.y.text = FALSE,
              legend.labs = c("High", "Low"),
              #pval = TRUE,
              conf.int = TRUE,
              palette = "Set1",
              ggtheme = theme_bw())
print(p)
#dev.off()
library(export)
graph2pdf(file=paste("TIMER2.CBX2.CEP55/Stemscore.OS.KM.pdf",sep = ''),height = 6,width = 5)

## stem score => COX
m = coxph(Surv(OS.time, OS) ~ Score, data =  Mid.dataset)
#‰πüÂèØ‰ΩøÁî®‰∫åÂàÜÁ±ªÂèò???
#meta$group=ifelse(gene>median(gene),'high','low') 
#m=coxph(Surv(time, event) ~ Score, data =  Mid.dataset)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se

summary(m)

## 
Mid.dataset$StemScore = ifelse(Mid.dataset$Score >= median(Mid.dataset$Score),"High","Low")

tab <- table(Mid.dataset$StemScore,Mid.dataset$RiskScore)
chisq.test(tab)
p<-ggboxplot(Mid.dataset, "StemScore", "total_risk_score", 
             fill = "StemScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("Risk score",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/Stemscore.Risk.pdf",height = 2.5,width = 2.5)




##################### TMB ###########################
#### Data download ####
library(TCGAbiolinks)
library(tidyverse)
maf <- GDCquery_Maf(
  tumor = "LIHC",
  save.csv = TRUE,
  directory = "./TMB",
  pipelines = "mutect2"
)

library(maftools)
LIHC.maf <- read.maf("./TMB/TCGA-LIHC/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/a630f0a0-39b3-4aab-8181-89c1dde8d3e2/TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf.gz",isTCGA = TRUE)

plotmafSummary(maf = LIHC.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = LIHC.maf, colors = vc_cols, top = 15)


LIHC.maf@data[["TumorVAF"]] <- LIHC.maf@data$t_alt_count / LIHC.maf@data$t_depth
plotVaf(maf = LIHC.maf, vafCol = 'TumorVAF')

#### TMB calculate ####
get_TMB <- function(file) {
  # ÈúÄË¶ÅÁî®Âà∞ÁöÑÂà?
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  # Âà†Èô§Ëøô‰∫õÁ™ÅÂèòÁ±ªÂûã
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  # ËØªÂèñÊñá‰ª∂
  df <- read_csv(file, col_select = use_cols)
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    # ËÆ°ÁÆó VAF
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}

#Â§ñÊòæÂ≠êÁöÑÊÄªÈïøÂ∫¶Êàë‰ΩøÁî®ÁöÑÊòØ 30MÔºåÊàëÁúãÂèñ 38M„Ä?40M ÁöÑÈÉΩÊúâÔºå‰∏çÊòØÂæàÁªü‰∏Ä„ÄÇÂ¶ÇÊûúÊàë‰ª¨‰ª• TMB ‰∏≠‰ΩçÂÄºÊù•Âå∫ÂàÜÈ´ò‰Ωé TMB ÁªÑÁöÑËØùÔºåÂπ∂‰∏çÊòØÈÄâÁî®‰∏Ä‰∏™Âõ∫ÂÆöÁöÑÂÄºÊù•ÂàÜÔºåÂÖ∂ÂÆûÊ≤°‰ªÄ‰πàÂΩ±Âìç„Ä?
LIHC.csv <- get_TMB('./TMB/TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf.csv')
load("Inter.miRNA.mRNA.Samples.361.Rdata")
LIHC.TMB <- LIHC.csv %>% mutate(SampleID = substr(Tumor_Sample_Barcode,1,15),Sample=substr(Tumor_Sample_Barcode,1,16)) %>% filter(TMB != max(TMB)) %>%
  filter(Sample %in% Inter.Sample.361) %>% merge.data.frame(.,Mid.dataset,by.x="SampleID",by.y = "Sample")

p<-ggboxplot(LIHC.TMB, "RiskScore", "TMB", 
             fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("TMB",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/TMB.Riskgroup.pdf",height = 2.5,width = 2.5)

write.csv(LIHC.csv,"TCGA-LIHC.TMB.csv",row.names = F)

load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2$Sample.ID = rownames(risk_score_table_multi_cox2)

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
colnames(mRNA.Exp.log2.TPM) = substr(colnames(mRNA.Exp.log2.TPM),1,15)

LIHC.TMB <- read.csv("TCGA-LIHC.TMB.csv")
LIHC.TMB$Sample.ID <- substr(LIHC.TMB$Tumor_Sample_Barcode,1,15)

mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000173894","ENSG00000138180"),
                                              rownames(risk_score_table_multi_cox2)] %>% t() %>%
  data.frame(check.names = F) %>% rownames_to_column("Sample.ID") %>% merge(.,risk_score_table_multi_cox2,by="Sample.ID") %>%
  merge(.,LIHC.TMB,by="Sample.ID") %>% mutate(CBX2Group = ifelse(ENSG00000173894 >= median(ENSG00000173894),"CBX2high","CBX2low")) %>%
  mutate(CEP55Group = ifelse(ENSG00000138180 >= median(ENSG00000138180),"CEP55high","CEP55low"))

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "RiskScore", "TMB", 
             fill = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("TMB",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/TMB.Riskgroup.pdf",height = 2.5,width = 2.5)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CBX2Group", "TMB", 
             fill = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("TMB",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/TMB.CBX2.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CEP55Group", "TMB", 
             fill = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("TMB",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T)
ggsave(p,filename = "TIMER2.CBX2.CEP55/TMB.CEP55.pdf",height = 2.5,width = 2.5)

#### TMB survival ####
library(survival)
library(survminer)
group=ifelse(mRNA.Exp.log2.TPM.Target$TMB > median(mRNA.Exp.log2.TPM.Target$TMB),"high","low")
if(length(table(group))==1) return(NULL) 
diff=survdiff(Surv(OS.time, OS) ~group,data = mRNA.Exp.log2.TPM.Target)

rt2 <- mRNA.Exp.log2.TPM.Target %>% dplyr::select(OS.time, OS) %>% mutate(group=group)

sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","group",sep = '')), data=rt2)

pValue=1-pchisq(diff$chisq,df=1)

p<-ggsurvplot(sfit1,pval =TRUE, data = rt2, 
              surv.median.line = "hv",
              legend.title = i,
              conf.int.style = "step",
              xlab = "Time in days",
              #break.time.by = 500,
              #risk.table = "abs_pct",
              #risk.table.y.text.col = T,
              #risk.table.y.text = FALSE,
              legend.labs = c("High", "Low"),
              #pval = TRUE,
              conf.int = TRUE,
              palette = "Set1",
              ggtheme = theme_bw())
print(p)
#dev.off()
library(export)
graph2pdf(file=paste("TIMER2.CBX2.CEP55/TMB.OS.KM.pdf",sep = ''),height = 2.5,width = 2.5)

## stem score => COX
m = coxph(Surv(OS.time, OS) ~ TMB, data =  mRNA.Exp.log2.TPM.Target)
#‰πüÂèØ‰ΩøÁî®‰∫åÂàÜÁ±ªÂèò???
#meta$group=ifelse(gene>median(gene),'high','low') 
#m=coxph(Surv(time, event) ~ Score, data =  Mid.dataset)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se



#"ENSG00000173894", "ENSG00000138180"




#####################  CBX2 CEP55 Expression level at Pancancer level   ############################
files <- list.files("../Pancancer/",pattern = "[.DEG.csv]",full.names = T)
files <- files[str_detect(files,"DEG.csv")]
library(RColorBrewer)
Data.set <- data.frame()
for (TCGA in files) {
  mid.dataset <- read.csv(TCGA,row.names = 1) #%>% remove_rownames() %>% column_to_rownames("x")
  Data.set <- mid.dataset[c("ENSG00000173894","ENSG00000138180"),] %>% data.frame() %>% dplyr::select(log2FoldChange,padj,Threshold) %>% 
    mutate(Threshold = ifelse(padj <= 0.01 & log2FoldChange > 1,"Up",if_else(padj <= 0.01 & log2FoldChange < -1,"Down","No"))) %>%
    mutate(PROJECT=(str_remove(TCGA,".*\\/") %>% str_remove("_.*"))) %>% mutate(SYMBOL=c("CBX2","CEP55")) %>% rbind.data.frame(Data.set)
  
}
#GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
#GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180
Data.set2 <- Data.set %>% remove_rownames() %>%
  mutate(Log2FC = cut(log2FoldChange,c(-Inf,-2,-1,0,1,2,Inf),
                      right=T,
                      labels=c("<-2","-2~-1","-1~0","0~1","1~2",">2"))) %>%
  mutate(Log2FC=as.character(Log2FC)) %>%
  mutate(Log2FC2=if_else(padj>0.05,">0.05","<=0.05")) %>%
  mutate(PROJECT=str_remove(PROJECT,"TCGA-"))

Data.set2$Log2FC[which(Data.set2$Log2FC2 == ">0.05")] <- "P>0.05"
Data.set2$Log2FC <- factor(Data.set2$Log2FC,levels = rev(c("P>0.05","<-2","-2~-1","-1~0","0~1","1~2",">2")))
write.csv(Data.set2,file = "CBX2.CEP55.Pancancer/Expression.csv",row.names = F)

Data.set2 <- read.csv(file = "CBX2.CEP55.Pancancer/Expression.csv")
Data.set2$Log2FC <- factor(Data.set2$Log2FC,levels = rev(c("P>0.05","<-2","-2~-1","-1~0","0~1","1~2",">2")))
colors <- c("#9E0142","#F46D43","#FEE08B","#E6F598","#66C2A5","#5E4FA2","#F7F7F7")
names(colors) = rev(c("P>0.05","<-2","-2~-1","-1~0","0~1","1~2",">2"))
library(tidyverse)
p0929.1 <- ggplot()+
  geom_point(data=Data.set2,aes(SYMBOL,PROJECT,size=-log10(padj),
                                fill=Log2FC),
             color="white",shape=21)+
  geom_point(data=Data.set2%>%filter(padj<=0.05),
             aes(SYMBOL,PROJECT,size=-log10(padj),fill=Log2FC),
             color="black",stroke=1.2,shape=21)+
  scale_size_continuous("-log10(P.adj)",range = c(4,6))+
  scale_shape_manual("Symbol",values = c(21,22))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_fill_nejm()+
  scale_fill_manual("log2FC",values = colors)+
  labs(x="",y="")

Survival.Data <- read.csv("CBX2.CEP55.Pancancer/Pancancer.CBX2.CEP55.Survival.csv")
Survival.Data2 <- Survival.Data %>% filter(Cancer %in% Data.set2$PROJECT) %>%
  mutate(Survival=factor(Survival,levels=c("OS","DSS","DFI","PFI"))) %>%
  mutate(Prognosis=if_else(KM>0.05,"P>0.05",if_else(HR>1,"Risky","Protective"))) %>%
  mutate(Prognosis=factor(Prognosis,levels = c("Risky","Protective","P>0.05")))

p0929.2 <- ggplot()+
  geom_point(data=Survival.Data2,aes(Survival,Cancer,size=-log10(KM),
                                     fill=Prognosis,group=Id),
             color="white",shape=22)+
  geom_point(data=Survival.Data2%>%filter(KM<=0.05),
             aes(Survival,Cancer,size=-log10(KM),fill=Prognosis,group=Id),
             color="black",stroke=1.2,shape=22)+
  scale_size_continuous("-log10(KM P)",range = c(4,6))+
  scale_shape_manual("Symbol",values = c(21,22))+
  facet_wrap(~Id)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
    strip.text = element_text(size=15))+
  #ggsci::scale_fill_nejm()+
  scale_fill_manual(values = c(brewer.pal(8,"Dark2")[c(4,1)],"#F7F7F7"))+
  #scale_fill_gradient2(low = "#053061",mid = "#F7F7F7",high = "#67001F")+
  labs(x="",y="")

p0929.3 <- p0929.2 %>% aplot::insert_left(p0929.1,width = .5)
ggsave(p0929.3,filename = "CBX2.CEP55.Pancancer/Pancancer.CBX2.CEP55.Survival.pdf",height = 10,width = 6)

## 
library(readxl)
library(tidyverse)
Pathway.Score <- read_xlsx("CBX2.CEP55.Pancancer/ExpressionAndRPPATable.xlsx")
Pathway.Score2 <- Pathway.Score %>% filter(fdr <= 0.05) %>%
  group_by(symbol,pathway,class) %>%
  summarise(Count=n()) %>%
  ungroup() %>%
  mutate(Percent=round(Count*100/32,2)) %>%
  dplyr::select(-Count) %>%
  spread(class,Percent,fill=0) %>%
  mutate(None=100-Activation-Inhibition) %>%
  gather(class,Percent,Activation:None) %>%
  data.frame() %>%
  mutate(class=factor(class,levels = rev(c("Activation","Inhibition","None"))))
  
p0930.1<-ggplot(Pathway.Score2,aes(Percent,pathway,fill=class))+
  geom_col()+
  facet_wrap(~symbol)+
  ggthemes::theme_few()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"),
        strip.text = element_text(size=12))+
  labs(x="Percent in 32 cancers",y="")+
  scale_fill_manual("",values = rev(c("#E64B35FF","#4DBBD5FF","#7E6148FF")))
  
ggsave(p0930.1,filename = "CBX2.CEP55.Pancancer/CBX2.CEP55.PathwayScore.pdf",height = 4,width = 4)


p<-ggplot(Data.set %>% filter(PROJECT != "TCGA-LIHC"),aes(Threshold,SYMBOL,color=log2FoldChange,size = -log10(padj))) +
  geom_tile(color="white",fill="white",size=0)+
  geom_point(shape=15) +
  facet_wrap(~PROJECT)+
  #geom_text(aes(label=signif),size=6,color="white",hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL) + 
  #scale_color_viridis_c()+
  theme_few()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=1, linetype="solid")) +
  scale_size(range=c(1,5))+
  theme(axis.text = element_text(colour = "black"))+
  scale_color_gradientn(colors = c(brewer.pal(8,"Set2")[8],brewer.pal(8,"Set1")[2], brewer.pal(8,"Set2")[1],brewer.pal(8,"Set1")[5],brewer.pal(8,"Paired")[5],brewer.pal(8,"Set1")[4],brewer.pal(8,"Set1")[1],brewer.pal(8,"Set1")[1]))

ggsave(p,filename = "CBX2.CEP55.Pancancer.pdf",height = 6,width = 8)


#####################  Pancancer => KM COX => CODE => Long time ####################
dirs <- list.dirs("../../Clinical.XENA/",full.names = T)
dirs <- dirs[str_detect(dirs,"//")]
#dir.create("../../KM.Cox.Pancancer")
Pancancer.CBX2.CEP55 <- data.frame()
###
Study.Gene.ENSEMBL <- c("ENSG00000173894","ENSG00000138180","ENSG00000073111",
                        "ENSG00000232995","ENSG00000215424","ENSG00000206195",
                        "ENSG00000240498",
                        "ENSG00000091831","ENSG00000107864")
for (dir in dirs) {
  project <- str_remove(dir,".*/")
  print(project)
  phenotype <- read.csv(paste(dir,"/TCGA-",project,".GDC_phenotype.tsv",sep = ""),sep = "\t",header = T) %>% 
    dplyr::filter(sample_type.samples == "Primary Tumor" | sample_type.samples == "Primary Blood Derived Cancer - Peripheral Blood")
  TPM <- read.csv(paste("../../mRNA.PanCancer.Exp/TCGA-",project,".mRNA.TPM.csv",sep = ""),check.names = F,row.names = 1)
  TPM.log <- log2(TPM+1)
  Survival.D <- read.csv(paste("../../Clinical.XENA/XENA.TCGA.",project,"_survival.txt",sep = ""),
                         row.names = 1,sep = "\t",header = T) %>%
    rownames_to_column("sample")
  Inter.sample <- intersect(intersect(phenotype$submitter_id.samples,colnames(TPM.log)) %>% substr(1,15),
                            Survival.D$sample)
  TPM.log2 <- TPM.log %>% magrittr::set_colnames(substr(colnames(.),1,15))
  mRNA.Exp.Target <- TPM.log2[,Inter.sample] %>% data.frame(check.names = F)
  mRNA.Exp.Target.Survival <- mRNA.Exp.Target[c("ENSG00000173894","ENSG00000138180"),] %>% t() %>% data.frame(check.names = F) %>% mutate(sample = rownames(.)) %>% 
    merge(.,Survival.D,by="sample") %>% dplyr::select(-Redaction)
  
  OS.Tab=data.frame()
  for(i in c("ENSG00000173894","ENSG00000138180")){
    if (i=="ENSG00000173894") {
      symbol="CBX2"
    }else{symbol="CEP55"}
    print(i)
    #symbol <- SNHG.DE[SNHG.DE$ENSEMBL == i,]$SYMBOL %>% as.character()
    #cox
    for (sur in c("OS","DSS","DFI","PFI")) {
      middata <- mRNA.Exp.Target.Survival %>%
        dplyr::select(c(i,sur,paste(sur,".time",sep = ""))) %>%
        na.omit()
      if (nrow(middata) >= 10) {
        Formula <- as.formula(paste("Surv(",sur,".time,",sur,") ~ ",i,sep = ""))
        cox <- coxph(Formula, data = middata)
        coxSummary = summary(cox)
        #KM
        group=ifelse(middata[[i]] > median(middata[[i]]),"high","low")
        if(length(table(group))==1) {
          OS.Tab=rbind(OS.Tab,
                       cbind(Id=symbol,
                             Survival=sur,
                             KM=NA,
                             HR=NA,
                             HR.95L=NA,
                             HR.95H=NA,
                             pvalue=NA))
        } else {
          new.formula <- as.formula(paste("Surv(",sur,".time,",sur,") ~ group",sep = ""))
          diff=survdiff(new.formula,data = middata)
          
          #rt2 <- mRNA.Exp.Target.Survival %>% dplyr::select(OS.time, OS) %>% mutate(group=group)
          #sfit1=survfit(as.formula(paste("Surv(OS.time, OS)~","group",sep = '')), data=rt2)
          
          pValue=1-pchisq(diff$chisq,df=1)
          OS.Tab=rbind(OS.Tab,
                       cbind(Id=symbol,
                             Survival=sur,
                             KM=pValue,
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR.95L=coxSummary$conf.int[,"lower .95"],
                             HR.95H=coxSummary$conf.int[,"upper .95"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
        }
      }
      Pancancer.CBX2.CEP55 <- OS.Tab %>% mutate(Cancer = project) %>%
        rbind(Pancancer.CBX2.CEP55)
      }
    }
}

write.csv(Pancancer.CBX2.CEP55,file = "CBX2.CEP55.Pancancer/Pancancer.CBX2.CEP55.Survival.csv",row.names = F)

Pancancer.CBX2.CEP55 <- read.csv("../LIHC/Pancancer.CBX2.CEP55.csv")
Pancancer.CBX2.CEP55 <- Pancancer.CBX2.CEP55 %>% mutate(Type = if_else(KM <= 0.05 & pvalue <= 0.05,if_else(HR > 1,"Risky","Protective"),"p>0.05"))
write.csv(Pancancer.CBX2.CEP55,file = "NHG.KM.Pancancer.csv",row.names = F)

Pancancer.CBX2.CEP55.Figure <- Pancancer.CBX2.CEP55 %>% dplyr::select(id,Cancer,Type,Gene.name) %>% spread(Cancer,Type) %>% 
  mutate(id = factor(id, levels = Study.Gene.ENSEMBL)) %>% dplyr::arrange(id) %>% 
  remove_rownames() %>% column_to_rownames("id") %>% remove_rownames() %>% column_to_rownames("Gene.name")

colors <- structure(c(brewer.pal(9,"Set1")[c(1,2)],brewer.pal(8,"Pastel2")[8]),names=c('Risky','Protective','p>0.05'))
Heatmap(Pancancer.CBX2.CEP55.Figure %>% as.matrix(), col = colors,border = "black",na_col = "grey60",rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns = T,
        heatmap_legend_param = list(title = "",title_position = "topcenter"))

graph2pdf(file="../LIHC/Pancancer.CBX2.CEP55.DE.KM.COX.Pancancer.pdf",height=4.5,width=8.5)

#####################  CBX2 CEP55 => Group => DEG ##############
library(TCGAbiolinks)
mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor")
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
#GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
#GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180

#371
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[,intersect(mRNA.Phenodata$Sample.ID,colnames(mRNA.Exp.log2.TPM))] %>%
  data.frame(check.names = F)

G_list <- read.csv("../LIHC/biomaRt.GeneID.tansfer.csv")
dat.2 <- mRNA.Exp.log2.TPM.Target %>% data.frame(check.names = F) %>%
  rownames_to_column("ensembl_gene_id") %>% data.frame(check.names = F) %>%
  merge(.,G_list[,1:2],by="ensembl_gene_id") %>% 
  data.frame(check.names = F) %>%
  na.omit() %>% dplyr::select(-ensembl_gene_id)
expr_mean=aggregate(.~hgnc_symbol,mean,data=dat.2) #calculate mean for same symbol
dat.3 <- expr_mean %>% remove_rownames() %>% column_to_rownames("hgnc_symbol")
#dat.3$Description = "na"

write.table(dat.3,"LIHC.Symbol.LogTPM.txt",sep = "\t",quote = F)


#### CBX2 => DEG => GO + KEGG ####
gene <- "ENSG00000173894"
gene.exp <- mRNA.Exp.log2.TPM.Target[gene,]
label <- if_else(gene.exp < median(as.numeric(gene.exp)), 0, 1)
#write.table(t(data.frame(label)),"LIHC.Symbol.LogTPM.CBX2.txt",sep = "\t")

group.low <- mRNA.Exp.log2.TPM.Target[,label == 0]
group.high <- mRNA.Exp.log2.TPM.Target[,label == 1]

DEGs <- TCGAanalyze_DEA(
  mat1 = group.low,
  mat2 = group.high,
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "CBX2_Low",
  Cond2type = "CBX2_High",
  fdr.cut = 0.01,
  logFC.cut = 0,
)

dim(DEGs) #1338 DEG
DEGs.CBX2 <- DEGs %>% rownames_to_column("ENSEMBL") %>% merge.data.frame(.,GeneType.Unique,by.x = "ENSEMBL",by.y="Gene.stable.ID")
write.csv(DEGs.CBX2,file = "CBX2.DEG.csv",row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

DEGs.CBX2 <- read.csv("CBX2.DEG.csv")

Up.CBX2 <- DEGs.CBX2 %>% filter(logFC >= 1 & adj.P.Val <= 0.01)
Down.CBX2 <- DEGs.CBX2 %>% filter(logFC <= -1 & adj.P.Val <= 0.01)

Up.go <- enrichGO(
  gene = Up.CBX2$ENSEMBL,keyType="ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)

write.csv(Up.go,"CBX2.DEG.UP.GO.2.csv")

Down.go <- enrichGO(
  gene = Down.CBX2$ENSEMBL,keyType="ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)

write.csv(Down.go,"CBX2.DEG.Down.GO.2.csv")

library(clusterProfiler)
DEGs.CBX2 <- read.csv("CBX2.DEG.csv")

Up.CBX2 <- DEGs.CBX2 %>% filter(logFC >= 1 & adj.P.Val <= 0.01)
Down.CBX2 <- DEGs.CBX2 %>% filter(logFC <= -1 & adj.P.Val <= 0.01)

Up.KEGG <- enrichKEGG(
  gene= Up.CBX2$NCBI.gene..formerly.Entrezgene..ID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)

save(Up.KEGG,file = "CBX2.Up.KEGG.Rdata")

Down.KEGG <- enrichKEGG(
  gene= Down.CBX2$NCBI.gene..formerly.Entrezgene..ID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)

save(Down.KEGG,file = "CBX2.Down.KEGG.Rdata")

DEGs.CEP55 <- read.csv("CEP55.DEG.csv")

Up.CEP55 <- DEGs.CEP55 %>% filter(logFC >= 1 & adj.P.Val <= 0.01)
Down.CEP55 <- DEGs.CEP55 %>% filter(logFC <= -1 & adj.P.Val <= 0.01)

Up.KEGG <- enrichKEGG(
  gene= Up.CEP55$NCBI.gene..formerly.Entrezgene..ID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)

save(Up.KEGG,file = "CEP55.Up.KEGG.Rdata")

Down.KEGG <- enrichKEGG(
  gene= Down.CEP55$NCBI.gene..formerly.Entrezgene..ID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)

save(Down.KEGG,file = "CEP55.Down.KEGG.Rdata")

#### KEGG Enrichment => Figure ####
load(file="CBX2.Up.KEGG.Rdata")
Up.KEGG@result %>% data.frame() %>% filter(p.adjust <= 0.05) %>% write.csv("CBX2.Up.KEGG.csv",row.names = F)
load(file="CBX2.Down.KEGG.Rdata")
Down.KEGG@result %>% data.frame() %>% filter(p.adjust <= 0.05) %>% write.csv("CBX2.Down.KEGG.csv",row.names = F)

load(file="CEP55.Up.KEGG.Rdata")
Up.KEGG@result %>% data.frame() %>% filter(p.adjust <= 0.05) %>% write.csv("CEP55.Up.KEGG.csv",row.names = F)
load(file="CEP55.Down.KEGG.Rdata")
Down.KEGG@result %>% data.frame() %>% filter(p.adjust <= 0.05) %>% write.csv("CEP55.Down.KEGG.csv",row.names = F)

CBX2.Up.KEGG <- read.csv("CBX2.Up.KEGG.csv") %>% arrange(desc(p.adjust)) %>%
  mutate(Description=factor(Description,levels = .$Description))

p12 <- ggplot(CBX2.Up.KEGG,aes(x=-log10(p.adjust),y=Description,size=Count,color=-log10(p.adjust)))+
  #geom_bar(stat="identity",position = position_dodge(0.9)) +
  geom_point()+
  #coord_flip()+
  ggthemes::theme_few() +
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="")+
  scale_color_gradientn(colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_x_continuous(limits = c(0,13))

ggsave(p12,filename = "CBX2.Up.KEGG.pdf",width = 5,height = 3)


CBX2.Down.KEGG <- read.csv("CBX2.Down.KEGG.csv") %>% arrange(desc(p.adjust)) %>%
  mutate(Description=factor(Description,levels = .$Description))

p13 <- ggplot(CBX2.Down.KEGG,aes(x=-log10(p.adjust),y=Description,size=Count,color=-log10(p.adjust)))+
  #geom_bar(stat="identity",position = position_dodge(0.9)) +
  geom_point()+
  #coord_flip()+
  ggthemes::theme_few() +
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="")+
  scale_color_gradientn(colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_x_continuous(limits = c(0,15))

ggsave(p13,filename = "CBX2.Down.KEGG.pdf",width = 7,height = 6)


CEP55.Up.KEGG <- read.csv("CEP55.Up.KEGG.csv") %>% arrange(desc(p.adjust)) %>%
  mutate(Description=factor(Description,levels = .$Description))

p16 <- ggplot(CEP55.Up.KEGG,aes(x=-log10(p.adjust),y=Description,size=Count,color=-log10(p.adjust)))+
  #geom_bar(stat="identity",position = position_dodge(0.9)) +
  geom_point()+
  #coord_flip()+
  ggthemes::theme_few() +
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="")+
  scale_color_gradientn(colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_x_continuous(limits = c(0,13))

ggsave(p16,filename = "CEP55.Up.KEGG.pdf",width = 7,height = 4)


CEP55.Down.KEGG <- read.csv("CEP55.Down.KEGG.csv") %>% arrange(desc(p.adjust)) %>%
  mutate(Description=factor(Description,levels = .$Description))

p15 <- ggplot(CEP55.Down.KEGG,aes(x=-log10(p.adjust),y=Description,size=Count,color=-log10(p.adjust)))+
  #geom_bar(stat="identity",position = position_dodge(0.9)) +
  geom_point()+
  #coord_flip()+
  ggthemes::theme_few() +
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(y="")+
  scale_color_gradientn(colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_x_continuous(limits = c(0,25))

ggsave(p15,filename = "CEP55.Down.KEGG.pdf",width = 7,height = 7)

#### CBX2 => GSEA (c2+c5) ####
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
#devtools::install_github("xjsun1221/tinyarray")
library(tinyarray)

#human <- msigdbr(species = "Homo sapiens")
#KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
#  dplyr::select(gs_exact_source,gene_symbol)

#GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
#  dplyr::select(gene_symbol,gs_exact_source,gs_subcat)
#GO_df = GO_df[GO_df$gs_subcat!="HPO",]
#GO_df = GO_df[,c(2,1)]

KEGG_df <- read.gmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO_df <- read.gmt("../../msigdb/c5.go.v7.5.1.symbols.gmt")
Hallmark <- read.gmt("../../msigdb/h.all.v7.5.1.symbols.gmt")

C2 <- read.gmt("../../msigdb/c2.all.v7.5.1.symbols.gmt")

DEGs.CBX2 <- read.csv(file = "CBX2.DEG.csv")

#dim(DEGs.CBX2) #1338 DEG
#DEGs.CBX2 <- DEGs %>% rownames_to_column("ENSEMBL") %>% merge.data.frame(.,GeneType.Unique,by.x = "ENSEMBL",by.y="Gene.stable.ID")
DEGs.CBX2.Mid <- aggregate(.~Gene.name,mean,data=(DEGs.CBX2 %>% dplyr::select(logFC,Gene.name))) %>% na.omit()

ge = DEGs.CBX2.Mid$logFC
names(ge) = DEGs.CBX2.Mid$Gene.name
ge = sort(ge,decreasing = T)

C2 <- read.gmt("../../msigdb/c2.all.v7.5.1.symbols.gmt")
em.C2 <- GSEA(ge, TERM2GENE = C2,eps = 1e-100)
save(em.C2,file = "GSEA/CBX2/GSEA.C2.Rdata")
write.csv(em.C2@result,file = "GSEA/CBX2/GSEA.C2.csv")

em.C2@result %>% filter(str_detect(ID,"REACTOME")) %>% View()

TERM <- C("CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP",
          "FISCHER_G2_M_CELL_CYCLE",
          "FISCHER_G1_S_CELL_CYCLE",
          "KONG_E2F3_TARGETS",
          "ISHIDA_E2F_TARGETS",
          "WHITFIELD_CELL_CYCLE_G2"
          "REACTOME_MITOTIC_PROMETAPHASE",
          "REACTOME_M_PHASE",
          "REACTOME_CELL_CYCLE_CHECKPOINTS")



library(fgsea)
GO_df$term <- unfactor(GO_df$term)
GO_df2 <- GO_df %>% split(x = .$gene, f = .$term)
em <- GSEA(ge, TERM2GENE = GO_df)
#Agsea_res <- fgsea(pathways = GO_df2, 
#                   stats = ge,
#                   minSize=5,
#                   maxSize=500,
#                   nperm=1000)
write.csv(em,file = "GSEA/CBX2/GSEA.GO.csv")
#write.csv(as.data.frame(Agsea_res),file = "GSEA/CBX2/GSEA.fgsea.GO.csv")
#gseaplot2(gsea,1,color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") 
check <- em@result$ID %in% c("GOBP_CHROMOSOME_SEGREGATION",
                             "GOBP_DNA_REPLICATION",
                             "GOBP_MEIOTIC_CELL_CYCLE",
                             "GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
                             "GOBP_MITOTIC_NUCLEAR_DIVISION",
                             "GOCC_SPINDLE",
                             "GOMF_TUBULIN_BINDING",
                             "GOBP_SPINDLE_ORGANIZATION",
                             "GOCC_MITOTIC_SPINDLE",
                             "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION")
index <- (1:nrow(em@result))[check]
gseaplot2(em,index,color="red",pvalue_table = T,base_size=10,ES_geom="line")
library(export)
graph2pdf(file="CBX2.GSEA.GO.gsea.pdf",height = 7,width = 7)
graph2pdf(file="CBX2.GSEA.GO.gsea2.pdf",height = 16,width = 10)
# (em@result[c(2,3,9,12:14,23,30,32,57,116),])$ID
#p2<-plotGseaTable(GO_df2[c("GOBP_CHROMOSOME_SEGREGATION",
#                       "GOBP_DNA_REPLICATION",
#                       "GOBP_MEIOTIC_CELL_CYCLE",
#                       "GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
#                       "GOBP_MITOTIC_NUCLEAR_DIVISION",
#                       "GOCC_SPINDLE",
#                       "GOMF_TUBULIN_BINDING",
#                       "GOBP_SPINDLE_ORGANIZATION",
#                       "GOCC_MITOTIC_SPINDLE",
#                       "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
#                       "GOMF_ARACHIDONIC_ACID_MONOOXYGENASE_ACTIVITY",
#                       "GOBP_XENOBIOTIC_CATABOLIC_PROCESS",
#                       "GOMF_AROMATASE_ACTIVITY")], ge, Agsea_res, 
#              gseaParam = 0.5)
#print(p2)
#ggsave(p2,filename = "CBX2.GSEA.GO.ggsea.pdf",height = 5,width = 8)
#graph2pdf(file="CBX2.GSEA.GO.ggsea.pdf",height = 5,width = 8)

## geneset_GO = msigdbr(species = "Homo sapiens",
#                     category = "C5", 
#                    subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
#geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#»•≥˝«∞◊∫KEGG_
#geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#Ω´¥Û–¥ªªŒ™–°–¥
#geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#Ω´_◊™ªØŒ™ø’∏Ò
## gsea_geneset_GO <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)

for (i in 1:500) {
  p<-gseaplot2(em, geneSetID = i, title = em$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CBX2/",em$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}

#intersect(rownames(Up.go@result),em$Description) %>% data.frame() %>% write.csv("CBX2.DEG.Up.GO.GSEA.Term.csv")
#intersect(rownames(Down.go@result),em$Description) %>% data.frame() %>% write.csv("CBX2.DEG.Down.GO.GSEA.Term.csv")

#### => KEGG 
em.kegg <- GSEA(ge, TERM2GENE = KEGG_df)
write.csv(em.kegg@result,file = "GSEA/CBX2/GSEA.KEGG.csv")

for (i in 68:500) {
  p<-gseaplot2(em.kegg, geneSetID = i, title = em.kegg$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CBX2/",em.kegg$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}

#### => Hallmarks
em.hallmark <- GSEA(ge, TERM2GENE = Hallmark)
write.csv(em.hallmark,file = "GSEA/CBX2/CBX2.GSEA.hallmark.csv")

for (i in 1:50) {
  p<-gseaplot2(em.hallmark, geneSetID = i, title = em.hallmark$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CBX2/",em.hallmark$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}

#### CEP55 => DEG => GO + KEGG ####
gene <- "ENSG00000138180"
gene.exp <- mRNA.Exp.log2.TPM.Target[gene,]
label <- if_else(gene.exp < median(as.numeric(gene.exp)), 0, 1)
group.low <- mRNA.Exp.log2.TPM.Target[,label == 0]
group.high <- mRNA.Exp.log2.TPM.Target[,label == 1]

DEGs <- TCGAanalyze_DEA(
  mat1 = group.low,
  mat2 = group.high,
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "CEP55_Low",
  Cond2type = "CEP55_High",
  fdr.cut = 0.01,
  logFC.cut = 1,
)

dim(DEGs) #3185 DEG
DEGs.CEP55 <- DEGs %>% rownames_to_column("ENSEMBL") %>% merge.data.frame(.,GeneType.Unique,by.x = "ENSEMBL",by.y="Gene.stable.ID")
write.csv(DEGs.CEP55,file = "CEP55.DEG.csv",row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
DEGs.CEP55 <- read.csv(file = "CEP55.DEG.csv")
Up.CEP55 <- DEGs.CEP55 %>% filter(logFC >= 1)
Down.CEP55 <- DEGs.CEP55 %>% filter(logFC <= -1)

Up.CEP55.go <- enrichGO(
  gene = Up.CEP55$ENSEMBL,keyType="ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)

write.csv(Up.CEP55.go,"CEP55.DEG.UP.GO.2.csv")

Down.CEP55.go <- enrichGO(
  gene = Down.CEP55$ENSEMBL,keyType="ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)

write.csv(Down.CEP55.go,"CEP55.DEG.Down.GO.2.csv")

## GSEA
DEGs <- TCGAanalyze_DEA(
  mat1 = group.low,
  mat2 = group.high,
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "CBX2_Low",
  Cond2type = "CBX2_High",
  fdr.cut = 1,
  logFC.cut = 0,
)

dim(DEGs) #56215 DEG
DEGs.CEP55 <- DEGs %>% rownames_to_column("ENSEMBL") %>% merge.data.frame(.,GeneType.Unique,by.x = "ENSEMBL",by.y="Gene.stable.ID")
write.csv(DEGs.CEP55,file = "GSEA/DEGs.CEP55.csv")

DEGs.CEP55.Mid <- aggregate(.~Gene.name,mean,data=(DEGs.CEP55 %>% dplyr::select(logFC,Gene.name))) %>% na.omit()
ge = DEGs.CEP55.Mid$logFC
names(ge) = DEGs.CEP55.Mid$Gene.name
ge = sort(ge,decreasing = T)
em2 <- GSEA(ge, TERM2GENE = GO_df)
dir.create("GSEA/CEP55/")
write.csv(em2,file = "GSEA/CEP55/GSEA.GO.csv")

intersect(rownames(Up.CEP55.go@result),em2$Description) %>% data.frame() %>% write.csv("CEP55.DEG.Up.GO.GSEA.Term.csv")
intersect(rownames(Down.CEP55.go@result),em2$Description) %>% data.frame() %>% write.csv("CEP55.DEG.Down.GO.GSEA.Term.csv")

for (i in 1:1000) {
  print(i)
  p<-gseaplot2(em2, geneSetID = 1, title = em$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CEP55/",em2$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}






#### CEP55 => GSEA (c2+c5) ####
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
#devtools::install_github("xjsun1221/tinyarray")
library(tinyarray)

#human <- msigdbr(species = "Homo sapiens")
#KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
#  dplyr::select(gs_exact_source,gene_symbol)

#GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
#  dplyr::select(gene_symbol,gs_exact_source,gs_subcat)
#GO_df = GO_df[GO_df$gs_subcat!="HPO",]
#GO_df = GO_df[,c(2,1)]

KEGG_df <- read.gmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO_df <- read.gmt("../../msigdb/c5.go.v7.5.1.symbols.gmt")
Hallmark <- read.gmt("../../msigdb/h.all.v7.5.1.symbols.gmt")

#DEGs.CEP55 <- read.csv(file = "CEP55.DEG.csv")

#dim(DEGs.CEP55) #1338 DEG
#DEGs.CEP55 <- DEGs %>% rownames_to_column("ENSEMBL") %>% merge.data.frame(.,GeneType.Unique,by.x = "ENSEMBL",by.y="Gene.stable.ID")
#DEGs.CEP55.Mid <- aggregate(.~Gene.name,mean,data=(DEGs.CEP55 %>% dplyr::select(logFC,Gene.name))) %>% na.omit()

#ge = DEGs.CEP55.Mid$logFC
#names(ge) = DEGs.CEP55.Mid$Gene.name
#ge = sort(ge,decreasing = T)

library(fgsea)
GO_df$term <- unfactor(GO_df$term)
GO_df2 <- GO_df %>% split(x = .$gene, f = .$term)
em <- GSEA(ge, TERM2GENE = GO_df)
#Agsea_res <- fgsea(pathways = GO_df2, 
#                   stats = ge,
#                   minSize=5,
#                   maxSize=500,
#                   nperm=1000)
write.csv(em,file = "GSEA/CEP55/GSEA.GO.csv")
#write.csv(as.data.frame(Agsea_res),file = "GSEA/CEP55/GSEA.fgsea.GO.csv")
#gseaplot2(gsea,1,color="red",pvalue_table = T,title="",base_size=10,ES_geom="line") 
check <- em@result$ID %in% c("GOBP_CHROMOSOME_SEGREGATION",
                             "GOBP_DNA_REPLICATION",
                             "GOBP_MEIOTIC_CELL_CYCLE",
                             "GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
                             "GOBP_MITOTIC_NUCLEAR_DIVISION",
                             "GOCC_SPINDLE",
                             "GOMF_TUBULIN_BINDING",
                             "GOBP_SPINDLE_ORGANIZATION",
                             "GOCC_MITOTIC_SPINDLE",
                             "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION")
index <- (1:nrow(em@result))[check]
gseaplot2(em,index,color="red",pvalue_table = T,base_size=10,ES_geom="line")
library(export)
graph2pdf(file="CEP55.GSEA.GO.gsea.pdf",height = 7,width = 7)
graph2pdf(file="CEP55.GSEA.GO.gsea2.pdf",height = 16,width = 10)
# (em@result[c(2,3,9,12:14,23,30,32,57,116),])$ID
#p2<-plotGseaTable(GO_df2[c("GOBP_CHROMOSOME_SEGREGATION",
#                       "GOBP_DNA_REPLICATION",
#                       "GOBP_MEIOTIC_CELL_CYCLE",
#                       "GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION",
#                       "GOBP_MITOTIC_NUCLEAR_DIVISION",
#                       "GOCC_SPINDLE",
#                       "GOMF_TUBULIN_BINDING",
#                       "GOBP_SPINDLE_ORGANIZATION",
#                       "GOCC_MITOTIC_SPINDLE",
#                       "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
#                       "GOMF_ARACHIDONIC_ACID_MONOOXYGENASE_ACTIVITY",
#                       "GOBP_XENOBIOTIC_CATABOLIC_PROCESS",
#                       "GOMF_AROMATASE_ACTIVITY")], ge, Agsea_res, 
#              gseaParam = 0.5)
#print(p2)
#ggsave(p2,filename = "CEP55.GSEA.GO.ggsea.pdf",height = 5,width = 8)
#graph2pdf(file="CEP55.GSEA.GO.ggsea.pdf",height = 5,width = 8)

## geneset_GO = msigdbr(species = "Homo sapiens",
#                     category = "C5", 
#                    subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
#geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#»•≥˝«∞◊∫KEGG_
#geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#Ω´¥Û–¥ªªŒ™–°–¥
#geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#Ω´_◊™ªØŒ™ø’∏Ò
## gsea_geneset_GO <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)

for (i in 1:500) {
  p<-gseaplot2(em, geneSetID = i, title = em$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CEP55/",em$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}

#intersect(rownames(Up.go@result),em$Description) %>% data.frame() %>% write.csv("CEP55.DEG.Up.GO.GSEA.Term.csv")
#intersect(rownames(Down.go@result),em$Description) %>% data.frame() %>% write.csv("CEP55.DEG.Down.GO.GSEA.Term.csv")

#### => KEGG 
em.kegg <- GSEA(ge, TERM2GENE = KEGG_df)
write.csv(em.kegg@result,file = "GSEA/CEP55/GSEA.KEGG.csv")

for (i in 1:100) {
  p<-gseaplot2(em.kegg, geneSetID = i, title = em.kegg$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CEP55/",em.kegg$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}

#### => Hallmarks
em.hallmark <- GSEA(ge, TERM2GENE = Hallmark)
write.csv(em.hallmark,file = "GSEA/CEP55/CEP55.GSEA.hallmark.csv")

for (i in 1:50) {
  p<-gseaplot2(em.hallmark, geneSetID = i, title = em.hallmark$Description[i])
  print(p)
  graph2pdf(file=paste("GSEA/CEP55/",em.hallmark$Description[i] %>% str_replace_all(":","_"),".pdf",sep = ""),width=2.5,height=3)
}


#### GO -> CBX2 ####
CBX2.Up <- read.csv("CBX2.DEG.UP.GO.2.csv")
CBX2.Up2 <- CBX2.Up %>% group_by(ONTOLOGY) %>% top_n(-10,p.adjust)
CBX2.Up2$ID <- factor(CBX2.Up2$ID,levels = rev(CBX2.Up2$ID))
p1<-ggplot(CBX2.Up2,aes(x=ID,y=-log10(p.adjust),fill=ONTOLOGY))+
  geom_bar(stat="identity",position = position_dodge(0.9)) +
  coord_flip()+
  theme_few() +
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=Description,y=0.1),hjust = 0, #vjust = 1.5,
            position = position_dodge(0.5),fontface = 'bold')+
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))+
  theme(legend.position = "top",axis.text.y = element_text(face = "bold"))+
  labs(x="")
ggsave(p1,filename = "CBX2.Up.pdf",height = 10,width = 6)


CBX2.Down <- read.csv("CBX2.DEG.Down.GO.2.csv")
CBX2.Down2 <- CBX2.Down %>% group_by(ONTOLOGY) %>% top_n(-10,p.adjust)
CBX2.Down2$ID <- factor(CBX2.Down2$ID,levels = rev(CBX2.Down2$ID))
p1<-ggplot(CBX2.Down2,aes(x=ID,y=-log10(p.adjust),fill=ONTOLOGY))+
  geom_bar(stat="identity",position = position_dodge(0.9)) +
  coord_flip()+
  theme_few() +
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=Description,y=0.1),hjust = 0, #vjust = 1.5,
            position = position_dodge(0.5),fontface = 'bold')+
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))+
  theme(legend.position = "top",axis.text.y = element_text(face = "bold"))+
  labs(x="")
ggsave(p1,filename = "CBX2.Down.pdf",height = 10,width = 6)






#### GO -> CEP55 ####
CEP55.Up <- read.csv("CEP55.DEG.UP.GO.2.csv")
CEP55.Up2 <- CEP55.Up %>% group_by(ONTOLOGY) %>% top_n(-10,p.adjust)
CEP55.Up2$ID <- factor(CEP55.Up2$ID,levels = rev(CEP55.Up2$ID))
p1<-ggplot(CEP55.Up2,aes(x=ID,y=-log10(p.adjust),fill=ONTOLOGY))+
  geom_bar(stat="identity",position = position_dodge(0.9)) +
  coord_flip()+
  theme_few() +
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=Description,y=0.1),hjust = 0, #vjust = 1.5,
            position = position_dodge(0.5),fontface = 'bold')+
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))+
  theme(legend.position = "top",axis.text.y = element_text(face = "bold"))+
  labs(x="")
ggsave(p1,filename = "CEP55.Up.pdf",height = 10,width = 6)


CEP55.Down <- read.csv("CEP55.DEG.Down.GO.2.csv")
CEP55.Down2 <- CEP55.Down %>% group_by(ONTOLOGY) %>% top_n(-10,p.adjust)
CEP55.Down2$ID <- factor(CEP55.Down2$ID,levels = rev(CEP55.Down2$ID))
p1<-ggplot(CEP55.Down2,aes(x=ID,y=-log10(p.adjust),fill=ONTOLOGY))+
  geom_bar(stat="identity",position = position_dodge(0.9)) +
  coord_flip()+
  theme_few() +
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=Description,y=0.1),hjust = 0, #vjust = 1.5,
            position = position_dodge(0.5),fontface = 'bold')+
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))+
  theme(legend.position = "top",axis.text.y = element_text(face = "bold"))+
  labs(x="")
ggsave(p1,filename = "CEP55.Down.pdf",height = 10,width = 6)
#### Figure => HALLMARK => GSEA ####
TERM <- c("HALLMARK_E2F_TARGETS",
          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
          "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_MITOTIC_SPINDLE",
          "HALLMARK_MYC_TARGETS_V1",
          "HALLMARK_MYC_TARGETS_V2",
          "HALLMARK_APOPTOSIS",
          "HALLMARK_ANGIOGENESIS",
          "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
CBX2 <- read.csv("GSEA/CBX2/CBX2.GSEA.hallmark.csv") %>%
  filter(ID %in% TERM) %>%
  dplyr::select(ID,NES,p.adjust) %>%
  mutate(SYMBOL="CBX2")

CEP55 <- read.csv("GSEA/CEP55/CEP55.GSEA.hallmark.csv")%>%
  filter(ID %in% TERM) %>%
  dplyr::select(ID,NES,p.adjust) %>%
  mutate(SYMBOL="CEP55")

HALLMARK <- rbind.data.frame(CBX2,CEP55) %>%
  arrange(p.adjust) %>%
  mutate(ID=str_remove_all(ID,"HALLMARK_")) %>%
  mutate(ID=factor(ID,levels = rev(unique(.$ID))))

p <- ggplot(HALLMARK,aes(SYMBOL,ID,color=NES,size=-log10(p.adjust)))+
  geom_point()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_color_gradientn(colors = brewer.pal(11,"RdBu")[4:1])

ggsave(p,filename = "CBX2.CEP55.GSEA.HALLMARK.pdf",height = 4,width = 7.5)

#list.files("GSEA/CBX2/",pattern = "csv$")

##################### CBX2 CEP55 => Group => Stemscore ######################
StemIndex <- read.csv("LIHC.361.csv")
write.csv(StemIndex,file = "LIHC.StemIndex.csv",row.names = F)
mRNA.Phenodata <- read.csv("../../mRNA.Phenotype.csv") %>% filter(Sample.Type == "Primary Tumor")
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000173894","ENSG00000138180"),StemIndex$SampleID] %>%
  t() %>% data.frame(check.names = F) %>% cbind(StemIndex)

mRNA.Exp.log2.TPM.Target$CBX2Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000173894 > median(mRNA.Exp.log2.TPM.Target$ENSG00000173894),"CBX2.High","CBX2.Low")
mRNA.Exp.log2.TPM.Target$CBX2Group <- factor(mRNA.Exp.log2.TPM.Target$CBX2Group,levels = c("CBX2.High","CBX2.Low"))
library(ggpubr)
library(ggthemes)

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CBX2Group", "Score", 
             fill = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("Stem score",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CBX2.stratified.StemScore.pdf",height = 2.5,width = 2.5)

mRNA.Exp.log2.TPM.Target$CEP55Group <- if_else(mRNA.Exp.log2.TPM.Target$ENSG00000138180 > median(mRNA.Exp.log2.TPM.Target$ENSG00000138180),"CEP55.High","CEP55.Low")
mRNA.Exp.log2.TPM.Target$CEP55Group <- factor(mRNA.Exp.log2.TPM.Target$CEP55Group,levels = c("CEP55.High","CEP55.Low"))

p<-ggboxplot(mRNA.Exp.log2.TPM.Target, "CEP55Group", "Score", 
             fill = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  labs(x="",y=paste("Stem score",sep = ""))+
  stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
ggsave(p,filename = "stratified/CEP55.stratified.StemScore.pdf",height = 2.5,width = 2.5)

write.csv(mRNA.Exp.log2.TPM.Target,"stratified/Stratified.StemScore.csv",row.names = F)

##################### HALLMARK SCORE => CEP55 + CBX2 ##################
HALLMARK.SCORE.LIHC <- read.csv("../LIHC-SNHG.Family/HALLMARK.ssGSEA.Pancancer/TCGA-LIHC.ssGSEA.csv",check.names = F,row.names = 1) %>% t() %>% data.frame(check.names = F)
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))

dir.create("HALLMARKS.Stratified")
HALLMARK.SCORE.LIHC.Target <- HALLMARK.SCORE.LIHC[rownames(mRNA.Exp.log2.TPM.Target),] %>% cbind(mRNA.Exp.log2.TPM.Target) %>%
  mutate(Sample=rownames(.)) %>%
  mutate(Sample.ID = substr(rownames(.),1,15)) %>% merge(.,risk_score_table_multi_cox2,by="Sample.ID")
write.csv(HALLMARK.SCORE.LIHC.Target,file = "HALLMARKS.Stratified/HALLMARK.SCORE.LIHC.Target.CBX2.CEP55.csv")


HALLMARK.SCORE.LIHC.Target$RiskScore <- factor(HALLMARK.SCORE.LIHC.Target$RiskScore,levels = c("High","Low"))
HALLMARK.SCORE.LIHC.Target$CBX2Group <- factor(HALLMARK.SCORE.LIHC.Target$CBX2Group,levels = c("CBX2.High","CBX2.Low"))
HALLMARK.SCORE.LIHC.Target$CEP55Group <- factor(HALLMARK.SCORE.LIHC.Target$CEP55Group,levels=c("CEP55.High","CEP55.Low"))

HALLMARK.SCORE.LIHC.Target <- read.csv("HALLMARKS.Stratified/HALLMARK.SCORE.LIHC.Target.CBX2.CEP55.csv")
HALLMARK.SCORE.LIHC.Target$RiskScore <- factor(HALLMARK.SCORE.LIHC.Target$RiskScore,levels = c("High","Low"))
HALLMARK.SCORE.LIHC.Target$CBX2Group <- factor(HALLMARK.SCORE.LIHC.Target$CBX2Group,levels = c("CBX2.High","CBX2.Low"))
HALLMARK.SCORE.LIHC.Target$CEP55Group <- factor(HALLMARK.SCORE.LIHC.Target$CEP55Group,levels=c("CEP55.High","CEP55.Low"))

HALLMARK.SCORE.LIHC.Target.Mid <- HALLMARK.SCORE.LIHC.Target %>% gather(Hallmark,Score,HALLMARK_TNFA_SIGNALING_VIA_NFKB:HALLMARK_PANCREAS_BETA_CELLS)
HALLMARK.SCORE.LIHC.Target.Mid$Hallmark <- HALLMARK.SCORE.LIHC.Target.Mid$Hallmark %>% str_remove("HALLMARK_")
Des.Hallmark <- c("MITOTIC_SPINDLE","WNT_BETA_CATENIN_SIGNALING",
                  "G2M_CHECKPOINT",
                  "PI3K_AKT_MTOR_SIGNALING","MTORC1_SIGNALING",
                  "E2F_TARGETS","MYC_TARGETS_V1",
                  "MYC_TARGETS_V2","EPITHELIAL_MESENCHYMAL_TRANSITION",
                  "GLYCOLYSIS","ANGIOGENESIS",
                  "APOPTOSIS",
                  "PEROXISOME",
                  "FATTY_ACID_METABOLISM","OXIDATIVE_PHOSPHORYLATION")

HALLMARK.SCORE.LIHC.Target.Mid <- HALLMARK.SCORE.LIHC.Target.Mid %>%
  filter(Hallmark %in% Des.Hallmark) %>%
  mutate(Hallmark=factor(Hallmark,levels = rev(Des.Hallmark)))

HALLMARK.SCORE.LIHC.Target.Mid2 <- HALLMARK.SCORE.LIHC.Target.Mid %>%
  dplyr::select(Hallmark,CBX2Group,CEP55Group,RiskScore,Score) %>%
  gather(Attribute,Group,CBX2Group:RiskScore) %>%
  mutate(Attribute=str_remove(Attribute,"Group")) %>%
  mutate(Group=str_remove(Group,".*\\.")) %>%
  mutate(Hallmark=factor(Hallmark,levels = rev(Des.Hallmark)))

T.W.test <- data.frame()
for (hallmark in unique(HALLMARK.SCORE.LIHC.Target.Mid2$Hallmark)) {
  for (gene in unique(HALLMARK.SCORE.LIHC.Target.Mid2$Attribute)) {
    High.g <- HALLMARK.SCORE.LIHC.Target.Mid2 %>%
      filter(Hallmark==hallmark) %>%
      filter(Attribute==gene) %>%
      filter(Group=="High")
    
    Low.g <- HALLMARK.SCORE.LIHC.Target.Mid2 %>%
      filter(Hallmark==hallmark) %>%
      filter(Attribute==gene) %>%
      filter(Group=="Low")
    
    Ttest <- t.test(as.numeric(High.g$Score),as.numeric(Low.g$Score),paired = F,alternative = "two.sided")
    Wtest <- wilcox.test(as.numeric(High.g$Score),as.numeric(Low.g$Score),paired = F,alternative = "two.sided")
    T.W.test <- data.frame(Hallmark=hallmark,
                           Attribute=gene,
                           Ttest.Pvalue=Ttest$p.value,
                           Wtest.Pvalue=Wtest$p.value
                           ) %>%
      rbind.data.frame(T.W.test)
  }
}

HALLMARK.SCORE.LIHC.Target.Mid3 <- HALLMARK.SCORE.LIHC.Target.Mid2 %>%
  group_by(Hallmark,Attribute,Group) %>%
  mutate(Score.Mean=mean(Score)) %>%
  mutate(Score.Median=median(Score)) %>%
  mutate(Score.SD=sd(Score)) %>%
  mutate(Score.CV=round(Score.SD/Score.Mean,5)) %>%
  ungroup() %>%
  dplyr::select(Hallmark,Attribute,Group,Score.Mean,Score.CV,Score.Median) %>%
  unique() %>%
  merge(T.W.test,by=c("Hallmark","Attribute"))


p4 <- ggplot()+
  geom_point(data = HALLMARK.SCORE.LIHC.Target.Mid3,
             aes(x=Group,
                 y=Hallmark,
                 size=Score.Median,
                 fill=Score.CV),
             color="white", shape=21)+
  scale_size_continuous(range = c(1,10))+
  #geom_hline(aes(yintercept = 1))+
  facet_wrap(~Attribute)+
  scale_fill_gradientn("ssGSEA score\nCV",colors = rev(c("white", "#E41A1C")))+
  #scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=HALLMARK.SCORE.LIHC.Target.Mid3%>%filter(Wtest.Pvalue<=0.05),
             aes(x=Group,
                 y=Hallmark,
                 size=Score.Median,
                 fill=Score.CV),
             color="black",shape=21,alpha=2, stroke=2)+
  scale_size_continuous("ssGSEA score\nmean",range = c(1,8))+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_boxplot()+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(#legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #ggsci::scale_fill_nejm()+
  labs(x="",y="")+
  #scale_x_log10() +
  theme(strip.text = element_text(size=15,color="black"))
  #geom_hline(aes(yintercept = 1))
  #geom_vline(aes(xintercept = 1))

ggsave(p4,filename = "HALLMARKS.Stratified/ssGSEA.CBX2.CEP55.Risk.PointPlot.pdf",height = 6,width = 10)

ggboxplot(HALLMARK.SCORE.LIHC.Target.Mid2, x = "Score", y = "Hallmark",
          color = "Group", palette = "jco",
          facet.by = "Attribute", short.panel.labs = FALSE)

p2 <- ggplot(HALLMARK.SCORE.LIHC.Target.Mid2,aes(x=Score, # %>% filter(Attribute=="CBX2")
                                           y=Hallmark,
                                           fill=Group))+
  #geom_violin()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  geom_boxplot()+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~Attribute,scales = "free_x")+
  #coord_flip()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="ssGSEA score",y="")+
  #scale_x_log10() +
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p2,filename = "HALLMARKS.Stratified/ssGSEA.CBX2.CEP55.Risk.pdf",height = 6,width = 10)




p<-ggboxplot(HALLMARK.SCORE.LIHC.Target.Mid, "Hallmark","Score", 
             color = "RiskScore",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "top")+
  labs(x="",y=paste("ssGSEA score",sep = ""))+
  stat_compare_means(aes(group = RiskScore),method = "wilcox.test",label="p.signif",hide.ns = T,label.y = 1.15,size=6)+
  coord_flip()+
  scale_y_sqrt()+ #limits = c(0,0.8)
  theme(axis.text.y = element_text(size = 12,face = "bold"),axis.text = element_text(colour = "black"),
        axis.line = element_line(size = 1))

ggsave(p,filename = "HALLMARK.SCORE.RiskScore.pdf",height = 12,width = 7)


p<-ggboxplot(HALLMARK.SCORE.LIHC.Target.Mid, "Hallmark","Score", color = "CBX2Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "top")+
  labs(x="",y=paste("ssGSEA score",sep = ""))+
  stat_compare_means(aes(group = CBX2Group),method = "wilcox.test",label="p.signif",hide.ns = T,label.y = 1.15,size=6)+
  coord_flip()+
  scale_y_sqrt()+ #limits = c(0,0.8)
  theme(axis.text.y = element_text(size = 12,face = "bold"),axis.text = element_text(colour = "black"))

ggsave(p,filename = "HALLMARK.SCORE.CBX2Group.pdf",height = 12,width = 7)

p<-ggboxplot(HALLMARK.SCORE.LIHC.Target.Mid, "Hallmark","Score", color = "CEP55Group",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "top")+
  labs(x="",y=paste("ssGSEA score",sep = ""))+
  stat_compare_means(aes(group = CEP55Group),method = "wilcox.test",label="p.signif",hide.ns = T,label.y = 1.15,size=6)+
  coord_flip()+
  scale_y_sqrt()+ #limits = c(0,0.8)
  theme(axis.text.y = element_text(size = 12,face = "bold"),axis.text = element_text(colour = "black"))

ggsave(p,filename = "HALLMARK.SCORE.CEP55Group.pdf",height = 12,width = 7)


for (name in colnames(HALLMARK.SCORE.LIHC)) {
  for (group in c("CBX2Group","CEP55Group","RiskScore")) {
    p<-ggboxplot(HALLMARK.SCORE.LIHC.Target, group, name, 
                 fill = group,palette = "Set1",add = "jitter")+
      theme_few()+
      #facet_wrap(~Immune.Cell)+
      theme(legend.position = "none")+
      labs(x="",y=name)+
      stat_compare_means(method = "t.test",label="p.format",hide.ns = T)
    ggsave(p,filename = paste("HALLMARKS.Stratified/",group,".",name,".pdf",sep = ""),height = 2.5,width = 2.5)
  }
}





##################### PD-1 PD-L1 => CBX2 + Risk score + CEP55 ###################
#### Risk score => Correlation #####
setwd("K:/TCGA/Anlysis/LIHC")
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))
Exp.data <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
Exp.data <- log2(Exp.data+1)
# PD-1 => PDCD1 => ENSG00000188389
# PD-L1 => CD274 => ENSG00000120217
# PD-L2 => CD273 => ENSG00000197646
# CTLA-4 => CTLA4 => ENSG00000163599
# TIGIT => TIGIT => ENSG00000181847
# LAG-3 => LAG3 => ENSG00000089692
# TIM-3 => HAVCR2 => ENSG00000135077
# CBX2 => CBX2 => ENSG00000173894
# CEP55 => CEP55 => ENSG00000138180
ICI.P.G.E <- data.frame(Protein = c("PD-1","PD-L1","PD-L2","CTLA-4","TIGIT","LAG-3","TIM-3","CBX2","CEP55"),
                        G.SYMBOL = c("PDCD1","CD274","CD273","CTLA4","TIGIT","LAG3","HAVCR2","CBX2","CEP55"),
                        G.ENSEMBL = c("ENSG00000188389","ENSG00000120217","ENSG00000197646",
                                      "ENSG00000163599","ENSG00000181847","ENSG00000089692",
                                      "ENSG00000135077","ENSG00000173894","ENSG00000138180"))

Exp.data2 <- Exp.data %>% rownames_to_column("ENSEMBL") %>%
  #filter(ENSEMBL %in% ICI.P.G.E$G.ENSEMBL) %>%
  merge(ICI.P.G.E,by.x="ENSEMBL",by.y="G.ENSEMBL") %>%
  dplyr::select(-ENSEMBL,-Protein) %>%
  remove_rownames() %>%
  column_to_rownames("G.SYMBOL") %>% t() %>% data.frame(check.names = F) %>%
  rownames_to_column("Sample") %>%
  mutate(Sample.ID = substr(Sample,1,15)) %>%
  merge(risk_score_table_multi_cox2,by="Sample.ID") %>% unique()
write.csv(Exp.data2,file = "ici.cbx2.cep55.expression.csv")
Exp.data2 <- Exp.data %>% rownames_to_column("ENSEMBL") %>%
  filter(ENSEMBL %in% ICI.P.G.E$G.ENSEMBL) %>% remove_rownames() %>%
  column_to_rownames("ENSEMBL") %>% t() %>% data.frame(check.names = F) %>%
  rownames_to_column("Sample") %>%
  mutate(Sample.ID = substr(Sample,1,15)) %>%
  merge(risk_score_table_multi_cox2,by="Sample.ID") %>% unique()



### Correlation => Risk Score
Cor.data <- Exp.data2 %>% gather(ENSEMBl,Expression,ENSG00000089692:ENSG00000197646) %>%
  data.frame(check.names = F) %>% merge(.,ICI.P.G.E,by.x="ENSEMBl",by.y="G.ENSEMBL") 

library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1")[1:5], 
                    brewer.pal(9, "Set1")[7:9],
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500",
                    "#FF6A6A", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00",
                    brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(8, "Set2"))

Cor.data2 <- Cor.data %>% filter(! Protein %in% c("CBX2","CEP55"))
library(ggpubr)
p<-ggplot(data = Cor.data2,aes(x=total_risk_score,y=Expression,color=Protein))+
  geom_point()+
  geom_smooth(formula = y~x,method = "lm")+
  facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Cor.data2,
           method = "pearson",label.x = 0.5,label.y = 8,size=4)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_color_manual(values = cell_type_cols)+
  labs(x="Risk score",y="log2(TPM+1)")+
  theme(strip.text.x = element_text(size=15))+
  scale_y_continuous(limits = c(0,8.1))

ggsave(p,filename = "ICI.RiskScore.Pearson.pdf",height = 2.7,width = 13)


p2<-ggplot(data = Cor.data2,aes(x=RiskScore,y=Expression,fill=RiskScore))+
  #geom_point()+
  geom_violin(cex=1.2)+           
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  geom_boxplot(width=0.1,cex=1.2,color="grey")+
  #geom_smooth(formula = y~x,method = "lm")+
  facet_wrap(~Protein,ncol = 7)+
  #stat_cor(data = Cor.data2,method = "pearson",label.x = 0.5,label.y = 8,size=4)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_fill_manual(values = cell_type_cols)+
  labs(x="",y="log2(TPM+1)")+
  stat_compare_means(comparisons = list(c("High","Low")),label = "p.value",method = "wilcox",label.y = 6.5,size=5)+
  scale_y_continuous(limits = c(0,8.1))+
  theme(strip.text.x = element_text(size=15))+
  scale_x_discrete(labels = c("High risk","Low risk"))

ggsave(p2,filename = "ICI.RiskScore.Boxplot.pdf",height = 2.7,width = 13)

######  Survival => Risk score + ICI ####
#Exp.data2
Exp.data3 <- Exp.data2 %>% dplyr::select(ICI.P.G.E$G.ENSEMBL)
colnames(Exp.data3) <- ICI.P.G.E$Protein %>% str_remove("-")
mid.data <- Exp.data2 %>% dplyr::select(total_risk_score,OS,OS.time,Sample)
colnames(mid.data)[1] = "RS"
Exp.data4 <- cbind.data.frame(Exp.data3,mid.data)

#### RS + PD1 ####
library(survminer)
library(survival)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("PD1","RS")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.PD1.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ RS + PD1,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### RS + PDL1 ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("PDL1","RS")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.PDL1.Survival.csv",row.names = F)


library(survival)
fit <- survfit(Surv(OS.time, OS) ~ RS + PDL1,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### RS + PDL2 ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("PDL2","RS")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.PDL2.Survival.csv",row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ RS + PDL2,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### RS + CTLA4 ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("RS","CTLA4")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.CTLA4.Survival.csv",row.names = F)


library(survival)
fit <- survfit(Surv(OS.time, OS) ~ RS + CTLA4,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### RS + TIGIT ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("TIGIT","RS")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.TIGIT.Survival.csv",row.names = F)


library(survival)
fit <- survfit(Surv(OS.time, OS) ~ RS + TIGIT,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### RS + LAG3 ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("LAG3","RS")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.LAG3.Survival.csv",row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ RS + LAG3,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### RS + TIM3 ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("TIM3","RS")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "ICI.RiskScore.TIM3.Survival.csv",row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ RS + TIM3,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### Survival => CBX2 + CEP55 ####
library(survminer)
cut <- surv_cutpoint(
  Exp.data4,
  time = "OS.time",
  event = "OS",
  variables = c("CBX2","CEP55")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data4$Sample)
write.csv(dat,file = "Survival.CEP55.CBX2.csv",row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ CBX2 + CEP55,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

###### Correlation => CBX2 ######
setwd("K:/TCGA/Anlysis/LIHC")
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))
Exp.data <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
Exp.data <- log2(Exp.data+1)

ICI.P.G.E <- data.frame(Protein = c("PD-1","PD-L1","CTLA-4","CBX2","CEP55"),
                        G.SYMBOL = c("PDCD1","CD274","CTLA4","CBX2","CEP55"),
                        G.ENSEMBL = c("ENSG00000188389","ENSG00000120217",
                                      "ENSG00000163599",
                                      "ENSG00000173894","ENSG00000138180"))

Exp.data2 <- Exp.data %>% rownames_to_column("ENSEMBL") %>%
  filter(ENSEMBL %in% ICI.P.G.E$G.ENSEMBL) %>% remove_rownames() %>%
  column_to_rownames("ENSEMBL") %>% t() %>% data.frame(check.names = F) %>%
  rownames_to_column("Sample") %>%
  mutate(Sample.ID = substr(Sample,1,15)) %>%
  merge(risk_score_table_multi_cox2,by="Sample.ID") %>% unique()

#### Pearson => CBX2 ####
P<- ggplot(data = Exp.data2,aes(x=ENSG00000173894,y=ENSG00000188389,color=cell_type_cols[1]))+
  geom_point()+
  geom_smooth(formula = y~x,method = "lm")+
  #facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Exp.data2,color="black",
           method = "pearson",label.x = 1,label.y = 5,size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) CD274")
  #theme(strip.text.x = element_text(size=15))+
  #scale_y_continuous(limits = c(0,8.1))
ggsave(P,filename = "ICI.CBX2.PDCD1.Pearson.pdf",height = 2,width = 2)

P2 <- ggplot(data = Exp.data2,aes(x=ENSG00000173894,y=ENSG00000120217))+
  geom_point(shape=19,color=cell_type_cols[2])+
  geom_smooth(formula = y~x,method = "lm")+
  #facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Exp.data2,color="black",
           method = "pearson",label.x = 1,label.y = 4,size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) CD274")
#theme(strip.text.x = element_text(size=15))+
#scale_y_continuous(limits = c(0,8.1))
ggsave(P2,filename = "ICI.CBX2.CD274.Pearson.pdf",height = 2,width = 2)

P3 <- ggplot(data = Exp.data2,aes(x=ENSG00000173894,y=ENSG00000163599))+
  geom_point(shape=19,color=cell_type_cols[4])+
  geom_smooth(formula = y~x,method = "lm")+
  #facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Exp.data2,color="black",
           method = "pearson",label.x = 1,label.y = 4,size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) CTLA4")
#theme(strip.text.x = element_text(size=15))+
#scale_y_continuous(limits = c(0,8.1))
ggsave(P3,filename = "ICI.CBX2.CTLA4.Pearson.pdf",height = 2,width = 2)

#### Pearson => CEP55 ####
P<- ggplot(data = Exp.data2,aes(x=ENSG00000138180,y=ENSG00000188389,
                                color=cell_type_cols[1]))+
  geom_point()+
  geom_smooth(formula = y~x,method = "lm")+
  #facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Exp.data2,color="black",
           method = "pearson",label.x = 1,label.y = 5,size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) PDCD1")
#theme(strip.text.x = element_text(size=15))+
#scale_y_continuous(limits = c(0,8.1))
ggsave(P,filename = "ICI.CEP55.PDCD1.Pearson.pdf",height = 2,width = 2)

P2 <- ggplot(data = Exp.data2,aes(x=ENSG00000138180,y=ENSG00000120217))+
  geom_point(shape=19,color=cell_type_cols[2])+
  geom_smooth(formula = y~x,method = "lm")+
  #facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Exp.data2,color="black",
           method = "pearson",label.x = 1,label.y = 4,size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) CD274")
#theme(strip.text.x = element_text(size=15))+
#scale_y_continuous(limits = c(0,8.1))
ggsave(P2,filename = "ICI.CEP55.CD274.Pearson.pdf",height = 2,width = 2)

P3 <- ggplot(data = Exp.data2,aes(x=ENSG00000138180,y=ENSG00000163599))+
  geom_point(shape=19,color=cell_type_cols[4])+
  geom_smooth(formula = y~x,method = "lm")+
  #facet_wrap(~Protein,ncol = 7)+
  stat_cor(data = Exp.data2,color="black",
           method = "pearson",label.x = 1,label.y = 4,size=3)+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="log2(TPM+1) CEP55",y="log2(TPM+1) CTLA4")
#theme(strip.text.x = element_text(size=15))+
#scale_y_continuous(limits = c(0,8.1))
ggsave(P3,filename = "ICI.CEP55.CTLA4.Pearson.pdf",height = 2,width = 2)


#### ICI CBX2 => Survival ####
library(survminer)
library(survival)
ICI.P.G.E <- data.frame(Protein = c("PD-1","PD-L1","CTLA-4","CBX2","CEP55"),
                        G.SYMBOL = c("PDCD1","CD274","CTLA4","CBX2","CEP55"),
                        G.ENSEMBL = c("ENSG00000188389","ENSG00000120217",
                                      "ENSG00000163599",
                                      "ENSG00000173894","ENSG00000138180"))
colnames(Exp.data2)[3:7] <- c("CD274","CEP55","CTLA4","CBX2","PDCD1")
cut <- surv_cutpoint(
  Exp.data2,
  time = "OS.time",
  event = "OS",
  variables = c("PDCD1","CBX2")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data2$Sample)
write.csv(dat,file = "ICI.CBX2.PD1.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ CBX2 + PDCD1,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

##
cut <- surv_cutpoint(
  Exp.data2,
  time = "OS.time",
  event = "OS",
  variables = c("CD274","CBX2")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data2$Sample)
write.csv(dat,file = "ICI.CBX2.CD274.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ CBX2 + CD274,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

##
cut <- surv_cutpoint(
  Exp.data2,
  time = "OS.time",
  event = "OS",
  variables = c("CTLA4","CBX2")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data2$Sample)
write.csv(dat,file = "ICI.CBX2.CTLA4.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ CBX2 + CTLA4,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

#### ICI CEP55 => Survival ####
library(survminer)
library(survival)
ICI.P.G.E <- data.frame(Protein = c("PD-1","PD-L1","CTLA-4","CBX2","CEP55"),
                        G.SYMBOL = c("PDCD1","CD274","CTLA4","CBX2","CEP55"),
                        G.ENSEMBL = c("ENSG00000188389","ENSG00000120217",
                                      "ENSG00000163599",
                                      "ENSG00000173894","ENSG00000138180"))
colnames(Exp.data2)[3:7] <- c("CD274","CEP55","CTLA4","CBX2","PDCD1") 
cut <- surv_cutpoint(
  Exp.data2,
  time = "OS.time",
  event = "OS",
  variables = c("PDCD1","CEP55")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data2$Sample)
write.csv(dat,file = "ICI.CEP55.PD1.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ CEP55 + PDCD1,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

##
cut <- surv_cutpoint(
  Exp.data2,
  time = "OS.time",
  event = "OS",
  variables = c("CD274","CEP55")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data2$Sample)
write.csv(dat,file = "ICI.CEP55.CD274.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ CEP55 + CD274,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)

##
cut <- surv_cutpoint(
  Exp.data2,
  time = "OS.time",
  event = "OS",
  variables = c("CTLA4","CEP55")
)
dat <- surv_categorize(cut) %>% data.frame() %>% mutate(Sample=Exp.data2$Sample)
write.csv(dat,file = "ICI.CEP55.CTLA4.Survival.csv",row.names = F)
theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

library(survival)
library(RColorBrewer)
fit <- survfit(Surv(OS.time, OS) ~ CEP55 + CTLA4,
               data = dat)
ggsurvplot(
  fit,
  risk.table = F,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0,4000),
  break.time.by = 1000,
  xlab = "Time (days)",
  palette = brewer.pal(8,"Set1"),
  ggtheme = theme.sur,
  risk.table.y.text.col = T,
  risk.table.y.text = F
)


########

##################### ICB => Response + Non-response ###########
setwd("K:/TCGA/Anlysis/LIHC")
data1 <- read.csv("TIGER.ICB.csv")
data2 <- read.csv("TIGER.ICB.PrePost.csv")

intersect(data1$Dataset,data2$Dataset)

setwd("K:/TCGA/TIGER.Database")
#intersect(data1$Dataset,data2$Dataset)
#[1] "Melanoma-GSE106128_DCs_treated"            "Melanoma-GSE115821_ALL"                   
#[3] "Melanoma-GSE115821_anti-CTLA-4"            "Melanoma-GSE91061_anti-PD-1"              
#[5] "Melanoma-Nathanson_2017_anti-CTLA-4"       "Melanoma-PRJEB23709_ALL"                  
#[7] "Melanoma-PRJEB23709_anti-CTLA-4+anti-PD-1" "Melanoma-PRJEB23709_anti-PD-1"
#### dataset1 => Melanoma-GSE115821_ALL #### 
data <- read.csv("Melanoma-GSE115821.Response.tsv",sep = "\t")
Clinical <- read.csv("Melanoma-GSE115821.Response (1).tsv",sep = "\t")
Clinical.Pre <- Clinical %>% filter(Treatment=="PRE")
#### dataset2 => Melanoma-GSE91061_anti-PD-1 #### 
data <- read.csv("Melanoma-GSE91061.Response.tsv",sep = "\t")
Clinical <- read.csv("Melanoma-GSE91061.Response (1).tsv",sep = "\t")
Clinical.Pre <- Clinical %>% filter(Treatment=="PRE") %>%
  filter(response != "UNK")

CBX2.CEP55.Exp <- data %>% filter(GENE_SYMBOL %in% c("CBX2","CEP55")) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical.Pre,by="sample_id") %>%
  mutate(Response=if_else(response %in% c("PR","CR"),"Yes","No")) %>%
  mutate(CBX2.Group = if_else(CBX2 > median(CBX2),"High","Low")) %>%
  mutate(CEP55.Group = if_else(CEP55 > median(CEP55),"High","Low"))

chisq.test(table(CBX2.CEP55.Exp$CBX2.Group,CBX2.CEP55.Exp$Response))

CBX2.Data <- CBX2.CEP55.Exp %>%
  dplyr::select(CBX2.Group,Response) %>%
  group_by(CBX2.Group,Response) %>%
  summarise(Count=n())

p <- ggplot(CBX2.Data,aes(x=CBX2.Group,y=Count,fill=Response))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_x_discrete(breaks=c("High","Low"),labels=c("CBX2 high","CBX2 low"))+
  ggsci::scale_fill_npg()+
  labs(x="",y="No. of samples")
ggsave(p,filename = "TIGER.CBX2.Melanoma-GSE91061_anti-PD-1.pdf",width = 4,height = 4)

CEP55.Data <- CBX2.CEP55.Exp %>%
  dplyr::select(CEP55.Group,Response) %>%
  group_by(CEP55.Group,Response) %>%
  summarise(Count=n())

p2 <- ggplot(CEP55.Data,aes(x=CEP55.Group,y=Count,fill=Response))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_x_discrete(breaks=c("High","Low"),labels=c("CEP55 high","CEP55 low"))+
  ggsci::scale_fill_npg()+
  labs(x="",y="No. of samples")
ggsave(p2,filename = "TIGER.CEP55.Melanoma-GSE91061_anti-PD-1.pdf",width = 4,height = 4)


chisq.test(table(CBX2.CEP55.Exp$CEP55.Group,CBX2.CEP55.Exp$Response))


#### dataset3 => Melanoma-PRJEB23709_anti-PD-1 #### 
data <- read.csv("Melanoma-PRJEB23709.Response.tsv",sep = "\t")
Clinical <- read.csv("Melanoma-PRJEB23709.Response (1).tsv",sep = "\t")
Clinical.Pre <- Clinical %>% filter(Treatment=="PRE") %>%
  filter(response != "UNK") %>% filter(Therapy=="anti-PD-1")

CBX2.CEP55.Exp <- data %>% filter(GENE_SYMBOL %in% c("CBX2","CEP55")) %>%
  remove_rownames() %>% column_to_rownames("GENE_SYMBOL") %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("sample_id") %>%
  merge(Clinical.Pre,by="sample_id") %>%
  mutate(Response=if_else(response %in% c("PR","CR"),"Yes","No")) %>%
  mutate(CBX2.Group = if_else(CBX2 > median(CBX2),"High","Low")) %>%
  mutate(CEP55.Group = if_else(CEP55 > median(CEP55),"High","Low"))

chisq.test(table(CBX2.CEP55.Exp$CBX2.Group,CBX2.CEP55.Exp$Response))

CBX2.Data <- CBX2.CEP55.Exp %>%
  dplyr::select(CBX2.Group,Response) %>%
  group_by(CBX2.Group,Response) %>%
  summarise(Count=n())

p <- ggplot(CBX2.Data,aes(x=CBX2.Group,y=Count,fill=Response))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_x_discrete(breaks=c("High","Low"),labels=c("CBX2 high","CBX2 low"))+
  ggsci::scale_fill_npg()+
  labs(x="",y="No. of samples")
ggsave(p,filename = "TIGER.CBX2.Melanoma-PRJEB23709_anti-PD-1.pdf",width = 4,height = 4)

CEP55.Data <- CBX2.CEP55.Exp %>%
  dplyr::select(CEP55.Group,Response) %>%
  group_by(CEP55.Group,Response) %>%
  summarise(Count=n())

p2 <- ggplot(CEP55.Data,aes(x=CEP55.Group,y=Count,fill=Response))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_x_discrete(breaks=c("High","Low"),labels=c("CEP55 high","CEP55 low"))+
  ggsci::scale_fill_npg()+
  labs(x="",y="No. of samples")
ggsave(p2,filename = "TIGER.CEP55.Melanoma-PRJEB23709_anti-PD-1.pdf",width = 4,height = 4)


chisq.test(table(CBX2.CEP55.Exp$CEP55.Group,CBX2.CEP55.Exp$Response))


##################### TCGA CBX2 CEP55 => DNA methylation => MEXPRESS ##########################
setwd("K:/TCGA/Clinical.XENA/LIHC")
list.files(".")
GPL <- read.csv("../../../HepG2-CBX2/Validate/GPL13534-11288.txt",sep = "\t") %>%
  filter(str_detect(UCSC_RefGene_Name,"CBX2") | str_detect(UCSC_RefGene_Name,"CEP55")) %>%
  dplyr::select(ID,UCSC_RefGene_Name,UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island)

for (i in 1:dim(GPL)[1]) {
  GPL$UCSC_RefGene_Name[i] = (GPL$UCSC_RefGene_Name[i] %>% str_split(";"))[[1]] %>% unique()
  GPL$UCSC_RefGene_Group[i] = (GPL$UCSC_RefGene_Group[i] %>% str_split(";"))[[1]] %>% unique() %>%
    paste0(collapse = ";")
  GPL$Relation_to_UCSC_CpG_Island[i] = (GPL$Relation_to_UCSC_CpG_Island[i] %>% str_split(";"))[[1]] %>% unique() %>%
    paste0(collapse = ";")
}
write.table(GPL,file = "CBX2.CEP55.MethylationProbe.Anno.txt",sep = "\t",row.names = F)

# data <- read.csv("TCGA-LIHC.methylation450.tsv",sep = "\t")

setwd("K:/HepG2-CBX2/Validate")

data <- read.csv("CBX2.CEP55.Probe.TCGA.txt",sep = "\t",check.names = F,row.names = 1) %>%
  na.omit()

Clinical <- read.csv("../../../TCGA/Clinical.XENA/XENA.TCGA.LIHC_survival.txt",sep = "\t") %>%
  dplyr::select(-sample,-Redaction) %>% unique(); dim(Clinical)

data2 <- data %>% t() %>% data.frame(check.names = F) %>% rownames_to_column("SampleID") %>%
  mutate(X_PATIENT=substr(SampleID,1,12)) %>%
  merge(Clinical,by="X_PATIENT")
write.csv(data2,file = "CBX2.CEP55.LIHC.Methylation.Survival.Data.csv",row.names = F)

data2.Meth <- data2[,2:43] %>% data.frame(check.names = F) %>%
  mutate(Group=if_else(str_detect(SampleID,"-11A"),"NTL","HCC"))

Meth.T.W.test <- data.frame()
for (ProbeID in Anno.data$ID) {
  NTL.data <- (data2.Meth %>% filter(Group=="NTL"))[[ProbeID]]
  HCC.data <- (data2.Meth %>% filter(Group=="HCC"))[[ProbeID]]
  TTest <- t.test(NTL.data,HCC.data,alternative = "two.sided",paired = F)
  WTest <- wilcox.test(NTL.data,HCC.data,alternative = "two.sided",paired = F)
  Meth.T.W.test <- data.frame(ProbeID=ProbeID,
                              NTL.mean=mean(NTL.data),
                              HCC.mean=mean(HCC.data),
                              T.pvalue=TTest$p.value,
                              W.pvalue=WTest$p.value) %>%
    rbind.data.frame(Meth.T.W.test)
}
write.csv(Meth.T.W.test,file = "TCGA.CBX2.CEP55.Methylation.DE.csv",row.names = F)

Meth.T.W.test$ProbeID <- factor(Meth.T.W.test$ProbeID,levels = Anno.data$ID)
pDE <- ggplot()+
  geom_segment(data=Meth.T.W.test,aes(x=HCC.mean, xend=NTL.mean, y=ProbeID, yend=ProbeID), 
               size=0.7,linetype=1,color="grey")+
  geom_point(data=Meth.T.W.test,aes(x=HCC.mean,y=ProbeID),size = 2, pch = 21, fill=brewer.pal(8,"Set1")[1],color="white")+
  geom_point(data=Meth.T.W.test %>% filter(T.pvalue <= 0.05),aes(x=HCC.mean,y=ProbeID),size = 2, pch = 21, fill=brewer.pal(8,"Set1")[1],color="black",stroke=1)+
  geom_point(data=Meth.T.W.test,aes(x=NTL.mean,y=ProbeID),size = 2, pch = 21, fill=brewer.pal(8,"Set1")[2],color="white")+
  geom_point(data=Meth.T.W.test %>% filter(T.pvalue <= 0.05),aes(x=NTL.mean,y=ProbeID),size = 2, pch = 21, fill=brewer.pal(8,"Set1")[2],color="black",stroke=1)+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Methylation level",y="")+
  scale_x_continuous(breaks = c(0,0.5,1))

Meth.T.W.test2 <- Meth.T.W.test %>% mutate(Label=if_else(T.pvalue>=0.05,"",
                                                         if_else(T.pvalue>=0.01,"*",
                                                                 if_else(T.pvalue>=0.001,"**","***"))))
pDE.P<-ggplot(data=Meth.T.W.test2,aes(x=T.pvalue,y=ProbeID))+
  geom_text(aes(label=Label),hjust = 0,size=5)+
  scale_x_continuous(expand = c(0,0))+
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank())

p117 <- pDE %>% aplot::insert_right(pDE.P,width = 0.1)

#### KM ####

log_rank_p <- apply(lncRNA.Down.Data , 2 , function(gene){
  # gene=exprSet[1,]
  mid.Survival$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(Surv(OS.time, OS)~group,data=mid.Survival)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

data2.2 <- data2 %>% filter(!str_detect(SampleID,"-11A"))

CBX2.CEP55.Methylation.Survival <- data.frame()
for (probe in colnames(data2.2)[3:43]) {
  for (index in c("OS","DFI","PFI","DSS")) {
    data3 <- data2.2 %>% dplyr::select(c(probe,index,paste(index,".time",sep = ""))) %>% na.omit()
    data3$group <- ifelse(data3[[probe]] > median(data3[[probe]]),'high','low')
    formula <- as.formula(paste("Surv(",paste(index,".time",sep = ""),",",index,")~group",sep = ""))
    data.survdiff=survdiff(formula,data=data3)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    
    formula2 <- as.formula(paste("Surv(",paste(index,".time",sep = ""),",",index,")~",probe,sep = ""))
    m = coxph(formula2, data =  data3)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    
    #summary(m)
    tmp <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)), 6)
    
    CBX2.CEP55.Methylation.Survival <- data.frame(Probe_id=probe,
                                                  SurvivalIndex=index,
                                                  KM.Pvalue=p.val,
                                                  HR = HR,
                                                  COX.Pvalue=1 - pchisq((beta/se)^2, 1),
                                                  HRz = (HR - 1) / HRse, 
                                                  HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                                  HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                                                  HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)) %>%
      rbind.data.frame(CBX2.CEP55.Methylation.Survival)
  }
}

write.csv(CBX2.CEP55.Methylation.Survival,file = "CBX2.CEP55.Methylation.Survival.Results.csv",row.names = F)
setwd("K:/TCGA/Clinical.XENA/LIHC")
Anno.data <- read.csv("CBX2.CEP55.MethylationProbe.Anno.txt",sep = "\t") %>%
  mutate(ID=factor(ID,levels = .$ID))
setwd("K:/TCGA/Anlysis/LIHC")
CBX2.CEP55.Methylation.Survival <- read.csv("CBX2.CEP55.Methylation.Survival.Results.csv")
Anno.data <- Anno.data %>% filter(ID %in% CBX2.CEP55.Methylation.Survival$Probe_id)
CBX2.CEP55.Methylation.Survival <- CBX2.CEP55.Methylation.Survival %>%
  mutate(Probe_id=factor(Probe_id,levels = Anno.data$ID)) %>%
  mutate(SurvivalIndex=factor(SurvivalIndex,levels = c("OS","DSS","PFI","DFI")))

p111 <- ggplot()+
  geom_point(data=CBX2.CEP55.Methylation.Survival,aes(x=SurvivalIndex,y=Probe_id,
                                                      size=-log10(KM.Pvalue),
                                                      fill=log10(HR)),shape=21,color="white") +
  geom_point(data=CBX2.CEP55.Methylation.Survival %>% filter(KM.Pvalue<=0.05),
             aes(x=SurvivalIndex,y=Probe_id,
                 size=-log10(KM.Pvalue),
                 fill=log10(HR)),shape=21,color="black")+
  scale_size_continuous(range = c(1,6))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  #scale_fill_gradientn(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0,limits=c(-15,15))
  scale_fill_gradientn(colours = c(brewer.pal(9,"YlGn")[c(5,6,7,8,9,9)],"white",
                                   brewer.pal(9,"YlOrRd")[c(4,4,5,5,6,6,7,7,8,8,8,9,9,9,9,9,9,9,9,9)]))#
ggsave(p111,filename = "TCGA.CBX2.CEP55.Methylation.Survival.Figure.pdf",height = 9,width = 5)

p112 <- ggplot(Anno.data,aes(1,ID,fill=UCSC_RefGene_Group))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank())+
  labs(x="",y="")+
  ggsci::scale_fill_nejm()

Anno.data$Relation_to_UCSC_CpG_Island <- as.character(Anno.data$Relation_to_UCSC_CpG_Island)

p113 <- ggplot(Anno.data,aes(1,ID,fill=Relation_to_UCSC_CpG_Island))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank())+
  labs(x="",y="")+
  ggsci::scale_fill_jama()

p110 <- ggplot(Anno.data,aes(1,ID,fill=UCSC_RefGene_Name))+
  geom_col()+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank())+
  labs(x="",y="")+
  ggsci::scale_fill_npg()

#### Expression correlation ####
Exp <- read.csv("../../../TCGA/mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
Exp.Log <- log2(Exp+1)
Phenodata <- read.csv("../../../TCGA/mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC")

mid.phenodata <- Phenodata %>% filter(Sample.ID %in% colnames(Exp)) %>% 
  filter(Sample.Type == "Solid Tissue Normal" | Sample.Type == "Primary Tumor") %>%
  mutate(Group = case_when(Sample.Type == "Solid Tissue Normal" ~ "NTL",
                           Sample.Type == "Primary Tumor" ~ "HCC")) %>%
  #mutate(Group = factor(.$Group,levels = c("NTL","HCC"))) %>%
  dplyr::select(Case.ID:Group)
#GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
#GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180
Exp.CBX2.CEP55 <- Exp.Log[c("ENSG00000173894","ENSG00000138180"),] %>%
  t() %>% data.frame(check.names = F) %>%
  rownames_to_column("SampleID") %>%
  merge(data2,by="SampleID")
write.csv(Exp.CBX2.CEP55,file = "TCGA.CBX2.CEP55.Methylation.Expression.Correlation.csv",row.names = F)

Cor.data <- data.frame()
for(gene in c("ENSG00000173894","ENSG00000138180")){
  if (gene == "ENSG00000173894") {
    symbol="CBX2"
  }else{
    symbol="CEP55"
  }
  
  Symbol.Probe <- Anno.data %>% filter(UCSC_RefGene_Name == symbol)
  for (probe in Symbol.Probe$ID) {
    test <- cor.test(Exp.CBX2.CEP55[[gene]],Exp.CBX2.CEP55[[probe]],
                     method = "spearman",exact = F)
    Cor.data <- data.frame(Symbol = symbol,
                           ProbeID=probe,
                           Rho=test$estimate,
                           Pvalue=test$p.value) %>%
      rbind.data.frame(Cor.data)
  }
}
write.csv(Cor.data,file = "TCGA.CBX2.CEP55.Methylation.Expression.Correlation.Results.csv",row.names = F)

Cor.data2 <- Cor.data %>% remove_rownames() %>% 
  mutate(ProbeID=factor(ProbeID,levels = Anno.data$ID))

p114 <- ggplot()+
  geom_point(data=Cor.data2,aes(x=Symbol,y=ProbeID,
                                size=-log10(Pvalue),
                                fill=Rho),shape=21,color="white") +
  geom_point(data=Cor.data2 %>% filter(Pvalue<=0.05),
             aes(x=Symbol,y=ProbeID,
                 size=-log10(Pvalue),
                 fill=Rho),shape=21,color="black")+
  scale_size_continuous(range = c(1,6))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0,limits=c(-.5,.5))
#scale_fill_gradientn(colours =rainbow(6))

### CBX2 CEP55 survival
Exp.CBX2.CEP55.2 <- Exp.CBX2.CEP55 %>% 
  filter(!str_detect(SampleID,"-11A"))
CBX2.CEP55.Expression.Survival <- data.frame()
for (gene in c("ENSG00000173894","ENSG00000138180")){
  if (gene == "ENSG00000173894") {
    symbol="CBX2"
  }else{
    symbol="CEP55"
  }
  for (index in c("OS","DFI","PFI","DSS")) {
    data3 <- Exp.CBX2.CEP55.2 %>% dplyr::select(c(gene,index,paste(index,".time",sep = ""))) %>% na.omit()
    data3$group <- ifelse(data3[[gene]] > median(data3[[gene]]),'high','low')
    formula <- as.formula(paste("Surv(",paste(index,".time",sep = ""),",",index,")~group",sep = ""))
    data.survdiff=survdiff(formula,data=data3)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    
    formula2 <- as.formula(paste("Surv(",paste(index,".time",sep = ""),",",index,")~",gene,sep = ""))
    m = coxph(formula2, data =  data3)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    
    #summary(m)
    tmp <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)), 6)
    
    CBX2.CEP55.Expression.Survival <- data.frame(Symbol=symbol,
                                                 SurvivalIndex=index,
                                                 KM.Pvalue=p.val,
                                                 HR = HR,
                                                 COX.Pvalue=1 - pchisq((beta/se)^2, 1),
                                                 HRz = (HR - 1) / HRse, 
                                                 HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                                                 HRCI.LL = exp(beta - qnorm(.975, 0, 1) * se),
                                                 HRCI.UL = exp(beta + qnorm(.975, 0, 1) * se)) %>%
      rbind.data.frame(CBX2.CEP55.Expression.Survival)
  }
}
write.csv(CBX2.CEP55.Expression.Survival,file = "TCGA.CBX2.CEP55.Expression.Survival.Data.csv",row.names = F)

p115 <- ggplot()+
  geom_point(data=CBX2.CEP55.Expression.Survival,aes(x=SurvivalIndex,y=Symbol,
                                                     size=-log10(KM.Pvalue),
                                                     fill=HR),shape=21,color="white") +
  geom_point(data=CBX2.CEP55.Expression.Survival %>% filter(KM.Pvalue<=0.05),
             aes(x=SurvivalIndex,y=Symbol,
                 size=-log10(KM.Pvalue),
                 fill=HR),shape=21,color="black")+
  scale_size_continuous(range = c(4,6))+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colors = brewer.pal(9,"YlOrRd")[5:9])



p116 <- p111 %>% 
  aplot::insert_right(pDE,width=.7) %>%
  aplot::insert_right(pDE.P,width=0.1) %>%
  aplot::insert_right(p110,width=.1) %>%
  aplot::insert_right(p112,width=.1) %>% 
  aplot::insert_right(p113,width=.1) %>%
  aplot::insert_right(p114,width=.4) %>% 
  aplot::insert_bottom(p115,height = 0.07)

ggsave(p116,filename = "TCGA.CBX2.CEP55.Methylation.Expression.Correlation.Sruvival2.pdf",height = 20,width = 10)

# 371+50
HCC.Samples <- mid.phenodata %>% filter(Group=="HCC")
NTL.Samples <- mid.phenodata %>% filter(Group=="NTL")

HCC.M1 <- data %>% dplyr::select(intersect(HCC.Samples$Sample.ID,colnames(.))) %>%
  magrittr::set_colnames(substr(colnames(.),1,12)); dim(HCC.M1)
NTL.M1 <- data %>% dplyr::select(intersect(NTL.Samples$Sample.ID,colnames(.))); dim(NTL.M1)

meta <- read.csv("../../../TCGA/Anlysis/LIHC/clinical.cart.2021-08-19/clinical.tsv",sep = "\t",header = T)
#meta <- column_to_rownames(meta,var = "case_submitter_id")

Meta.data=meta[,colnames(meta) %in% c("case_submitter_id",
                                      "race",
                                      "gender",
                                      "ajcc_pathologic_t",
                                      "ajcc_pathologic_n",
                                      "ajcc_pathologic_m",
                                      "age","Stage")] %>% unique() %>%
  mutate(Stage=str_remove(Stage,"Stage ")) %>%
  magrittr::set_colnames(c("PATIENT","Age","Gender","Stage","AJCC_M","AJCC_N","AJCC_T")) %>% 
  merge(Clinical,by.x="PATIENT",by.y="X_PATIENT")

HCC.M1 <- HCC.M1 %>% dplyr::select(Meta.data$PATIENT)




##################### TCGA CBX2 CEP55 => CNV => Expression => Survival #######
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
library(readxl)
library(ggplot2)

CBX2CNV <- read_xlsx("CNV.CBX2.CEP55/CBX2.CnvAndExpressionTable.xlsx")
CEP55CNV <- read_xlsx("CNV.CBX2.CEP55/CEP55.CnvAndExpressionTable.xlsx")

CNV <- rbind(CBX2CNV,CEP55CNV) %>% data.frame()
p1 <- ggplot()+
  geom_point(data=CNV,aes(symbol,symbol,fill=spm,size=-log10(fdr)),
             shape=21,color="white")+
  geom_point(data=CNV %>% filter(fdr<=0.05),aes(symbol,symbol,fill=spm,size=-log10(fdr)),
             shape=21,color="black",stroke=1.5)+
  scale_size_continuous(range = c(4,6))+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_fill_gradient2("Spearman'CC",low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2))
  labs(y="CNV",x="Expression")

CBX2.CNV.S <- read_xlsx("CNV.CBX2.CEP55/CBX2.CnvAndSurvivalTable.xlsx")
CEP55.CNV.S <- read_xlsx("CNV.CBX2.CEP55/CEP55.CnvAndSurvivalTable.xlsx")

CNV.S <- rbind.data.frame(CBX2.CNV.S,CEP55.CNV.S) %>%
  mutate(sur_type=factor(sur_type,levels=c("OS","DSS","DFI","PFS")))

p2 <- ggplot()+
  geom_point(data=CNV.S,
             aes(sur_type,symbol,fill=-log10(log_rank_p),size=-log10(log_rank_p)),
             shape=21,color="white")+
  geom_point(data=CNV.S %>% filter(log_rank_p<=0.05),
             aes(sur_type,symbol,fill=-log10(log_rank_p),size=-log10(log_rank_p)),
             shape=21,color="black",stroke=1.5)+
  scale_size_continuous(range = c(4,6))+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_fill_gradient2("-log10(KM P)",low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)
  
CBX2.CNV.T <- read_xlsx("CNV.CBX2.CEP55/CBX2.CnvSummaryTable.xlsx")
CEP55.CNV.T <- read_xlsx("CNV.CBX2.CEP55/CEP55.CnvSummaryTable.xlsx")
CNV.T <- rbind.data.frame(CBX2.CNV.T,CEP55.CNV.T) %>%
  mutate(None=100-a_total-d_total) %>% 
  magrittr::set_colnames(c("Cancer","Symbol","Amp","Del","Hete.Amp",
                           "Hete.Del","Homo.Amp","Homo.Del","Enterz","None"))

CNV.T1 <- CNV.T %>% dplyr::select(-Enterz,-Amp.Total,-Del.Total) %>%
  gather(Type,Percent,Hete.Amp:None) %>%
  mutate(Type=factor(Type,levels = rev(c("Hete.Amp","Hete.Del",
                                     "Homo.Amp","Homo.Del","None"))))
write.csv(CNV.T,file="CNV.CBX2.CEP55/CBX2.CEP55.CNV.Type.csv",row.names = F)
p3 <- ggplot(CNV.T1,aes(Percent,Symbol,fill=Type))+
  geom_col()+
  ggthemes::theme_few()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="CNV percent",y="")+
  scale_fill_manual(values = brewer.pal(8,"Set1")[c(9,2,4,5,1)])+
  scale_x_continuous(expand = c(0,0))
 
CNV.T2 <- CNV.T %>% dplyr::select(Symbol,Amp,Del,None) %>%
  gather(Type,Percent,Amp:None) %>%
  mutate(Type=factor(Type,levels =c("None","Del","Amp")))

p4 <- ggplot(CNV.T2,aes(Percent,Symbol,fill=Type))+
  geom_col()+
  ggthemes::theme_few()+
  theme(legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="CNV percent",y="")+
  scale_fill_manual(values = brewer.pal(8,"Set1")[c(9,7,8)])+
  scale_x_continuous(expand = c(0,0))

p6 <- p2 %>% aplot::insert_left(p1,width = 0.5) %>%
  aplot::insert_right(p3,width = 0.9) %>%
  aplot::insert_right(p4,width = 0.9)
ggsave(p6,filename = "CNV.CBX2.CEP55/CBX2.CEP.CNV2.pdf",height = 12,width = 14)

#### CNV => Expression ####
library(data.table)
CNV.D <- fread("../../GSCA/LIHC.cnv.tsv.gz")
CBX2.CNV <- CNV.D %>% filter(symbol == "CBX2") %>%
  dplyr::select(-cancer_type,-entrez) %>%
  remove_rownames() %>%
  column_to_rownames("symbol") %>%
  t() %>%
  data.frame(check.names = F) %>%
  magrittr::set_rownames(substring(rownames(.),1,16)) %>%
  mutate(Type=if_else(CBX2==1,"Amp.Hete",
                      if_else(CBX2==2,"Amp.Homo",
                              if_else(CBX2==-1,"Del.Hete",
                                      if_else(CBX2==-2,"Del.Home","None"))))) %>%
  rownames_to_column("SampleID")

TPM.data <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",
                     row.names = 1,check.names = F)
TPM.data <- log2(TPM.data+1)
CBX2.TPM <- TPM.data["ENSG00000173894",] %>% t() %>%
  data.frame(check.names = F) %>%
  rownames_to_column("SampleID") %>%
  merge(CBX2.CNV,by="SampleID")

library(ggpubr)

p12 <- ggplot(CBX2.TPM,aes(Type,ENSG00000173894,fill=Type))+
  geom_violin(cex=1.2)+           
  geom_boxplot(width=0.3,cex=1.2,color="white")+
  #ggbeeswarm::geom_beeswarm(priority = "descending")
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        #axis.ticks.y = element_blank(),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="log2(TPM+1) CBX2")+
  #stat_compare_means(label = "p.signif")+
  stat_compare_means(label.x = 2)+
  ggsci::scale_fill_npg()

ggsave(p12,filename = "CBX2.CNV.Expression.pdf",height = 3,width = 4)


##################### TCGA CBX2 CEP55 => Expression => Cytokine #######
mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

Pheno <- read.csv("../../mRNA.Phenotype.csv") %>% 
  filter(Project.ID=="TCGA-LIHC") %>% filter(Sample.Type == "Primary Tumor")

Cor.data <- mRNA.Exp.log2.TPM %>% dplyr::select(Pheno$Sample.ID) %>% 
  rownames_to_column("ENSEMBL")

GeneType <- read.table("../../../../Metastasis/TissueMicrobiome/Human.GRC38.GeneType.txt",header = T,sep = ",")
GeneType.Unique <- GeneType %>% dplyr::select(Gene.stable.ID,Gene.name) %>% unique()

Immune.regulators <- read_xlsx("../../Immunomodulators.xlsx")

Immune.regulators.New <- data.frame()
for (i in 1:nrow(Immune.regulators)) {
  genelist <- Immune.regulators$Gene[i] %>% str_split(",") %>% unlist()
  for (gene in genelist) {
    Immune.regulators.New <- data.frame(Type=Immune.regulators$Category[i],
                                        Gene=gene) %>%
      rbind.data.frame(Immune.regulators.New)
  }
}
Immune.regulators.New$Gene <- Immune.regulators.New$Gene %>% str_remove(" ")
write.csv(Immune.regulators.New,file = "../../Immunomodulators.csv",row.names = F)

Des.Gene <- GeneType.Unique %>% data.frame() %>%
  filter(Gene.name %in% c(Immune.regulators.New$Gene,"CBX2","CEP55"))

table(Des.Gene$Gene.name)

Cor.data2 <- Cor.data %>% merge(Des.Gene,by.x="ENSEMBL",by.y="Gene.stable.ID")

Pearson.data <- data.frame()
CBX2.CEP55.DATA <- Cor.data2 %>% filter(Gene.name %in% c("CBX2","CEP55"))
Other.DATA <- Cor.data2 %>% filter(!Gene.name %in% c("CBX2","CEP55"))
for (gene in c("CBX2","CEP55")) {
  for (gene2 in Other.DATA$Gene.name) {
    test <- cor.test(as.numeric((CBX2.CEP55.DATA %>% filter(Gene.name == gene))[2:372]),
                     as.numeric((Other.DATA %>% filter(Gene.name == gene2))[2:372]),
                     exact = F,method = "pearson")
    Pearson.data <- data.frame(Gene1=gene,
                               Gene2=gene2,
                               Rho=test$estimate,
                               Pvalue=test$p.value) %>%
      rbind.data.frame(Pearson.data)
  }
}
Pearson.data2 <- Pearson.data %>% na.omit() %>% 
  merge(Immune.regulators.New,by.x="Gene2",by.y="Gene")
write.csv(Pearson.data2,file = "CBX2.CEP55.Chemokine.Pearson.csv",row.names = F)

Pearson.data2 <- read.csv(file = "CBX2.CEP55.Chemokine.Pearson.csv")
library(tidyverse)
library(pheatmap)
Pearson.data.Pvalue <- Pearson.data2 %>% arrange(Type) %>%
  dplyr::select(-Rho) %>% 
  spread(Gene1,Pvalue) %>%
  arrange(Type,Gene2) %>% 
  remove_rownames() %>%
  column_to_rownames("Gene2") %>% dplyr::select(-Type)

Pearson.data.Rho <- Pearson.data2 %>% arrange(Type) %>%
  dplyr::select(-Pvalue) %>% 
  spread(Gene1,Rho)%>%
  arrange(Type,Gene2) %>% 
  remove_rownames() %>%
  column_to_rownames("Gene2")
library(RColorBrewer)
annotation_row = data.frame(Class = factor(Pearson.data.Rho$Type))
Pearson.data.Rho2 <- Pearson.data.Rho %>% dplyr::select(-Type)
rownames(annotation_row) = rownames(Pearson.data.Rho2)

colors <- brewer.pal(8,"Set1")[1:7]
names(colors) <- unique(Pearson.data.Rho$Type)

anno_color <- list(Class = c(`Antigen-presenting molecules`=brewer.pal(8,"Set1")[1],
                             `Chemokines and receptors`=brewer.pal(8,"Set1")[2],
                             `Immunoinhibitors`=brewer.pal(8,"Set1")[3],
                             `Immunostimulators`=brewer.pal(8,"Set1")[4],
                             `Interferon and receptors`=brewer.pal(8,"Set1")[5],
                             `Interleukin and receptors`=brewer.pal(8,"Set1")[6],
                             `Other immunomodulators`=brewer.pal(8,"Set1")[7]))

pheatmap(Pearson.data.Rho2, 
         display_numbers = matrix(ifelse(Pearson.data.Pvalue > 0.05, "",
                                         if_else(Pearson.data.Pvalue > 0.01,"*",
                                                 if_else(Pearson.data.Pvalue > 0.001,"**","***"))), 
                                  nrow(Pearson.data.Rho2)), 
         annotation_row = annotation_row, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)[11:50],
         border_color = "black",
         show_rownames = F,
         annotation_colors = anno_color,
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         #cellwidth = 30,
         #cellheight = 3,
         #fontsize_row = 5,
         fontsize_col = 7)

#pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

##################### TCGA CBX2 CEP55 => Expression => Drugs #######
#### Gdsc ####
library(readxl)
data <- read_xlsx("Drugs/CBX2.CEP55.GdscIC50AndExprTable.xlsx")
data2 <- data %>% filter(fdr<=0.001) %>% arrange(cor)
data3 <- data %>% filter(drug %in% data2$drug) %>%
  arrange(cor) %>% mutate(drug=factor(drug,data2$drug))
p123 <- ggplot()+
  geom_point(data=data3,aes(symbol,drug,fill=cor,size=-log10(fdr)),shape=21,color="white")+
  geom_point(data=data3 %>% filter(fdr<=0.05),aes(symbol,drug,fill=cor,size=-log10(fdr)),
             shape=21,color="black",stroke=1.5)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(#legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9,angle = 90,hjust=0.5),
        axis.text.y = element_text(size=9),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_fill_gradient2("Spearman'CC",low = "#2166AC", mid = "white", high = "#B2182B",midpoint = 0)+ #,,limits = c(-2, 2))
  labs(y="",x="")
ggsave(p123,filename = "Drugs/CBX2.CEP55.Spearman.pdf",height = 10,width = 4)
#### Ctrp ####
setwd("K:/TCGA/Anlysis/LIHC")
data <- read_xlsx("Drugs/CtrpDrugIC50AndExpTable.xlsx")
data2 <- data %>% filter(fdr<=0.001) %>% arrange(cor) %>% top_n(50,cor)
data3 <- data %>% filter(drug %in% data2$drug) %>%
  arrange(cor) %>% 
  mutate(drug=factor(drug,unique(data2$drug)))

p123 <- ggplot()+
  geom_point(data=data3,aes(symbol,drug,fill=cor,size=-log10(fdr)),shape=21,color="white")+
  geom_point(data=data3 %>% filter(fdr<=0.05),aes(symbol,drug,fill=cor,size=-log10(fdr)),
             shape=21,color="black",stroke=1.5)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_fill_gradient2("Spearman'CC",low = "#2166AC", mid = "white", high = "#B2182B",midpoint = 0)+ #,,limits = c(-2, 2))
  labs(y="",x="")
ggsave(p123,filename = "Drugs/CBX2.CEP55.CTRP.Spearman.pdf",height = 10,width = 4.5)

#### Overlap ####
Ctrp <- read_xlsx("Drugs/CtrpDrugIC50AndExpTable.xlsx") #%>%
  filter(fdr<=0.01)
Gdsc <- read_xlsx("Drugs/CBX2.CEP55.GdscIC50AndExprTable.xlsx") #%>%
  filter(fdr<=0.01)

Inter.Drugs <- intersect(Ctrp$drug,Gdsc$drug)
Ctrp.Mid <- Ctrp %>% filter(drug %in% Inter.Drugs) %>% mutate(Database="CTRP")
Gdsc.Mid <- Gdsc %>% filter(drug %in% Inter.Drugs) %>% mutate(Database="GDSC")

Drugs <- rbind.data.frame(Ctrp.Mid,Gdsc.Mid)

p123 <- ggplot()+
  geom_point(data=Drugs,aes(symbol,drug,fill=cor,size=-log10(fdr)),
             shape=21,color="white")+
  geom_point(data=Drugs %>% filter(fdr<=0.05),aes(symbol,drug,fill=cor,size=-log10(fdr)),
             shape=21,color="black",stroke=1)+
  scale_size_continuous(range = c(1,4))+
  ggthemes::theme_few()+
  facet_wrap(~Database)+
  #coord_flip()+
  theme(#legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  scale_fill_gradient2("Spearman'CC",low = "#053061", mid = "white", high = "#67001F",midpoint = 0)+ #,,limits = c(-2, 2))
  labs(y="",x="")
ggsave(p123,filename = "Drugs/CBX2.CEP55.Intersect.Spearman.pdf",height = 4,width = 4)


##################### TCGA CBX2 CEP55 => Expression => TIGER #############
setwd("K:/TCGA/Anlysis/LIHC")
library(tidyverse)
data <- read.csv("TIGER.ICB.csv")
colnames(data)
p1 <- ggplot(data,aes(Log2FC,X.LOG10.Pvalue.,shape=Symbol,color=CancerType))+
  geom_point(size=3)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Log2FC",y="-log10(P)")+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2)+
  geom_vline(aes(xintercept=1),linetype=2)+
  geom_vline(aes(xintercept=-1),linetype=2)+
  scale_x_continuous(limits = c(-1.6,1.6))+
  ggsci::scale_color_nejm()

ggsave(p1,filename = "TIGER.ICB.pdf",width = 5.5,height = 3.2)

###
data2 <- read.csv("TIGER.ICB.PrePost.csv")
colnames(data2)
p2 <- ggplot(data2,aes(Log2FC,X.NAME.,shape=Symbol,color=CancerType))+
  geom_point(size=3)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Log2FC",y="-log10(P)")+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2)+
  geom_vline(aes(xintercept=1),linetype=2)+
  geom_vline(aes(xintercept=-1),linetype=2)+
  scale_x_continuous(limits = c(-1.2,1.2))+
  ggsci::scale_color_nejm()

ggsave(p2,filename = "TIGER.ICB.PrePost.pdf",width = 4,height = 2.5)

###
data3 <- read.csv("TIGER.ICB.Survival.csv")
colnames(data3)
p3 <- ggplot(data3,aes(HR,X.NAME.,shape=Symbol,color=CancerType))+
  geom_point(size=3)+
  ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title = element_text(size=15),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="HR",y="-log10(P)")+
  geom_hline(aes(yintercept=-log10(0.05)),linetype=2)+
  geom_vline(aes(xintercept=1),linetype=2)+
  #geom_vline(aes(xintercept=-1),linetype=2)+
  #scale_x_continuous(limits = c(1.2))+
  ggsci::scale_color_nejm()

ggsave(p3,filename = "TIGER.ICB.Survival.pdf",width = 5,height = 3)

############################### CBX2 AML #######################
#### ATAC ####
library(tidyverse)
library(magrittr)
setwd("K:/HepG2-CBX2")
Data <- read.csv("AML.CBX2/GSE193477.U937/U937.Overlapped.Readcount.txt",header = F,sep = "\t")
CountData <- Data[,c(1:4,11,18,25)] %>%
  data.frame(check.names = F) %>%
  set_colnames(c("Chr","Start","End","WT1","WT2","KO1","KO2"))
write.csv(CountData,file = "AML.CBX2/GSE193477.U937/U937.Overlapped.Readcount.Merge.txt",row.names = F)

dim(CountData) # 32212
library(DESeq2)
countData <- CountData[,c(4:7)] %>% data.frame(check.names = F) %>%
  set_rownames(paste("Peaks",1:nrow(CountData),sep = "_"))

condition <- factor(c(rep("WT",2),rep("KO",2)))
colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","KO","WT"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  cbind.data.frame(CountData)
write.csv(res1,file = "AML.CBX2/GSE193477.U937/U937.Overlapped.Readcount.DESeq2.csv",row.names = F)
#### 
Sig.Peak <- res1 %>% filter(padj <= 0.05)
KO.Up.Peak <- res1 %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1)
KO.Down.Peak <- res1 %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1)

write.csv(KO.Up.Peak,file = "AML.CBX2/GSE193477.U937/U937.Overlapped.Readcount.DESeq2.KO.Up.csv",row.names = F)
write.csv(KO.Down.Peak,file = "AML.CBX2/GSE193477.U937/U937.Overlapped.Readcount.DESeq2.KO.Down.csv",row.names = F)

write.table(KO.Up.Peak[,c("Chr","Start","End")],file = "AML.CBX2/GSE193477.U937/U937.Overlapped.KO.UP.bed",row.names = F,sep = "\t",quote = F)
write.table(KO.Down.Peak[,c("Chr","Start","End")],file = "AML.CBX2/GSE193477.U937/U937.Overlapped.KO.Down.bed",row.names = F,sep = "\t",quote = F)

#### ChIPseeker => Peaks Up ####
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)

KO.Up.Peak <- readPeakFile("AML.CBX2/GSE193477.U937/U937.Overlapped.KO.UP.bed")
peakAnno <- annotatePeak(KO.Up.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
peakAnno.df <- as.data.frame(peakAnno)
gene <- bitr(peakAnno.df$geneId,fromType = 'ENTREZID',
             toType = c('SYMBOL'),
             OrgDb="org.Hs.eg.db")
# bitr_kegg()
peakAnno.df2 <- merge(peakAnno.df,gene,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(peakAnno.df2,file = "AML.CBX2/GSE193477.U937/U937.Overlapped.KO.UP.Peaks.ChIPseeker.csv",row.names = F)

peakAnno.GO <- enrichGO(
  gene = peakAnno.df2$geneId,keyType="ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(peakAnno.GO,file = "AML.CBX2/GSE193477.U937/U937.Overlapped.KO.UP.Peaks.GOAnnotation.csv",row.names = F)

#### ChIPseeker => Peaks Down ####
KO.Down.Peak <- readPeakFile("U937.Overlapped.KO.Down.bed")
peakAnno.Down <- annotatePeak(KO.Down.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
peakAnno.df.Down <- as.data.frame(peakAnno.Down)
gene.Down <- bitr(peakAnno.df.Down$geneId,fromType = 'ENTREZID',
             toType = c('SYMBOL'),
             OrgDb="org.Hs.eg.db")
# bitr_kegg()
peakAnno.df2.Down <- merge(peakAnno.df.Down,gene.Down,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(peakAnno.df2.Down,file = "U937.Overlapped.KO.Down.Peaks.ChIPseeker.csv",row.names = F)

#### ChIPseeker => WT Specific peaks ####
library(clusterProfiler)
library(ChIPseeker)
setwd("K:/HepG2-CBX2/AML.CBX2/GSE193477.U937")
WT.Spec.Peak <- readPeakFile("U937.CBX2_AML.WT.ATAC.IDR.bed2")
KO.Spec.Peak <- readPeakFile("U937.CBX2_AML.CBX2KO.ATAC.IDR.bed2")

WT.Spec.Peak.Anno <- annotatePeak(WT.Spec.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
WT.Spec.Peak.Anno <- as.data.frame(WT.Spec.Peak.Anno)
gene.Anno <- bitr(WT.Spec.Peak.Anno$geneId,fromType = 'ENTREZID',
                  toType = c('SYMBOL'),
                  OrgDb="org.Hs.eg.db")
WT.Spec.Peak.Anno.M <- merge(WT.Spec.Peak.Anno,gene.Anno,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(WT.Spec.Peak.Anno.M,file = "U937.WT.Specific.Peaks.ChIPseeker.csv",row.names = F)

WT.Spec.Peak.Anno.M.GO <- enrichGO(
  gene = WT.Spec.Peak.Anno.M$geneId,keyType="ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(WT.Spec.Peak.Anno.M.GO,file = "U937.WT.Specific.Peaks.GOAnnotation.csv",row.names = F)

WT.Spec.Peak.Anno.M.KEGG <- enrichKEGG(gene = na.omit(WT.Spec.Peak.Anno$geneId),
                                       organism = "hsa",
                                       #keyType = "kegg",
                                       pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       minGSSize = 3,use_internal_data = FALSE)

write.csv(WT.Spec.Peak.Anno.M.KEGG@result,file = "U937.WT.Specific.Peaks.KEGGAnnotation.csv",row.names = F)


#### ChIPseeker => KO Specific peaks ####
KO.Spec.Peak.Anno <- annotatePeak(KO.Spec.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
KO.Spec.Peak.Anno <- as.data.frame(KO.Spec.Peak.Anno)
gene.Anno <- bitr(KO.Spec.Peak.Anno$geneId,fromType = 'ENTREZID',
                  toType = c('SYMBOL'),
                  OrgDb="org.Hs.eg.db")
KO.Spec.Peak.Anno.M <- merge(KO.Spec.Peak.Anno,gene.Anno,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(KO.Spec.Peak.Anno.M,file = "U937.KO.Specific.Peaks.ChIPseeker.csv",row.names = F)

KO.Spec.Peak.Anno.M.GO <- enrichGO(
  gene = KO.Spec.Peak.Anno.M$geneId,keyType="ENTREZID",
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(KO.Spec.Peak.Anno.M.GO,file = "U937.KO.Specific.Peaks.GOAnnotation.csv",row.names = F)


KO.Spec.Peak.Anno.M.KEGG <- enrichKEGG(gene = na.omit(KO.Spec.Peak.Anno$geneId),
                                       organism = "hsa",
                                       #keyType = "kegg",
                                       pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       minGSSize = 3,use_internal_data = FALSE)

write.csv(KO.Spec.Peak.Anno.M.KEGG@result,file = "U937.KO.Specific.Peaks.KEGGAnnotation.csv",row.names = F)





###### RNA -> GSEA ####
setwd("K:/HepG2-CBX2")
RNA.Data <- read.csv("AML.CBX2/GSE193477.U937/shCBX2VSControl.diff",sep = "\t")
colnames(RNA.Data)
## GSEA
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)
dir.create("AML.CBX2/GSE193477.U937/Hallmark.GSEA")
ge = RNA.Data$GFOLD.0.01.
names(ge) = RNA.Data$GeneSymbol
ge = sort(ge,decreasing = T)
Hallmark <- read.gmt("../../TCGA/msigdb/h.all.v7.5.1.symbols.gmt")

em <- GSEA(ge, TERM2GENE = Hallmark,pvalueCutoff = 1,pAdjustMethod = "BH")
write.csv(em@result,file = "AML.CBX2/GSE193477.U937/Hallmark.GSEA/RNA.Hallmark.GSEA.csv",row.names = F)
save(em,file = "AML.CBX2/GSE193477.U937/Hallmark.GSEA/Hallmark.GSEA.Rdata")

KEGG <- read.gmt("../../TCGA/msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
em.KEGG <- GSEA(ge, TERM2GENE = KEGG,pvalueCutoff = 0.05,pAdjustMethod = "BH")
save(em.KEGG,file = "AML.CBX2/GSE193477.U937/Hallmark.GSEA/KEGG.GSEA.Rdata")
write.csv(em.KEGG@result,file = "RNA.KEGG.GSEA.csv",row.names = F)
library(GseaVis)
setwd("K:/HepG2-CBX2/AML.CBX2/GSE193477.U937/Hallmark.GSEA")
gseaNb(object = em.KEGG,
       geneSetID = "KEGG_SPLICEOSOME",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

gseaNb(object = em.KEGG,
       geneSetID = "KEGG_CELL_CYCLE",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

gseaNb(object = em.KEGG,
       geneSetID = "KEGG_RIBOSOME",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)


gseaNb(object = em.KEGG,
       geneSetID = "KEGG_DNA_REPLICATION",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

### GO
GO <- read.gmt("../../../../../TCGA/msigdb/c5.go.v7.5.1.symbols.gmt")
em.GO <- GSEA(ge, TERM2GENE = GO,pvalueCutoff = 0.05,pAdjustMethod = "BH")
save(em.GO,file = "GO.GSEA.Rdata")
write.csv(em.GO@result,file = "RNA.GO.GSEA.csv",row.names = F)
em.GO.res <- em.GO@result %>% arrange(NES)

terms1 <- c('GOBP_CHROMATIN_REMODELING',
           'GOBP_CHROMOSOME_SEGREGATION',
           'GOBP_DNA_REPLICATION',
           'GOBP_MITOTIC_NUCLEAR_DIVISION',
           'GOBP_MRNA_PROCESSING',
           'GOBP_NUCLEAR_CHROMOSOME_SEGREGATION',
           'GOBP_RNA_SPLICING',
           'GOBP_RIBOSOME_BIOGENESIS',
           'GOBP_MEIOTIC_CELL_CYCLE',
           'GOBP_TELOMERE_ORGANIZATION',
           'GOCC_SPINDLE',
           'GOCC_PROTEIN_DNA_COMPLEX')

terms2 <- c("GOMF_PHOSPHOLIPASE_ACTIVITY",
            "GOCC_MHC_PROTEIN_COMPLEX")

lapply(terms1, function(x){
  gseaNb(object = em.GO,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,subPlot = 2)
}) -> gseaList3
cowplot::plot_grid(plotlist = gseaList3,ncol = 4,align = 'hv')

lapply(terms2, function(x){
  gseaNb(object = em.GO,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,subPlot = 2)
}) -> gseaList4
cowplot::plot_grid(plotlist = gseaList4,ncol = 2,align = 'hv')

# Hallmark
library(GseaVis)
terms1 <- c('HALLMARK_MYC_TARGETS_V1',
            'HALLMARK_E2F_TARGETS',
            'HALLMARK_G2M_CHECKPOINT',
            'HALLMARK_MYC_TARGETS_V2',
            'HALLMARK_MITOTIC_SPINDLE',
            'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
            'HALLMARK_FATTY_ACID_METABOLISM',
            'HALLMARK_DNA_REPAIR')

terms2 <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_KRAS_SIGNALING_DN",
            "HALLMARK_P53_PATHWAY")

lapply(terms1, function(x){
  gseaNb(object = em,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,subPlot = 2)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 4,align = 'hv')

lapply(terms2, function(x){
  gseaNb(object = em,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,subPlot = 2)
}) -> gseaList2

# combine
cowplot::plot_grid(plotlist = gseaList2,ncol = 3,align = 'hv')

gseaNb(object = em,
       geneSetID = 'HALLMARK_MYC_TARGETS_V1',
       subPlot = 2,
       termWidth = 30,addPval = T,newGsea = T,
       pvalX = 0.75,pvalY = 0.8,
       pCol = 'black',pHjust = 0)

if (F) {
  library(fgsea)
  Hallmark$term <- unfactor(Hallmark$term)
  Hallmark2 <- Hallmark %>% split(x = .$gene, f = .$term)
  em <- fgsea(Hallmark2, stats = ge, nperm = 10000,minSize=3, maxSize=100000)
  #write.csv(unlist(em),file = "AML.CBX2/GSE193477.U937/Hallmark.GSEA/RNA.Hallmark.GSEA.csv",row.names = F)
  data.table::fwrite(em, file="AML.CBX2/GSE193477.U937/Hallmark.GSEA/RNA.Hallmark.GSEA.csv", sep=",", sep2=c("", " ", ""))
  
  ## NES
  em2 <- as.data.frame(em) %>% arrange(desc(NES)) %>% 
    mutate(pathway=str_remove(pathway,"HALLMARK_")) %>% 
    filter(NES >= 1 | NES <= -1) %>%
    mutate(mycolor=if_else(NES >= 1,"NES>=1","NES<=-1"))
  
  p <- ggplot(em2, aes(x=reorder(pathway,NES), y=NES,fill = mycolor)) +
    geom_segment( aes(x=reorder(pathway,NES), xend=reorder(pathway,-NES), 
                      y=0, yend=NES,color=mycolor),
                  size=0.5,linetype=1)+# π”√reorder()≈≈–Ú±‰¡ø
    geom_point(size = 5, pch = 21, color="black") + 
    #‘⁄…¢µ„…œœ‘ æ ˝÷µ≤¢±£¡Ù¡ΩŒª–° ˝
    geom_text(aes(label =sprintf("%.2f",NES)), color = "black", size = 3,hjust=-0.5)+ 
    #scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set2")[2:1])+
    ggsci::scale_fill_lancet()+
    ggsci::scale_color_lancet()+
    coord_flip() +
    ggthemes::theme_few()+
    theme(legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=15),
          panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
    scale_y_continuous("NES",breaks  =c(-3,-2,-1,0,1,2,3,4),limits = c(-3.5,4.3))+
    labs(x="")
  
  ggsave(p,filename = "AML.CBX2/GSE193477.U937/Hallmark.GSEA/RNA.Hallmark.GSEA.NES.pdf",height = 6,width = 6)
  
  plotEnrichment(em[["HALLMARK_MYC_TARGETS_V1"]],
                 ge)
  
}





###### RNA -> GO + KEGG ####
setwd("K:/HepG2-CBX2")
library(tidyverse)
library(magrittr)
RNA.Data <- read.csv("AML.CBX2/GSE193477.U937/shCBX2VSControl.diff",sep = "\t",skip = "#")
RNA.Data.Meta <- read.csv("AML.CBX2/GSE193477.U937/shCBX2VSControl.diff.ext",sep = "\t",skip = "#")
dim(RNA.Data)
dim(RNA.Data.Meta)
RNA.Data.All <- cbind(RNA.Data,RNA.Data.Meta)
Sig.RNA.Down <- RNA.Data.All %>% filter(E.FDR == 1) %>% filter(GFOLD.0.01. >= 1)
Sig.RNA.Up <- RNA.Data.All %>% filter(E.FDR == 1) %>% filter(GFOLD.0.01. <= -1)

# Sig.RNA.Down => GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

Sig.RNA.Down.GO <- enrichGO(
  gene = Sig.RNA.Down$GeneSymbol,
  keyType="SYMBOL",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write.csv(Sig.RNA.Down.GO@result,file = "AML.CBX2/GSE193477.U937/shCBX2VSControl.RNA.Down.GO.csv",row.names = F)

gene.Anno.K <- bitr(Sig.RNA.Down$GeneSymbol,fromType = 'SYMBOL',
                  toType = c('ENTREZID'),
                  OrgDb="org.Hs.eg.db")

gene.Anno.K2 <- merge(gene.Anno.K,Sig.RNA.Down,by.x="SYMBOL",by.y="GeneSymbol")
write.csv(gene.Anno.K2,file = "AML.CBX2/GSE193477.U937/shCBX2VSControl.RNA.Down.Gene.csv",row.names = F)
library(clusterProfiler)
gene.Anno.K.A <- enrichKEGG(gene = na.omit(gene.Anno.K$ENTREZID),
                            organism = "hsa",
                            #keyType = "kegg",
                            pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize = 3,use_internal_data = FALSE)
save(gene.Anno.K.A,file="shCBX2VSControl.RNA.Down.KEGG.Rdata")

###
Sig.RNA.Up.GO <- enrichGO(
  gene = Sig.RNA.Up$GeneSymbol,
  keyType="SYMBOL",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write.csv(Sig.RNA.Up.GO@result,file = "AML.CBX2/GSE193477.U937/shCBX2VSControl.RNA.Up.GO.csv",row.names = F)

gene.Anno.U <- bitr(Sig.RNA.Up$GeneSymbol,fromType = 'SYMBOL',
                    toType = c('ENTREZID'),
                    OrgDb="org.Hs.eg.db")

gene.Anno.U.A <- enrichKEGG(gene = na.omit(gene.Anno.U$ENTREZID),
                            organism = "hsa",
                            #keyType = "kegg",
                            pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize = 3,use_internal_data = FALSE)
save(gene.Anno.U.A,file="shCBX2VSControl.RNA.Up.KEGG.Rdata")


###### ATAC + RNA ######
#### ATAC KO Up => RNA KO Up ####


#### ATAC KO Down => RNA KO Down ####
RNA.Data <- read.csv("AML.CBX2/GSE193477.U937/shCBX2VSControl.diff",sep = "\t")
RNA.Data.Meta <- read.csv("AML.CBX2/GSE193477.U937/shCBX2VSControl.diff.ext",sep = "\t")
# dim(RNA.Data)
# dim(RNA.Data.Meta)
RNA.Data.All <- cbind(RNA.Data,RNA.Data.Meta)
Sig.RNA.Down <- RNA.Data.All %>% filter(E.FDR == 1) %>% filter(GFOLD.0.01. >= 1) %>%
  arrange(desc(GFOLD.0.01.))
Sig.RNA.Up <- RNA.Data.All %>% filter(E.FDR == 1) %>% filter(GFOLD.0.01. <= -1)%>%
  arrange(GFOLD.0.01.)

ATAC.KO.Down <- read.csv("AML.CBX2/GSE193477.U937/U937.Overlapped.KO.Down.Peaks.ChIPseeker.csv") #%>%
  filter(str_detect(annotation,"Promoter"))
ATAC.KO.Up <- read.csv("AML.CBX2/GSE193477.U937/U937.Overlapped.KO.Up.Peaks.ChIPseeker.csv") %>%
  filter(str_detect(annotation,"Promoter"))
#ATAC.WT.Spec <- read.csv("AML.CBX2/GSE193477.U937/U937.WT.Specific.Peaks.ChIPseeker.csv") %>%
#  filter(str_detect(annotation,"Promoter"))
#ATAC.KO.Spec <- read.csv("AML.CBX2/GSE193477.U937/U937.KO.Specific.Peaks.ChIPseeker.csv") %>%
#  filter(str_detect(annotation,"Promoter"))



KO.Down <- intersect(ATAC.KO.Down$SYMBOL,
                   Sig.RNA.Down$GeneSymbol)
KO.Up <- intersect(ATAC.KO.Up$SYMBOL,
                     Sig.RNA.Up$GeneSymbol)

Hallmark <- read.gmt("../../TCGA/msigdb/h.all.v7.5.1.symbols.gmt")
Inter.UP <- Hallmark %>% filter(gene %in% KO.Up)

Inter.Down <- Hallmark %>% filter(term %in% terms1)
Down.Sig <- merge(Sig.RNA.Down,Inter.Down,by.x="GeneSymbol",by.y="gene") %>%
  data.frame() %>%
  arrange(desc(GFOLD.0.01.))


KO.Up.GO <- enrichGO(
  gene = KO.Up,
  keyType="SYMBOL",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)



############################## Mouse -> MEF CBX2 ##################################
setwd("K:/HepG2-CBX2/GSE156413.MEFs.CBX2")
###### RNA ####
Countdata <- read.table("MEF.Cbx2.count.txt",header = T,sep = "\t")
Countdata2 <- Countdata[,c(1,7:10)] %>% data.frame() %>%
  remove_rownames() %>% 
  column_to_rownames("Geneid") %>%
  set_colnames(str_remove(colnames(.),".sorted.bam"))

library(DESeq2)
countData <- Countdata2 %>% data.frame(check.names = F) %>%
  set_colnames(str_remove(colnames(.),"RNA_"))

condition <- factor(c(rep("KO",2),rep("WT",2)))
colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","KO","WT"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  na.omit() %>% rownames_to_column("geneid") %>%
  merge(countData%>% rownames_to_column("geneid"),by="geneid")

write.csv(res1,file = "MEF.Cbx2.DESeq2.csv",row.names = F)

x <- countData/Countdata$Length
expMatrix_tpm <- t( t(x) / colSums(x) ) * 1e6
write.csv(expMatrix_tpm,file = "MEF.Cbx2.TPM.csv")

#### GSEA
ge = -res1$log2FoldChange
names(ge) = res1$geneid %>% toupper()
ge = sort(ge,decreasing = T)
Hallmark <- read.gmt("../../../../TCGA/msigdb/h.all.v7.5.1.symbols.gmt")

em2 <- GSEA(ge, TERM2GENE = Hallmark,pvalueCutoff = 1,pAdjustMethod = "BH")
write.csv(em2@result,file = "MEF.Hallmark.GSEA.csv",row.names = F)
save(em2,file = "MEF.Hallmark.GSEA.Rdata")

library(GseaVis)
gseaNb(object = em2,
       geneSetID = "HALLMARK_E2F_TARGETS",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

gseaNb(object = em2,
       geneSetID = "HALLMARK_G2M_CHECKPOINT",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)


MEF.KOUp <- res1 %>% filter(log2FoldChange >= 1) %>%
  filter(padj <= 0.05)

MEF.KODown <- res1 %>% filter(log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>% arrange(log2FoldChange)


library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
setwd("K:/HepG2-CBX2/GSE156413.MEFs.CBX2")

MEF.KOUp.GO <- enrichGO(
  gene = MEF.KOUp$geneid,
  keyType="SYMBOL",
  OrgDb = "org.Mm.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write.csv(MEF.KOUp.GO@result,file = "RNA.Up.GO.csv",row.names = F)

MEF.KODown.GO <- enrichGO(
  gene = MEF.KODown$geneid,
  keyType="SYMBOL",
  OrgDb = "org.Mm.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write.csv(MEF.KODown.GO@result,file = "RNA.Down.GO.csv",row.names = F)

###### ATAC ####
setwd("K:/HepG2-CBX2/GSE156413.MEFs.CBX2")
library(clusterProfiler)
library(ChIPseeker)
library(tidyverse)
library(magrittr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(org.Mm.eg.db)
#### WT-Specific ####
WT.Spec.Peak <- readPeakFile("MEF.Cbx2_MEF.WT.ATAC.IDR.bed2")
WT.Spec.Peak.Anno <- annotatePeak(WT.Spec.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
WT.Spec.Peak.Anno <- as.data.frame(WT.Spec.Peak.Anno)
gene.Anno <- bitr(WT.Spec.Peak.Anno$geneId,fromType = 'ENTREZID',
                  toType = c('SYMBOL'),
                  OrgDb="org.Mm.eg.db")
WT.Spec.Peak.Anno.M <- merge(WT.Spec.Peak.Anno,gene.Anno,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(WT.Spec.Peak.Anno.M,file = "MEF.Cbx2_WT.Specific.Peaks.ChIPseeker.csv",row.names = F)

WT.Spec.Peak.Anno.M.GO <- enrichGO(
  gene = WT.Spec.Peak.Anno.M$geneId,keyType="ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(WT.Spec.Peak.Anno.M.GO,file = "MEF.Cbx2_WT.Specific.Peaks.GOAnnotation.csv",row.names = F)

WT.Spec.Peak.Anno.M.KEGG <- enrichKEGG(gene = na.omit(WT.Spec.Peak.Anno$geneId),
                                       organism = "mouse",
                                       #keyType = "kegg",
                                       pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       minGSSize = 3,use_internal_data = FALSE)

write.csv(WT.Spec.Peak.Anno.M.KEGG@result,file = "MEF.Cbx2_WT.Specific.Peaks.KEGGAnnotation.csv",row.names = F)
#### KO-Specific ####
KO.Spec.Peak <- readPeakFile("MEF.Cbx2_MEF.Cbx2KO.ATAC.IDR.bed2")
KO.Spec.Peak.Anno <- annotatePeak(KO.Spec.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
KO.Spec.Peak.Anno <- as.data.frame(KO.Spec.Peak.Anno)
gene.Anno <- bitr(KO.Spec.Peak.Anno$geneId,fromType = 'ENTREZID',
                  toType = c('SYMBOL'),
                  OrgDb="org.Mm.eg.db")
KO.Spec.Peak.Anno.M <- merge(KO.Spec.Peak.Anno,gene.Anno,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(KO.Spec.Peak.Anno.M,file = "MEF.Cbx2_KO.Specific.Peaks.ChIPseeker.csv",row.names = F)

KO.Spec.Peak.Anno.M.GO <- enrichGO(
  gene = KO.Spec.Peak.Anno.M$geneId,keyType="ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(KO.Spec.Peak.Anno.M.GO,file = "MEF.Cbx2_KO.Specific.Peaks.GOAnnotation.csv",row.names = F)

KO.Spec.Peak.Anno.M.KEGG <- enrichKEGG(gene = na.omit(KO.Spec.Peak.Anno$geneId),
                                       organism = "mouse",
                                       #keyType = "kegg",
                                       pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       minGSSize = 3,use_internal_data = FALSE)

write.csv(KO.Spec.Peak.Anno.M.KEGG@result,file = "MEF.Cbx2_KO.Specific.Peaks.KEGGAnnotation.csv",row.names = F)

#### DEseq2 ####
setwd("K:/HepG2-CBX2/GSE156413.MEFs.CBX2")
data <- read.table("Overlapped.Readcount.txt",sep = "\t",header = F)
countData <- data[,c(1:3,4,11,18,25)] %>% data.frame() %>%
  mutate(Peak=paste(V1,V2,V3,sep = "_")) %>%
  dplyr::select(-V1,-V2,-V3) %>%
  remove_rownames() %>% 
  column_to_rownames("Peak") %>%
  set_colnames(c("KO1","KO2","WT1","WT2"))

condition <- factor(c(rep("KO",2),rep("WT",2)))
colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","KO","WT"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  cbind.data.frame(countData)
write.csv(res1,file = "ATAC.DESeq2.csv")

#### DAP -> WT #### 
data <- read.table("Overlapped.Readcount.txt",sep = "\t",header = F)
Mid <- data[,c(1:3,4,11,18,25)] %>% data.frame() %>%
  mutate(Peak=paste(V1,V2,V3,sep = "_"))
WT.DAP <- res1 %>% filter(padj <= 0.05) %>%
  filter(log2FoldChange <= -1) %>%
  rownames_to_column("Peak") %>%
  merge(Mid,by="Peak")

write.csv(WT.DAP,file = "ATAC.DESeq2.Down.csv",row.names = F)
write.table(WT.DAP[,c("V1","V2","V3")],file = "ATAC.DESeq2.Down.bed",sep = "\t",row.names = F,quote = F,col.names = F)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(org.Mm.eg.db)
library(ChIPseeker)
DAP.Down.Peak <- readPeakFile("ATAC.DESeq2.Down.bed")
DAP.Down.Peak.Anno <- annotatePeak(DAP.Down.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
DAP.Down.Peak.Anno <- as.data.frame(DAP.Down.Peak.Anno)
gene.Anno <- bitr(DAP.Down.Peak.Anno$geneId,fromType = 'ENTREZID',
                  toType = c('SYMBOL'),
                  OrgDb="org.Mm.eg.db")
DAP.Down.Peak.Anno.M <- merge(DAP.Down.Peak.Anno,gene.Anno,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(DAP.Down.Peak.Anno.M,file = "MEF.Cbx2_DAP.Down.Specific.Peaks.ChIPseeker.csv",row.names = F)

DAP.Down.Peak.Anno.M.GO <- enrichGO(
  gene = DAP.Down.Peak.Anno.M$geneId,keyType="ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(DAP.Down.Peak.Anno.M.GO,file = "MEF.Cbx2_DAP.Down.Specific.Peaks.GOAnnotation.csv",row.names = F)

DAP.Down.Peak.Anno.M.KEGG <- enrichKEGG(gene = na.omit(DAP.Down.Peak.Anno$geneId),
                                       organism = "mouse",
                                       #keyType = "kegg",
                                       pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       minGSSize = 3,use_internal_data = FALSE)

write.csv(DAP.Down.Peak.Anno.M.KEGG@result,file = "MEF.Cbx2_DAP.Down.Specific.Peaks.KEGGAnnotation.csv",row.names = F)

#### DAP -> KO ####
data <- read.table("Overlapped.Readcount.txt",sep = "\t",header = F)
Mid <- data[,c(1:3,4,11,18,25)] %>% data.frame() %>%
  mutate(Peak=paste(V1,V2,V3,sep = "_"))
KO.DAP <- res1 %>% filter(padj <= 0.05) %>%
  filter(log2FoldChange >= 1) %>%
  rownames_to_column("Peak") %>%
  merge(Mid,by="Peak")

write.csv(KO.DAP,file = "ATAC.DESeq2.Up.csv",row.names = F)
write.table(KO.DAP[,c("V1","V2","V3")],file = "ATAC.DESeq2.Up.bed",sep = "\t",row.names = F,quote = F,col.names = F)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(org.Mm.eg.db)

DAP.Up.Peak <- readPeakFile("ATAC.DESeq2.Up.bed")
DAP.Up.Peak.Anno <- annotatePeak(DAP.Up.Peak, tssRegion=c(-3000, 3000), TxDb=txdb)
DAP.Up.Peak.Anno <- as.data.frame(DAP.Up.Peak.Anno)
gene.Anno <- bitr(DAP.Up.Peak.Anno$geneId,fromType = 'ENTREZID',
                  toType = c('SYMBOL'),
                  OrgDb="org.Mm.eg.db")
DAP.Up.Peak.Anno.M <- merge(DAP.Up.Peak.Anno,gene.Anno,by.y="ENTREZID",by.x="geneId",all.x=T)
write.csv(DAP.Up.Peak.Anno.M,file = "MEF.Cbx2_DAP.Up.Specific.Peaks.ChIPseeker.csv",row.names = F)

DAP.Up.Peak.Anno.M.GO <- enrichGO(
  gene = DAP.Up.Peak.Anno.M$geneId,keyType="ENTREZID",
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)
write.csv(DAP.Up.Peak.Anno.M.GO,file = "MEF.Cbx2_DAP.Up.Specific.Peaks.GOAnnotation.csv",row.names = F)

DAP.Up.Peak.Anno.M.KEGG <- enrichKEGG(gene = na.omit(DAP.Up.Peak.Anno$geneId),
                                        organism = "mouse",
                                        #keyType = "kegg",
                                        pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                        pAdjustMethod = "BH",
                                        minGSSize = 3,use_internal_data = FALSE)

write.csv(DAP.Up.Peak.Anno.M.KEGG@result,file = "MEF.Cbx2_DAP.Up.Specific.Peaks.KEGGAnnotation.csv",row.names = F)



#### Test Overlap ####
Inter.Up <- intersect(toupper(MEF.KOUp$geneid),Sig.RNA.Up$GeneSymbol)
##  [1] "ADGRE1"   "ARHGAP8"  "DSP"      "FCER1G"   "LRRC25"   "NEFM"     "PPL"      "RASGEF1B"
##  [9] "RASL12"   "SLC15A3"  "SYT7"     "SYTL1"    "THY1"     "TMEM8B" 
Inter.Down <- intersect(toupper(MEF.KODown$geneid),Sig.RNA.Down$GeneSymbol)
## WT1
save(Sig.RNA.Up,Sig.RNA.Down,MEF.KODown,MEF.KOUp,file = "RNAseq.DEseq2.Rdata")


##################### GSE112227.M.Bone.RNA ######################
setwd("K:/HepG2-CBX2/GSE112227.M.Bone.RNA")
library(tidyverse)
library(magrittr)
library(DESeq2)
library(GseaVis)
library(ChIPseeker)
library(clusterProfiler)
data <- read.table("MouseBone.CBX2.count.txt",header = T,sep = "\t")
countData <- data[,c(1,7:12)] %>% data.frame() %>%
  set_colnames(str_remove(colnames(.),".sorted.bam")) %>%
  remove_rownames() %>%
  column_to_rownames("Geneid")

condition <- factor(c(rep("KO",3),rep("WT",3)))
colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","KO","WT"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  na.omit() %>% rownames_to_column("geneid") %>%
  merge(countData%>% rownames_to_column("geneid"),by="geneid")
write.csv(res1,file = "DESeq2.csv",row.names = F)

#### GSEA => Hallmark+KEGG+GO ####
ge = -res1$log2FoldChange
names(ge) = toupper(res1$geneid)
ge = sort(ge,decreasing = T)
library(ChIPseeker)
library(clusterProfiler)
Hallmark <- read.gmt("../../../../TCGA/msigdb/h.all.v7.5.1.symbols.gmt")

em <- GSEA(ge, TERM2GENE = Hallmark,pvalueCutoff = 0.05,pAdjustMethod = "BH")
write.csv(em@result,file = "RNA.Hallmark.GSEA.csv",row.names = F)
save(em,file = "Hallmark.GSEA.Rdata")

library(GseaVis)
gseaNb(object = em,
       geneSetID = "HALLMARK_E2F_TARGETS",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

gseaNb(object = em,
       geneSetID = "HALLMARK_G2M_CHECKPOINT",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)


KEGG <- read.gmt("../../../../TCGA/msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
em.KEGG <- GSEA(ge, TERM2GENE = KEGG,pvalueCutoff = 0.05,pAdjustMethod = "BH")
save(em.KEGG,file = "KEGG.GSEA.Rdata")
write.csv(em.KEGG@result,file = "RNA.KEGG.GSEA.csv",row.names = F)

gseaNb(object = em.KEGG,
       geneSetID = "KEGG_DNA_REPLICATION",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

gseaNb(object = em.KEGG,
       geneSetID = "KEGG_CELL_CYCLE",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

gseaNb(object = em.KEGG,
       geneSetID = "KEGG_SPLICEOSOME",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)


GO <- read.gmt("../../../../../TCGA/msigdb/c5.go.v7.5.1.symbols.gmt")
em.GO <- GSEA(ge, TERM2GENE = GO,pvalueCutoff = 0.05,pAdjustMethod = "BH")
save(em.GO,file = "GO.GSEA.Rdata")
write.csv(em.GO@result,file = "RNA.GO.GSEA.csv",row.names = F)

load(file = "GO.GSEA.Rdata")
terms1 <- c("GOBP_CHROMOSOME_SEGREGATION",
           "GOBP_NUCLEAR_CHROMOSOME_SEGREGATION",
           "GOBP_DNA_RECOMBINATION",
           "GOBP_DNA_REPLICATION",
           "GOBP_DNA_DEPENDENT_DNA_REPLICATION",
           "GOBP_DNA_REPAIR",
           "GOCC_PROTEIN_DNA_COMPLEX",
           "GOCC_SPINDLE")

terms2 <-c("GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS",
           "GOBP_ENERGY_RESERVE_METABOLIC_PROCESS")

lapply(terms1, function(x){
  gseaNb(object = em.GO,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,subPlot = 2)
}) -> gseaList5
cowplot::plot_grid(plotlist = gseaList5,ncol = 2,align = 'hv')

lapply(terms2, function(x){
  gseaNb(object = em.GO,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,subPlot = 2)
}) -> gseaList6
cowplot::plot_grid(plotlist = gseaList6,ncol = 2,align = 'hv')

#### KEGG ####
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Mm.eg.db))
setwd("K:/HepG2-CBX2/GSE112227.M.Bone.RNA")
res1 <- read.csv("DESeq2.csv")

Up.Gene <- res1 %>% filter(padj <= 0.05) %>%
  filter(log2FoldChange >= 1)
Down.Gene <- res1 %>% filter(padj <= 0.05) %>%
  filter(log2FoldChange <= -1)


UpGene = bitr(Up.Gene$geneid, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

UpGene.KEGG <- enrichKEGG(gene= UpGene$ENTREZID,
                          organism     = 'mouse',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 3,use_internal_data = FALSE)
write.csv(UpGene.KEGG@result,file = "DESeq2.UpGene.KEGG.csv",row.names = F)

DownGene = bitr(Down.Gene$geneid, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

DownGene.KEGG <- enrichKEGG(gene= DownGene$ENTREZID,
                          organism     = 'mouse',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          minGSSize = 3,use_internal_data = FALSE)
write.csv(DownGene.KEGG@result,file = "DESeq2.DownGene.KEGG.csv",row.names = F)

#### Test ####
load(file = "../GSE156413.MEFs.CBX2/RNAseq.DEseq2.Rdata")

Up.inter <- intersect(Up.Gene$geneid,MEF.KOUp$geneid)
# [1] "Adam8"     "Als2cl"    "Cdh5"      "Csf2rb2"   "Esm1"      "Gm6093"    "Ifitm1"   
# [8] "Kcnj8"     "Lamc3"     "Lilr4b"    "Lilrb4a"   "Map3k9"    "Ms4a6d"    "Msr1"     
# [15] "Perp"      "Prss16"    "Serpinb1a" "Shtn1"     "Syt7"      "Vsig2"
Down.inter <- intersect(Down.Gene$geneid,MEF.KODown$geneid)
#[1] "Acan"    "Col24a1" "Dcn"     "Fmn2"    "Fxyd6"   "Gfra1"   "Gm42362" "Lef1"   
#[9] "Lum"     "Map2"    "Mndal"   "Pdgfrl"  "Postn"   "Shox2"   "Slc14a1" "Tnn"
intersect(toupper(Down.inter),Sig.RNA.Down$GeneSymbol)
intersect(toupper(Up.inter),Sig.RNA.Up$GeneSymbol)




##################### GSE179160.M.Cbx2KO.K4me3 ##########################
setwd("K:/HepG2-CBX2/GSE179160.M.Cbx2KO.K4me3")
library(clusterProfiler)
library(ChIPseeker)
library(tidyverse)
library(magrittr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#### DAP ####
library(tidyverse)
data <- read.csv("HSPC.Count.txt",sep = "\t",header = T,skip=1) %>%
  mutate(Peak=paste(Chr,Start,End,sep = "_"))
data2 <- data %>% mutate(Peak=paste(Chr,Start,End,sep = "_")) %>%
  dplyr::select(-Chr,-Start,-End,-Geneid,-Strand,-Length) %>%
  remove_rownames() %>%
  column_to_rownames("Peak") %>%
  magrittr::set_colnames(str_remove_all(colnames(.),".sorted.bam"))

data3 <- data2 %>% mutate(KO1=HSPC.Cbx2KO.H3K4me3.Rep1,
                          KO2=HSPC.Cbx2KO.H3K4me3.Rep2,
                          WT1=HSPC.WT.H3K4me3.Rep1,
                          WT2=HSPC.WT.H3K4me3.Rep2) %>%
  dplyr::select(KO1,KO2,WT1,WT2)
  
data3[data3<0]=0
library(DESeq2)
condition <- factor(c(rep("KO",2),rep("WT",2)))
colData <- data.frame(row.names=colnames(data3), condition)

dds <- DESeqDataSetFromMatrix(countData = data3, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","KO","WT"))

res2 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  na.omit() %>% arrange(padj) %>%
  rownames_to_column("Peak") %>%
  merge(data,by="Peak")%>% arrange(padj)

res2 %>% filter(padj <= 0.05)


##################### GSE86084.M.ESC.CBX2 => RNA DOX ####################
setwd("K:/HepG2-CBX2/GSE86084.M.ESC.CBX2")
data <- read.csv("GSE86084.MouseESC.CBX2.count.txt",sep="\t",skip = 1,header = T)
data2 <- data %>% remove_rownames() %>% column_to_rownames("Geneid") %>%
  dplyr::select(ends_with("bam")) %>%
  magrittr::set_colnames(str_remove_all(colnames(.),".sorted.bam"))

countData <- data2[,1:4]
condition <- factor(c(rep(c("WT","OE"),times=2)))
colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","OE","WT"))

res2 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  na.omit() %>% arrange(padj) %>%
  rownames_to_column("geneid") %>%
  #rownames_to_column("Peak") %>%
  merge(countData%>% rownames_to_column("geneid"),by="geneid")%>% 
  arrange(padj)
write.csv(res2,file = "GSE86084.M.ESC.CBX2.WT1.DESeq2.csv",row.names = F)

ge = res2$log2FoldChange
names(ge) = toupper(res2$geneid)
ge = sort(ge,decreasing = T)
library(ChIPseeker)
library(clusterProfiler)
Hallmark <- read.gmt("../../../../TCGA/msigdb/h.all.v7.5.1.symbols.gmt")
KEGG <- read.gmt("../../../../TCGA/msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO <- read.gmt("../../../../TCGA/msigdb/c5.go.v7.5.1.symbols.gmt")

em <- GSEA(ge, TERM2GENE = Hallmark,pvalueCutoff = 0.05,pAdjustMethod = "BH",eps = 1e-100)
write.csv(em@result,file = "GSE86084.M.ESC.CBX2.WT1.Hallmark.GSEA.csv",row.names = F)

library(GseaVis)
gseaNb(object = em,
       geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)


em <- GSEA(ge, TERM2GENE = KEGG,pvalueCutoff = 0.05,pAdjustMethod = "BH",eps = 1e-100)
write.csv(em@result,file = "GSE86084.M.ESC.CBX2.WT1.KEGG.GSEA.csv",row.names = F)

em <- GSEA(ge, TERM2GENE = GO,pvalueCutoff = 0.05,pAdjustMethod = "BH",eps = 1e-100)
write.csv(em@result,file = "GSE86084.M.ESC.CBX2.WT1.GO.GSEA.csv",row.names = F)

gseaNb(object = em,
       geneSetID = "GOCC_ACTIN_CYTOSKELETON",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)

library(org.Mm.eg.db)
geneList = res2$log2FoldChange
names(geneList) = res2$geneid
geneList = sort(geneList,decreasing = T)
s2g=select(org.Mm.eg.db,names(geneList),
           'ENTREZID','SYMBOL')
s2g=s2g[!is.na(s2g$ENTREZID),]
s2g=s2g[s2g$SYMBOL %in% names(geneList),]
geneList=geneList[s2g$SYMBOL] 
names(geneList)=s2g$ENTREZID

em <- gseGO(geneList = geneList, 
      OrgDb='org.Mm.eg.db',
      ont = 'ALL',
      #nPerm = 1000,
      #minGSSize = 10,
      eps=1e-50,
      pvalueCutoff = 0.5,
      by="fgsea",
      verbose = FALSE)
em <- gseKEGG(geneList = geneList, 
            #OrgDb='org.Mm.eg.db',
            organism = "mmu",
            #ont = 'ALL',
            #nPerm = 1000,
            #minGSSize = 10,
            pvalueCutoff = 0.05,
            verbose = FALSE)

countData <- data2[,5:10]
condition <- factor(c(rep(c("WT","OE"),times=3)))
colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, parallel = FALSE)
res <- results(dds1,contrast = c("condition","OE","WT"))

res2 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE) %>%
  na.omit() %>% arrange(padj) %>%
  rownames_to_column("geneid") %>%
  #rownames_to_column("Peak") %>%
  merge(countData%>% rownames_to_column("geneid"),by="geneid")%>% 
  arrange(padj)
write.csv(res2,file = "GSE86084.M.ESC.CBX2.WT2.DESeq2.csv",row.names = F)

ge = res2$log2FoldChange
names(ge) = toupper(res2$geneid)
ge = sort(ge,decreasing = T)
library(ChIPseeker)
library(clusterProfiler)
Hallmark <- read.gmt("../../../../TCGA/msigdb/h.all.v7.5.1.symbols.gmt")
KEGG <- read.gmt("../../../../TCGA/msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
GO <- read.gmt("../../../../TCGA/msigdb/c5.go.v7.5.1.symbols.gmt")

em <- GSEA(ge, TERM2GENE = Hallmark,pvalueCutoff = 1,pAdjustMethod = "BH",eps = 1e-100)
write.csv(em@result,file = "GSE86084.M.ESC.CBX2.WT2.Hallmark.GSEA.csv",row.names = F)

library(GseaVis)
gseaNb(object = em,
       geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
       addPval = T,
       pvalX = 0.6,pvalY = 0.75,
       pCol = 'black',
       pHjust = 0,subPlot = 2)


em <- GSEA(ge, TERM2GENE = KEGG,pvalueCutoff = 1,pAdjustMethod = "BH",eps = 1e-100)
write.csv(em@result,file = "GSE86084.M.ESC.CBX2.WT2.KEGG.GSEA.csv",row.names = F)

em <- GSEA(ge, TERM2GENE = GO,pvalueCutoff = 0.5,pAdjustMethod = "BH",eps = 1e-100)
write.csv(em@result,file = "GSE86084.M.ESC.CBX2.WT2.GO.GSEA.csv",row.names = F)

library(org.Mm.eg.db)
geneList = res2$log2FoldChange
names(geneList) = res2$geneid
geneList = sort(geneList,decreasing = T)
s2g=select(org.Mm.eg.db,names(geneList),
           'ENTREZID','SYMBOL')
s2g=s2g[!is.na(s2g$ENTREZID),]
s2g=s2g[s2g$SYMBOL %in% names(geneList),]
geneList=geneList[s2g$SYMBOL] 
names(geneList)=s2g$ENTREZID

em <- gseGO(geneList = geneList, 
            OrgDb='org.Mm.eg.db',
            ont = 'ALL',
            #nPerm = 1000,
            #minGSSize = 10,
            eps=1e-50,
            pvalueCutoff = 0.5,
            by="fgsea",
            verbose = FALSE)
em <- gseKEGG(geneList = geneList, 
              #OrgDb='org.Mm.eg.db',
              organism = "mmu",
              #ont = 'ALL',
              #nPerm = 1000,
              #minGSSize = 10,
              pvalueCutoff = 0.05,
              verbose = FALSE)





##########################################

##################### Cell cycle + MAPK pathway => gene list ##################
setwd("K:/TCGA/Anlysis/LIHC")

library(clusterProfiler)
KEGG <- read.gmt("../../msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
KEGG.MAPK <- KEGG %>% filter(str_detect(term,"MAPK"))
write.table(KEGG.MAPK,file = "KEGG.MAPK.Pathway.txt",sep = "\t",row.names = F,quote = F)

KEGG.CellCycle <- KEGG %>% filter(str_detect(term,"CELL")) %>%
  filter(str_detect(term,"CYCLE"))
write.table(KEGG.CellCycle,file = "KEGG.CellCycle.Pathway.txt",sep = "\t",row.names = F,quote = F)


##################### TCGA CBX2 CEP55 => pathway correlation ##########################
setwd("K:/HepG2-CBX2/Validate")
library(tidyverse)
Exp <- read.csv("../../../TCGA/mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
Exp.Log <- log2(Exp+1)
Phenodata <- read.csv("../../../TCGA/mRNA.Phenotype.csv") %>% filter(Project.ID == "TCGA-LIHC")

mid.phenodata <- Phenodata %>% filter(Sample.ID %in% colnames(Exp)) %>% 
  filter(Sample.Type == "Solid Tissue Normal" | Sample.Type == "Primary Tumor") %>%
  mutate(Group = case_when(Sample.Type == "Solid Tissue Normal" ~ "NTL",
                           Sample.Type == "Primary Tumor" ~ "HCC")) %>%
  #mutate(Group = factor(.$Group,levels = c("NTL","HCC"))) %>%
  dplyr::select(Case.ID:Group)

Pathway.Score <- read.csv("../../../TCGA/GSCA/LIHC.PAS.tsv",sep = "\t")
Pathway.Score.D <- Pathway.Score %>% dplyr::select(-cancer_type) %>%
  spread(pathway,score) 

#GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
#GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180
CBX2.CEP55.Exp <- Exp.Log[c("ENSG00000173894","ENSG00000138180"),] %>%
  t() %>% data.frame(check.names = F) %>%
  magrittr::set_colnames(c("CBX2","CEP55")) %>%
  rownames_to_column("Sample") %>%
  mutate(PatientID=substr(.$Sample,1,12)) %>%
  merge(Pathway.Score.D,by.x="PatientID",by.y="barcode") %>%
  filter(!str_detect(Sample,"-11A"))
Cor.data <- data.frame()
for (i in c("CBX2","CEP55")) {
  for (j in colnames(CBX2.CEP55.Exp)[5:14]) {
    test <- cor.test(CBX2.CEP55.Exp[[i]],CBX2.CEP55.Exp[[j]],method = "spearman",exact = F)
    Cor.data <- data.frame(Gene=i,
                           Pathway=j,
                           Pvalue=test$p.value,
                           Rho=test$estimate) %>%
      rbind(Cor.data)
  }
}

write.table(Cor.data,"CBX2.CEP55.PathwayScore.txt",sep = "\t",row.names = F)
Cor.data <- read.csv("CBX2.CEP55.PathwayScore.txt",sep = "\t")
Cor.data2 <- Cor.data %>% mutate(Size=abs(Rho))
Cor.data3 <- Cor.data2%>%filter(Pvalue<=0.05)
P43<-ggplot()+
  geom_point(data=Cor.data2,aes(x=Gene,y=Pathway,fill=Rho,size=Size),color="white",shape=21)+
  #geom_tile(data=dat_all,aes(SNHGs,Pathway),fill="white",color=NA)+
  #geom_point(data=dat_all,aes(SNHGs,Pathway,size=-log10(p.adjust),color=NES),shape=20)+
  scale_size_continuous(range = c(1,8))+
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",midpoint = 0)+ #,,limits = c(-2, 2)
  geom_point(data=Cor.data3,aes(Gene,Pathway,size=Size,alpha=5,fill=Rho,stroke=2),color="black",shape=21)+
  scale_size_continuous(range = c(1,8))+
  ggthemes::theme_few()+
  theme(#legend.position = "top",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")
ggsave(P43,filename = "CBX2.PathwayScore.Correlation.pdf",width = 4,height = 5)



##################### Tumor-ATAC-BW => Figure plot ###############
setwd("K:/HepG2-CBX2/HepG2.ATAC")
library(transPlotR)
library(rtracklayer)
file <- list.files(pattern = '.bw$')
file <- file[!str_detect(file,"SRR")]
mybw <- loadBigWig(file)
gtf <- import('hg19.ncbiRefSeq.gtf',format = "gtf") %>%
  data.frame()
## https://mp.weixin.qq.com/s/M7iMdi43vcgrINDoympUdA
RAD23A <-
  trackVis(bWData = mybw,
           gtf.file = gtf,
           gene.name = "RAD23A",
           extend.up = 1000,
           extend.dn = 1000,
           xAxis.info = F,
           theme = "bw",
           yAxis.info = F,new.yaxis = T,
           pos.ratio = c(0.06,0.8),
           color = jjAnno::useMyCol(platte = "stallion",n = 6))

# peak
bd <- bedVis(bdFile = file,
             chr = "chr19",
             track.width = 0.3,
             show.legend = T)

# combine
RAD23A %>% insert_bottom(p,height = 0.3) %>%
  insert_bottom(bd,height = 0.15)
############# CBX2
CBX2 <-
  trackVis(bWData = mybw,
           gtf.file = gtf,
           gene.name = "CBX2",
           extend.up = 500,
           extend.dn = 500,
           xAxis.info = F,
           theme = "bw",
           scales = "fixed",
           yAxis.info = F,new.yaxis = T,
           pos.ratio = c(0.06,0.8),
           yinfo.text.size = 3,
           color = rep((jjAnno::useMyCol(platte = "stallion",n = 2))[2:1],c(9,9)))

p <-
  trancriptVis(gtfFile = gtf,
               gene = "CBX2",
               relTextDist = -0.5,
               exonWidth = 0.5,
               exonFill = (jjAnno::useMyCol(platte = "stallion",n = 5))[5],
               arrowCol = (jjAnno::useMyCol(platte = "stallion",n = 5))[5],
               textLabelSize = 2,
               addNormalArrow = FALSE,
               newStyleArrow = TRUE,
               xAxis.info = FALSE,
               textLabel = 'gene_name') +
  xlab('')

p2 <- CBX2 %>% aplot::insert_bottom(p,height = 0.1)
ggsave(p2,file="CBX2.ATAC.Tumor.pdf",height=4,width=3)

######### CEP55
CEP55 <-
  trackVis(bWData = mybw,
           gtf.file = gtf,
           gene.name = "CEP55",
           extend.up = 2000,
           extend.dn = 2000,
           xAxis.info = F,
           theme = "bw",
           scales = "fixed",
           yAxis.info = F,new.yaxis = T,
           pos.ratio = c(0.06,0.8),
           yinfo.text.size = 3,
           color = rep((jjAnno::useMyCol(platte = "stallion",n = 2))[2:1],c(9,9)))

p.1 <-
  trancriptVis(gtfFile = gtf,
               gene = "CEP55",
               relTextDist = -0.5,
               exonWidth = 0.5,
               exonFill = (jjAnno::useMyCol(platte = "stallion",n = 5))[5],
               arrowCol = (jjAnno::useMyCol(platte = "stallion",n = 5))[5],
               textLabelSize = 2,
               addNormalArrow = FALSE,
               newStyleArrow = TRUE,
               xAxis.info = FALSE,
               textLabel = 'gene_name') +
  xlab('')

p3 <- CEP55 %>% aplot::insert_bottom(p.1,height = 0.1)
ggsave(p3,file="CEP55.ATAC.Tumor.pdf",height=4,width=3)

##################### Mouse => HCC => GSE94583 #######################
setwd("K:/HepG2-CBX2/Validate/GSE94583_RAW")
library(R.utils)
library(tidyverse)
library(magrittr)
#fileNames <- list.files("./")
#sapply (fileNames, gunzip)
Count.Data <- data.frame()
files <- list.files("./")
for(i in 1:length(files)){
  filename <- str_remove(files[i],"_.*")
  Count.Data <- read.table(file = files[i]) %>%
    set_colnames(c("GeneName",filename)) %>%
    arrange(GeneName) %>% remove_rownames() %>%
    column_to_rownames("GeneName") %>% t() %>%
    rbind(Count.Data)
  print(dim(Count.Data))
}

Count.Data <- Count.Data %>% t()
Count.Data2 <- Count.Data[-c(1:5),]
write.csv(Count.Data2,file = "GSE94583.Count.csv")

Counts <- apply(Count.Data2,2,sum)

library(readxl)
Metadata <- read_xlsx("Metadata.xlsx")

CPM.Data <- (Count.Data2[c("Cbx2","Cep55"),]/Counts) %>% t() %>%
  data.frame() %>% rownames_to_column("Sample") %>% 
  merge(Metadata,by="Sample")

library(ggpubr)
library(ggplot2)
CPM.Data$Attri <- factor(CPM.Data$Attri,levels = c("Normal","PNR","LCS","LCS_Metastasis",
                                                   "PPTR","PPTR_Metastasis"))
p1<-ggplot(CPM.Data %>% filter(Tissue != "Organoids"),aes(Attri,Cbx2,fill=Attri))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="CPM of Cbx2")
  
ggsave(p1,filename = "GSE94583.Cbx2.pdf",height = 4,width = 4)
  

p2<-ggplot(CPM.Data %>% filter(Tissue != "Organoids"),aes(Attri,Cep55,fill=Attri))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  ggthemes::theme_few()+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="CPM of Cep55")+
  scale_y_sqrt()

ggsave(p2,filename = "GSE94583.Cep55.pdf",height = 4,width = 4)

##################### Mouse => HCC => GSE116463 ##################
setwd("K:/HepG2-CBX2/Validate/GSE116463_RAW")
#fileNames <- list.files("./")
#sapply (fileNames, untar)
files <-  list.files(".",recursive = T)
files2 <- files[str_detect(files,"genes.fpkm_tracking")]
RPKM.data <- data.frame()
for (file in files2) {
  filename=str_remove(file,"_.*")
  RPKM.data <- read.table(file,header = T,sep = "\t") %>%
    dplyr::select(gene_short_name,FPKM) %>%
    filter(gene_short_name %in% c("Cbx2","Cep55")) %>%
    mutate(Sample=filename) %>%
    rbind(RPKM.data)
  print(dim(RPKM.data))
}

RPKM.data$Group <- rep(c("Normal_MC","Tumor_MC","Tumor_LC","Normal_LC","Tumor_Chow","Normal_Chow"),
                       c(10,10,10,10,8,10))
RPKM.data$Group <- factor(RPKM.data$Group,levels = c("Normal_Chow","Tumor_Chow","Normal_MC","Tumor_MC","Normal_LC","Tumor_LC"))
write.csv(RPKM.data,file = "GSE116463.RPKM.Cbx2.csv",row.names = F)

p3 <- ggplot(RPKM.data,aes(Group,FPKM,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  ggthemes::theme_few()+
  facet_wrap(~gene_short_name,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="FPKM")+
  theme(strip.text = element_text(size=15,color="black"))
  
ggsave(p3,filename = "Mouse.GSE116463.FPKM.Cbx2.pdf",height = 2.5,width = 4)

##################### Human => HCC => GSE140462 ##################
setwd("K:/HepG2-CBX2/Validate/GSE140462_RAW")
library(R.utils)
#fileNames <- list.files("./")
#sapply (fileNames, gunzip)
files <- list.files(".")
count.data <- data.frame()
for (file in files) {
  data <- read.table(file,header = T)
  filename <- str_remove(file,"_.*")
  count.data <- data[,c(1,7)] %>% data.frame() %>%
    remove_rownames() %>% 
    column_to_rownames("Geneid") %>%
    set_colnames(filename) %>% t() %>%
    rbind(count.data)
  print(dim(count.data))
}

count.data <- count.data %>% t()
length.gene.kb <- data$Length/1000
rpk <- count.data/length.gene.kb
fpkm <- t(t(rpk)/colSums(count.data)* 10^6)
write.csv(fpkm,file = "GSE140462.FPKM.csv")

fpkm.cbx2 <- fpkm[c("CBX2","CEP55"),] %>% t() %>%
  data.frame() %>% rownames_to_column("Sample") %>%
  gather(gene,FPKM,-Sample) %>%
  mutate(Group=rep(c("HCC","Normal","HCC","Normal"),c(7,7,7,7))) %>%
  mutate(Group=factor(Group,levels = c("Normal","HCC"))) %>%
  mutate(Line=rep(1:7,length=28))


p4 <- ggplot(fpkm.cbx2,aes(Group,FPKM,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="FPKM")+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p4,filename = "Human.GSE140462.FPKM.CBX2.pdf",height = 2.5,width = 4)

##################### Mouse => HCC => GSE153580 ##############
setwd("K:/HepG2-CBX2/Validate/GSE153580")
files <- list.files(".",recursive = T)
files <- files[str_detect(files,"output")]
Count.data <- data.frame()
for (file in files) {
  filename=str_remove(file,"/.*")
  Count.data <- read.table(file,header = F,sep = "\t") %>%
    set_colnames(c("ENSEMBL","Count")) %>%
    mutate(ENSEMBL=str_remove(ENSEMBL,"\\..*")) %>%
    remove_rownames() %>%
    column_to_rownames("ENSEMBL") %>%
    set_colnames(filename) %>%
    t() %>%
    rbind(Count.data)
  print(dim(Count.data))
}

write.csv(Count.data,file = "Mouse.GSE153580.Count.csv")

Count.data <- Count.data %>% t()
Counts.sum <- apply(Count.data, 2, sum)
Count.data.CPM <- Count.data/Counts.sum
Count.data2 <- Count.data.CPM[c("ENSMUSG00000025577","ENSMUSG00000024989"),] %>% 
  t() %>% set_colnames(c("Cbx2","Cep55")) %>% data.frame(check.names = F) %>%
  rownames_to_column("Sample") %>%
  gather(gene,FPKM,-Sample) %>%
  mutate(Group=rep(c("STZ_12w","STZ_10w","STZ_6w","Ctrl_12w","Ctrl_10w","Ctrl_6w"),each=3,times=2)) %>%
  mutate(Group=factor(Group,levels = c("Ctrl_6w","STZ_6w","Ctrl_10w","STZ_10w","Ctrl_12w","STZ_12w")))

p5 <- ggplot(Count.data2,aes(Group,FPKM,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="CPM")+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p5,filename = "Human.GSE140462.CPM.CBX2.pdf",height = 2.5,width = 4) 
  
  

##################### Mouse => HCC ÷Œ¡∆ => GSE179352 ############
setwd("K:/HepG2-CBX2/Validate/GSE179352")

data <- read.csv("GSE179352_FPKM.txt",header = T,sep = "\t")
data2 <- data %>% filter(gene_id %in% c("Cbx2","Cep55")) %>%
  gather(Sample,FPKM,-gene_id) %>%
  mutate(Group = rep(c("Ctrl","IFN-¶¡","IFN-¶¡+anti-PD-1","anti-PD-1"),c(6,6,6,6))) %>%
  mutate(Group=factor(Group,levels = c("Ctrl","IFN-¶¡","anti-PD-1","IFN-¶¡+anti-PD-1")))

p6 <- ggplot(data2,aes(Group,FPKM,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene_id,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="FPkm")+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p6,filename = "Mouse.GSE179352.FPKM.CBX2.pdf",height = 2.5,width = 4) 


##################### Human => HCC => GSE201868 => No ########
setwd("K:/HepG2-CBX2/Validate/GSE201868")
data <- read.table("Tongji_cohort_Metadata.txt")

##################### Human => HCC => GSE62905 = CD133 #############
setwd("K:/HepG2-CBX2/Validate/GSE62905")
files <- list.files(".",pattern = "txt$")
data <- data.frame()
for (file in files) {
  filename <- paste((str_remove(file,"\\..*") %>% str_split("_"))[[1]][2:3],collapse = "_")
  data <- read.table(file,sep = "\t",header = T) %>%
    filter(gene_id %in% c("CBX2","CEP55")) %>%
    set_colnames(c("gene","locus","FPKM")) %>%
    mutate(Sample=filename) %>%
    rbind(data)
}
write.csv(data,"Human.GSE62905.FPKM.CBX2.csv",row.names = F)

p7 <- ggplot(data,aes(Sample,FPKM,fill=Sample,color=Sample))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=4)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="FPkm")+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p7,filename = "Human.GSE62905.FPKM.CBX2.pdf",height = 2,width = 4) 



##################### Human => HCC => GSE23450 = CD133 ########
setwd("K:/HepG2-CBX2/Validate/HCC.CSC")
GPL <- read.csv("../GPL570-55999.csv") %>%
  filter(str_detect(Gene.Symbol,"CBX2") | str_detect(Gene.Symbol,"CEP55"))
data <- read.table("GSE23450_series_matrix.txt",header = T,sep = "\t") %>%
  filter(ID_REF %in% GPL$ID)
data2 <- data %>% mutate(Symbol=GPL$Gene.Symbol) %>%
  gather(Sample,Expression,GSM575346:GSM575349) %>%
  mutate(Group=rep(c("Huh7_CD133+_MAS5","Huh7_CD133-_MAS5","Huh7_CD133+_RMA","Huh7_CD133-_RMA"),each=5)) %>%
  mutate(Group=factor(Group,levels=c("Huh7_CD133+_MAS5","Huh7_CD133-_MAS5","Huh7_CD133+_RMA","Huh7_CD133-_RMA"))) %>%
  mutate(ID_REF=paste(Symbol,ID_REF,sep="."))
  

p8 <- ggplot(data2,aes(Group,Expression,fill=Group,color=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=3)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~ID_REF,scales = "free_x",ncol = 5)+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="Expression")+
  scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p8,filename = "Human.GSE23450.FPKM.CBX2.pdf",height = 3,width = 12) 

##################### Human => HCC => GSE23451 = CD133 ########
setwd("K:/HepG2-CBX2/Validate/HCC.CSC")
GPL <- read.csv("../GPL570-55999.csv") %>%
  filter(str_detect(Gene.Symbol,"CBX2") | str_detect(Gene.Symbol,"CEP55"))
data <- read.table("GSE23451_series_matrix.txt",header = T,sep = "\t") %>%
  filter(ID_REF %in% GPL$ID)
data2 <- data %>% mutate(Symbol=GPL$Gene.Symbol) %>%
  gather(Sample,Expression,GSM575351:GSM575355) %>%
  mutate(Group=rep(c("PLC8024_CD133+_MAS5","PLC8024_CD133-_MAS5","PLC8024_CD133+_RMA","PLC8024_CD133-_RMA"),each=5)) %>%
  mutate(Group=factor(Group,levels=c("PLC8024_CD133+_MAS5","PLC8024_CD133-_MAS5","PLC8024_CD133+_RMA","PLC8024_CD133-_RMA"))) %>%
  mutate(ID_REF=paste(Symbol,ID_REF,sep="."))


p9 <- ggplot(data2,aes(Group,Expression,fill=Group,color=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=3)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~ID_REF,scales = "free_x",ncol = 5)+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="Expression")+
  scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p9,filename = "Human.GSE23451.FPKM.CBX2.pdf",height = 3,width = 12) 

##################### Human => HCC => GSE56771 = CD133 ##########
setwd("K:/HepG2-CBX2/Validate/HCC.CSC")
GPL <- read.csv("../GPL6244-17930.csv") %>% 
  filter(str_detect(gene_assignment,"CBX2") | str_detect(gene_assignment,"CEP55"))
data <- read.table("GSE56771_series_matrix.txt",header = T,sep = "\t") %>%
  filter(ID_REF %in% GPL$ID)
data2 <- data %>% mutate(Symbol=c("CEP55","CBX2")) %>%
  gather(Sample,Expression,GSM1368818:GSM1368821) %>%
  mutate(Group=rep(c("Huh7_CD133+","Huh7_CD133-","Huh7_ALDH+","Huh7_ALDH-"),each=2)) %>%
  mutate(Group=factor(Group,levels=c("Huh7_CD133+","Huh7_CD133-","Huh7_ALDH+","Huh7_ALDH-")))

p9 <- ggplot(data2,aes(Group,Expression,fill=Group,color=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=3)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~Symbol,scales = "free_x",ncol = 5)+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="Expression")+
  scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p9,filename = "Human.GSE56771.FPKM.CBX2.pdf",height = 2.5,width = 5)

##################### Human => HCC => GSE84223 = CD44 ####################
files <- list.files("./GSE84223_RAW/",pattern = "xls$")
library(readxl)
data.d <- data.frame()
for (file in files) {
  filename=str_remove(file,"_.*")
  data.d <- read_xls(file.path("./GSE84223_RAW/",file)) %>% filter(symbol %in% c("CBX2","CEP55")) %>%
    dplyr::select(symbol,`global normalization`) %>%
    set_colnames(c("gene","Expression")) %>%
    mutate(Sample=filename) %>%
    rbind(data.d)
}
data.d2 <- data.d %>% mutate(Group=rep(c("Huh7_CD44+","Huh7_CD44-","Normal_Hepatocytes"),each=2))
data.d2$Group <- factor(data.d2$Group,levels = c("Normal_Hepatocytes","Huh7_CD44+","Huh7_CD44-"))
p10 <- ggplot(data.d2,aes(Group,Expression,fill=Group,color=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=3)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_x",ncol = 5)+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="Expression")+
  scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p10,filename = "Human.GSE84223.FPKM.CBX2.pdf",height = 2.5,width = 5)




##################### Human => HCC => GSE112705 RNA+Ribo #######################
setwd("K:/HepG2-CBX2/Validate/GSE112705")
data <- read.table("GSE112705_RPF_RNA_readCounts.txt",header = T,sep = "\t",fill = 0)
data %>% filter(geneSymbol %in% c("CBX2","CEP55"))
data2 <- data %>% dplyr::select(-geneSymbol) %>%
  remove_rownames() %>% column_to_rownames("geneID")
data2[is.na.data.frame(data2)]=0
Count.s <- apply(data2, 2, sum)
data2.CPM <- data2/Count.s
data2.CPM2 <- data2.CPM[c("ENSG00000138180.15","ENSG00000173894.10"),] %>%
  mutate(gene=c("CBX2","CEP55")) %>%
  gather(Sample,FPKM,-gene) %>%
  mutate(Patient=str_remove(Sample,"\\..*"),
         Seq=str_remove(Sample,".*\\."))
Str <- unlist(str_split(data2.CPM2$Sample,"\\.")) 
data2.CPM2$Group=R.utils::capitalize(Str[seq(2,length(Str),3)])
data2.CPM2$Group <- factor(data2.CPM2$Group,levels = c("Tumor","Normal"))
p11 <- ggplot(data2.CPM2,aes(Group,FPKM,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Patient))+
  ggthemes::theme_few()+
  facet_grid(Seq~gene,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="CMP")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p11,filename = "Human.GSE112705.FPKM.CBX2.pdf",height = 4,width = 4)

##################### Human => HCC => GSE212604 iRFA ÷Œ¡∆ ############
library(tidyverse)
library(magrittr)
library(ggplot2)
setwd("K:/HepG2-CBX2/Validate/GSE212604")
data <- read.table("GSE212604.iRFA.count.txt",header=T,row.names = 1)
gene.length <- data$Length
data2 <- data %>% dplyr::select(contains("bam"))
colnames(data2) <- str_remove(colnames(data2),"\\..*")

kb<-gene.length/1000
rpk <- data2/kb
fpkm <- t(t(rpk)/colSums(data2) * 10^6)
write.csv(fpkm,file="Human.GSE212604.iRFA.FPKM.csv")

fpkm.cbx2 <- fpkm[c("CBX2","CEP55"),] %>% data.frame() %>%
  rownames_to_column("gene") %>%
  gather(Sample,FPKM,SRR21403484:SRR21403503) %>%
  mutate(Group=rep(c("iRFA","HCC"),c(20,20)))
  
p12 <- ggplot(fpkm.cbx2,aes(Group,FPKM,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Patient))+
  ggthemes::theme_few()+
  facet_grid(~gene,scales = "free_x")+
  coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  ggsci::scale_color_nejm()+
  labs(x="",y="FPKM")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))

ggsave(p12,filename = "Human.GSE212604.FPKM.iRFA.CBX2.pdf",height = 2.5,width = 4)


##################### Human => HCC => GSE104580 TACE ÷Œ¡∆ ################
setwd("K:/HepG2-CBX2/Validate/GSE104580")
GPL <- read.csv("../GPL570-55999.csv") %>%
  filter(str_detect(Gene.Symbol,"CBX2") | str_detect(Gene.Symbol,"CEP55"))

data <- read.csv("GSE104580_series_matrix.txt",header = T,sep = "\t") %>%
  filter(ID_REF %in% GPL$ID)
data2 <- data %>% mutate(Symbol=GPL$Gene.Symbol) %>%
  gather(Sample,Expression,GSM2803655:GSM2803801) %>%
  mutate(Group=rep(c("Response","Non-Response"),c(405,330))) %>%
  mutate(ID_REF=paste(Symbol,ID_REF,sep=".")) %>%
  group_by(Group,Sample,Symbol) %>%
  mutate(Expression = mean(Expression))

p13 <- ggplot(data2,aes(Group,Expression,fill=Group))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12,angle = 90),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Expression")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("Response","Non-Response")),
                             method = "t.test")

R.cbx <- data2 %>% filter(Group=="Response")
N.R.cbx <- data2 %>% filter(Group=="Non-Response")
ggsave(p13,filename = "Human.GSE104580.TACE.FPKM.CBX2.pdf",height = 4,width = 5) 


##################### HCC celline CBX2+CEP55 ####################
setwd("K:/HepG2-CBX2/Validate")
data <- read.table("RNAseqFPKM.txt",sep = "\t",header = T,check.names = F,row.names = 1)
data2 <- data %>% filter(Gene_name %in% c("CBX2","CEP55","AFP")) %>%
  remove_rownames() %>% 
  column_to_rownames("Gene_name") %>%
  t() %>% data.frame() %>%
  rownames_to_column("Cellline")

library(ggpubr)
library(tidyverse)
ggplot(data = data2,aes(x=log2(AFP+1),y=log2(CEP55+1),color=Cellline))+ #
  geom_point()+
  geom_smooth(formula = y~x,method = "lm")+ #
  #facet_wrap(~Project,nrow = 3,ncol = 7)+
  stat_cor(method = "pearson")+ #,label.x = 200,label.y = 2,size=10
  ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  #scale_color_manual(values = cell_type_cols)+
  labs(x="AFP log2(FPKM+1)",y="CEP55 log2(FPKM+1)")

cor.test(log2(data2$AFP+1),log2(data2$CBX2+1),method = "spearman")
cor.test(log2(data2$AFP+1),log2(data2$CEP55+1),method = "spearman")

p14<-ggdotchart(data2, x = "Cellline", y = "CBX2",
           color = "Cellline",                               
           sorting = "ascending",                        
           add = "segments",                             
           ylab="CBX2 FPKM", xlab="Cell line", 
           rotate = TRUE,
           dot.size = 3 
           
)+ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
ggsave(p14,filename = "HCC.CellLine.CBX2.pdf",height = 7,width = 3)

p15<-ggdotchart(data2, x = "Cellline", y = "CEP55",
                color = "Cellline",                               
                sorting = "ascending",                        
                add = "segments",                             
                ylab="CEP55 FPKM", xlab="Cell line", 
                rotate = TRUE,
                dot.size = 3 
                
)+ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))
ggsave(p15,filename = "HCC.CellLine.CEP55.pdf",height = 7,width = 3)



##################### Human => HCC => GSE112221 RRBS => No use ##################
data <- read.csv("GSE112221_RRBS_TAB.txt",sep = "\t")
data2 <- data %>% filter(Gene %in% c("CBX2","CEP55"))

CEP55.M <- data2 %>% filter(Gene=="CEP55") %>%
  dplyr::select(X756:X674T)
P16 <- apply(CEP55.M,2, function(x){sum(as.numeric(x))}) %>% data.frame() %>%
  set_colnames("Methylation level") %>%
  rownames_to_column("Sample") %>%
  mutate(Group=c("Normal","Normal","Cirrhosis","HCC","Cirrhosis","HCC","Cirrhosis","HCC","Cirrhosis","HCC")) %>%
  mutate(Group=factor(Group,levels = c("Normal","Cirrhosis","HCC"))) %>%
  ggplot(aes(Group,`Methylation level`,color=Group)) +
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level")
ggsave(P16, filename = "RRBS.CEP55.pdf",height = 2.5,width = 2.5)
  
CBX2.M <- data2 %>% filter(Gene=="CBX2") %>%
  dplyr::select(X756:X674T)
P17 <- apply(CBX2.M,2, function(x){sum(na.omit(as.numeric(x)))}) %>% data.frame() %>%
  set_colnames("Methylation level") %>%
  rownames_to_column("Sample") %>%
  mutate(Group=c("Normal","Normal","Cirrhosis","HCC","Cirrhosis","HCC","Cirrhosis","HCC","Cirrhosis","HCC")) %>%
  mutate(Group=factor(Group,levels = c("Normal","Cirrhosis","HCC"))) %>%
  ggplot(aes(Group,`Methylation level`,color=Group)) +
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level")
ggsave(P17, filename = "RRBS.CBX2.pdf",height = 2.5,width = 2.5)


##################### Mouse => HCC => GSE77725 5mC Probe => No use ##########################
library(tidyverse)
setwd("K:/HepG2-CBX2/Validate/GSE77725")
GPL <- read.csv("GPL16743-30696.txt",sep = "\t",skip = "#") %>%
  filter(str_detect(GENE_SYMBOL,"Cbx2") | str_detect(GENE_SYMBOL,"Cep55"))

data <- read.csv("GSE77725_MeDIP_NASH_study.txt",sep = "\t") %>%
  filter(probe_id %in% GPL$ID)

data2 <- data %>% merge(.,GPL %>% dplyr::select(ID,GENE_SYMBOL),by.x="probe_id",by.y="ID") %>%
  group_by(GENE_SYMBOL) %>% 
  mutate(WT1=sum(WT1_MC_NASHstudy),
         WT2=sum(WT2_MC_NASHstudy),
         Tumor1=sum(Tum1_MC_NASHstudy),
         Tumor2=sum(Tum2_MC_NASHstudy)) %>%
  dplyr::select(WT1:Tumor2) %>% unique()



##################### Human => HCC => GSE63775 MCTA => No use #####################
setwd("K:/HepG2-CBX2/Validate/GSE63775")
files <- list.files(".",pattern = "bed$")
Count.data <- data.frame()
for (file in files) {
  filename <- str_remove(file,"\\..*")
  Count.data <- read.table(file,sep = "\t") %>% mutate(Pos=paste(V1,V2,V3,sep="_")) %>%
    remove_rownames() %>% column_to_rownames("Pos") %>%
    dplyr::select(V4) %>% mutate(MePM=(V4*1e6/sum(V4))) %>%
    dplyr::select(MePM) %>%
    set_colnames(filename) %>% 
    t() %>% rbind(Count.data)
  #print(dim(data))
}

Count.data2 <- Count.data %>% t() %>% data.frame()

Count.data2[c("chr17_77751378_77751673","chr17_77752051_77752772","chr10_95256133_95256562"),] %>%
  t() %>% data.frame() %>% mutate(Group=rep(c("Normal","Adjacent","HCC"),c(5,27,27)))



##################### º◊ª˘ªØ–æ∆¨ =>  GSE18081 ∏Œ∞© ∏Œ”≤ªØ => No use #################
##################### º◊ª˘ªØ–æ∆¨ =>  GSE113017 => #################
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL13534-11288.txt",sep = "\t",skip = "#") %>% 
  filter(str_detect(UCSC_RefGene_Name,"CBX2") | str_detect(UCSC_RefGene_Name,"CEP55"))
GPL.ID <- GPL %>% filter(!str_detect(UCSC_RefGene_Group,"Body"))

data <- read.csv("GSE113017_series_matrix2.txt",sep = "\t")
# metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>% set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(UCSC_RefGene_Name,"CBX2"))
CEP55.ID <- GPL.ID %>% filter(str_detect(UCSC_RefGene_Name,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p33.data <- apply(data2.CBX2,2,mean) %>% data.frame() %>% set_colnames("Methylation_level") %>%
  rownames_to_column("Sample") %>%
  mutate(Group2=rep(c("NTL","HCC"),times=30))

p33 <- ggstatsplot::ggbetweenstats(
  data = p33.data,
  x = Group2,
  y = Methylation_level,
  pairwise.comparisons = F,
  p.adjust.method = "fdr",
  type = "robust",
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Attri,scales = "free_y",ncol = 3)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)

ggsave(p33,filename = "Human.GSE113019.CBX2.Methylation.pdf",height = 3.1,width = 2.8)


data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p34.data <- apply(data2.CEP55,2,function(x){mean(na.omit(x))}) %>% data.frame() %>% set_colnames("Methylation_level") %>%
  rownames_to_column("Sample") %>%
  mutate(Group2=rep(c("NTL","HCC"),times=30))

p34 <- ggstatsplot::ggbetweenstats(
  data = p34.data,
  x = Group2,
  y = Methylation_level,
  pairwise.comparisons = F,
  p.adjust.method = "fdr",
  type = "robust",
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Attri,scales = "free_y",ncol = 3)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)

ggsave(p34,filename = "Human.GSE113019.CEP55.Methylation.pdf",height = 3.1,width = 2.8)


##################### º◊ª˘ªØ–æ∆¨ =>  GSE113019 => Right #################
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL13534-11288.txt",sep = "\t",skip = "#") %>% 
  filter(str_detect(UCSC_RefGene_Name,"CBX2") | str_detect(UCSC_RefGene_Name,"CEP55"))
GPL.ID <- GPL %>% filter(!str_detect(UCSC_RefGene_Group,"Body"))

data <- read.csv("GSE113019_series_matrix.txt",sep = "\t")
metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>%
  set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(UCSC_RefGene_Name,"CBX2"))
CEP55.ID <- GPL.ID %>% filter(str_detect(UCSC_RefGene_Name,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p17 <- apply(data2.CBX2,2,mean) %>% data.frame() %>% set_colnames("Methylation_level") %>%
  rownames_to_column("Sample") %>% merge(metadata,by="Sample") %>%
  mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",
                 if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%
  ggplot(aes(Group2,Methylation_level,fill=Group2))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("Adjacent","PrimaryT"),c("Adjacent","ReccurentT")),
                             method = "t.test")

ggsave(p17,filename = "Human.GSE113019.CBX2.Methylation.pdf",height = 3.2,width = 3.5)

data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p18 <- apply(data2.CEP55,2,mean) %>% data.frame() %>% set_colnames("Methylation_level") %>%
  rownames_to_column("Sample") %>% merge(metadata,by="Sample") %>%
  mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",
                         if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%
  ggplot(aes(Group2,Methylation_level,fill=Group2))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("Adjacent","PrimaryT"),c("Adjacent","ReccurentT")),
                             method = "t.test")

ggsave(p18,filename = "Human.GSE113019.CEP55.Methylation.pdf",height = 3.2,width = 3.5)

##################### º◊ª˘ªØ–æ∆¨ Human GSE136583 => GPL21145 => No sig ###############
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL21145_MethylationEPIC_15073387_v-1-0.csv")

GPL.ID <- GPL %>% filter(str_detect(UCSC_RefGene_Name,"CBX2") | str_detect(UCSC_RefGene_Name,"CEP55")) %>%
  filter(!str_detect(UCSC_RefGene_Group,"Body"))

data <- read.csv("GSE136583_series_matrix2.txt",sep = "\t")
#metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>% set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(UCSC_RefGene_Name,"CBX2")) #%>% filter(str_detect(Regulatory_Feature_Group,"Promoter"))
CEP55.ID <- GPL.ID %>% filter(str_detect(UCSC_RefGene_Name,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$IlmnID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p19 <- apply(data2.CBX2,2,mean) %>% data.frame() %>% set_colnames("Methylation_level") %>%
  rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  mutate(Group2 =rep(c("HCC","NTL"),c(31,31))) %>%
  ggplot(aes(Group2,Methylation_level,fill=Group2))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "t.test")

ggsave(p19,filename = "Human.GSE136583.CBX2.Methylation.NO.pdf",height = 3.2,width = 3.5)

data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$IlmnID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p20 <- apply(data2.CEP55,2,mean) %>% data.frame() %>% set_colnames("Methylation_level") %>%
  rownames_to_column("Sample") %>% 
  mutate(Group2 =rep(c("HCC","NTL"),c(31,31))) %>%
  #merge(metadata,by="Sample") %>%
  #mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%
  ggplot(aes(Group2,Methylation_level,fill=Group2))+
  geom_violin()+
  #geom_boxplot()+
  ggbeeswarm::geom_beeswarm(priority = "descending",size=2)+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~Symbol,scales = "free_y",ncol = 5)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "t.test")

ggsave(p20,filename = "Human.GSE136583.CEP55.Methylation.NO.pdf",height = 3.2,width = 3.5)

##################### º◊ª˘ªØ–æ∆¨ Human GSE37988 => GPL8490 =>  ###############
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL8490-65.txt",sep = "\t")

GPL.ID <- GPL %>% filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) 

data <- read.csv("GSE37988_series_matrix2.txt",sep = "\t")
#metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>% set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(Symbol,"CBX2")) #%>% filter(str_detect(Regulatory_Feature_Group,"Promoter"))
CEP55.ID <- GPL.ID %>% filter(str_detect(Symbol,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p21.data <- data2.CBX2 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =rep(c("HCC","NTL"),c(62,62))) #%>% group_by(Group2) %>% summarise(Mean=mean(Methylation_level))

library(ggstatsplot)

p21<-ggstatsplot::ggbetweenstats(
  data = p21.data,
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)

ggsave(p21,filename = "Human.GSE37988.CBX2.Methylation.pdf",height = 3.2,width = 2.8)

data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p22.data <- data2.CEP55 %>% data.frame() %>% #set_colnames("Methylation_level") %>%
  data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #rownames_to_column("Sample") %>% 
  mutate(Group2 =rep(c("HCC","NTL"),c(62,62))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  filter(Methylation_level>0)
  #merge(metadata,by="Sample") %>%
  #mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%

p22<-ggstatsplot::ggbetweenstats(
  data = p22.data,
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)
ggsave(p22,filename = "Human.GSE37988.CEP55.Methylation.NO.pdf",height = 3.2,width = 2.8)

##################### º◊ª˘ªØ–æ∆¨ Human GSE44909 => GPL8490 =>  ###############
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL8490-65.txt",sep = "\t")

GPL.ID <- GPL %>% filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) 

data <- read.csv("GSE44909_series_matrix2.txt",sep = "\t")
#metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>% set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(Symbol,"CBX2")) #%>% filter(str_detect(Regulatory_Feature_Group,"Promoter"))
CEP55.ID <- GPL.ID %>% filter(str_detect(Symbol,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p23.data <- data2.CBX2 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  #mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =c(rep(c("HCC","NTL"),times=12),rep("Normal",8))) #%>% group_by(Group2) %>% summarise(Mean=mean(Methylation_level))

library(ggstatsplot)

p23<-ggstatsplot::ggbetweenstats(
  data = p23.data %>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)

ggsave(p23,filename = "Human.GSE44909.CBX2.Methylation.pdf",height = 3.2,width = 2.8)

data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p24.data <- data2.CEP55 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  #mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =c(rep(c("HCC","NTL"),times=12),rep("Normal",8)))
#merge(metadata,by="Sample") %>%
#mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%

p24<-ggstatsplot::ggbetweenstats(
  data = p24.data%>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)
ggsave(p24,filename = "Human.GSE44909.CEP55.Methylation.pdf",height = 3.2,width = 2.8)

##################### º◊ª˘ªØ–æ∆¨ Human GSE57956 => GPL8490 =>  ###############
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL8490-65.txt",sep = "\t")

GPL.ID <- GPL %>% filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) 

data <- read.csv("GSE57956_series_matrix2.txt",sep = "\t") %>%
  dplyr::select(-GSM1398525,-GSM1398606,-GSM1398607,-GSM1398608)
#metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>% set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(Symbol,"CBX2")) #%>% filter(str_detect(Regulatory_Feature_Group,"Promoter"))
CEP55.ID <- GPL.ID %>% filter(str_detect(Symbol,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p25.data <- data2.CBX2 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  #mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =rep(c("NTL","HCC"),times=58)) #%>% group_by(Group2) %>% summarise(Mean=mean(Methylation_level))

library(ggstatsplot)

p25<-ggstatsplot::ggbetweenstats(
  data = p25.data %>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)

ggsave(p25,filename = "Human.GSE57956.CBX2.Methylation.pdf",height = 3.2,width = 2.8)

data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p26.data <- data2.CEP55 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  #mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =rep(c("NTL","HCC"),times=58))
#merge(metadata,by="Sample") %>%
#mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%

p26<-ggstatsplot::ggbetweenstats(
  data = p26.data%>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)
ggsave(p26,filename = "Human.GSE57956.CEP55.Methylation.NO.pdf",height = 3.2,width = 2.8)


##################### º◊ª˘ªØ–æ∆¨ Human GSE73003 => GPL8490 =>  ###############
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL8490-65.txt",sep = "\t")

GPL.ID <- GPL %>% filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) 

data <- read.csv("GSE73003_series_matrix2.txt",sep = "\t")
#metadata <- read.csv("GSE113019.Sample.txt",sep = "\t",header = F) %>% set_colnames(c("Sample","Group"))
CBX2.ID <- GPL.ID %>% filter(str_detect(Symbol,"CBX2")) #%>% filter(str_detect(Regulatory_Feature_Group,"Promoter"))
CEP55.ID <- GPL.ID %>% filter(str_detect(Symbol,"CEP55"))

data2.CBX2 <- data %>% filter(ID_REF %in% CBX2.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p27.data <- data2.CBX2 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  #mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =rep(c("HCC","NTL"),times=20)) #%>% group_by(Group2) %>% summarise(Mean=mean(Methylation_level))

library(ggstatsplot)

p27<-ggstatsplot::ggbetweenstats(
  data = p27.data %>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)

ggsave(p27,filename = "Human.GSE73003.CBX2.Methylation.pdf",height = 3.2,width = 2.8)

data2.CEP55 <- data %>% filter(ID_REF %in% CEP55.ID$ID) %>% remove_rownames() %>%
  column_to_rownames("ID_REF")

p28.data <- data2.CEP55 %>% data.frame() %>% rownames_to_column("Probeid") %>%
  gather(Sample,Methylation_level,-Probeid) %>%
  #set_colnames("Methylation_level") %>%
  #rownames_to_column("Sample") %>% #merge(metadata,by="Sample") %>%
  #mutate(Group2 =rep(c("HCC","NTL"),c(124,124))) %>%
  mutate(Methylation_level=as.numeric(as.character(Methylation_level))) %>%
  group_by(Sample) %>%
  mutate(Methylation_level=mean(Methylation_level)) %>%
  dplyr::select(Sample,Methylation_level) %>%
  unique() %>% data.frame() %>% mutate(Group2 =rep(c("HCC","NTL"),times=20))
#merge(metadata,by="Sample") %>%
#mutate(Group2 =if_else(str_detect(Group,"N1"),"Adjacent",if_else(str_detect(Group,"T1"),"PrimaryT","ReccurentT"))) %>%

p28<-ggstatsplot::ggbetweenstats(
  data = p28.data%>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = T)
ggsave(p28,filename = "Human.GSE73003.CEP55.Methylation.NO.pdf",height = 3.2,width = 2.8)

##################### º◊ª˘ªØ GSE75041∏¥∑¢ =°∑ No use ####################
##################### GSE82176 HELP => hg19 #######################
library(tidyverse)
# CBX2 TSS 77751978 chr17
# CEP55 TSS 95256389 chr10
data <- read.csv("GSE82176_Combined_processed_HELP_tagging_data.txt",sep="\t")
CBX2.data <- data %>% filter(chr=="chr17") %>%
  filter(pos >= 77751978-2000) %>%
  filter(pos <= 77751978+2000)
write.csv(CBX2.data,file = "GSE82176.HELP.CBX2.csv",row.names = F)

CEP55.data <- data %>% filter(chr=="chr10") %>%
  filter(pos >= 95256389-2000) %>%
  filter(pos <= 95256389+2000)
write.csv(CEP55.data,file = "GSE82176.HELP.CEP55.csv",row.names = F)

CBX2.data <- read.csv("GSE82176.HELP.CBX2.csv")
p29.data <- CBX2.data %>% dplyr::select(contains("_")) %>%
  gather(Sample,Methylation_level) %>% group_by(Sample) %>%
  summarise(Methylation_level=mean(Methylation_level)) %>%
  mutate(Group2=if_else(str_detect(Sample,"_C"),"Normal",if_else(str_detect(Sample,"NT"),"NTL","HCC")))

p29 <- ggstatsplot::ggbetweenstats(
  data = p29.data %>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CBX2")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)

ggsave(p29,filename = "Human.GSE82176.CBX2.Methylation.pdf",height = 3.2,width = 2.8)

CEP55.data <- read.csv("GSE82176.HELP.CEP55.csv")
p30.data <- CEP55.data %>% dplyr::select(contains("_")) %>%
  gather(Sample,Methylation_level) %>% group_by(Sample) %>%
  summarise(Methylation_level=mean(Methylation_level)) %>%
  mutate(Group2=if_else(str_detect(Sample,"_C"),"Normal",if_else(str_detect(Sample,"NT"),"NTL","HCC")))

p30 <- ggstatsplot::ggbetweenstats(
  data = p30.data %>% filter(Group2!="Normal"),
  x = Group2,
  y = Methylation_level,
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Probeid,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level CEP55")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)

ggsave(p30,filename = "Human.GSE82176.CEP55.Methylation.pdf",height = 3.2,width = 2.8)

##################### Human º◊ª˘ªØ GSE55752 ###################
setwd("K:/HepG2-CBX2/Validate/GSE55752")
files <- list.files(".",pattern = "txt$")
CBX2.files <- files[str_detect(files,"CBX2")]
CEP55.files <- files[str_detect(files,"CEP55")]
data <- data.frame()
for(file in files){
  gene <- ((str_split(file,"\\."))[[1]])[2]
  data <- read.csv(file,sep = " ",header = T) %>% mutate(gene=gene) %>%
    rbind(data)
}

data2 <- data %>% mutate(Group2=if_else(str_detect(Sample,"_P"),"NTL","HCC")) %>%
  mutate(Line=rep(1:4,each=4,length=32))
CBX2.data <- data2 %>% filter(gene=="CBX2") %>% mutate(Group2=if_else(str_detect(Sample,"_P"),"NTL","HCC"))
p31 <- ggstatsplot::ggbetweenstats(
  data = CBX2.data,
  x = Group2,
  y = S_Methylation_Level,
  pairwise.comparisons = F,
  p.adjust.method = "fdr",
  type = "robust",
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~Attri,scales = "free_y",ncol = 3)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)

ggsave(p31,filename = "Human.GSE55752.CBX2.Methylation.pdf",height = 3.1,width = 2.8)

p32 <- ggstatsplot::ggbetweenstats(
  data = data2 %>% filter(gene=="CEP55"),
  x = Group2,
  y = S_Methylation_Level,
  #pairwise.comparisons = TRUE,
  p.adjust.method = "fdr",
  type = "robust",
  messages = FALSE
)+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y",ncol = 2)+ #
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  #ggsci::scale_color_nejm()+
  labs(x="",y="Methylation level")+
  #scale_y_log10()+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),
                             method = "wilcox",paired = F)
ggsave(p32,filename = "Human.GSE55752.CEP55.Methylation.pdf",height = 3.1,width = 2.8)

##################### 

##################### Array => CBX2 + CEP55 GSE112790 #################
setwd("K:/HepG2-CBX2/Validate")
GPL <- read.csv("GPL570-55999.txt",sep = "\t") %>%
  filter(str_detect(Gene.Symbol,"CBX2") | str_detect(Gene.Symbol,"CEP55")) %>%
  dplyr::select(ID,Gene.Symbol)
data <- read.csv("GSE112790_series_matrix2.txt",sep = "\t") %>%
  merge(GPL,by.x="ID_REF",by.y="ID")

data2 <- data %>% gather(Sample,Expression,GSM3083512:GSM3083709) %>%
  group_by(Sample,Gene.Symbol) %>%
  summarise(Expression=mean(Expression)) %>% data.frame() %>%
  mutate(Group2=rep(c("NTL","HCC"),c(30,366)))

p35 <- ggplot(data2,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~Gene.Symbol,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p35,filename = "Human.GSE112790.CBX2.CEP55.Array.Expression.pdf",height = 3.1,width = 5)

##################### Array => CBX2 + CEP55 GSE25097 HCCDB3 ###################
setwd("K:/HepG2-CBX2/Validate/HCCDB3-GSE25097")
data.expression <- read.csv("GSE25097.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB3.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  filter(TYPE!="Cirrhotic") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC",
                        if_else(str_detect(TYPE,"Healthy"),"HL","NTL"))) %>%
  gather(gene,Expression,CBX2:CEP55)

p45 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  scale_y_log10()+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL"),
                                                c("HCC","HL")),method = "t.test")

ggsave(p45,filename = "../Human.GSE25097.CBX2.CEP55.Array.Expression.pdf",height = 3.1,width = 5.5)


##################### Array => CBX2 + CEP55 GSE22058 HCCDB1 #################
setwd("K:/HepG2-CBX2/Validate/HCCDB1-GSE22058")
data.expression <- read.csv("GSE22058-GPL6793.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB1.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

p44 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p44,filename = "../Human.GSE22058.CBX2.CEP55.Array.Expression.pdf",height = 3.1,width = 5)


##################### Array => CBX2 + CEP55 => Survival => GSE14520 ##########################
setwd("K:/HepG2-CBX2/Validate/GSE14520.Survival")
library(tidyverse)
library(readxl)
data.expression <- read.csv("GSE14520-GPL3921.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB6.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB6.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

p36 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p36,filename = "Human.GSE14520.CBX2.CEP55.Array.Expression.NO.pdf",height = 3.1,width = 5)

CBX2.CEP55.Survival <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>% merge(data.sample,by="Sample") %>%
  dplyr::select(Sample,CBX2,CEP55,PATIENT_ID) %>%
  merge(data.clinical,by.x="PATIENT_ID",by.y="Sample") %>%
  gather(Gene,Expression,CBX2:CEP55) %>%
  group_by(PATIENT_ID,Gene) %>%
  mutate(Expression=mean(Expression)) %>%
  ungroup() %>% dplyr::select(-Sample) %>%
  unique() %>% 
  mutate(Status=if_else(STATUS=="Dead",1,0))

CBX2.Survival <- CBX2.CEP55.Survival %>% 
  filter(Gene=="CBX2") %>%
  mutate(Exp=if_else(Expression>median(Expression),"High","Low")) %>%
  mutate(SURVIVAL_TIME=as.numeric(as.character(SURVIVAL_TIME))) %>%
  mutate(Status=as.numeric(as.character(Status)))
library(survival)
library(survminer)
# RµƒCatPredi∞¸°¢cutoff∞¸
# survminer ∞¸µƒsurv_cutpoint()∫Ø ˝
# https://zhuanlan.zhihu.com/p/424346991
res.cut <- surv_cutpoint(CBX2.Survival, # ˝æ›ºØ
                         time = "SURVIVAL_TIME", #…˙¥Ê ±º‰
                         event = "Status", #…˙¥Ê◊¥Ã¨
                         variables = "Expression" #–Ë“™º∆À„µƒ ˝æ›¡–√˚
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(SURVIVAL_TIME, Status)~Expression, data=res.cat)

p37 <- ggsurvplot(fit,pval =TRUE, data = res.cat, #risk.table = TRUE,
           surv.median.line = "hv", #fun = "event" cumulative events
           legend.title = "CBX2",
           legend.labs = c("High", "Low"),
           conf.int.style = "ribbon",# "step",
           xlab = "Time in days",
           risk.table = F, # "abs_pct",
           risk.table.y.text.col = F,
           risk.table.y.text = FALSE,
           conf.int = TRUE,
           palette = "Set1",
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           font.legend = c(12, "plain", "black"),
           ggtheme = ggthemes::theme_few())
p37
#ggsave(p37,filename = "Survival.GSE14520.CBX2.pdf",height = 3,width = 3)

CEP55.Survival <- CBX2.CEP55.Survival %>% 
  filter(Gene=="CEP55") %>%
  mutate(Exp=if_else(Expression>median(Expression),"High","Low")) %>%
  mutate(SURVIVAL_TIME=as.numeric(as.character(SURVIVAL_TIME))) %>%
  mutate(Status=as.numeric(as.character(Status)))

res.cut <- surv_cutpoint(CEP55.Survival, # ˝æ›ºØ
                         time = "SURVIVAL_TIME", #…˙¥Ê ±º‰
                         event = "Status", #…˙¥Ê◊¥Ã¨
                         variables = "Expression" #–Ë“™º∆À„µƒ ˝æ›¡–√˚
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(SURVIVAL_TIME, Status)~Expression, data=res.cat)

p38 <- ggsurvplot(fit,pval =TRUE, data = res.cat, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event" cumulative events
                  legend.title = "CEP55",
                  legend.labs = c("High", "Low"),
                  conf.int.style = "ribbon",# "step",
                  xlab = "Time in days",
                  risk.table = F, # "abs_pct",
                  risk.table.y.text.col = F,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "plain", "black"),
                  font.legend = c(12, "plain", "black"),
                  ggtheme = ggthemes::theme_few())
p38
#ggsave(p38,filename = "Survival.GSE14520.CEP55.pdf",height = 3,width = 3)

##################### Array => CBX2 + CEP55 => Survival => GSE10143 => No use ##########################
setwd("K:/HepG2-CBX2/Validate/GSE10143.Survival")
library(tidyverse)
library(readxl)
data.expression <- read.csv("GSE10143.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB7.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB7.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

##################### Array => CBX2 + CEP55 => Survival => GSE9843 => No use √ª”–…˙¥Ê ˝æ› ##########################
setwd("K:/HepG2-CBX2/Validate/GSE9843.Survival")
data.expression <- read.csv("GSE9843.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB8.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB8.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

##################### Array => CBX2 + CEP55 => Survival => GSE19977 => No use √ª”–…˙¥Ê ˝æ› ##########################
setwd("K:/HepG2-CBX2/Validate/GSE19977.Survival")
data.expression <- read.csv("GSE19977.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB9.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB9.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

##################### Array => CBX2 + CEP55 => Survival => GSE54236 => No use ±Ì¥Ô¡ø∑÷Œˆø…“‘+√ª”–…˙¥Ê ˝æ› ##########################
setwd("K:/HepG2-CBX2/Validate/GSE54236.Survival")
data.expression <- read.csv("GSE54236.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB12.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB12.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

p37 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p37,filename = "../Human.GSE54236.CBX2.CEP55.Array.Expression.pdf",height = 3.1,width = 5)

##################### Array => CBX2 + CEP55 => Survival => GSE63898 => No use ±Ì¥Ô¡ø∑÷Œˆø…“‘+√ª”–…˙¥Ê ˝æ› ##########################
setwd("K:/HepG2-CBX2/Validate/GSE63898.Survival")
data.expression <- read.csv("GSE63898.gene.txt",sep = "\t")
data.sample <- read.csv("HCCDB13.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB13.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

p38 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p38,filename = "../Human.GSE63898.CBX2.CEP55.Array.Expression.pdf",height = 3.1,width = 5)

##################### HGT => CBX2 + CEP55 => Survival => TCGA-LIHC => ##########################
##################### Array => CBX2 + CEP55 => Survival => GSE76427 => ±Ì¥Ô¡ø+…˙¥Ê ##########################
setwd("K:/HepG2-CBX2/Validate/GSE76427.Survival")
data.expression <- read.csv("HCCDB17_mRNA_level3.txt",sep = "\t")
data.sample <- read.csv("HCCDB17.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB17.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

p39 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p39,filename = "../Human.GSE76427.CBX2.CEP55.Array.Expression.pdf",height = 3.1,width = 5)
#######
CBX2.CEP55.Survival <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>% merge(data.sample,by="Sample") %>%
  dplyr::select(Sample,CBX2,CEP55,PATIENT_ID) %>%
  merge(data.clinical,by.x="PATIENT_ID",by.y="PATIENT_ID") %>%
  gather(Gene,Expression,CBX2:CEP55) %>%
  #group_by(PATIENT_ID,Gene) %>%
  #mutate(Expression=mean(Expression)) %>%
  #ungroup() %>% dplyr::select(-Sample) %>%
  #unique() %>% 
  mutate(Status=if_else(OS_STATUS=="Dead",1,0))

CBX2.Survival <- CBX2.CEP55.Survival %>% 
  filter(Gene=="CBX2") %>%
  #mutate(Exp=if_else(Expression>median(Expression),"High","Low")) %>%
  mutate(OS_YEAR=as.numeric(as.character(OS_YEAR))) %>%
  mutate(Status=as.numeric(as.character(Status)))
library(survival)
library(survminer)
# RµƒCatPredi∞¸°¢cutoff∞¸
# survminer ∞¸µƒsurv_cutpoint()∫Ø ˝
# https://zhuanlan.zhihu.com/p/424346991
res.cut <- surv_cutpoint(CBX2.Survival, # ˝æ›ºØ
                         time = "OS_YEAR", #…˙¥Ê ±º‰
                         event = "Status", #…˙¥Ê◊¥Ã¨
                         variables = "Expression" #–Ë“™º∆À„µƒ ˝æ›¡–√˚
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS_YEAR, Status)~Expression, data=res.cat)

p40 <- ggsurvplot(fit,pval =TRUE, data = res.cat, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event" cumulative events
                  legend.title = "CBX2",
                  legend.labs = c("High", "Low"),
                  conf.int.style = "ribbon",# "step",
                  xlab = "Time in years",
                  risk.table = F, # "abs_pct",
                  risk.table.y.text.col = F,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "plain", "black"),
                  font.legend = c(12, "plain", "black"),
                  ggtheme = ggthemes::theme_few())
p40
#ggsave(p40,filename = "Survival.GSE76427.CBX2.pdf",height = 3,width = 3)

CEP55.Survival <- CBX2.CEP55.Survival %>% 
  filter(Gene=="CEP55") %>%
  #mutate(Exp=if_else(Expression>median(Expression),"High","Low")) %>%
  mutate(OS_YEAR=as.numeric(as.character(OS_YEAR))) %>%
  mutate(Status=as.numeric(as.character(Status)))

res.cut <- surv_cutpoint(CEP55.Survival, # ˝æ›ºØ
                         time = "OS_YEAR", #…˙¥Ê ±º‰
                         event = "Status", #…˙¥Ê◊¥Ã¨
                         variables = "Expression" #–Ë“™º∆À„µƒ ˝æ›¡–√˚
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS_YEAR, Status)~Expression, data=res.cat)

p41 <- ggsurvplot(fit,pval =TRUE, data = res.cat, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event" cumulative events
                  legend.title = "CEP55",
                  legend.labs = c("High", "Low"),
                  conf.int.style = "ribbon",# "step",
                  xlab = "Time in years",
                  risk.table = F, # "abs_pct",
                  risk.table.y.text.col = F,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "plain", "black"),
                  font.legend = c(12, "plain", "black"),
                  ggtheme = ggthemes::theme_few())
p41
#ggsave(p41,filename = "Survival.GSE76427.CEP55.pdf",height = 3,width = 3)

##################### HGT => CBX2 + CEP55 => Survival => ICGC-LIRI-JP =>±Ì¥Ô¡ø∑÷Œˆø…“‘+…˙¥Ê∑÷Œˆ ##########################
setwd("K:/HepG2-CBX2/Validate/ICGC-LIRI-JP")
data.expression <- read.csv("HCCDB18_mRNA_level3.txt",sep = "\t")
data.sample <- read.csv("HCCDB18.sample.txt",sep = "\t",row.names = 1) %>%
  t() %>% data.frame() %>% rownames_to_column("Sample")
data.clinical <- read.csv("HCCDB18.patient.txt",sep = "\t",row.names = 1,check.names = F)%>%
  t() %>% data.frame() %>% rownames_to_column("Sample") %>%
  na.omit()

CBX2.CEP55.Exp <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>%
  merge(data.sample,by="Sample") %>%
  mutate(Group2=if_else(str_detect(TYPE,"HCC"),"HCC","NTL")) %>%
  gather(gene,Expression,CBX2:CEP55)

p41 <- ggplot(CBX2.CEP55.Exp,aes(Group2,Expression,fill=Group2))+
  geom_violin()+
  geom_boxplot(color="white")+
  ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(x="",y="Expression")+
  theme(strip.text = element_text(size=15,color="black"))+
  ggpubr::stat_compare_means(comparisons = list(c("HCC","NTL")),method = "t.test")

ggsave(p41,filename = "../Human.ICGC-LIRI-JP.CBX2.CEP55.HGT.Expression.pdf",height = 3.1,width = 5)

##### Survival #####
CBX2.CEP55.Survival <- data.expression %>% 
  filter(str_detect(Symbol,"CBX2") | str_detect(Symbol,"CEP55")) %>%
  remove_rownames() %>% column_to_rownames("Symbol") %>%
  dplyr::select(-Entrez_ID) %>% 
  t() %>% data.frame() %>% 
  rownames_to_column("Sample") %>% merge(data.sample,by="Sample") %>%
  dplyr::select(Sample,CBX2,CEP55,PATIENT_ID) %>%
  merge(data.clinical,by.x="PATIENT_ID",by.y="Sample") %>%
  gather(Gene,Expression,CBX2:CEP55) %>%
  group_by(PATIENT_ID,Gene) %>%
  mutate(Expression=mean(Expression)) %>%
  ungroup() %>% dplyr::select(-Sample) %>%
  unique() %>% 
  mutate(Status=if_else(STATUS=="Dead",1,0)) %>% data.frame()

CBX2.Survival <- CBX2.CEP55.Survival %>% 
  filter(Gene=="CBX2") %>%
  #mutate(Exp=if_else(Expression>median(Expression),"High","Low")) %>%
  mutate(SURVIVAL_TIME=as.numeric(as.character(SUR))) %>%
  mutate(Status=as.numeric(as.character(Status)))
library(survival)
library(survminer)
# RµƒCatPredi∞¸°¢cutoff∞¸
# survminer ∞¸µƒsurv_cutpoint()∫Ø ˝
# https://zhuanlan.zhihu.com/p/424346991
res.cut <- surv_cutpoint(CBX2.Survival, # ˝æ›ºØ
                         time = "SURVIVAL_TIME", #…˙¥Ê ±º‰
                         event = "Status", #…˙¥Ê◊¥Ã¨
                         variables = "Expression" #–Ë“™º∆À„µƒ ˝æ›¡–√˚
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(SURVIVAL_TIME, Status)~Expression, data=res.cat)

p42 <- ggsurvplot(fit,pval =TRUE, data = res.cat, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event" cumulative events
                  legend.title = "CBX2",
                  legend.labs = c("High", "Low"),
                  conf.int.style = "ribbon",# "step",
                  xlab = "Time in months",
                  risk.table = F, # "abs_pct",
                  risk.table.y.text.col = F,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "plain", "black"),
                  font.legend = c(12, "plain", "black"),
                  ggtheme = ggthemes::theme_few())
p42
#ggsave(p42,filename = "Survival.ICGC-LIRI-JP.CBX2.pdf",height = 3,width = 3)

CEP55.Survival <- CBX2.CEP55.Survival %>% 
  filter(Gene=="CEP55") %>%
  #mutate(Exp=if_else(Expression>median(Expression),"High","Low")) %>%
  mutate(SURVIVAL_TIME=as.numeric(as.character(SUR))) %>%
  mutate(Status=as.numeric(as.character(Status)))

res.cut <- surv_cutpoint(CEP55.Survival, # ˝æ›ºØ
                         time = "SURVIVAL_TIME", #…˙¥Ê ±º‰
                         event = "Status", #…˙¥Ê◊¥Ã¨
                         variables = "Expression" #–Ë“™º∆À„µƒ ˝æ›¡–√˚
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(SURVIVAL_TIME, Status)~Expression, data=res.cat)

p43 <- ggsurvplot(fit,pval =TRUE, data = res.cat, #risk.table = TRUE,
                  surv.median.line = "hv", #fun = "event" cumulative events
                  legend.title = "CEP55",
                  legend.labs = c("High", "Low"),
                  conf.int.style = "ribbon",# "step",
                  xlab = "Time in months",
                  risk.table = F, # "abs_pct",
                  risk.table.y.text.col = F,
                  risk.table.y.text = FALSE,
                  conf.int = TRUE,
                  palette = "Set1",
                  font.x = c(14, "bold", "black"),
                  font.y = c(14, "bold", "black"),
                  font.tickslab = c(12, "plain", "black"),
                  font.legend = c(12, "plain", "black"),
                  ggtheme = ggthemes::theme_few())
p43
#ggsave(p43,filename = "Survival.Survival.ICGC-LIRI-JP.CEP55.pdf",height = 3,width = 3)





#####################  CBX2 CEP55 EXPRESSION LEVEL ACROSS STAGE  #####################
load("Survival.Signature.RiskScore.Rdata") #risk_score_table_multi_cox2,multi_variate_cox_2,ph_hypo_table

meta <- read.csv("clinical.cart.2021-08-19/clinical.tsv",sep = "\t",header = T)
#meta <- column_to_rownames(meta,var = "case_submitter_id")
meta=meta[,colnames(meta) %in% c("case_submitter_id",
                                 "race",
                                 "gender",
                                 "ajcc_pathologic_t",
                                 "ajcc_pathologic_n",
                                 "ajcc_pathologic_m",
                                 "age","Stage")] %>% unique()

mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
colnames(mRNA.Exp.log2.TPM) = substr(colnames(mRNA.Exp.log2.TPM),1,12)
mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[c("ENSG00000173894","ENSG00000138180"),intersect(colnames(mRNA.Exp.log2.TPM),meta$case_submitter_id)] %>% t() %>% data.frame() %>%
  rownames_to_column("case_submitter_id") %>% merge(.,meta,by="case_submitter_id")
mRNA.Exp.log2.TPM.Target$Stage <- str_remove(mRNA.Exp.log2.TPM.Target$Stage,"Stage ")
mRNA.Exp.log2.TPM.Target$Stage <- factor(mRNA.Exp.log2.TPM.Target$Stage,levels = c("I","II","III","IV"))

p<-ggboxplot(mRNA.Exp.log2.TPM.Target %>% select(ENSG00000173894,Stage) %>% na.omit(), "Stage", "ENSG00000173894", 
             fill = "Stage",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("log2(TPM+1) CBX2",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T,comparisons = list(c("I","II"),c("I","III"),c("II","III")))
ggsave(p,filename = "TIMER2.CBX2.CEP55/Stage.CBX2.pdf",height = 2.5,width = 2.5)


p<-ggboxplot(mRNA.Exp.log2.TPM.Target %>% select(ENSG00000138180,Stage) %>% na.omit(), "Stage", "ENSG00000138180", 
             fill = "Stage",palette = "Set1",add = "jitter")+
  theme_few()+
  #facet_wrap(~Immune.Cell)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(colour = "black"))+
  labs(x="",y=paste("log2(TPM+1) CEP55",sep = ""))+
  stat_compare_means(method = "wilcox.test",label="p.format",hide.ns = T,comparisons = list(c("I","II"),c("I","III"),c("II","III")))
ggsave(p,filename = "TIMER2.CBX2.CEP55/Stage.CEP55.pdf",height = 2.5,width = 2.5)


#####################  CBX2 ChIP  ############################
setwd("K:/HepG2-CBX2")

cor <- data.frame()
for (i in 1:nrow(CBX2.CEP55.data)) {
  for (j in 1:nrow(Other.data)) {
    print(CBX2.CEP55.data$Symbol[i])
    print(Other.data$Symbol[j])
    test <- cor.test(as.numeric(as.character(CBX2.CEP55.data[i,2:213])),
                     as.numeric(as.character(Other.data[j,2:213])),method="pearson",exact = F,alternative = "two.sided")
    cor <- data.frame(Gene1=CBX2.CEP55.data$Symbol[i],
                      Gene2=Other.data$Symbol[j],
                      Rho=test$estimate,
                      Pvalue=test$p.value) %>%
      rbind.data.frame(cor)
  }
}










######### EZH2 CBX2 H3K27 Peak position barplot #########
library(tidyverse)
setwd("K:/HepG2-CBX2")
library(org.Hs.eg.db)
library(ChIPseeker)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak <- readPeakFile("HepG2.CBX2.bed")
CBX2.peakAnno <- annotatePeak(peak,tssRegion = c(-3000, 3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
#write.table(as.data.frame(peakAnno),"HepG2.CBX2.peak.annotation.tsv",sep="\t",row.names = F,quote = F)
#plotAnnoPie(peakAnno)
CBX2.peakAnno@annoStat %>% mutate(Type="CBX2") %>%
  write.csv("SourceData.F5A.csv")
peak <- readPeakFile("EZH2.IDR.bed")
EZH2.peakAnno <- annotatePeak(peak,tssRegion = c(-3000, 3000),TxDb = txdb,annoDb = "org.Hs.eg.db")

peak <- readPeakFile("K27me3.IDR.bed")
K27me3.peakAnno <- annotatePeak(peak,tssRegion = c(-3000, 3000),TxDb = txdb,annoDb = "org.Hs.eg.db")

middata.annoPeak <- CBX2.peakAnno@annoStat %>% mutate(Type="CBX2") %>% rbind.data.frame(EZH2.peakAnno@annoStat %>% mutate(Type="EZH2")) %>% 
  rbind.data.frame(K27me3.peakAnno@annoStat %>% mutate(Type="H3K27me3"))

middata.annoPeak$Type <- factor(middata.annoPeak$Type,levels = c("CBX2","EZH2","H3K27me3"))

library(ggthemes)
library(ggsci)
library(RColorBrewer)
p<-ggplot(middata.annoPeak,aes(x=Type,y=Frequency,fill=Feature))+geom_bar(stat = "identity",position = "stack")+
  theme_few()+labs(x="")+
  scale_fill_manual(values = c(brewer.pal(8,"Dark2"),brewer.pal(12,"Paired")[c(2,5,6,10)]))
  #scale_fill_npg()
ggsave(p,filename = "PeakAnnotation.pdf",height = 4,width = 5)

## Peak annotation
CBX2.peakAnno.Data <- as.data.frame(CBX2.peakAnno)
CBX2.peakAnno.Data.Promoter <- CBX2.peakAnno.Data %>% subset(str_detect(annotation,"Promoter"))
#GO
library(clusterProfiler)
library(org.Hs.eg.db)
ego <- enrichGO(gene          = CBX2.peakAnno.Data.Promoter$SYMBOL,
                keyType       = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

write.csv(as.data.frame(ego),file = "CBX2.Promoter.GO.csv")
setwd("K:/HepG2-CBX2")
CBX2.Peak.GO <- read.csv("CBX2.Promoter.GO.csv")
CBX2.Peak.GO2 <- CBX2.Peak.GO %>% group_by(ONTOLOGY) %>% top_n(-10,p.adjust) %>%
  filter(ONTOLOGY != "BP") %>%
  mutate(ID=factor(ID,levels = .$ID))
p1<-ggplot(CBX2.Peak.GO2,aes(x=ID,y=-log10(p.adjust),fill=ONTOLOGY))+
  geom_bar(stat="identity",position = position_dodge(0.9)) +
  coord_flip()+
  ggthemes::theme_few() +
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=Description,y=0.1),hjust = 0, #vjust = 1.5,
            position = position_dodge(0.5),fontface = 'bold')+
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))+
  #scale_fill_d3()+
  theme(legend.position = "top",axis.text.y = element_text(face = "bold",colour = "black"))+
  labs(x="")
ggsave(p1,filename = "CBX2.ChIP.GO2.pdf",height = 8,width = 6)

### CBX2 EZH2 H3K27me3 CHIP + RNA-SEQ OVERLAP
CBX2.Statified.Expression.Down <- read.csv("../TCGA/Anlysis/LIHC/CBX2.DEG.csv") %>% filter(adj.P.Val <= 0.01) %>% filter(abs(logFC) >= 1)

library(org.Hs.eg.db)
library(ChIPseeker)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

Overlap.peak <- readPeakFile("Overlap.Peak.bed")
Overlap.peakAnno <- annotatePeak(Overlap.peak,tssRegion = c(-3000, 3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
Overlap.peakAnno.Data <- as.data.frame(Overlap.peakAnno)
write.csv(Overlap.peakAnno.Data,"Overlap.peakAnno.Data.csv")

library(ggVennDiagram)
library(ggplot2)

Overlap.peakAnno.Data <- read.csv("Overlap.peakAnno.Data.csv")
dim(Overlap.peakAnno.Data)
inter.peaks <- intersect(CBX2.Statified.Expression.Down$ENSEMBL,Overlap.peakAnno.Data$ENSEMBL)
Overlap.peakAnno.Data <- Overlap.peakAnno.Data %>% filter(str_detect(annotation,"Promoter"))
inter.peaks2 <- Overlap.peakAnno.Data[Overlap.peakAnno.Data$ENSEMBL %in% inter.peaks,] %>% dim()
  data.frame() %>% filter(str_detect(annotation,"Promoter"))

CBX2.Statified.Expression.Down %>% filter(ENSEMBL %in% inter.peaks)
  

Overlap.peakAnno.Data$SYMBOL

library(ggVennDiagram)
library(RColorBrewer)
color <- brewer.pal(3, "Set3")
ggVennDiagram(list(CBX2.Statified.Expression.Down$ENSEMBL,Overlap.peakAnno.Data$ENSEMBL), 
              category.names = c("CBX2 stratified DEGs","Peaks annotated genes"))+
  theme(legend.position = "none")+
  scale_fill_distiller(palette = "RdBu")


#library(ggVennDiagram)
#x = list("RNA DEGs"=na.omit(CBX2.Statified.Expression.Down$ENSEMBL), "Peaks-related gene"=na.omit(Overlap.peakAnno.Data$ENSEMBL))
#ggve

CBX2.Statified.Expression.Down <- read.csv("../TCGA/Anlysis/LIHC/CBX2.DEG.csv") %>% filter(adj.P.Val <= 0.01) %>% filter(logFC <= -1)
Inter <- intersect(CBX2.Statified.Expression.Down$ENSEMBL,Overlap.peakAnno.Data$ENSEMBL)
Overlap.peakAnno.Data %>% filter(ENSEMBL %in% Inter)

CBX2.Statified.Expression.Up <- read.csv("../TCGA/Anlysis/LIHC/CBX2.DEG.csv") %>% filter(adj.P.Val <= 0.01) %>% filter(logFC >= 1)
Inter <- intersect(CBX2.Statified.Expression.Up$ENSEMBL,Overlap.peakAnno.Data$ENSEMBL)
Overlap.peakAnno.Data %>% filter(ENSEMBL %in% Inter)


LabelGene <- c("MYB","CDO1","LHX4","DMBX1","CNR1","STMN1","ZIC5","DLX4","HOXC6")
DEG <- read.csv("../TCGA/Anlysis/LIHC/CBX2.DEG.csv")

DEG <- Overlap.peakAnno.Data %>% subset(SYMBOL %in% LabelGene) %>% dplyr::select(ENSEMBL,SYMBOL) %>% merge.data.frame(DEG,by="ENSEMBL",all = T)

library(RColorBrewer)
library(ggrepel) #±Í«©”√

#»∑∂® «…œµ˜ªπ «œ¬µ˜£¨”√”⁄∏¯Õº÷–µ„…œ…´
DEG$threshold = factor(ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) >= 1, ifelse(DEG$logFC >= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
#df$gene <- row.names(df) #ÃÌº”“ª¡–ª˘“Ú√˚£¨“‘±„±∏◊¢
p<-ggplot(df,aes(x=logFC,y= -log10(adj.P.Val),color=threshold))+
  geom_point(data = DEG[! is.na(DEG$SYMBOL),],size = 3)+ 
  geom_point(data = DEG[is.na(DEG$SYMBOL),],size = 1,alpha=0.3)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+#»∑∂®µ„µƒ—’…´
  geom_text_repel(data=DEG[! is.na(DEG$SYMBOL),],aes(label = SYMBOL),
    size = 4,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#ÃÌº”πÿ◊¢µƒµ„µƒª˘“Ú√˚
  ylab('-log10 (P.adjust)')+#–ﬁ∏ƒy÷·√˚≥∆
  xlab('log2 (FoldChange)')+#–ﬁ∏ƒx÷·√˚≥∆
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#ÃÌº”∫·œﬂ|logFoldChange|>0.25
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5) +#ÃÌº” ˙œﬂpadj<0.05
  theme_few()+
  theme(axis.title.x = element_text(size = 15, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 15,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # …Ë÷√Õº¿˝±Í«©Œƒ◊÷
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # –ﬁ∏ƒX÷·…œ◊÷ÃÂ¥Û–°£¨
                                   color = "black", # —’…´
                                   face = "bold", #  face»°÷µ£∫plain∆’Õ®£¨boldº”¥÷£¨italic–±ÃÂ£¨bold.italic–±ÃÂº”¥÷
                                   vjust = 0.5, # Œª÷√
                                   hjust = 0.5, 
                                   angle = 0), #Ω«∂»
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )

ggsave(p,filename = 'Overlap.Gene.Representative.pdf',height = 5,width = 6)

### CBX2 CEP55
setwd("K:/HepG2-CBX2")
LabelGene <- c("CEP55")
DEG <- read.csv("../TCGA/Anlysis/LIHC/CBX2.DEG.csv")
DEG$Label = ifelse(DEG$Gene.name =="CEP55","CEP55",NA)

DEG$threshold = factor(ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) >= 1, ifelse(DEG$logFC >= 1 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
#df$gene <- row.names(df) #ÃÌº”“ª¡–ª˘“Ú√˚£¨“‘±„±∏◊¢
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(ggrepel)
p<-ggplot(df,aes(x=logFC,y= -log10(adj.P.Val),color=threshold))+
  geom_point(data = DEG[! is.na(DEG$Label),],size = 2)+ 
  geom_point(data = DEG[is.na(DEG$Label),],size = 1,alpha=0.3)+
  scale_color_manual(values=c("#4393C3","#00000033","#FC4E2A"))+#»∑∂®µ„µƒ—’…´
  geom_text_repel(data=DEG[! is.na(DEG$Label),],aes(label = Label),
                  size = 4,
                  color = "black",
                  segment.color = "black", show.legend = FALSE )+#ÃÌº”πÿ◊¢µƒµ„µƒª˘“Ú√˚
  ylab('-log10 (P.adjust)')+#–ﬁ∏ƒy÷·√˚≥∆
  xlab('log2 (FoldChange)')+#–ﬁ∏ƒx÷·√˚≥∆
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#ÃÌº”∫·œﬂ|logFoldChange|>0.25
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5) +#ÃÌº” ˙œﬂpadj<0.05
  theme_few()+
  theme(axis.title.x = element_text(size = 15, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 15,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # …Ë÷√Õº¿˝±Í«©Œƒ◊÷
                                   size = 10, 
                                   face = "bold"),
        axis.text.x = element_text(size = 13,  # –ﬁ∏ƒX÷·…œ◊÷ÃÂ¥Û–°£¨
                                   color = "black", # —’…´
                                   face = "bold", #  face»°÷µ£∫plain∆’Õ®£¨boldº”¥÷£¨italic–±ÃÂ£¨bold.italic–±ÃÂº”¥÷
                                   vjust = 0.5, # Œª÷√
                                   hjust = 0.5, 
                                   angle = 0), #Ω«∂»
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0) 
  )

ggsave(p,filename = 'Overlap.Gene.CEP55.Representative.pdf',height = 5,width = 6)



## ≤Èø¥UP DOWN
mRNA.Exp <- read.csv("../TCGA/mRNA.PanCancer.Exp/PanCancer.TCGA-LIHC.mRNA.Exp.csv",check.names = F,row.names = 1)
mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)

mRNA.Exp[c("ENSG00000129596","ENSG00000173894"),] %>% t() %>% data.frame() %>% mutate(CBX2 = ifelse(ENSG00000173894 >= median(ENSG00000173894),"High","Low")) %>%
  group_by(CBX2) %>% summarize(mean=mean(ENSG00000129596))


library(clusterProfiler)
library(org.Hs.eg.db)
ego <- enrichGO(gene          = Overlap.peakAnno.Data$SYMBOL,
                keyType       = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

write.csv(as.data.frame(ego),file = "MergePeak.Overlap.GO.csv")
MergePeakGO <- read.csv("MergePeak.Overlap.GO.csv")
MergePeakGO2 <- MergePeakGO %>% group_by(ONTOLOGY) %>% top_n(-10,p.adjust) %>%
  mutate(ID=factor(ID,levels = .$ID))
p1<-ggplot(MergePeakGO2,aes(x=ID,y=-log10(p.adjust),fill=ONTOLOGY))+
  geom_bar(stat="identity",position = position_dodge(0.9)) +
  coord_flip()+
  theme_few() +
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label=Description,y=0.1),hjust = 0, #vjust = 1.5,
            position = position_dodge(0.5),fontface = 'bold')+
  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5"))+
  #scale_fill_d3()+
  theme(legend.position = "top",axis.text.y = element_text(face = "bold",colour = "black"))+
  labs(x="")

ggsave(p1,filename = "Overlap.Peak.GO.pdf",height = 10,width = 6)




##################################### 
if (F) {
  mRNA.Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",check.names = F,row.names = 1)
  mRNA.Exp.log2.TPM <- log2(mRNA.Exp+1)
  load("Inter.miRNA.mRNA.Samples.361.Rdata")
  mRNA.Exp.log2.TPM.Target <- mRNA.Exp.log2.TPM[,Inter.Sample.361] %>%
    data.frame(check.names = F)
  ### CBX2/CEP55 => EXPRESSION + GROUP
  #saveRDS(mRNA.Exp.log2.TPM.Target,file = "Figure6.CBX2.CEP55.Expression.Group.Rds")
  Figure.D1 <- readRDS("Figure6.CBX2.CEP55.Expression.Group.Rds")
  ### Age Gender
  #save(risk.score.meta,file = "risk.score.meta.Rdata")
  load(file="risk.score.meta.Rdata")
  risk.score.meta$Stage <- risk.score.meta$Stage %>% str_remove_all(".* ")
  risk.score.meta$Stage <- factor(risk.score.meta$Stage,levels = c("I","II","III","IV"))
  risk.score.meta$gender <- factor(risk.score.meta$gender,levels = c("male","female"))
  colnames(risk.score.meta)[10:12] = c("M","N","T")
  risk.score.meta$T <- factor(risk.score.meta$T,levels =paste("T",c(1,2,3,4,"X"),sep = ""))
  risk.score.meta$M <- factor(risk.score.meta$M,levels = paste("M", c("0","1","X"),sep = ""))
  risk.score.meta$N <- factor(risk.score.meta$N,levels = paste("N", c("0","1","X"),sep = ""))
  saveRDS(risk.score.meta,"Figure6.Age.Stage.Rds")
  Figure.D2 <- readRDS("Figure6.Age.Stage.Rds")
  ### CTNNB1 TP53
  library(data.table)
  mutation.D <- fread("../../Clinical.XENA/LIHC/TCGA-LIHC.mutect2_snv.tsv") %>%
    filter(gene %in% c("CTNNB1","TP53")) %>% 
    dplyr::select(Sample_ID,gene) %>%
    mutate(Mutation="Yes") %>%
    unique() %>%
    spread(gene,Mutation,fill="No")
  saveRDS(mutation.D,"Figure6.CTNNB1.TP53.Rds")
  Figure.D3 <- readRDS("Figure6.CTNNB1.TP53.Rds")
  ### Immune Gene + AFP
  library(readxl)
  Immune <- read_xlsx("../../Immune.Representative.xlsx")
  Immune.Gene.Type <- mid.Gene.Type %>% filter(HGNC.symbol %in% Immune$SYMBOL) %>%
    filter(Gene.stable.ID %in% colnames(mRNA.Exp.log2.TPM.Target)) %>%
    mutate(HGNC.symbol=factor(HGNC.symbol,levels=Immune$SYMBOL)) %>%
    arrange(HGNC.symbol)
  Immune.E <- mRNA.Exp.log2.TPM.Target[c(Immune.Gene.Type$Gene.stable.ID),] %>%
    t() %>%
    data.frame(check.names = F) %>%
    magrittr::set_colnames(Immune.Gene.Type$HGNC.symbol) %>%
    rownames_to_column("Sample.ID")
  saveRDS(Immune.E,file = "Figure6.Immune.Exp.Rds")
  Figure.D4 <- readRDS("Figure6.Immune.Exp.Rds")
  
  AFP.E <- mRNA.Exp.log2.TPM.Target["ENSG00000081051",] %>%
    t() %>%
    data.frame(check.names = F) %>%
    magrittr::set_colnames("AFP") %>%
    rownames_to_column("Sample.ID")
  saveRDS(AFP.E,file = "Figure6.AFP.Exp.Rds")
  Figure.D5 <- readRDS("Figure6.AFP.Exp.Rds")
  ### TMB IPS ...
  ESTIMATE.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 2) %>%
    filter(str_detect(Type,"LIHC")) %>% dplyr::select(-Type)
  IPS.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 3) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  TMB.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 4) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  MSI.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 6) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  NEO.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 7) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  Purity.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 8) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  Ploidy.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 9) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  HRD.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 10) %>%
    filter(str_detect(CODE,"LIHC")) %>% dplyr::select(-CODE)
  TIDE.D <- read.csv("../../TIDE/TCGA-LIHC.TIDE.csv") %>% 
    mutate(Sample.ID = substr(Patient,1,15)) %>%
    dplyr::select(Dysfunction,Exclusion,TIDE,Sample.ID)
  Grade.D <- read_xlsx("../../PanCancer.DiverseIndex.xlsx",sheet = 18) %>%
    filter(str_detect(...1,"LIHC")) %>% dplyr::select(-...1)
  Index.D <- merge(ESTIMATE.D,IPS.D,by="ID",all = T) %>%
    merge(TMB.D,by.x="ID",by.y="SampleName",all=T) %>%
    merge(MSI.D,by.x="ID",by.y="SampleName",all=T) %>%
    merge(NEO.D,by.x="ID",by.y="SampleName",all=T) %>%
    merge(Purity.D,by.x="ID",by.y="SampleName",all=T) %>%
    merge(Ploidy.D,by.x="ID",by.y="SampleName",all=T) %>%
    merge(HRD.D,by.x="ID",by.y="SampleName",all=T) %>%
    merge(TIDE.D,by.x="ID",by.y="Sample.ID",all=T) %>%
    merge(Grade.D,by.x="ID",by.y="SampleName",all=T) %>%
    filter(!str_detect(ID,"-11"))
  saveRDS(Index.D,file = "Figure6.TMB.MSI...Rds")
  Figure.D6 <- readRDS(file = "Figure6.TMB.MSI...Rds")
  ### Immune Index
  MidData <- read_xlsx("../../PanCancer.Immunity.Index.xlsx") %>%
    filter(`TCGA Study`=="LIHC")
  MidData2 <- MidData[,c(1,5:7,9:14,26:32,37:58)]
  saveRDS(MidData2,file = "Figure6.ImmuneCell...Rds")
  Figure.D7 <- readRDS(file = "Figure6.ImmuneCell...Rds")
  ### Signature => CancerSEA
  MIDDATA <- read.csv("../HCC.Recurrence/Pathway.CancerSEA.LIHC.csv",check.names = F,
                      row.names = 1)
  saveRDS(MIDDATA,file = "Figure6.CancerSEA...Rds")
  Figure.D8 <- readRDS(file = "Figure6.CancerSEA...Rds")
  ### Signature 
  data <- fread("../../TCGA-LIHC.SYMBOLS.TPM.csv") %>%
    remove_rownames() %>% column_to_rownames("V1") %>%
    remove_rownames() %>% column_to_rownames("hgnc_symbol")
  data <- log2(data+1)
  
  Signature <- read_xlsx("../../ICB.Related.Signature.xlsx") %>%
    gather(Signature,SYMBOL) %>%
    na.omit()
  Signature.S <- split(Signature$SYMBOL,Signature$Signature)
  library(GSVA)
  LIHC.GSEA <- gsva(expr=as.matrix(data),
                    gset.idx.list=Signature.S, 
                    method="ssgsea",
                    kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                    verbose=T)
  write.csv(LIHC.GSEA,file = "../ICB-Related.LIHC.csv")
  saveRDS(t(LIHC.GSEA),file = "Figure6.ICB.Signature...Rds")
  Figure.D9 <- readRDS(file = "Figure6.ICB.Signature...Rds")
  ####
  Tms <- read.csv("../../Tms.DeMixt/LIHC_TmS_summary.csv") %>%
    dplyr::select(Sample.ID,TmS)
  colnames(Figure.D1)[1:2]=c("CBX2","CEP55")
  Figure6.Data <- Figure.D1 %>% mutate(sampleid=substr(SampleID,1,15)) %>%
    merge(Figure.D2,by="sampleid",all = T) %>%
    dplyr::select(-OS.time,-OS,-total_risk_score,-RiskScore,-Score) %>%
    merge(Figure.D3,by.x = "SampleID",by.y = "Sample_ID",all=T) %>%
    merge(Figure.D4,by.x="SampleID",by.y="Sample.ID",all=T) %>%
    merge(Figure.D5,by.x="SampleID",by.y="Sample.ID",all=T) %>%
    merge(Figure.D6,by.x="sampleid",by.y="ID",all=T) %>%
    merge(Figure.D7,by.x="case_submitter_id",by.y="TCGA Participant Barcode",all=T) %>%
    merge(Figure.D8 %>% t() %>% data.frame(check.names = F) %>%
            rownames_to_column("SampleID"),by="SampleID",all=T) %>%
    merge(Figure.D9 %>% data.frame(check.names = F) %>%
            rownames_to_column("SampleID"),by="SampleID",all=T) %>%
    merge(Tms,by.x="case_submitter_id",by.y="Sample.ID",all=T)
  saveRDS(Figure6.Data,file = "Figure6.Data.Test.Rds")
  write.csv(Figure6.Data,file = "Figure6.Data.Test.csv")
}

#### Heatmap => CBX2 ####
Figure6.Data <- read.csv("Figure6.Data.csv",check.names = F)
Figure6.Data2 <- Figure6.Data %>% filter(!is.na(CBX2))
Figure6.Data2$age <- if_else(Figure6.Data2$age>60,">60","<=60")
Figure6.Data2$Stage <- as.character(Figure6.Data2$Stage)
Figure6.Data2$Stage <- if_else(is.na(Figure6.Data2$Stage),"Unknown",Figure6.Data2$Stage)
Figure6.Data2$gender <- as.character(Figure6.Data2$gender)
Figure6.Data2$gender <- if_else(Figure6.Data2$gender=="female","Female","Male")
Figure6.Data2$T <- as.character(Figure6.Data2$T)
Figure6.Data2$T <- if_else(is.na(Figure6.Data2$T)|Figure6.Data2$T=="TX","Unknown",Figure6.Data2$T)
Figure6.Data2$M <- as.character(Figure6.Data2$M)
Figure6.Data2$M <- if_else(is.na(Figure6.Data2$M)|Figure6.Data2$M=="MX","Unknown",Figure6.Data2$M)
Figure6.Data2$N <- as.character(Figure6.Data2$N)
Figure6.Data2$N <- if_else(is.na(Figure6.Data2$N) | Figure6.Data2$N=="NX","Unknown",Figure6.Data2$N)
Figure6.Data2$CTNNB1 = if_else(is.na(Figure6.Data2$CTNNB1),"No",Figure6.Data2$CTNNB1)
Figure6.Data2$TP53 = if_else(is.na(Figure6.Data2$TP53),"No",Figure6.Data2$TP53)
Figure6.Data2$Grade <- as.character(Figure6.Data2$Grade)
Figure6.Data2$Grade <- if_else(is.na(Figure6.Data2$Grade),"Unknown",Figure6.Data2$Grade)

fisher.test(table(Figure6.Data2$CTNNB1,Figure6.Data2$CBX2Group))
# p-value = 0.000444 odds ratio 0.4263922 
fisher.test(table(Figure6.Data2$TP53,Figure6.Data2$CBX2Group))

Figure6.HD <- Figure6.Data2 %>% dplyr::select(-Grade,-age,-gender,-T,-N,-M,-TP53,
                                              -CTNNB1,-Stage,-CEP55Group,-CBX2Group) %>%
  dplyr::select(-case_submitter_id,-SampleID) %>%
  #dplyr::select(-Proliferation) %>% 
  remove_rownames() %>%
  column_to_rownames("sampleid")
for (name in colnames(Figure6.HD)) {
  Figure6.HD[[name]] = as.numeric(as.character(Figure6.HD[[name]]))
}

#Figure6.HD.S %>% filter(p <= 0.05) %>% View()

Figure6.Scale <- scale(Figure6.HD)

Figure6S.S <- Figure6.Data2 %>% dplyr::select(Grade,age,gender,T,N,M,TP53,
                                              CTNNB1,Stage,CEP55Group,CBX2Group) %>%
  cbind(Figure6.Scale)
write.csv(Figure6S.S,file = "SourceData.F6D.csv")

Figure6.HD.T <- Figure6.HD %>% mutate(group=Figure6.Data2$CBX2Group) %>%
  gather(Signature,Score,-group) %>%
  group_by(Signature) %>%
  rstatix::t_test(Score~group,detailed = T)%>%
  mutate(Signature = factor(Signature,levels = colnames(Figure6.Scale))) %>%
  arrange(Signature) %>%
  rstatix::add_significance("p") %>%
  mutate(p.signif = if_else(p.signif=="ns","",p.signif))

Figure6.HD.W <- Figure6.HD %>% mutate(group=Figure6.Data2$CBX2Group) %>%
  gather(Signature,Score,-group) %>%
  group_by(Signature) %>%
  rstatix::wilcox_test(Score~group,detailed = T) %>%
  mutate(Signature = factor(Signature,levels = colnames(Figure6.Scale))) %>%
  arrange(Signature) %>%
  rstatix::add_significance("p")

library(ComplexHeatmap)
library(circlize)
nejm <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
npg <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF")
Subtype = Figure6.Data2$CBX2Group
Age1 = Figure6.Data2$age
Gender1 = Figure6.Data2$gender
Stage1 = Figure6.Data2$Stage
Grade1 = Figure6.Data2$Grade
CEP55G = Figure6.Data2$CEP55Group
T1 = Figure6.Data2$T
M1 = Figure6.Data2$M
N1 = Figure6.Data2$N
TP53 = Figure6.Data2$TP53
CTNNB1 = Figure6.Data2$CTNNB1

harow = rowAnnotation(foo = anno_text(Figure6.HD.T$p.signif, just='left'))

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = nejm[2:1]),
                                        labels = c("CBX2 low","CBX2 high"), 
                                        labels_gp = gpar(col = "white", fontsize = 12)),
                       Age= factor(Age1,levels = c(">60","<=60")),
                       Gender=factor(Gender1,levels = c("Female","Male")),
                       Stage=factor(Stage1,levels = c("I","II","III","IV","Unknown")),
                       #T_stage=factor(T1,levels = c("T1","T2","T3","T4","Unknown")),
                       #M_stage=factor(M1,levels = c("M0","M1","Unknown")),
                       #N_stage=factor(N1,levels = c("N0","N1","Unknown")),
                       Grade=factor(Grade1,levels = c("G1","G2","G3","G4","Unknown")),
                       CEP55= factor(CEP55G,levels = c("CEP55high","CEP55low","Unknown")),
                       TP53 = factor(TP53,levels = c("Yes","No","Unknown")),
                       CTNNB1 = factor(CTNNB1,levels = c("Yes","No","Unknown")),
                       col = list(Age = c(">60" = "#BC3C2999", "<=60" = "#FFDC9199"),
                                  Gender = c("Female" = nejm[3],"Male"="#ABDDA4"),
                                  #ISUP = c("1" = "#FFFFBF" , "2-3" = "#ABDDA4", "4" = "#9970AB",'5'="#D6604D",'Unknown'='#FFFFFF'),
                                  #TMPRSS2ERG = c("Positive" = "#BC3C2999", "Negative" = "#FFDC9199",'Unknown'='#FFFFFF'),
                                  Stage = c('I'="#FFFFBF" ,'II' ="#ABDDA4","III" = "#9970AB",'IV'="#D6604D" , 'Unknown'='#FFFFFF'),
                                  #T_stage =  c('T1'="#FFFFBF" ,'T2' ="#ABDDA4","T3" = "#9970AB",'T4'="#D6604D" , 'Unknown'='#FFFFFF'),
                                  #M_stage =  c('M0' ="#ABDDA4","M1" = "#9970AB",'Unknown'='#FFFFFF'),
                                  #N_stage =  c('N0' ="#ABDDA4","N1" = "#9970AB",'Unknown'='#FFFFFF'),
                                  Grade = c('G1'="#FFFFBF" ,'G2' ="#ABDDA4","G3" = "#9970AB",'G4'="#D6604D",'Unknown'='#FFFFFF'),
                                  CEP55 = c("CEP55high" = "#BC3C29FF","CEP55low"="#ABDDA4"),
                                  TP53 = c("Yes" = "#BC3C29FF", "No"="white"),
                                  CTNNB1 = c("Yes" = "#BC3C29FF", "No"="white")),#∏˜∏ˆ¡Ÿ¥≤–≈œ¢µƒ—’…´£¨Anatomic_location—’…´ ‰»Î“≤Õ¨∆‰À˚£¨∏¸∫√—’…´¥Ó≈‰–≈œ¢«Î≤ŒøºFigureYa28color
                       show_legend = rep(TRUE, 7),# «∑Ò“™œ‘ æannotation legend
                       annotation_height = unit(rep(5,8), "mm"),#¡Ÿ¥≤annotationµƒ∏ﬂ∂»
                       annotation_legend_param = list(
                         Age = list(title = "Age"),
                         Gender = list(title = "Gender"),
                         Stage = list(title = "Stage"),
                         #T_stage = list(title = "T_stage"),
                         #M_stage = list(title = "M_stage"),
                         #N_stage = list(title = "N_stage"),
                         Grade = list(title = "Grade"),
                         #CEP55= CEP55G,
                         CEP55 = list(title = "CEP55 expression",
                                      at=c('CEP55high','CEP55low'),
                                      labels=c('High','Low')),
                         TP53 = list(title = "TP53"),
                         CTNNB1 = list(title = "CTNNB1")),
                       simple_anno_size = unit(0.3, "cm"),
                       na_col = "white",
                       annotation_name_side = "left"
                       )

ht_opt(legend_border='black',
       heatmap_border = TRUE,
       annotation_border = TRUE) 
col_fun <- colorRamp2(c(-2, 0,2), c("#377EB8", "white", "#E41A1C"))
subtype <- Subtype
ht <- Heatmap(t(Figure6.Scale), col = col_fun, 
              name = "Score",
              cluster_rows = F, 
              cluster_columns = T,         
              show_row_names = TRUE, 
              row_names_side = "left", #––√˚Œª÷√
              show_column_names = FALSE,
              column_split = subtype, #∏˘æ› ˝æ›–ﬁ∏ƒ
              row_split = c(rep(1,3),rep(2,6),rep(3,15),rep(4,1),rep(5,6),rep(6,3),rep(7,6)),#∏˘æ› ˝æ›–ﬁ∏ƒ
              row_gap = unit(c(1,1,1,1,1,1) ,"mm"),
              column_gap = unit(c(2), "mm"),
              row_title = NULL,column_title = NULL,
              right_annotation = harow,
              top_annotation = ha)
pdf("Figure6.F6.pdf",width = 16,height = 10)
draw(ht,padding = unit(c(0.1, 4, 0.1, 0.5), "cm"),
     annotation_legend_side = "right", 
     heatmap_legend_side = "right")
dev.off()

#### Heatmap => CEP55 ####
Figure6.Data <- read.csv("Figure6.Data.csv",check.names = F)
Figure6.Data2 <- Figure6.Data %>% filter(!is.na(CBX2)) %>%
  arrange(CEP55Group)
Figure6.Data2$age <- if_else(Figure6.Data2$age>60,">60","<=60")
Figure6.Data2$Stage <- as.character(Figure6.Data2$Stage)
Figure6.Data2$Stage <- if_else(is.na(Figure6.Data2$Stage),"Unknown",Figure6.Data2$Stage)
Figure6.Data2$gender <- as.character(Figure6.Data2$gender)
Figure6.Data2$gender <- if_else(Figure6.Data2$gender=="female","Female","Male")
Figure6.Data2$T <- as.character(Figure6.Data2$T)
Figure6.Data2$T <- if_else(is.na(Figure6.Data2$T)|Figure6.Data2$T=="TX","Unknown",Figure6.Data2$T)
Figure6.Data2$M <- as.character(Figure6.Data2$M)
Figure6.Data2$M <- if_else(is.na(Figure6.Data2$M)|Figure6.Data2$M=="MX","Unknown",Figure6.Data2$M)
Figure6.Data2$N <- as.character(Figure6.Data2$N)
Figure6.Data2$N <- if_else(is.na(Figure6.Data2$N) | Figure6.Data2$N=="NX","Unknown",Figure6.Data2$N)
#Figure6.Data2$CTNNB1 = if_else(is.na(Figure6.Data2$CTNNB1),"No",Figure6.Data2$CTNNB1)
#Figure6.Data2$TP53 = if_else(is.na(Figure6.Data2$TP53),"No",Figure6.Data2$TP53)
Figure6.Data2$Grade <- as.character(Figure6.Data2$Grade)
Figure6.Data2$Grade <- if_else(is.na(Figure6.Data2$Grade),"Unknown",Figure6.Data2$Grade)

Figure6.HD <- Figure6.Data2 %>% dplyr::select(-Grade,-age,-gender,-T,-N,-M,-TP53,
                                              -CTNNB1,-Stage,-CEP55Group,-CBX2Group) %>%
  dplyr::select(-case_submitter_id,-SampleID) %>%
  #dplyr::select(-Proliferation) %>% 
  remove_rownames() %>%
  column_to_rownames("sampleid")
for (name in colnames(Figure6.HD)) {
  Figure6.HD[[name]] = as.numeric(as.character(Figure6.HD[[name]]))
}

#Figure6.HD.S %>% filter(p <= 0.05) %>% View()

Figure6.Scale <- scale(Figure6.HD)

Figure6.HD.T <- Figure6.HD %>% mutate(group=Figure6.Data2$CEP55Group) %>%
  gather(Signature,Score,-group) %>%
  group_by(Signature) %>%
  rstatix::t_test(Score~group,detailed = T)%>%
  mutate(Signature = factor(Signature,levels = colnames(Figure6.Scale))) %>%
  arrange(Signature) %>%
  rstatix::add_significance("p") %>%
  mutate(p.signif = if_else(p.signif=="ns","",p.signif))
Figure6.HD.W <- Figure6.HD %>% mutate(group=Figure6.Data2$CEP55Group) %>%
  gather(Signature,Score,-group) %>%
  group_by(Signature) %>%
  rstatix::wilcox_test(Score~group,detailed = T) %>%
  mutate(Signature = factor(Signature,levels = colnames(Figure6.Scale))) %>%
  arrange(Signature) %>%
  rstatix::add_significance("p")

library(ComplexHeatmap)
library(circlize)
nejm <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
npg <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF")
Subtype = Figure6.Data2$CEP55Group
Age1 = Figure6.Data2$age
Gender1 = Figure6.Data2$gender
Stage1 = Figure6.Data2$Stage
Grade1 = Figure6.Data2$Grade
#CEP55G = Figure6.Data2$CEP55Group
#T1 = Figure6.Data2$T
#M1 = Figure6.Data2$M
#N1 = Figure6.Data2$N
#TP53 = Figure6.Data2$TP53
#CTNNB1 = Figure6.Data2$CTNNB1

harow = rowAnnotation(foo = anno_text(Figure6.HD.T$p.signif, just='left'))

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = nejm[2:1]),
                                        labels = c("CEP55 low","CEP55 high"), 
                                        labels_gp = gpar(col = "white", fontsize = 12)),
                       Age= factor(Age1,levels = c(">60","<=60")),
                       Gender=factor(Gender1,levels = c("Female","Male")),
                       Stage=factor(Stage1,levels = c("I","II","III","IV","Unknown")),
                       #T_stage=factor(T1,levels = c("T1","T2","T3","T4","Unknown")),
                       #M_stage=factor(M1,levels = c("M0","M1","Unknown")),
                       #N_stage=factor(N1,levels = c("N0","N1","Unknown")),
                       Grade=factor(Grade1,levels = c("G1","G2","G3","G4","Unknown")),
                       #CEP55= factor(CEP55G,levels = c("CEP55high","CEP55low","Unknown")),
                       #TP53 = factor(TP53,levels = c("Yes","No","Unknown")),
                       #CTNNB1 = factor(CTNNB1,levels = c("Yes","No","Unknown")),
                       col = list(Age = c(">60" = "#BC3C2999", "<=60" = "#FFDC9199"),
                                  Gender = c("Female" = nejm[3],"Male"="#ABDDA4"),
                                  #ISUP = c("1" = "#FFFFBF" , "2-3" = "#ABDDA4", "4" = "#9970AB",'5'="#D6604D",'Unknown'='#FFFFFF'),
                                  #TMPRSS2ERG = c("Positive" = "#BC3C2999", "Negative" = "#FFDC9199",'Unknown'='#FFFFFF'),
                                  Stage = c('I'="#FFFFBF" ,'II' ="#ABDDA4","III" = "#9970AB",'IV'="#D6604D" , 'Unknown'='#FFFFFF'),
                                  #T_stage =  c('T1'="#FFFFBF" ,'T2' ="#ABDDA4","T3" = "#9970AB",'T4'="#D6604D" , 'Unknown'='#FFFFFF'),
                                  #M_stage =  c('M0' ="#ABDDA4","M1" = "#9970AB",'Unknown'='#FFFFFF'),
                                  #N_stage =  c('N0' ="#ABDDA4","N1" = "#9970AB",'Unknown'='#FFFFFF'),
                                  Grade = c('G1'="#FFFFBF" ,'G2' ="#ABDDA4","G3" = "#9970AB",'G4'="#D6604D",'Unknown'='#FFFFFF')
                                  #CEP55 = c("CEP55high" = "#BC3C29FF","CEP55low"="#ABDDA4"),
                                  #TP53 = c("Yes" = "#BC3C29FF", "No"="white"),
                                  #CTNNB1 = c("Yes" = "#BC3C29FF", "No"="white")
                                  ),#∏˜∏ˆ¡Ÿ¥≤–≈œ¢µƒ—’…´£¨Anatomic_location—’…´ ‰»Î“≤Õ¨∆‰À˚£¨∏¸∫√—’…´¥Ó≈‰–≈œ¢«Î≤ŒøºFigureYa28color
                       show_legend = rep(TRUE, 4),# «∑Ò“™œ‘ æannotation legend
                       annotation_height = unit(rep(5,5), "mm"),#¡Ÿ¥≤annotationµƒ∏ﬂ∂»
                       annotation_legend_param = list(
                         Age = list(title = "Age"),
                         Gender = list(title = "Gender"),
                         Stage = list(title = "Stage"),
                         #T_stage = list(title = "T_stage"),
                         #M_stage = list(title = "M_stage"),
                         #N_stage = list(title = "N_stage"),
                         Grade = list(title = "Grade")
                         #CEP55= CEP55G,
                         #CEP55 = list(title = "CEP55 expression",at=c('CEP55high','CEP55low'),labels=c('High','Low')),
                         #TP53 = list(title = "TP53"),
                         #CTNNB1 = list(title = "CTNNB1")
                         ),
                       simple_anno_size = unit(0.3, "cm"),
                       na_col = "white",
                       annotation_name_side = "left"
)

ht_opt(legend_border='black',
       heatmap_border = TRUE,
       annotation_border = TRUE) 
col_fun <- colorRamp2(c(-2, 0,2), c("#377EB8", "white", "#E41A1C"))
subtype <- Subtype
ht <- Heatmap(t(Figure6.Scale), col = col_fun, 
              name = "Score",
              cluster_rows = F, 
              cluster_columns = T,         
              show_row_names = TRUE, 
              row_names_side = "left", #––√˚Œª÷√
              show_column_names = FALSE,
              column_split = subtype, #∏˘æ› ˝æ›–ﬁ∏ƒ
              row_split = c(rep(1,3),rep(2,6),rep(3,15),rep(4,1),rep(5,6),rep(6,3),rep(7,6)),#∏˘æ› ˝æ›–ﬁ∏ƒ
              row_gap = unit(c(1,1,1,1,1,1) ,"mm"),
              column_gap = unit(c(2), "mm"),
              row_title = NULL,column_title = NULL,
              right_annotation = harow,
              top_annotation = ha)
pdf("Figure6.F6.CEP55.pdf",width = 14,height = 10)
draw(ht,padding = unit(c(0.1, 4, 0.1, 0.5), "cm"),
     annotation_legend_side = "right", 
     heatmap_legend_side = "right")
dev.off()
#### Heatmap => Immune genes #####
Figure6.S <- read_xlsx("Figure6.S.xlsx")

Figure6.SD <- Figure6.S %>% dplyr::select(-Grade,-age,-gender,-T,-N,-M,-TP53,
                                              -CTNNB1,-Stage,-CEP55Group,-CBX2Group) %>%
  dplyr::select(-case_submitter_id,-SampleID) %>%
  #dplyr::select(-Proliferation) %>% 
  remove_rownames() %>%
  column_to_rownames("sampleid")
for (name in colnames(Figure6.SD)) {
  Figure6.SD[[name]] = as.numeric(as.character(Figure6.SD[[name]]))
}

Figure6.scale <- scale(Figure6.SD)

Figure6.SD.T <- Figure6.SD %>% 
  mutate(group=Figure6.S$CBX2Group) %>%
  gather(Signature,Score,-group) %>%
  group_by(Signature) %>%
  rstatix::t_test(Score~group,detailed = T) %>%
  mutate(Signature = factor(Signature,levels = colnames(Figure6.scale))) %>%
  arrange(Signature) %>%
  rstatix::add_significance("p")

CEP55G = Figure6.S$CEP55Group

harow = rowAnnotation(foo = anno_text(Figure6.SD.T$p.signif, just='left'))

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = nejm[2:1]),
                                        labels = c("CBX2 low","CBX2 high"), 
                                        labels_gp = gpar(col = "white", fontsize = 12)),
                       #Age= factor(Age1,levels = c(">60","<=60")),
                       #Gender=factor(Gender1,levels = c("Female","Male")),
                       #Stage=factor(Stage1,levels = c("I","II","III","IV","Unknown")),
                       #T_stage=factor(T1,levels = c("T1","T2","T3","T4","Unknown")),
                       #M_stage=factor(M1,levels = c("M0","M1","Unknown")),
                       #N_stage=factor(N1,levels = c("N0","N1","Unknown")),
                       #Grade=factor(Grade1,levels = c("G1","G2","G3","G4","Unknown")),
                       CEP55= factor(CEP55G,levels = c("CEP55high","CEP55low","Unknown")),
                       #TP53 = factor(TP53,levels = c("Yes","No","Unknown")),
                       #CTNNB1 = factor(CTNNB1,levels = c("Yes","No","Unknown")),
                       col = list(#Age = c(">60" = "#BC3C2999", "<=60" = "#FFDC9199"),
                                  #Gender = c("Female" = nejm[3],"Male"="#ABDDA4"),
                                  #ISUP = c("1" = "#FFFFBF" , "2-3" = "#ABDDA4", "4" = "#9970AB",'5'="#D6604D",'Unknown'='#FFFFFF'),
                                  #TMPRSS2ERG = c("Positive" = "#BC3C2999", "Negative" = "#FFDC9199",'Unknown'='#FFFFFF'),
                                  #Stage = c('I'="#FFFFBF" ,'II' ="#ABDDA4","III" = "#9970AB",'IV'="#D6604D" , 'Unknown'='#FFFFFF'),
                                  #T_stage =  c('T1'="#FFFFBF" ,'T2' ="#ABDDA4","T3" = "#9970AB",'T4'="#D6604D" , 'Unknown'='#FFFFFF'),
                                  #M_stage =  c('M0' ="#ABDDA4","M1" = "#9970AB",'Unknown'='#FFFFFF'),
                                  #N_stage =  c('N0' ="#ABDDA4","N1" = "#9970AB",'Unknown'='#FFFFFF'),
                                  #Grade = c('G1'="#FFFFBF" ,'G2' ="#ABDDA4","G3" = "#9970AB",'G4'="#D6604D",'Unknown'='#FFFFFF'),
                                  CEP55 = c("CEP55high" = "#BC3C29FF","CEP55low"="#ABDDA4")
                                  #TP53 = c("Yes" = "#BC3C29FF", "No"="white"),
                                  #CTNNB1 = c("Yes" = "#BC3C29FF", "No"="white")
                                  ),#∏˜∏ˆ¡Ÿ¥≤–≈œ¢µƒ—’…´£¨Anatomic_location—’…´ ‰»Î“≤Õ¨∆‰À˚£¨∏¸∫√—’…´¥Ó≈‰–≈œ¢«Î≤ŒøºFigureYa28color
                       show_legend = rep(TRUE, 3),# «∑Ò“™œ‘ æannotation legend
                       annotation_height = unit(rep(5,2), "mm"),#¡Ÿ¥≤annotationµƒ∏ﬂ∂»
                       annotation_legend_param = list(
                         #Age = list(title = "Age"),
                         #Gender = list(title = "Gender"),
                         #Stage = list(title = "Stage"),
                         #T_stage = list(title = "T_stage"),
                         #M_stage = list(title = "M_stage"),
                         #N_stage = list(title = "N_stage"),
                         #Grade = list(title = "Grade"),
                         #CEP55= CEP55G,
                         CEP55 = list(title = "CEP55 expression",
                                      at=c('CEP55high','CEP55low'),
                                      labels=c('High','Low'))
                         #TP53 = list(title = "TP53"),
                         #CTNNB1 = list(title = "CTNNB1")
                         ),
                       simple_anno_size = unit(0.3, "cm"),
                       na_col = "white",
                       annotation_name_side = "left"
)

ht_opt(legend_border='black',
       heatmap_border = TRUE,
       annotation_border = TRUE) 
col_fun <- colorRamp2(c(-2, 0,2), c("#377EB8", "white", "#E41A1C"))
subtype <- Subtype
ht <- Heatmap(t(Figure6.scale), col = col_fun, 
              name = "Score",
              cluster_rows = F, 
              cluster_columns = T,         
              show_row_names = TRUE, 
              row_names_side = "left", #––√˚Œª÷√
              show_column_names = FALSE,
              column_split = subtype, #∏˘æ› ˝æ›–ﬁ∏ƒ
              #row_split = c(rep(1,3),rep(2,6),rep(3,15),rep(4,1),rep(5,6),rep(6,3),rep(7,6)),#∏˘æ› ˝æ›–ﬁ∏ƒ
              #row_gap = unit(c(1,1,1,1,1,1) ,"mm"),
              column_gap = unit(c(2), "mm"),
              row_title = NULL,column_title = NULL,
              right_annotation = harow,
              top_annotation = ha)
pdf("Figure6.S111.pdf",width = 10,height = 6)
draw(ht,padding = unit(c(0.1, 4, 0.1, 0.5), "cm"),
     annotation_legend_side = "right", 
     heatmap_legend_side = "right")
dev.off()

Figure6.S <- read_xlsx("Figure6.S.xlsx")

Figure6.SD <- Figure6.S %>% dplyr::select(-Grade,-age,-gender,-T,-N,-M,-TP53,
                                          -CTNNB1,-Stage,-CEP55Group,-CBX2Group) %>%
  dplyr::select(-case_submitter_id,-SampleID) %>%
  #dplyr::select(-Proliferation) %>% 
  remove_rownames() %>%
  column_to_rownames("sampleid")
for (name in colnames(Figure6.SD)) {
  Figure6.SD[[name]] = as.numeric(as.character(Figure6.SD[[name]]))
}

Figure6.scale <- scale(Figure6.SD)

Figure6.SD.T <- Figure6.SD %>% 
  mutate(group=Figure6.S$CEP55Group) %>%
  gather(Signature,Score,-group) %>%
  group_by(Signature) %>%
  rstatix::t_test(Score~group,detailed = T) %>%
  mutate(Signature = factor(Signature,levels = colnames(Figure6.scale))) %>%
  arrange(Signature) %>%
  rstatix::add_significance("p")

Subtype = Figure6.S$CEP55Group

harow = rowAnnotation(foo = anno_text(Figure6.SD.T$p.signif, just='left'))

ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = nejm[1:2]),
                                        labels = c("CEP55 high","CEP55 low"), 
                                        labels_gp = gpar(col = "white", fontsize = 12)),
                       #Age= factor(Age1,levels = c(">60","<=60")),
                       #Gender=factor(Gender1,levels = c("Female","Male")),
                       #Stage=factor(Stage1,levels = c("I","II","III","IV","Unknown")),
                       #T_stage=factor(T1,levels = c("T1","T2","T3","T4","Unknown")),
                       #M_stage=factor(M1,levels = c("M0","M1","Unknown")),
                       #N_stage=factor(N1,levels = c("N0","N1","Unknown")),
                       #Grade=factor(Grade1,levels = c("G1","G2","G3","G4","Unknown")),
                       #CEP55= factor(CEP55G,levels = c("CEP55high","CEP55low","Unknown")),
                       #TP53 = factor(TP53,levels = c("Yes","No","Unknown")),
                       #CTNNB1 = factor(CTNNB1,levels = c("Yes","No","Unknown")),
                       #col = list(#Age = c(">60" = "#BC3C2999", "<=60" = "#FFDC9199"),
                         #Gender = c("Female" = nejm[3],"Male"="#ABDDA4"),
                         #ISUP = c("1" = "#FFFFBF" , "2-3" = "#ABDDA4", "4" = "#9970AB",'5'="#D6604D",'Unknown'='#FFFFFF'),
                         #TMPRSS2ERG = c("Positive" = "#BC3C2999", "Negative" = "#FFDC9199",'Unknown'='#FFFFFF'),
                         #Stage = c('I'="#FFFFBF" ,'II' ="#ABDDA4","III" = "#9970AB",'IV'="#D6604D" , 'Unknown'='#FFFFFF'),
                         #T_stage =  c('T1'="#FFFFBF" ,'T2' ="#ABDDA4","T3" = "#9970AB",'T4'="#D6604D" , 'Unknown'='#FFFFFF'),
                         #M_stage =  c('M0' ="#ABDDA4","M1" = "#9970AB",'Unknown'='#FFFFFF'),
                         #N_stage =  c('N0' ="#ABDDA4","N1" = "#9970AB",'Unknown'='#FFFFFF'),
                         #Grade = c('G1'="#FFFFBF" ,'G2' ="#ABDDA4","G3" = "#9970AB",'G4'="#D6604D",'Unknown'='#FFFFFF'),
                         #CEP55 = c("CEP55high" = "#BC3C29FF","CEP55low"="#ABDDA4")
                         #TP53 = c("Yes" = "#BC3C29FF", "No"="white"),
                         #CTNNB1 = c("Yes" = "#BC3C29FF", "No"="white")
                       #),#∏˜∏ˆ¡Ÿ¥≤–≈œ¢µƒ—’…´£¨Anatomic_location—’…´ ‰»Î“≤Õ¨∆‰À˚£¨∏¸∫√—’…´¥Ó≈‰–≈œ¢«Î≤ŒøºFigureYa28color
                       show_legend = rep(TRUE, 2),# «∑Ò“™œ‘ æannotation legend
                       annotation_height = unit(rep(5,1), "mm"),#¡Ÿ¥≤annotationµƒ∏ﬂ∂»
                       annotation_legend_param = list(
                         #Age = list(title = "Age"),
                         #Gender = list(title = "Gender"),
                         #Stage = list(title = "Stage"),
                         #T_stage = list(title = "T_stage"),
                         #M_stage = list(title = "M_stage"),
                         #N_stage = list(title = "N_stage"),
                         #Grade = list(title = "Grade"),
                         #CEP55= CEP55G,
                         #CEP55 = list(title = "CEP55 expression",
                         #             at=c('CEP55high','CEP55low'),
                         #             labels=c('High','Low'))
                         #TP53 = list(title = "TP53"),
                         #CTNNB1 = list(title = "CTNNB1")
                       ),
                       simple_anno_size = unit(0.3, "cm"),
                       na_col = "white",
                       annotation_name_side = "left"
)

ht_opt(legend_border='black',
       heatmap_border = TRUE,
       annotation_border = TRUE) 
col_fun <- colorRamp2(c(-2, 0,2), c("#377EB8", "white", "#E41A1C"))
subtype <- Subtype
ht <- Heatmap(t(Figure6.scale), col = col_fun, 
              name = "Score",
              cluster_rows = F, 
              cluster_columns = T,         
              show_row_names = TRUE, 
              row_names_side = "left", #––√˚Œª÷√
              show_column_names = FALSE,
              column_split = subtype, #∏˘æ› ˝æ›–ﬁ∏ƒ
              #row_split = c(rep(1,3),rep(2,6),rep(3,15),rep(4,1),rep(5,6),rep(6,3),rep(7,6)),#∏˘æ› ˝æ›–ﬁ∏ƒ
              #row_gap = unit(c(1,1,1,1,1,1) ,"mm"),
              column_gap = unit(c(2), "mm"),
              row_title = NULL,column_title = NULL,
              right_annotation = harow,
              top_annotation = ha)
pdf("Figure6.S112.pdf",width = 9,height = 6)
draw(ht,padding = unit(c(0.1, 4, 0.1, 0.5), "cm"),
     annotation_legend_side = "right", 
     heatmap_legend_side = "right")
dev.off()
#### IPS => CBX2 ####
library(readxl)
data <- read_xlsx("Figure6.S2.xlsx")
data$IPS <- as.numeric(data$IPS)
load("Survival.Signature.RiskScore.Rdata")
risk_score_table_multi_cox2 <- risk_score_table_multi_cox2 %>% mutate(Sample.ID = rownames(.))
TIDE.Data <- read.csv("../LIHC.FPKM.TIDE.csv") %>% 
  mutate(Sample.ID = substr(Patient,1,15)) %>%
  merge(.,data,by.x="Sample.ID",by.y="sampleid")
write.csv(TIDE.Data,file = "SourceData.F6GH.csv")
ImmunoPhenoScore <- read.csv("TCIA-ClinicalData.tsv",sep = "\t") %>%
  merge(data,by.x="barcode",by.y="case_submitter_id")

p12 <- ggplot(TIDE.Data,aes(CBX2Group,TIDE.x,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="TIDE",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.2,label.y = 0.8)
p12
ggsave(p12,filename = "CBX2.Group.TIDE.pdf",height = 2.5,width = 2)

p12 <- ggplot(TIDE.Data,aes(CBX2Group,MDSC,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="MDSC",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.2)
p12
ggsave(p12,filename = "CBX2.Group.MDSC.pdf",height = 2.5,width = 2)

p12 <- ggplot(TIDE.Data,aes(CEP55Group,MDSC,fill=CEP55Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CEP55Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CEP55Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CEP55Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="MDSC",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.2)
p12
ggsave(p12,filename = "CEP55.Group.MDSC.pdf",height = 2.5,width = 2)


p121 <- ggplot(TIDE.Data,aes(CBX2Group,Dysfunction.x,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="Dysfunction",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.2,label.y = 0.8)
p121
ggsave(p121,filename = "CBX2.Group.TIDE.Dysfunction.pdf",height = 2.5,width = 2)

p122 <- ggplot(TIDE.Data,aes(CBX2Group,Exclusion.x,fill=CBX2Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CBX2Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CBX2Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CBX2Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="Exclusion",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.1,label.y = 2.1)
p122
ggsave(p122,filename = "CBX2.Group.TIDE.Exclusion.pdf",height = 2.5,width = 2)

p12 <- ggplot(TIDE.Data,aes(CEP55Group,TIDE.x,fill=CEP55Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CEP55Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CEP55Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CEP55Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="TIDE",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.2,label.y = 0.8)
p12
ggsave(p12,filename = "CEP55.Group.TIDE.pdf",height = 2.5,width = 2)

p121 <- ggplot(TIDE.Data,aes(CEP55Group,Dysfunction.x,fill=CEP55Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CEP55Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CEP55Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CEP55Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="Dysfunction",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.2,label.y = 0.8)
p121
ggsave(p121,filename = "CEP55.Group.TIDE.Dysfunction.pdf",height = 2.5,width = 2)

p122 <- ggplot(TIDE.Data,aes(CEP55Group,Exclusion.x,fill=CEP55Group))+
  #geom_violin()+
  #geom_boxplot(color="white")+
  geom_boxplot(linetype="dashed",color="black")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..,
                   fill=CEP55Group),
               color="black")+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=CEP55Group),width=0.3,size=1)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=CEP55Group),width=0.3,size=1)+
  #geom_jitter()+
  #ggbeeswarm::geom_beeswarm(priority = "descending")+
  #geom_line(aes(group=Line))+
  ggthemes::theme_few()+
  #facet_wrap(~gene,scales = "free_y")+
  #coord_flip()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  ggsci::scale_fill_nejm()+
  labs(y="Exclusion",x="")+
  theme(strip.text = element_text(size=15,color="black"))+
  #scale_y_log10()+
  ggpubr::stat_compare_means(method = "wilcox.test",label.x=1.1,label.y = 2.1)
p122
ggsave(p122,filename = "CEP55.Group.TIDE.Exclusion.pdf",height = 2.5,width = 2)

######## scRNA => CBX2 CEP55 ########
##### GSE125449 ##### 
library(ggplot2)
metadata <- read.csv("LIHC_GSE125449_aPDL1aCTLA4_CellMetainfo_table.tsv",sep = "\t",
                     check.names = F)
data<-read.csv("CBX2.CEP55.GSE125449.csv",sep=",",header=T) %>%
  merge(metadata,by.x="X",by.y="Cell")
write.csv(data,"SourceData.F6AB.csv")
Theme2<-ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))  #
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
            "#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080",
            "#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`Celltype (major-lineage)`),size=0.1)+
  scale_color_manual("Cell type",values = allcolour)+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p,filename = "scRNA.GSE125449.UMAP.pdf",height = 2.7,width = 5)
library(RColorBrewer)
p1<-ggplot(data %>% arrange(CBX2) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CBX2`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Reds")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE125449.UMAP.CBX2.pdf",height = 2.7,width = 5)

##### GSE140228_10X #####
library(ggplot2)
metadata <- read.csv("LIHC_GSE140228_10X_CellMetainfo_table.tsv",sep = "\t",
                     check.names = F)
data<-read.csv("CBX2.CEP55.GSE140228_10X.csv",sep=",",header=T) %>%
  merge(metadata,by.x="X",by.y="Cell")
Theme2<-ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))  #
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
            "#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080",
            "#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`Celltype (major-lineage)`),size=0.1)+
  scale_color_manual("Cell type",values = allcolour)+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p,filename = "scRNA.GSE140228_10X.UMAP.pdf",height = 4,width = 5)
library(RColorBrewer)
p1<-ggplot(data %>% arrange(CEP55) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CEP55`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Reds")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE140228_10X.UMAP.CEP55.pdf",height = 4,width = 5)


##### GSE146115 #####
library(ggplot2)
metadata <- read.csv("LIHC_GSE146115_CellMetainfo_table.tsv",sep = "\t",
                     check.names = F)
data<-read.csv("CBX2.CEP55.GSE146115.csv",sep=",",header=T) %>%
  merge(metadata,by.x="X",by.y="Cell")
Theme2<-ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))  #
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
            "#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080",
            "#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`Celltype (major-lineage)`),size=0.1)+
  scale_color_manual("Cell type",values = allcolour)+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p,filename = "scRNA.GSE146115.UMAP.pdf",height = 4,width = 5)
library(RColorBrewer)
p1<-ggplot(data %>% arrange(CEP55) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CEP55`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Reds")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE146115.UMAP.CEP55.pdf",height = 4,width = 5)

p1<-ggplot(data %>% arrange(CBX2) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CBX2`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE146115.UMAP.CBX2.pdf",height = 4,width = 5)

##### GSE166635 #####
library(ggplot2)
metadata <- read.csv("LIHC_GSE166635_CellMetainfo_table.tsv",sep = "\t",
                     check.names = F)
data<-read.csv("CBX2.CEP55.GSE166635.csv",sep=",",header=T) %>%
  merge(metadata,by.x="X",by.y="Cell")
Theme2<-ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))  #
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
            "#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080",
            "#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`Celltype (major-lineage)`),size=0.1)+
  scale_color_manual("Cell type",values = allcolour)+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p,filename = "scRNA.GSE166635.UMAP.pdf",height = 4,width = 5)
library(RColorBrewer)
p1<-ggplot(data %>% arrange(CEP55) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CEP55`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Reds")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE166635.UMAP.CEP55.pdf",height = 4,width = 5)

p1<-ggplot(data %>% arrange(CBX2) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CBX2`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE166635.UMAP.CBX2.pdf",height = 4,width = 5)
##### GSE98638 #####
library(ggplot2)
metadata <- read.csv("LIHC_GSE98638_CellMetainfo_table.tsv",sep = "\t",
                     check.names = F)
data<-read.csv("CBX2.CEP55.GSE98638.csv",sep=",",header=T) %>%
  merge(metadata,by.x="X",by.y="Cell")
Theme2<-ggthemes::theme_few()+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))  #
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
            "#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080",
            "#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p<-ggplot(data,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`Celltype (major-lineage)`),size=0.1)+
  scale_color_manual("Cell type",values = allcolour)+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p,filename = "scRNA.GSE98638.UMAP.pdf",height = 4,width = 5)
library(RColorBrewer)
p1<-ggplot(data %>% arrange(CEP55) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CEP55`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Reds")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE98638.UMAP.CEP55.pdf",height = 4,width = 5)

p1<-ggplot(data %>% arrange(CBX2) %>%
             mutate(),
           aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=`CBX2`),size=0.5)+
  #scale_color_manual("Cell type",values = allcolour)+
  scale_color_gradientn(colors=brewer.pal(9,"Greens")[c(1,3,5,7,9)])+
  #geom_label(aes(label=`Celltype (major-lineage)`),
  #           nudge_x=0,alpha=.5,size=5)+
  labs(x="UMAP_1",y="UMAP_2")+
  theme_bw()+theme(text=element_text(size=18))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=4)))
ggsave(p1,filename = "scRNA.GSE98638.UMAP.CBX2.pdf",height = 4,width = 5)
#### HPA ####
CBX2.HPA <- read.csv("HPA.CBX2.txt",sep = "\t",header = F)
CEP55.HPA <- read.csv("HPA.CEP55.txt",sep = "\t",header = F)
HPA.D <- rbind(CBX2.HPA,CEP55.HPA)
write.csv(HPA.D,file = "SourceData.F6C.csv")
p1 <- ggplot(HPA.D %>% filter(V2=="CBX2") %>%
         mutate(V4=factor(V4,levels = rev(CBX2.HPA$V4))),aes(V2,V4))+
  geom_tile(aes(fill=V7))+
  geom_point(aes(size=V6),shape=1)+
  scale_size_continuous("no. of cells",range = c(0,5),breaks = c(0,4,8))+
  labs(x="",y="")+
  scale_fill_gradientn("Z score",colors=brewer.pal(9,"Reds"),breaks = c(0,2,4))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=0.5)))+
  theme(legend.position = "top")
p2 <- ggplot(HPA.D %>% filter(V2=="CEP55") %>%
               mutate(V4=factor(V4,levels = rev(CBX2.HPA$V4))),aes(V2,V4))+
  geom_tile(aes(fill=V7))+
  geom_point(aes(size=V6),shape=1)+
  scale_size_continuous("no. of cells",range = c(0,5),breaks = c(0,20,40))+
  labs(x="",y="")+
  scale_fill_gradientn("Z score",colors=brewer.pal(9,"Greens"),breaks = c(0,10,20))+
  Theme2+guides(colour = guide_legend(override.aes = list(size=0.5)))+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        legend.position = "top")
p1 %>% aplot::insert_right(p2,width = 1) %>%
  ggsave(filename="scRNA.HPA.pdf",height=4,width=4.3)
p1 %>% aplot::insert_right(p2,width = 1) %>%
  ggsave(filename="scRNA.HPA2.pdf",height=6,width=2.7)


##### CERES #####
HCC.CERES <- read.csv("../../Depmap/CBX2 Expression Public 22Q4.csv") #%>% 
  #filter(`Lineage.Subtype` == "Hepatocellular Carcinoma")
load("TCGA.ICGC.Pearson.CEP55.Inter0.5.Rdata")
Inter.Genes2 <- Inter.Genes %>% 
  dplyr::select(Gene.name,Rho.x,Rho.y,Pvalue.x,Pvalue.y) %>%
  magrittr::set_colnames(c("Gene","TCGA-LIHC.Rho","ICGC-LIRI-JP.Rho",
                           "TCGA-LIHC.Pvalue","ICGC-LIRI-JP.Pvalue"))

Inter.Genes.FigData <- data.frame(Gene=c(Inter.Genes2$Gene,Inter.Genes2$Gene),
                                  Rho=c(Inter.Genes2$`TCGA-LIHC.Rho`,
                                        Inter.Genes2$`ICGC-LIRI-JP.Rho`),
                                  Pvalue=c(Inter.Genes2$`TCGA-LIHC.Pvalue`,
                                           Inter.Genes2$`ICGC-LIRI-JP.Pvalue`)) %>%
  mutate(Dataset2=rep(c("TCGA-LIHC","ICGC-LIRI-JP"),c(362,362))) %>%
  arrange(Gene,Dataset2,Rho)
Select.Gene <- Inter.Genes.FigData %>%
  filter(Rho >= 0.85 | Rho <= -0.6)

Inter.Genes.FigData2 <- Inter.Genes.FigData %>% filter(Gene %in% Select.Gene$Gene)
write.csv(Inter.Genes.FigData2,file = "CEP55.Pearson.SYMBOLS.csv")

CERES.data <- read.csv("../../DepMap/CRISPR_(Project_Score,_Chronos).csv")
CERES.Target <- CERES.data %>% dplyr::select(X,CEP55,unique(Inter.Genes.FigData2$Gene)) %>%
  merge(HCC.CERES %>% dplyr::select(Depmap.ID,Primary.Disease),by.x="X",by.y="Depmap.ID")

CEP55.CERES <- CERES.Target %>%
  remove_rownames() %>% column_to_rownames("X") %>%
  gather(SYMBOL,Chronos,-Primary.Disease) %>%
  group_by(Primary.Disease,SYMBOL) %>%
  summarise(Chronos=median(Chronos))

p1234 <- ggplot(CEP55.CERES %>%
         mutate(SYMBOL=factor(SYMBOL,levels = c("CEP55",unique(Inter.Genes.FigData2$Gene)))),
       aes(SYMBOL,y=Primary.Disease))+
  geom_tile(aes(fill=Chronos))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle=90,vjust=0.5,hjust=1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colors = c(brewer.pal(9,"Reds")[9:1],"white","Blue"))
ggsave(p1234,filename = "CEP55.Chronos.pdf",
       height = 4,width = 10)

load("TCGA.ICGC.Pearson.CBX2.Inter.Rdata")
Inter.Genes2 <- Inter.Genes %>% 
  dplyr::select(Gene.name,Rho.x,Rho.y,Pvalue.x,Pvalue.y) %>%
  magrittr::set_colnames(c("Gene","TCGA-LIHC.Rho","ICGC-LIRI-JP.Rho",
                           "TCGA-LIHC.Pvalue","ICGC-LIRI-JP.Pvalue"))

Inter.Genes.FigData <- data.frame(Gene=c(Inter.Genes2$Gene,Inter.Genes2$Gene),
                                  Rho=c(Inter.Genes2$`TCGA-LIHC.Rho`,
                                        Inter.Genes2$`ICGC-LIRI-JP.Rho`),
                                  Pvalue=c(Inter.Genes2$`TCGA-LIHC.Pvalue`,
                                           Inter.Genes2$`ICGC-LIRI-JP.Pvalue`)) %>%
  mutate(Dataset2=rep(c("TCGA-LIHC","ICGC-LIRI-JP"),c(55,55)))
write.csv(Inter.Genes.FigData,file = "CBX2.Pearson.SYMBOLS.csv")

CERES.Target <- CERES.data %>% dplyr::select(X,CBX2,unique(Inter.Genes.FigData$Gene)) %>%
  merge(HCC.CERES %>% dplyr::select(Depmap.ID,Primary.Disease),by.x="X",by.y="Depmap.ID")

CBX2.CERES <- CERES.Target %>%
  remove_rownames() %>% column_to_rownames("X") %>%
  gather(SYMBOL,Chronos,-Primary.Disease) %>%
  group_by(Primary.Disease,SYMBOL) %>%
  summarise(Chronos=median(Chronos))

p12345 <- ggplot(CBX2.CERES %>%
         mutate(SYMBOL=factor(SYMBOL,levels = c("CBX2",unique(Inter.Genes.FigData$Gene)))),
       aes(SYMBOL,y=Primary.Disease))+
  geom_tile(aes(fill=Chronos))+
  theme(#legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9,angle=90,vjust=0.5,hjust=1),
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradientn(colors = c(brewer.pal(9,"Reds")[9:1],"white","Blue"))
ggsave(p12345,filename = "CBX2.Chronos.pdf",
       height = 4,width = 10)
### Boxplot ###
Chronos.Box <- CERES.data %>% 
  dplyr::select(X,POLD1,ORC1,MCM2,GINS1,CDT1,EFTUD2,ORC6,PLK1,
                CDK1,KIF11,KIF23,NDC80,TOP2A,NUF2) %>%
  merge(HCC.CERES %>% dplyr::select(Depmap.ID,Primary.Disease),by.x="X",by.y="Depmap.ID") %>%
  remove_rownames() %>% column_to_rownames("X") %>%
  gather(SYMBOL,Chronos,-Primary.Disease)

p123 <- ggplot(Chronos.Box,aes(reorder(SYMBOL,Chronos),Chronos))+
  geom_boxplot(aes(fill=SYMBOL))+
  ggthemes::theme_few()+
  theme(legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),#,angle=90,vjust=0.5,hjust=1
    axis.text.y = element_text(size=9),
    #axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="")+
  #ggthemes::scale_fill_stata()+
  scale_fill_manual(values = allcolour)+
  geom_hline(yintercept=-1,linetype=2)
ggsave(p123,filename = "CBX2.CEP55.Pearson.Chronos.pdf",
       height = 2,width=8)


##### CBX2 CEP55 => BRCA #####
data <- read.csv("../../mRNA.PanCancer.Exp/TCGA-BRCA.mRNA.TPM.csv",
                 row.names = 1,check.names = F)
library(data.table)
survival.data <- fread("../../GSCA/BRCA.survival.tsv.gz")
CBX2.D <- data["ENSG00000173894",] %>% t() %>%
  data.frame(check.names = F) %>%
  magrittr::set_colnames("CBX2") %>%
  mutate(CBX2=log2(CBX2+1)) %>%
  rownames_to_column("SampleID") %>%
  filter(str_detect(SampleID,"-01")) %>%
  mutate(sample_name=substr(SampleID,1,12)) %>%
  merge(survival.data,by="sample_name")
#CBX2.D$group <- if_else(CBX2.D$CBX2 >= median(CBX2.D$CBX2),"High","Low")
#diff=survdiff(Surv(os_days,os_status) ~ group,data = CBX2.D)
#pValue=1-pchisq(diff$chisq,df=1)
res.cut <- surv_cutpoint(CBX2.D,
                         time = "os_days",
                         event = "os_status",
                         variables = "CBX2")
dat <- surv_categorize(res.cut) #%>% data.frame() %>% mutate(Sample=SNHGs.ICI.data$Sample.ID)
write.csv(dat,"SourceData.F7C.csv")
fit <- survfit(Surv(os_days,os_status) ~ CBX2,
               data = dat)

cox <- coxph(Surv(os_days,os_status) ~ CBX2, data = CBX2.D)
summary(cox)
p <- ggsurvplot(fit,pval =TRUE, data = dat, 
           surv.median.line = "hv",
           legend.title = "CBX2",
           conf.int.style = "step",
           xlab = "Time (days)",
           #break.time.by = 500,
           risk.table = "abs_pct",
           #risk.table.y.text.col = T,
           #risk.table.y.text = FALSE,
           legend.labs = c("High", "Low"),
           #pval = TRUE,
           conf.int = TRUE,
           palette = "Set1",
           #ggtheme = theme.sur,
           risk.table.y.text.col = T,
           risk.table.y.text = F,
           ggtheme = theme.sur)
print(p)
graph2pdf(file="CBX2.BRCA.Survival.pdf",height = 3.5,width = 3.5)

theme.sur <- ggthemes::theme_few()+
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))

CBX2.CEP55.D <- data[c("ENSG00000173894",
                       "ENSG00000138180"),] %>% t() %>%
  data.frame(check.names = F) %>%
  magrittr::set_colnames(c("CBX2","CEP55")) %>%
  #mutate(CBX2=log2(CBX2+1)) %>%
  rownames_to_column("SampleID") %>%
  gather(SYMBOL,TPM,-SampleID) %>%
  mutate(TPM=log2(TPM+1)) %>%
  filter(str_detect(SampleID,"-01")) %>%
  spread(SYMBOL,TPM)
  #mutate(sample_name=substr(SampleID,1,12)) %>%
  #merge(survival.data,by="sample_name")

p21 <- ggscatter(CBX2.CEP55.D, x = "CBX2", y = "CEP55",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray")
          
)+labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) CEP55")+
  stat_cor(method = "pearson", 
           label.x = 1, label.y = 8)+Theme2
ggsave(p21,filename = "CBX2.CEP55.BRCA.Pearson.pdf",height = 2.5,width = 2.5)


##### CBX2 HIF1A MAP2K2 TOP2A => LIHC #####
data <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",
                 row.names = 1,check.names = F)

CBX2.CEP55.D <- data[c("ENSG00000173894",
                       "ENSG00000100644",
                       "ENSG00000126934",
                       "ENSG00000131747"),] %>% t() %>%
  data.frame(check.names = F) %>%
  magrittr::set_colnames(c("CBX2","HIF1A","MAP2K2","TOP2A")) %>%
  #mutate(CBX2=log2(CBX2+1)) %>%
  rownames_to_column("SampleID") %>%
  gather(SYMBOL,TPM,-SampleID) %>%
  mutate(TPM=log2(TPM+1)) %>%
  filter(str_detect(SampleID,"-01")) %>%
  spread(SYMBOL,TPM)
#mutate(sample_name=substr(SampleID,1,12)) %>%
#merge(survival.data,by="sample_name")

p21 <- ggscatter(CBX2.CEP55.D, x = "CBX2", y = "TOP2A",
                 add = "reg.line", conf.int = TRUE,    
                 add.params = list(fill = "lightgray")
                 
)+labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) TOP2A")+
  stat_cor(method = "pearson", 
           label.x = 1, label.y = 8)+Theme2
ggsave(p21,filename = "CBX2.TOP2A.LIHC.Pearson.pdf",height = 2.5,width = 2.5)

##### CBX2 CEP55 #####
files <- list.files("../../mRNA.PanCancer.Exp/")
Pearson.D <- data.frame()
for(file in files){
  data <- read.csv(file.path("../../mRNA.PanCancer.Exp/",file),
                   row.names = 1,check.names = F)
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  CBX2.CEP55.D <- data[c("ENSG00000173894",
                         "ENSG00000138180"),] %>% t() %>%
    data.frame(check.names = F) %>%
    magrittr::set_colnames(c("CBX2","CEP55")) %>%
    #mutate(CBX2=log2(CBX2+1)) %>%
    rownames_to_column("SampleID") %>%
    gather(SYMBOL,TPM,-SampleID) %>%
    mutate(TPM=log2(TPM+1)) %>%
    filter(str_detect(SampleID,"-01") | str_detect(SampleID,"-03")) %>%
    spread(SYMBOL,TPM)
  test <- cor.test(CBX2.CEP55.D$CBX2,CBX2.CEP55.D$CEP55,
                   method = "pearson",alternative = "two.sided",exact = F)
  Pearson.D <- data.frame(p=test$p.value,r=test$estimate,project=project,
                          n=nrow(CBX2.CEP55.D)) %>%
    rbind.data.frame(Pearson.D)
}
write.csv(Pearson.D,file="Pancancer.CBX2.CEP55.Pearson.csv",row.names = F)
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))

Pearson.D2 <- Pearson.D %>% mutate(xAxis=paste(project,"(n=",n,")",sep = "")) %>%
  arrange(r) %>% mutate(xAxis=factor(xAxis,levels = .$xAxis))

p131 <- ggplot(Pearson.D2,aes(xAxis,y=r))+
  geom_segment(aes(x=xAxis,xend=xAxis,y=0,yend=r))+
  geom_point(aes(fill=xAxis,size=-log10(p)),color="black",shape=21)+
  Theme2+
  theme(legend.position = "none",
        axis.text.x = element_text(size=9,angle=90,vjust=0.5,hjust=1))+
  labs(x="",y="Pearson'CC")+
  scale_fill_manual(values = color)+
  geom_hline(yintercept = c(-0.2,0.2,0.4,0.6,0.8),linetype=2)+
  geom_hline(yintercept = c(0))+
  scale_size_continuous(range = c(2,6))+
  ylim(-0.3,1)+
  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4,0.6,0.8))
ggsave(p131,filename = "PanCancer.CBX2.CEP55.Pearson.pdf",
       height = 3,width = 8)

p132 <- ggplot(Pearson.D2,aes(xAxis,y=r))+
  geom_segment(aes(x=xAxis,xend=xAxis,y=0,yend=r))+
  geom_point(aes(size=-log10(p)),color="black",shape=21)+
  Theme2+
  theme(#legend.position = "none",
        axis.text.x = element_text(size=9,angle=90,vjust=0.5,hjust=1))+
  labs(x="",y="Pearson'CC")+
  scale_fill_manual(values = color)+
  geom_hline(yintercept = c(-0.2,0.2,0.4,0.6,0.8),linetype=2)+
  geom_hline(yintercept = c(0))+
  scale_size_continuous(range = c(2,6))+
  ylim(-0.3,1)+
  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4,0.6,0.8))
ggsave(p132,filename = "PanCancer.CBX2.CEP55.Pearson2.pdf",
       height = 3,width = 8)

###### CBX2 CEP55 => P21 CDKN1A CDKN1B 	CDKN1C P16 CDKN2A CDKN2B CDKN2C #######
#### Mutation ####
library(data.table)
files <- list.files("../../mRNA.PanCancer.Exp/")
for (file in files) {
  project = file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  Mutation.D <- fread(paste("../../Clinical.XENA/",project,
                            "/TCGA-",project,".mutect2_snv.tsv",sep = "")) %>%
    filter(gene == "CDKN1A")
}

#### Expression ####
# ENSG00000124762 CDKN1A
# ENSG00000111276 CDKN1B
# ENSG00000129757 CDKN1C
# ENSG00000147889 CDKN2A
# ENSG00000147883 CDKN2B
# ENSG00000123080 CDKN2C
files <- list.files("../../mRNA.PanCancer.Exp/")
Pearson.D <- data.frame()
for(file in files){
  data <- read.csv(file.path("../../mRNA.PanCancer.Exp/",file),
                   row.names = 1,check.names = F)
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  CBX2.CEP55.D <- data[c("ENSG00000173894",
                         "ENSG00000138180",
                         "ENSG00000124762","ENSG00000111276","ENSG00000129757",
                         "ENSG00000147889","ENSG00000147883","ENSG00000123080"),] %>% t() %>%
    data.frame(check.names = F) %>%
    magrittr::set_colnames(c("CBX2","CEP55","CDKN1A","CDKN1B","CDKN1C","CDKN2A","CDKN2B","CDKN2C")) %>%
    #mutate(CBX2=log2(CBX2+1)) %>%
    rownames_to_column("SampleID") %>%
    gather(SYMBOL,TPM,-SampleID) %>%
    mutate(TPM=log2(TPM+1)) %>%
    filter(str_detect(SampleID,"-01") | str_detect(SampleID,"-03")) %>%
    spread(SYMBOL,TPM)
  for (gene in c("CDKN1A","CDKN1B","CDKN1C","CDKN2A","CDKN2B","CDKN2C")) {
    test <- cor.test(CBX2.CEP55.D$CBX2,CBX2.CEP55.D[[gene]],
                     method = "pearson",alternative = "two.sided",exact = F)
    Pearson.D <- data.frame(p=test$p.value,r=test$estimate,project=project,
                            n=nrow(CBX2.CEP55.D),source="CBX2",target=gene) %>%
      rbind.data.frame(Pearson.D)
    
    test <- cor.test(CBX2.CEP55.D$CEP55,CBX2.CEP55.D[[gene]],
                     method = "pearson",alternative = "two.sided",exact = F)
    Pearson.D <- data.frame(p=test$p.value,r=test$estimate,project=project,
                            n=nrow(CBX2.CEP55.D),source="CEP55",target=gene) %>%
      rbind.data.frame(Pearson.D)
  }
}
write.csv(Pearson.D,file="Pancancer.CBX2.CEP55.CellCycleInhibitor.Pearson.csv",row.names = F)
library(ggsci)
color = c(pal_d3("category20")(20),pal_d3("category20b")(20),pal_d3("category20c")(20),pal_d3("category10")(10))

#Pearson.D %>% filter(target=="CDKN1A")
Pearson.D2 <- Pearson.D %>% mutate(xAxis=paste(project,"(n=",n,")",sep = "")) #%>%
  #arrange(r) %>% mutate(xAxis=factor(xAxis,levels = .$xAxis))

data <- read.csv("../../PANCANCER_IC_Tue Jun 28 14_55_17 2022.csv")

ggplot(Pearson.D2 %>% filter(source=="CBX2" & target=="CDKN1A") %>%
         arrange(r) %>% mutate(xAxis=factor(xAxis,levels = .$xAxis)),
       aes(xAxis,y=r))+
  geom_segment(aes(x=xAxis,xend=xAxis,y=0,yend=r))+
  geom_point(aes(size=-log10(p)),color="black",shape=21)+
  Theme2+
  theme(#legend.position = "none",
    axis.text.x = element_text(size=9,angle=90,vjust=0.5,hjust=1))+
  labs(x="",y="Pearson'CC")+
  scale_fill_manual(values = color)+
  geom_hline(yintercept = c(-0.2,0.2,0.4,0.6,0.8),linetype=2)+
  geom_hline(yintercept = c(0))+
  scale_size_continuous(range = c(2,6))+
  ylim(-0.3,1)+
  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4,0.6,0.8))
write.csv(Pearson.D2,file = "CBX2.CEP55.CDKN1A.CDKN2A.Pancancer.Pearson.csv")

files <- list.files("../../mRNA.PanCancer.Exp/")
CBX2.CEP55.D <- data.frame()
for(file in files){
  data <- read.csv(file.path("../../mRNA.PanCancer.Exp/",file),
                   row.names = 1,check.names = F)
  project <- file %>% str_remove_all("\\..*") %>% str_remove_all("TCGA-")
  print(project)
  CBX2.CEP55.D <- data[c("ENSG00000173894",
                         "ENSG00000138180",
                         "ENSG00000124762","ENSG00000111276","ENSG00000129757",
                         "ENSG00000147889","ENSG00000147883","ENSG00000123080"),] %>% t() %>%
    data.frame(check.names = F) %>%
    magrittr::set_colnames(c("CBX2","CEP55","CDKN1A","CDKN1B","CDKN1C","CDKN2A","CDKN2B","CDKN2C")) %>%
    #mutate(CBX2=log2(CBX2+1)) %>%
    rownames_to_column("SampleID") %>%
    gather(SYMBOL,TPM,-SampleID) %>%
    mutate(TPM=log2(TPM+1)) %>%
    filter(str_detect(SampleID,"-01") | str_detect(SampleID,"-03")) %>%
    spread(SYMBOL,TPM) %>%
    mutate(Cancer=project) %>%
    rbind.data.frame(CBX2.CEP55.D)
}

cor.test(CBX2.CEP55.D$CBX2,
         CBX2.CEP55.D$CDKN1A%>% as.numeric(),method="pearson")


library(ggpubr)
p211 <- ggscatter(CBX2.CEP55.D, x = "CBX2", y = "CDKN1A",
                 #add = "reg.line", 
                 conf.int = TRUE, color = "Cancer",size = 1,
                 add.params = list(fill = "lightgray")
                   
  )+labs(x="log2(TPM+1) CBX2",y="log2(TPM+1) CDKN1A")+
  stat_cor(method = "pearson", 
           label.x = 1, label.y = 11)+Theme2+
  geom_smooth(method="lm",formula = "y~x",color="red")+
  scale_color_manual(values = color)+
  theme(legend.position = "none")
  
ggsave(p211,filename = "CBX2.CDKN1A.PanCancer1.Pearson.pdf",height = 2.5,width = 2.5)

data <- read.csv("../../Pancancer.Stemindex/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",
                 sep="\t",row.names = 1,header = T)
rownames(data) <- rownames(data) %>% str_remove_all("\\|.*")
#CBX2.Exp <- data["CBX2",] %>% data.frame() %>% dplyr::select(contains("\\.01A"))
#CDKN1A.Exp <- data["CDKN1A",] %>% data.frame() %>% dplyr::select(contains("\\.01A"))
# rownames(data)[str_detect(rownames(data),"CBX2")] # CBX2|84733
# rownames(data)[str_detect(rownames(data),"CEP55")] # CEP55|55165
# rownames(data)[str_detect(rownames(data),"CDKN1A")] # CDKN1A|1026
# rownames(data)[str_detect(rownames(data),"CDKN2A")] # CDKN2A|1029
cor.test(data["CBX2|84733",] %>% as.numeric(),
         data["CDKN1A|1026",]%>% as.numeric(),method="pearson")

middata <- data[c("CBX2|84733","CEP55|55165","CDKN1A|1026"),]



###### CBX2 CEP55 DRUGS ######
Exp <- read.csv("../../mRNA.PanCancer.Exp/TCGA-LIHC.mRNA.TPM.csv",
                row.names = 1,check.names = F)
Exp <- log2(Exp+1)
CBX2.CEP55.E <- Exp[c("ENSG00000173894",
                      "ENSG00000138180"),] %>% t() %>%
  data.frame(check.names = F) %>%
  magrittr::set_colnames(c("CBX2","CEP55")) %>%
  #mutate(CBX2=log2(CBX2+1)) %>%
  rownames_to_column("SampleID") #%>%
  #gather(SYMBOL,TPM,-SampleID)
data <- read.csv("../../DepMap/oncoPredict/PRISM.AUC.LIHC.csv",row.names=1) %>%
  merge(CBX2.CEP55.E,by.x="V1",by.y="SampleID") %>%
  mutate(CBX2G=if_else(CBX2 >= median(CBX2),"High","Low"),
         CEP55G=if_else(CEP55 >= median(CEP55),"High","Low"))

CBX2.Drugs <- data %>% dplyr::select(-CBX2,-CEP55,-CEP55G) %>%
  remove_rownames() %>% column_to_rownames("V1") %>%
  gather(Drugs,AUC,-CBX2G) %>%
  na.omit() %>%
  group_by(Drugs) %>%
  rstatix::t_test(AUC~CBX2G,detailed = T)
write.csv(CBX2.Drugs,file = "CBX2.PRISM.DRUGS.csv")

data1 <- read.csv("../../DepMap/oncoPredict/CTRP.AUC.LIHC.csv",row.names=1) %>%
  merge(CBX2.CEP55.E,by.x="V1",by.y="SampleID") %>%
  mutate(CBX2G=if_else(CBX2 >= median(CBX2),"High","Low"),
         CEP55G=if_else(CEP55 >= median(CEP55),"High","Low"))

CBX2.Drugs.CTRP <- data1 %>% dplyr::select(-CBX2,-CEP55,-CEP55G) %>%
  remove_rownames() %>% column_to_rownames("V1") %>%
  gather(Drugs,AUC,-CBX2G) %>%
  na.omit() %>%
  group_by(Drugs) %>%
  rstatix::t_test(AUC~CBX2G,detailed = T)
write.csv(CBX2.Drugs.CTRP,file = "CBX2.CTRP.DRUGS.csv")

# BAY.87.2243 »±—ı”’µº“Ú◊”-1 (HIF-1) “÷÷∆º¡°£
# 
# regorafenib
# Lenvatinib
# sorafenib
# Cabozantinib
# docetaxel and paclitaxel
PRISM.T <- data %>% dplyr::select(BAY.87.2243..BRD.BRD.K70463136.001.01.5.,
                                  binimetinib..BRD.BRD.K82244583.001.01.3.,
                                  #dinaciclib..BRD.BRD.K13662825.001.07.5.,
                                  voreloxin..BRD.BRD.K23677682.003.01.2.,
                                  floxuridine..BRD.BRD.K47832606.001.30.1.,
                                  triptolide..BRD.BRD.K39484304.001.16.5.,
                                  topotecan..BRD.BRD.K55696337.003.24.4.,
                                  CBX2G)
library(ggbeeswarm)
p1 <- ggplot(PRISM.T,aes(x=CBX2G,y=BAY.87.2243..BRD.BRD.K70463136.001.01.5.))+
  geom_boxplot(aes(fill=CBX2G))+
  geom_beeswarm()+
  Theme2+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5,size=12),
        axis.text = element_text(size=12))+
  labs(x="",y="BAY 87-2243")+
  ggtitle("HIF-1 inhibitor")+
  scale_fill_nejm()+
  scale_x_discrete(labels=c("CBX2 high","CBX2 low"))
ggsave(p1,filename = "CBX2.Drugs.BAY 87-2243.pdf",height = 2.5,width = 2.5)

p12 <- ggplot(PRISM.T,aes(x=CBX2G,y=floxuridine..BRD.BRD.K47832606.001.30.1.))+
  geom_boxplot(aes(fill=CBX2G))+
  geom_beeswarm()+
  Theme2+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5,size=12),
        axis.text = element_text(size=12))+
  labs(x="",y="Estimate AUC")+
  ggtitle("Floxuridine\nOncology antimetabolites")+
  scale_fill_nejm()+
  scale_x_discrete(labels=c("CBX2 high","CBX2 low"))
ggsave(p12,filename = "CBX2.Drugs.Floxuridine.pdf",height = 2.8,width = 2.5)

p12 <- ggplot(PRISM.T,aes(x=CBX2G,y=voreloxin..BRD.BRD.K23677682.003.01.2.))+
  geom_boxplot(aes(fill=CBX2G))+
  geom_beeswarm()+
  Theme2+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5,size=12),
        axis.text = element_text(size=12))+
  labs(x="",y="Estimate AUC")+
  ggtitle("Voreloxin\ntopoisomerase II inhibitor")+
  scale_fill_nejm()+
  scale_x_discrete(labels=c("CBX2 high","CBX2 low"))
ggsave(p12,filename = "CBX2.Drugs.Voreloxin.pdf",height = 2.8,width = 2.5)

p12 <- ggplot(PRISM.T,aes(x=CBX2G,y=voreloxin..BRD.BRD.K23677682.003.01.2.))+
  geom_boxplot(aes(fill=CBX2G))+
  geom_beeswarm()+
  Theme2+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5,size=12),
        axis.text = element_text(size=12))+
  labs(x="",y="Estimate AUC")+
  ggtitle("Voreloxin\ntopoisomerase II inhibitor")+
  scale_fill_nejm()+
  scale_x_discrete(labels=c("CBX2 high","CBX2 low"))
ggsave(p12,filename = "CBX2.Drugs.Voreloxin.pdf",height = 2.8,width = 2.5)

# regorafenib
# Lenvatinib
# sorafenib
CTRP.T <- data1 %>% dplyr::select(indisulam..CTRP.411874.,
                                  paclitaxel..CTRP.26956.,
                                  ruxolitinib..CTRP.639450.,
                                  olaparib..CTRP.411867.,
                                  regorafenib..CTRP.628613.,
                                  lenvatinib..CTRP.628639.,
                                  sorafenib..CTRP.349006.,
                                  CBX2G)
p121 <- ggplot(CTRP.T,aes(x=CBX2G,y=regorafenib..CTRP.628613.))+
  geom_boxplot(aes(fill=CBX2G))+
  geom_beeswarm()+
  Theme2+
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5,size=12),
        axis.text = element_text(size=12))+
  labs(x="",y="Estimate AUC")+
  ggtitle("Regorafenib")+
  scale_fill_nejm()+
  scale_x_discrete(labels=c("CBX2 high","CBX2 low"))
ggsave(p121,filename = "CBX2.Drugs.CTRP.Regorafenib.pdf",height = 2.8,width = 2.5)

######################################################
setwd("K:/HepG2-CBX2")
library(readxl)
library(RColorBrewer)
data <- read_xlsx("SnvSummaryTable.xlsx")
p1 <- ggplot(data,aes(cancertype,symbol,fill=EffectiveMut))+
  geom_tile(color="grey50")+
  geom_text(aes(label=EffectiveMut))+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size=12),
    panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_x_discrete(breaks="LIHC",labels="SNV")+
  scale_fill_gradientn("Mutation freq.",colors = brewer.pal(11,"RdBu")[5:2],
                       breaks=c(0,1,2),labels=c(0,1,2))

cnv.data <- read_xlsx("CnvSummaryTable.xlsx")
cnv.data2 <- cnv.data %>% dplyr::select(-entrez,-a_total,-d_total) %>%
  mutate(None=100-a_hete-d_hete-a_homo-d_homo) %>%
  gather(Type,Percent,-c("symbol","cancertype")) %>%
  mutate(Type=factor(Type,levels = c("None","d_hete","d_homo","a_hete","a_homo")))

p2 <- ggplot(cnv.data2,aes(Percent,symbol,fill=Type))+
  geom_col()+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size=9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="Percentage",y="")+
  scale_fill_manual("CNV type",values = c("grey60","#1B9E77","#7570B3","#D95F02","#E6AB02"))

data3 <- read_xlsx("CnvAndExpressionTable.xlsx")
p3 <- ggplot()+
  geom_point(data=data3,
             aes(cancertype,symbol,fill=spm,size=-log10(fdr)),
             shape=21,color="white",stroke=1)+
  geom_point(data=data3 %>% filter(fdr <= 0.05),
             aes(cancertype,symbol,fill=spm,size=-log10(fdr)),
             shape=21,color="black",stroke=1)+
  ggthemes::theme_few()+
  theme(legend.position = "bottom",
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        panel.border = element_rect(fill=NA,color="black", size=1.2, linetype="solid"))+
  labs(x="",y="")+
  scale_fill_gradient2("Spearman'CC",mid = "white",high = brewer.pal(11,"RdBu")[2],
                       low = brewer.pal(11,"RdBu")[11],limits=c(-0.2,0.2))
p4 <- p2 %>% aplot::insert_left(p1,width = 0.2) %>%
  aplot::insert_right(p3,width = 0.2)

ggsave(p4,filename = "CBX2.CEP55.CNV.SNV.pdf",height = 3,width = 10)

### cnv expression
library(data.table)
data <- fread("../../GSCA/LIHC.cnv.tsv.gz")
CBX2.CNV <- data %>% filter(symbol == "CBX2")


######################################################


######################################################


######################################################


###############  scCancer => Linux ####################
checkPkg <- function(pkg){
  return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")

library(devtools)
if(!checkPkg("harmony")) install_github("immunogenomics/harmony")

if(!checkPkg("RcppArmadillo")) install.packages("RcppArmadillo")
if(!checkPkg("RcppProgress")) install.packages("RcppProgress")
if(!checkPkg("NNLM")) install_github("linxihui/NNLM")

if(!checkPkg("liger")) install_github("MacoskoLab/liger")
devtools::install_github("wguo-research/scCancer")



############### ÷◊¡ˆ≤Ó“Ï∑÷Œˆ Wilcoxon ###############################

# read data

readCount<-read.table(file="examples/examples.countMatrix.tsv", 
                      header = T, row.names = 1)
conditions<-read.table(file="examples/examples.conditions.tsv", header = F)
conditions<-factor(t(conditions))

# edger TMM normalize

y <- DGEList(counts=readCount,group=conditions)
##Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)

# Run the Wilcoxon rank-sum test for each gene

pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")

# Calculate fold-change for each gene

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

# Output results based on FDR threshold

outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05
write.table(outRst[outRst$FDR<fdrThres,], file="examples/examples.WilcoxonTest.rst.tsv",sep="\t", quote=F,row.names = T,col.names = T)




######################################################

for (i in Inter.miRNA) {
  print(i)
}

clinical.Raw <- read.csv("..//../TCGA.clinical.csv")
clinical.Raw.2 <- clinical.Raw %>% mutate(PatientID = str_remove_all(xmls,".*clinical\\.") %>% str_remove_all("\\.xml")) %>%
  dplyr::select(PatientID,vital_status,tumor_tissue_site,histological_type,gender,
                days_to_birth,days_to_death,days_to_last_followup,race_list,bcr_patient_barcode,
                bcr_patient_uuid,icd_o_3_site,icd_o_3_histology,icd_10,stage_event,diagnosis,
                follow_ups,weight,height,venous_invasion,adjacent_hepatic_tissue_inflammation_extent_type,
                viral_hepatitis_serologies,alcohol_history_documented)


GeneType[GeneType$HGNC.symbol =="MCM3AP-AS1",] %>% unique() #ENSG00000215424
#GeneType[GeneType$HGNC.symbol =="MAP2",] %>% unique() #ENSG00000078018
GeneType[GeneType$HGNC.symbol =="MCM2",] %>% unique() #ENSG00000073111
GeneType[GeneType$HGNC.symbol =="LRRC1",] %>% unique() #ENSG00000137269

GeneType[GeneType$HGNC.symbol =="DUXAP8",] %>% unique() #ENSG00000206195

GeneType[GeneType$HGNC.symbol =="lnc-RGS5-1",] %>% unique() #ENSG00000232995

GeneType[GeneType$HGNC.symbol =="CKAP2L",] %>% unique() #ENSG00000169607
GeneType[GeneType$HGNC.symbol =="RMI2",] %>% unique() #ENSG00000175643


GeneType[GeneType$HGNC.symbol =="CBX2",] %>% unique() #ENSG00000173894
GeneType[GeneType$HGNC.symbol =="CEP55",] %>% unique() #ENSG00000138180

GeneType[GeneType$HGNC.symbol =="ESR1",] %>% unique() #ENSG00000091831
GeneType[GeneType$HGNC.symbol =="CPEB3",] %>% unique() #ENSG00000107864



#ENSG00000240498	CDKN2B-AS1	Up	CDKN2B-AS1









