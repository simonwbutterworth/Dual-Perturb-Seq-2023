##### GLOBAL PACKAGES AND FUNCTIONS FOR DUAL-PERTURB-SEQ-2023 #####

##### Libraries ##### 
library(Seurat)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggrepel)
library(Hotelling)
library(dichromat)
library(VISION)

##### Functions #####

ImportFeatureBC <- function(path.to.feature.bc.matrix, sample.name, sample.id.gex, sample.id.crispr) {
  #Import 10x data
  data <- Read10X(path.to.feature.bc.matrix)
  data.crispr <- data$`CRISPR Guide Capture`
  data.gex <- data$`Gene Expression`
  #Fix TGA4 gene names
  rownames(data.gex) <- sub("TGA4___", "TGME49", rownames(data.gex))
  #Split gene expresion data into GRCh38 and TGA4 genomes
  data.gex.grch38 <- data.gex[grepl("GRCh38",rownames(data.gex)),]
  rownames(data.gex.grch38) <- sub("GRCh38_", "", rownames(data.gex.grch38))
  data.gex.tga4 <- data.gex[grepl("TGME49_",rownames(data.gex)),]
  rownames(data.gex.tga4) <- sub("TGME49_", "", rownames(data.gex.tga4))
  #Create a Seurat object with the counts
  data <- CreateSeuratObject(counts=data.gex.grch38, assay="Hs", project=sample.name)
  data[["Tg"]] <- CreateAssayObject(data.gex.tga4)
  data[["sgRNA"]] <- CreateAssayObject(data.crispr)
  #Add metadata
  data$orig.ident <- sample.name
  data$sample.id.gex <- sample.id.gex
  data$sample.id.gex <- sample.id.crispr
  data$percent.mt <- PercentageFeatureSet(data, assay="Hs", pattern="^MT-")
  data$percent.tg <- data$nCount_Tg/(data$nCount_Tg+data$nCount_Hs)*100
  data
}

ImportFeatureBC_Mtx <- function(matrix, cells, features, sample.name, sample.id.gex, sample.id.crispr) {
  # Import data as matrix
  data <- ReadMtx(mtx = matrix, cells = cells, features = features)
  # Separate matrix into Hs gene, Tg gene, and sgRNA counts
  data.hs <- data[grepl("GRCh38_",rownames(data)),]
  data.tg <- data[grepl("TGA4___",rownames(data)),]
  data.sgrna <- data[rownames(data)[!rownames(data)%in%rownames(data.hs)&!rownames(data)%in%rownames(data.tg)],]
  # Fix gene names
  rownames(data.hs) <- gsub(x = rownames(data.hs), pattern = "GRCh38_", replacement = "")
  rownames(data.tg) <- gsub(x = rownames(data.tg), pattern = "TGA4___", replacement = "")
  # Create Seurat
  data <- CreateSeuratObject(counts=data.hs, project = sample.name, assay = "Hs")
  data[["Tg"]] <- CreateAssayObject(data.tg)
  data[["sgRNA"]] <- CreateAssayObject(data.sgrna)
  # Add metadata
  data$sample.id.gex <- sample.id.gex
  data$sample.id.gex <- sample.id.crispr
  data$percent.mt <- PercentageFeatureSet(data, assay="Hs", pattern="^MT-")
  data$percent.tg <- data$nCount_Tg/(data$nCount_Tg+data$nCount_Hs)*100
  # Return Seurat object
  data
}


Assign_sgRNAs <- function(object, reference.file) {
  # Identity the most highly expressed sgRNA in each cell and match to lookup reference to identify target gene
  # If nFeature_sgRNA != 1, returns ND

  # Create data frame for guide counts
  sgrna.counts <- data.frame(object@assays$sgRNA@counts) %>% t()
  # If nCount_sgrna==0, assign ND, else identify most highly expressed guide 
  object$guide <- ifelse(object$nFeature_sgRNA==1, yes = colnames(sgrna.counts)[max.col(sgrna.counts, ties.method ="random")], no = "ND")
  # Read in lookup file and select guide ID, guide target, and guide identity
  lookup <- read.csv(reference.file) %>% select(guide.r, gene, identity) %>% rename(guide = guide.r) 
  # Extract object metadata and left join to lookup
  lookup <- left_join(object[[]], lookup, by="guide")
  # Put guide target gene & identity back into object metadata
  object$guide.target.gene <- lookup$gene
  object$guide.target.gene[is.na(object$guide.target.gene)] <- "ND" 
  object$guide.target.identity <- lookup$identity 
  object$guide.target.identity[is.na(object$guide.target.identity)] <- "ND"
  # Return object as output
  object
}

CellsPerGene <- function(object) {
  data.frame(ident=object$guide.target.identity) %>% group_by(ident) %>% tally()
}

t2.test <- function(seurat.object, npcs=20) {
  results <- data.frame(ident="NA", t2.stat=NA, p.value=NA)
  idents <- data.frame(table(Idents(seurat.object))) %>% filter(Freq>=3) #Only test genes with at least 3 cells
  for (i in idents$Var1)  {
    hotelling.res <- Hotelling::hotelling.test(x = seurat.object@reductions$pca@cell.embeddings[which(Idents(seurat.object)==i),1:npcs],
                                               y = seurat.object@reductions$pca@cell.embeddings[which(Idents(seurat.object)!=i),1:npcs])
    results.test=data.frame(ident=i, t2.stat=hotelling.res$stats$statistic, p.value=hotelling.res$pval)
    results <- rbind(results, results.test)
    print(i)
  }
  results <- filter(results, ident!="NA")
  results$p.value.bh <- p.adjust(results$p.value, method = "BH")
  results
}

t2.test.excluding.controls <- function(seurat.object, controls, npcs=20) {
  results <- data.frame(ident="NA", t2.stat=NA, p.value=NA)
  idents <- data.frame(table(Idents(seurat.object))) %>% filter(Freq>=3) #Only test genes with at least 3 cells
  for (i in idents$Var1)  {
    hotelling.res <- Hotelling::hotelling.test(x = seurat.object@reductions$pca@cell.embeddings[which(Idents(seurat.object)==i),1:npcs],
                                               y = seurat.object@reductions$pca@cell.embeddings[which(Idents(seurat.object)!=i&!Idents(seurat.object)%in%controls),1:npcs])
    results.test=data.frame(ident=i, t2.stat=hotelling.res$stats$statistic, p.value=hotelling.res$pval)
    results <- rbind(results, results.test)
    print(i)
  }
  results <- filter(results, ident!="NA")
  results$p.value.bh <- p.adjust(results$p.value, method = "BH")
  results
}

