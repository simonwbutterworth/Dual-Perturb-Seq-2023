##### Identification of differentially expressed PID genes sets in dual perturb-seq screen #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
controls <- c("MYR1", "MYR2", "MYR3", "MYR4", "ROP17", "GRA45", "GRA16", "GRA24", "GRA28", "IST", "HCE1", "ROP16", "PPM3C", "MAF1B")

##### Vision signature scoring for PID gene sets #####

# Get PCA embeddings & metadata
pca.embeddings <- as.data.frame(Embeddings(screen, reduction="pca"))
screen$guide.target.identity <- factor(screen$guide.target.identity)
metadata <- screen[[]]

# Run Vision with Reactome gene sets
screen.vision <- Vision(data=screen@assays$SCT@data[,], 
                        latentSpace=pca.embeddings, 
                        meta=metadata, 
                        projection_methods=NULL, 
                        signatures=c("../inputs/gene_sets/c2.cp.pid.v7.4.symbols.gmt"))
screen.vision <- calcSignatureScores(object=screen.vision)

# Run Wilcoxon rank-sum test, 1 vs. all (excluding known effectors), for every target with at least 30 cells
cells.per.gene.30 <- CellsPerGene(screen) %>% filter(n>=30) %>% pull(ident)
gene.sets.wilcox <- data.frame(target=NA, gene.set=NA, p.value=NA, auc=NA)
gene.clock <- 1
for (i in cells.per.gene.30) {
  for (j in  colnames(screen.vision@SigScores)) {
    print(paste(i, "(", gene.clock, "of", length(unique(screen.vision@metaData$guide.target.identity)), ") - ", j))
    wilcox <- wilcox.test(x = screen.vision@SigScores[,j][screen.vision@metaData$guide.target.identity==i],
                          y = screen.vision@SigScores[,j][screen.vision@metaData$guide.target.identity!=i&!screen.vision@metaData$guide.target.identity%in%controls],
                          alternative = "two.sided", paired = F, exact = F)
    auc <- wilcox$statistic/(length(screen.vision@SigScores[,j][screen.vision@metaData$guide.target.identity==i])*length(screen.vision@SigScores[,j][screen.vision@metaData$guide.target.identity!=i&!screen.vision@metaData$guide.target.identity%in%controls]))
    gene.sets.wilcox.j <- data.frame(target=i,
                                     gene.set=j,
                                     p.value=wilcox$p.value,
                                     auc=unname(auc))
    gene.sets.wilcox <- rbind(gene.sets.wilcox, gene.sets.wilcox.j)
  }
  gene.clock <- gene.clock+1
}
gene.sets.wilcox <- gene.sets.wilcox %>%
  filter(!is.na(target)) %>%
  mutate(p.value.bh=p.adjust(p.value, method="BH")) %>%
  relocate(p.value, .before=p.value.bh)

# Get mean signature score for each ident and Z-transform
screen.vision.average <- data.frame(screen.vision@SigScores[,])
screen.vision.average$target <- screen.vision@metaData$guide.target.identity
screen.vision.average <- screen.vision.average %>%
  pivot_longer(cols = !target, names_to = "gene.set", values_to = "score") %>%
  filter(target%in%cells.per.gene.30) %>%
  group_by(target, gene.set) %>%
  dplyr::summarise(mean.score=mean(score), .groups = "drop") %>%
  group_by(gene.set) %>%
  mutate_at(vars(mean.score), list(mean.score.z=scale)) %>%
  ungroup()

# Join Wilcox test results and mean/z-scores and save
screen.vision.average <- left_join(x=screen.vision.average, y=gene.sets.wilcox, by=join_by(target, gene.set))
write.csv(screen.vision.average, file = "../outputs/supplementary_data/SUPPLEMENTARY_DATA_8_PID_GENE_SETS.csv", row.names = F)

##### Plot gene set signature score heatmap #####

# Select gene sets significant for at least one effector 
hits <- read.csv("../outputs/supplementary_data/SUPPLEMENTARY_DATA_4A_PERTURBATION_SCORES_COMBINED.csv") %>% filter(p.value.bh<=0.01) %>% pull(target)
screen.vision.heatmap <- filter(screen.vision.average, target%in%hits)
screen.vision.heatmap <- screen.vision.heatmap %>%
  group_by(gene.set) %>%
  filter(min(p.value.bh)<=0.01)

# Pivot wider for heatmap plotting
screen.vision.heatmap <- screen.vision.heatmap %>%
  select(target, gene.set, mean.score.z) %>%
  pivot_wider(names_from = gene.set, values_from = mean.score.z) %>%
  column_to_rownames(var = "target") %>%
  as.matrix() %>%
  t()

# Get number of differentially expressed gene sets
nrow(screen.vision.heatmap)

# Plot heatmap for main figure - from RStudio, export as PDF 10inx10in for a plot ~80mmx80mm at print scale
slateblue.goldenrod <- colorRampPalette(colors = c("slateblue", "grey90", "goldenrod"))
dev.off()
gplots::heatmap.2(screen.vision.heatmap, col=slateblue.goldenrod(99),  scale = "none", labRow = F,
                  trace = "none", density.info = "none", lhei=c(1, 10), lwid=c(1,10), cexRow=2, cexCol=2, margins = c(8,8))

# Plot heatmap for supplementary figure - from RStudio, export as PDF 30inx40in for a full-page figure at print scale
dev.off()
gplots::heatmap.2(screen.vision.heatmap, col=slateblue.goldenrod(99),  scale = "none",
                  trace = "none", density.info = "none", lhei=c(1, 20), lwid=c(1,10), cexRow=1.5, cexCol=1.5, margins = c(8,32))
dev.off()

