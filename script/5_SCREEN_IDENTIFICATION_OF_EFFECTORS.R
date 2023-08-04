##### Identification of effector proteins in dual perturb-seq screen #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
controls <- c("MYR1", "MYR2", "MYR3", "MYR4", "ROP17", "GRA45", "GRA16", "GRA24", "GRA28", "IST", "HCE1", "ROP16", "PPM3C", "MAF1B")

##### Illustration of perturbation score test #####

# Main panel
ggplot()+
  geom_point(aes(x=screen@reductions$pca@cell.embeddings[,"PC_1"],
                 y=screen@reductions$pca@cell.embeddings[,"PC_2"]),
             size=0.5, alpha=1, colour="grey80")+
  geom_point(aes(x=screen@reductions$pca@cell.embeddings[,"PC_1"][screen$guide.target.identity=="MYR1"],
                 y=screen@reductions$pca@cell.embeddings[,"PC_2"][screen$guide.target.identity=="MYR1"]),
             size=0.5, alpha=1, colour="slateblue")+
  theme_pubr(border=T)+
  theme(axis.text = element_blank(), axis.ticks.length = unit(0, "pt"), axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0))+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(expand=expansion(mult = c(0,0)))+
  scale_x_continuous(expand=expansion(mult = c(0,0)))
ggsave(file="../outputs/fig_s3/FIG_S3J_PCA_PLOT.png", width = 3.2, height = 3.2, scale = 1, dpi = 300)

# PC1 distribution
ggplot()+
  geom_density(aes(x=screen@reductions$pca@cell.embeddings[,"PC_1"]),
               size=1, colour="grey70", fill="grey70", alpha=0.5)+
  geom_density(aes(x=screen@reductions$pca@cell.embeddings[,"PC_1"][screen$guide.target.identity=="MYR1"]),
               size=1, colour="slateblue", fill="slateblue", alpha=0.5)+
  theme_pubr(border=T)+
  theme(axis.text = element_blank(), axis.ticks.length = unit(0, "pt"), axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0))+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(expand=expansion(mult = c(0,0.05)))+
  scale_x_continuous(expand=expansion(mult = c(0,0)))
ggsave(file="../outputs/fig_s3/FIG_S3J_PC1_DISTRIBUTION.png", width = 3.2, height = 0.5, scale = 1, dpi = 300)

# PC3 distribution
ggplot()+
  geom_density(aes(x=screen@reductions$pca@cell.embeddings[,"PC_2"]),
               size=1, colour="grey70", fill="grey70", alpha=0.5)+
  geom_density(aes(x=screen@reductions$pca@cell.embeddings[,"PC_2"][screen$guide.target.identity=="MYR1"]),
               size=1, colour="slateblue", fill="slateblue", alpha=0.5)+
  theme_pubr(border=T)+
  theme(axis.text = element_blank(), axis.ticks.length = unit(0, "pt"), axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0))+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(expand=expansion(mult = c(0,0.05)))+
  scale_x_continuous(expand=expansion(mult = c(0,0)))
ggsave(file="../outputs/fig_s3/FIG_S3J_PC2_DISTRIBUTION.png", width = 3.2, height = 0.5, scale = 1, dpi = 300)

##### Perturbation scores - all data combined #####

# Calculate perturbation scores
t2.results <- t2.test.excluding.controls(screen, controls = controls, npcs=20) %>%
  mutate(significance=ifelse(p.value.bh<=0.01, yes="Significant", no="Non-significant")) %>%
  mutate(significance=ifelse(ident%in%controls, yes="Control effector", no=significance)) %>%
  mutate(significance=ifelse(ident%in%c("222100", "GRA59"), yes="New effector", no=significance)) %>%
  mutate(significance=factor(significance, levels=c("Control effector", "New effector", "Significant", "Non-significant"))) %>%
  left_join(CellsPerGene(screen), by="ident") %>%
  rename(n.cells=n, target=ident) %>%
  relocate(significance, .after=n.cells)
write.csv(t2.results, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_4A_PERTURBATION_SCORES_COMBINED.csv", row.names = F)

# Plot perturbation scores
ggplot(t2.results, aes(x=t2.stat, y=pmin(-log10(p.value.bh),15), label=target, colour=significance))+
  geom_hline(yintercept = 2, colour="grey80", size=0.5, linetype="dashed")+
  geom_point(aes(size=pmax(pmin(log10(n.cells),3),1)), alpha=0.75)+
  scale_x_log10(limits=c(8,10000), oob=squish)+
  scale_y_continuous(limits=c(0,15.5), labels=c(0, 5, 10, ">15"), expand=expansion(add=c(1,0)))+
  scale_colour_manual(values=c("slateblue", "goldenrod", "grey50", "grey80"))+
  scale_radius(limits=c(1,3), breaks=c(1,2,3), labels = c(10,100,1000), range=c(2,10))+
  theme_pubr(border=T, legend="none")+
  labs(x="Perturbation score", y="-Log10(adjusted p-value)")+
  theme(text=element_text(size=15), 
        plot.margin = margin(t = 20,10,0,0), line=element_line(size = 1), panel.border = element_rect(size=1, fill = NA))
ggsave("../outputs/fig_2/FIG_2B_PERTURBATION_SCORES.png", width = 5.4, height = 5.9, scale = 1, dpi = 300)

#Identify hits
hits <- filter(t2.results, p.value.bh<=0.01)

##### Perturbation scores - unstimulated vs IFNG (shared PCA space) ##### 

# Perturbation scores - unstimulated data
t2.results.unstim <- t2.test.excluding.controls(subset(screen, ifng=="Unstimulated"), controls = controls, npcs=20) %>%
  left_join(CellsPerGene(subset(screen, ifng=="Unstimulated")), by="ident") %>%
  rename(n.cells=n, target=ident)   
write.csv(t2.results.unstim, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_4B_PERTURBATION_SCORES_UNSTIMULATED.csv", row.names = F)

# Perturbation scores - IFNG-stimulated data
t2.results.ifng   <- t2.test.excluding.controls(subset(screen, ifng=="IFNG"), controls = controls, npcs=20) %>%
  left_join(CellsPerGene(subset(screen, ifng=="IFNG")), by="ident") %>%
  rename(n.cells=n, target=ident) 
write.csv(t2.results.unstim, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_4C_PERTURBATION_SCORES_IFNG.csv", row.names=F)

# Join perturbation scores from unstimulated and IFNG
t2.results.unstim.ifng <- full_join(t2.results.unstim, t2.results.ifng, by="target", suffix=c(".unstim", ".ifng"))
t2.results.unstim.ifng <- left_join(t2.results.unstim.ifng, rename(CellsPerGene(screen), target=ident), by="target")
t2.results.unstim.ifng$controls <- ifelse(t2.results.unstim.ifng$target%in%hits$target, yes=2, no=3)
t2.results.unstim.ifng$controls <- ifelse(t2.results.unstim.ifng$target%in%controls, yes=1, no=t2.results.unstim.ifng$controls)
t2.results.unstim.ifng <- arrange(t2.results.unstim.ifng, desc(controls))

# Plot perturbation score comparison
ggplot(t2.results.unstim.ifng, aes(x=t2.stat.unstim, y=t2.stat.ifng, label=target, colour=factor(controls)))+
  geom_abline(intercept = 0, slope = 1,  colour="grey80", size=0.5)+
  geom_point(aes(size=pmax(pmin(log10(n),3),1)), alpha=0.75)+
  geom_text_repel(data=filter(t2.results.unstim.ifng, target%in%hits$target),
                  max.overlaps = Inf, size=5, box.padding = 0.5)+
  scale_x_log10(limits=c(10,1000))+scale_y_log10(limits=c(10,1000))+
  scale_colour_manual(values=c("slateblue", "grey50", "grey80"))+
  theme_pubr(border=T, legend="none", base_size = 15)+
  labs(x="Perturbation score (t\u00b2 statistic): Unstimulated", y="Perturbation score (t\u00b2 statistic): IFN\u03b3")+
  theme(plot.margin = margin(t = 0, r = 0, b = 5,l = 5), line=element_line(size = 1), panel.border = element_rect(size=1, fill = NA))
ggsave("../outputs/fig_s3/FIG_S3K_UNSTIMULATED_IFNG_PERTURBATION_SCORE_COMPARISON.png", width = 5.9, height = 5.9, scale = 1, dpi = 300)

# Plot adjusted p-value comaparison
ggplot(t2.results.unstim.ifng, aes(x=pmax(log10(p.value.bh.unstim),-15), y=pmax(log10(p.value.bh.ifng), -15),label=target, colour=factor(controls)))+
  geom_hline(yintercept = -2, colour="grey80", size=0.5, linetype="dashed")+
  geom_vline(xintercept = -2, colour="grey80", size=0.5, linetype="dashed")+
  geom_point(aes(size=pmax(pmin(log10(n),3),1)), alpha=0.75)+
  geom_text_repel(data=filter(t2.results.unstim.ifng, target%in%hits$target), max.overlaps = Inf,  size=5, , box.padding = 0.5)+
  scale_colour_manual(values=c("slateblue", "grey50", "grey80"))+
  theme_pubr(border=T, legend="none", , base_size = 15)+
  scale_x_reverse(labels=c("> -15", -10, -5, 0))+
  scale_y_reverse(labels=c("> -15", -10, -5, 0))+
  labs(x="Log10(adjusted p-value): Unstimulated", y="Log10(adjusted p-value): IFN\u03b3")+
  theme(plot.margin = margin(t = 0, r = 0, b = 5,l = 5), line=element_line(size = 1), panel.border = element_rect(size=1, fill = NA))
ggsave("../outputs/fig_s3/FIG_S3K_UNSTIMULATED_IFNG_PVALUE_COMPARISON.png", width = 5.9, height = 5.9, scale = 1, dpi = 300)

##### Pseudobulk PCA #####

# Identify targets with at least 30 cells, corresponding to 25th percentile
quantile(CellsPerGene(screen)$n)
cells.per.gene.30 <- filter(CellsPerGene(screen), n>=30)

# Average SCT scale.data for each target
screen.average <- AverageExpression(screen, assays="SCT", slot="scale.data", group.by = "ident", return.seurat = T)
screen.average@assays$SCT@data <- screen.average@assays$SCT@scale.data
screen.average$orig.ident <- names(screen.average$orig.ident)

# Filter to only retain targets with at least 30 cells
screen.average <- subset(screen.average, subset=orig.ident%in%cells.per.gene.30$ident)

# Run PCA
screen.average <- FindVariableFeatures(screen.average, assay = "SCT")
screen.average <- RunPCA(screen.average)

#Extract embeddings & join to perturbation score and cells per gene data
screen.average.pca <- data.frame(Embeddings(screen.average, reduction = "pca")) %>%
  rownames_to_column(var = "target") %>%
  left_join(t2.results, by="target") %>%
  arrange(desc(significance))

# Calculate percentage of total variance explained by each PC
mat <- GetAssayData(screen.average, assay = "SCT", slot = "scale.data")
total.variance <- sum(matrixStats::rowVars(mat))
eig.values = (screen.average@reductions$pca@stdev)^2
variance.explained = eig.values/total.variance*100
variance.explained[1:4]

# Plot pseudobulk PC1 and PC2 embeddings
ggplot(screen.average.pca, aes(x=PC_1, y=PC_2))+
  geom_point(aes(size=pmax(pmin(log10(n.cells),3),1), colour=significance), alpha=0.75)+
  geom_text_repel(data=filter(screen.average.pca, significance!="Non-significant"), aes(label=target, colour=significance), size=5, box.padding = 0.75, min.segment.length = 1, max.overlaps = 15)+
  scale_colour_manual(values=c("slateblue", "goldenrod", "grey50", "grey80"))+
  scale_radius(limits=c(1,3), breaks=c(1,2,3), labels = c(10,100,1000), range=c(2,10))+
  theme_pubr(border=T, legend="none")+
  labs(x=NULL, y=NULL)+
  theme(text=element_text(size=15), 
        axis.text = element_blank(), axis.ticks.length = unit(0, "pt"), axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0))
ggsave("../outputs/fig_2/FIG_2C_PSEUDOBULK_PC1_PC2.png", width = 5.1, height = 5.1, scale = 1, dpi = 300)

# Plot pseudobulk PC3 and PC4 embeddings
ggplot(screen.average.pca, aes(x=PC_3, y=PC_4))+
  geom_point(aes(size=pmax(pmin(log10(n.cells),3),1), colour=significance), alpha=0.75)+
  geom_text_repel(data=filter(screen.average.pca, significance!="Non-significant"), aes(label=target, colour=significance), size=5, box.padding = 0.75, min.segment.length = 1, max.overlaps = 15)+
  scale_colour_manual(values=c("slateblue", "goldenrod", "grey50", "grey80"))+
  scale_radius(limits=c(1,3), breaks=c(1,2,3), labels = c(10,100,1000), range=c(2,10))+
  theme_pubr(border=T, legend="none")+
  labs(x=NULL, y=NULL)+
  theme(text=element_text(size=15), 
        axis.text = element_blank(), axis.ticks.length = unit(0, "pt"), axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0))
ggsave("../outputs/fig_2/FIG_2D_PSEUDOBULK_PC3_PC4.png", width = 5.1, height = 5.1, scale = 1, dpi = 300)

##### Calculate correlation between pseudobulk transcriptomes and cluster hits ##### 

# Get SCT average scale.data as a matrix
screen.average.scale.data <- GetAssayData(screen.average, assay = "SCT", slot = "data")

# Calulate pearson correlation coefficients
screen.average.pearson <- cor(x = screen.average.scale.data, method = "pearson")

# Filter to retain only Hotelling significant hits
screen.average.pearson <- screen.average.pearson[hits$target, hits$target]

# Plot heatmap - from RStudio, export as PDF 10inx10in for a plot ~80mmx80mm at print scale
slateblue.goldenrod <- colorRampPalette(colors = c("slateblue", "grey90", "goldenrod"))
dev.off()
gplots::heatmap.2(screen.average.pearson[hits$target, hits$target],col=slateblue.goldenrod(99), scale="none", revC=T,
                  trace = "none", density.info = "none", lhei=c(1, 10), lwid=c(1,10), cexRow=2, cexCol=2, margins = c(8,8))
dev.off()