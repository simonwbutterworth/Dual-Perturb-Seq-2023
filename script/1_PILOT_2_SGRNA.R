##### Analysis of Pilot 2 sgRNA Experiment #####
rm(list = ls())
gc()
source("./FUNCTIONS.R")

##### Import data and filter cell barcodes #####

# Import data
pilot <- ImportFeatureBC_Mtx(matrix="../data/GSE229505_BUT290A1_BUT290A2_matrix.mtx.gz",
                             cells="../data/GSE229505_BUT290A1_BUT290A2_barcodes.tsv.gz",
                             features="../data/GSE229505_BUT290A1_BUT290A2_features.tsv.gz",
                             sample.name = "PILOT_2_SGRNA", sample.id.gex = "BUT290A1", sample.id.crispr = "BUT290A2")

# Assign sgRNAs
pilot <- Assign_sgRNAs(object = pilot, reference.file = "../inputs/seurat_sgrna_reference/PILOT_2_SGRNA_REFERENCE_221025.csv")

# Plot nCount_Hs and nCount_Tg
ggplot(mapping = aes(x=pilot$nCount_Tg, y=pilot$nCount_Hs))+
  geom_point(colour="grey50")+
  scale_x_log10(limits=c(10,100000), oob=squish, expand=c(0,0), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits=c(100,300000), oob=squish, expand=c(0,0), labels = trans_format("log10", math_format(10^.x)))+
  annotate(geom="rect", xmin=3000, xmax=50000, ymin=10000, ymax=200000, colour="grey30", fill=NA, size=1)+
  theme_pubr(legend = "none")+
  labs(x="Tg UMI", y="Hs UMI", colour="sgRNA UMI")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=15, b = 5, l = 5))
ggsave(filename = "../outputs/fig_s2/FIG_S2A_NCOUNT.png", width=4, height=4, dpi=300, scale=1)    

# Plot percent MT UMIs
ggplot(mapping=aes(x=pilot$orig.ident, y=pilot$percent.mt))+
  geom_jitter(colour="grey50")+
  scale_y_continuous(limits=c(0,100), expand=c(0,0))+
  theme_pubr(legend = "none")+
  labs(x="", y="% Mitochondrial UMIs")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=5, b = 5, l = 5))+
  annotate(geom="rect", xmin=0.5, xmax=1.5, ymin=0.5, ymax=10, colour="grey30", fill=NA, size=1)
ggsave(filename = "../outputs/fig_s2/FIG_S2A_PERCENT_MT.png", width=4, height=4, dpi=300, scale=1)    

# Filter by nCount and percent MT
pilot <- subset(pilot, subset=nCount_Hs>=10000&nCount_Tg>=3000&percent.mt<=10)

# Plot nFeature_sgRNA
ggplot(mapping = aes(x=pilot$nFeature_sgRNA))+
  geom_histogram(binwidth = 1, colour="grey50", fill="grey50")+
  theme_pubr()+
  scale_y_continuous(limits=c(0,650), expand=c(0,0))+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=5, b = 5, l = 5))+
  labs(x="Unique sgRNAs per cell", y="Count")+
  annotate(geom="rect", xmin=0.4, xmax=1.6, ymin=5, ymax=620, colour="grey30", fill=NA, size=1)
ggsave(filename = "../outputs/fig_s2/FIG_S2A_NFEATURE_SGRNA.png", width=4, height=4, dpi=300, scale=1)    

# Filter by nFeature_sgRNA
pilot <- subset(pilot, nFeature_sgRNA==1)

# Plot cells per target gene
cells.per.gene.pilot <- CellsPerGene(pilot)
cells.per.gene.pilot$ident <- factor(cells.per.gene.pilot$ident, levels=c("UPRT", "MYR1"))
ggplot(cells.per.gene.pilot, aes(x=ident, y=n, fill=ident))+
  geom_col()+
  scale_y_continuous(expand=c(0,0), limit=c(0,400))+
  theme_pubr(legend="none")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=0, b = 5, l = 5))+
  labs(x="", y="Cells per target gene")+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(MYR1)"))+
  scale_fill_manual(values=c("grey50", "slateblue"))
ggsave(filename = "../outputs/fig_s2/FIG_S2A_CELLS_PER_TARGET_GENE.png", width=4, height=4, dpi=300, scale=1)    

##### Hs UMAP plot & DE genes #####

# SCTransform, PCA, UMAP
pilot <- SCTransform(pilot, vst.flavor="V2", assay="Hs")
pilot <- RunPCA(pilot, assay="SCT")
pilot <- RunUMAP(pilot, dims=1:10)
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], 
                   y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], 
                   colour=pilot$guide.target.identity))+
       geom_point(size=2)+
       theme_pubr(border = T, legend="none")+
       theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
             plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
       scale_colour_manual(values = c("slateblue", "grey70"))+
       labs(x=NULL, y=NULL, title=NULL)+
       scale_x_reverse()+scale_y_reverse()
ggsave(filename = "../outputs/fig_1/FIG_1G_HS_UMAP_SGRNA.png", width=3.7, height=3.7, dpi=300, scale=1)

# Find DE genes for MYR1 vs UPRT
Idents(pilot) <- pilot$guide.target.identity
markers <- FindMarkers(pilot, ident.1 = "MYR1", ident.2 = "UPRT", assay = "SCT", logfc.threshold = -Inf, min.pct = 0)
markers$target <- "MYR1"
markers$p_val_bh <- p.adjust(markers$p_val, method="BH")
markers <- rownames_to_column(markers, var="gene")
head(markers %>% filter(avg_log2FC>0) %>% arrange(p_val_bh))
head(markers %>% filter(avg_log2FC<0)  %>% arrange(p_val_bh))
nrow(markers %>% filter(p_val_bh<=0.01))
write.csv(markers, file = "../outputs/supplementary_data/SUPPLEMENTARY_DATA_1A_HS_DE_GENES_MYR1_VS_UPRT.csv")

# Volcano plot of DE genes with top five up- and down-regulated genes labelled
ggplot(markers, aes(x=avg_log2FC, y=-log10(p_val_bh)))+
  geom_point(alpha=1, colour="grey50")+
  geom_text_repel(data = filter(markers, gene%in%c("BIRC3", "PTMA", "HSP90AA1", "PHLDA2", "ARNTL2",  "ANXA1", "FTL", "THBS1", "PDLIM2", "REXO2")),  aes(label=gene), min.segment.length = 0, size=5)+
  scale_x_continuous(limits=c(-1.5,1.5), expand = c(0,0))+
  scale_y_continuous(limits=c(0,95), expand=c(0,0))+
  theme_pubr()+
  theme(text=element_text(size=15), plot.margin = margin(t = 0, r = 10, b = 5.5, l = 5))+
  labs(x="Average Log2FC", y="-Log10(adjusted p-value)")
ggsave(filename = "../outputs/fig_s2/FIG_S2B_HS_DE_GENES_PLOT.png", width=8.2, height=3.7, dpi=300, scale=1)

# UMAP plots for top two up- and down-regulated genes: BIRC3, PTMA, FTL, ANXA1
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], colour=pilot@assays$SCT@data["BIRC3",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_x_reverse()+scale_y_reverse()
ggsave(filename = "../outputs/fig_1/FIG_1H_HS_UMAP_BIRC3.png", width=1.7, height=1.7, dpi=300, scale=1)

ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"],y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"],colour=pilot@assays$SCT@data["PTMA",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_x_reverse()+scale_y_reverse()
ggsave(filename = "../outputs/fig_1/FIG_1H_HS_UMAP_PTMA.png", width=1.7, height=1.7, dpi=300, scale=1)

ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"],y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"],colour=pilot@assays$SCT@data["FTL",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(), plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_x_reverse()+scale_y_reverse()
ggsave(filename = "../outputs/fig_1/FIG_1H_HS_UMAP_FTL.png", width=1.7, height=1.7, dpi=300, scale=1)

ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], colour=pilot@assays$SCT@data["ANXA1",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(), plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_x_reverse()+scale_y_reverse()
ggsave(filename = "../outputs/fig_1/Fig_1H_HS_UMAP_ANXA1.png", width=1.7, height=1.7, dpi=300, scale=1)

##### Tg UMAP plots & cell cycle markers #####

# SCTransform, PCA, UMAP
pilot <- SCTransform(pilot, vst.flavor="v2", assay="Tg")
pilot <- RunPCA(pilot)
pilot <- RunUMAP(pilot, dims=1:10)
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], colour=pilot$guide.target.identity))+
  geom_point(size=2)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_manual(values = c("slateblue", "grey70"))+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_1/FIG_1I_TG_UMAP_SGRNA.png", width=3.7, height=3.7, dpi=300, scale=1)

# Toxo cell cycle markers: 
# GAPDH2 G1a Phase TGME49-269190
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], colour=pilot@assays$SCT@data["TGME49-269190",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_1/FIG_1J_TG_UMAP_GAPDH2_G1A.png", width=1.7, height=1.7, dpi=300, scale=1)

#ACP G1b Phase TGME49-264080
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], colour=pilot@assays$SCT@data["TGME49-264080",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_1/FIG_1J_TG_UMAP_ACP_G1B.png", width=1.7, height=1.7, dpi=300, scale=1)

#ROP1 S Phase TGME49-309590
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"], colour=pilot@assays$SCT@data["TGME49-309590",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_1/FIG_1J_TG_UMAP_ROP1_S.png", width=1.7, height=1.7, dpi=300, scale=1)

#MIC1 M Phase TGME49-291890
ggplot(mapping=aes(x=pilot@reductions$umap@cell.embeddings[,"UMAP_1"], y=pilot@reductions$umap@cell.embeddings[,"UMAP_2"],colour=pilot@assays$SCT@data["TGME49-291890",]))+
  geom_point(size=1.5)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  scale_colour_gradient(low="grey80", high="darkorchid")+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_1/FIG_1J_TG_UMAP_MIC1_M.png", width=1.7, height=1.7, dpi=300, scale=1)

# Find DE Tg genes for MYR1 vs UPRT
markers.tg <- FindMarkers(pilot, ident.1 = "MYR1", ident.2 = "UPRT", assay = "SCT", logfc.threshold = -Inf, min.pct = 0)
markers.tg$target <- "MYR1"
markers.tg$p_val_bh <- p.adjust(markers.tg$p_val, method="BH")
markers.tg <- rownames_to_column(markers.tg, var="gene")
head(markers.tg %>% filter(avg_log2FC>0) %>% arrange(p_val_bh))
head(markers.tg %>% filter(avg_log2FC<0)  %>% arrange(p_val_bh))
write.csv(markers.tg, file = "../outputs/supplementary_data/SUPPLEMENTARY_DATA_1B_TG_DE_GENES_MYR1_VS_UPRT.csv")

ggplot(markers.tg, aes(x=avg_log2FC, y=-log10(p_val_bh)))+
  geom_point(alpha=1, colour="grey50")+
  geom_text_repel(data = filter(markers.tg, p_val_bh<=0.01),  aes(label=c("HSP90","HSP70","RPL31","UPRT")), min.segment.length = 0, size=5)+
  scale_x_continuous(limits=c(-1.5,1.5), expand = c(0,0))+
  scale_y_continuous(limits=c(0,95), expand=c(0,0))+
  theme_pubr()+
  theme(text=element_text(size=15), plot.margin = margin(t = 0, r = 10, b = 5.5, l = 5))+
  labs(x="Average Log2FC", y="-Log10(adjusted p-value)")
ggsave(filename = "../outputs/fig_s2/FIG_S2C_TG_DE_GENES_PLOT.png", width=8.2, height=3.7, dpi=300, scale=1)
