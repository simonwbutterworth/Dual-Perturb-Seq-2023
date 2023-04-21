##### Analysis of Pilot 24 sgRNA Experiment #####
rm(list = ls())
gc()
source("./FUNCTIONS.r")
gene.order <- c("UPRT", "MYR1", "HCE1", "IST", "GRA16", "GRA24", "GRA28", "ROP16", "GRA6", "GRA7", "ROP18", "GRA15",  "GRA18", "GRA25", "ROP1",  "ROP47", "GRA14", "PP2CHN", "202200", "202620", "217530", "262400", "304955", "GRA59")

##### Import data and filter cell barcodes #####

# Import data
pilot24 <- ImportFeatureBC_Mtx(matrix="../data/GSE229505_BUT290A19_BUT290A20_matrix.mtx.gz",
                             cells="../data/GSE229505_BUT290A19_BUT290A20_barcodes.tsv.gz",
                             features="../data/GSE229505_BUT290A19_BUT290A20_features.tsv.gz",
                             sample.name = "PILOT_24_SGRNA", sample.id.gex = "BUT290A19", sample.id.crispr = "BUT290A20")


# Filter cell barcodes
ggplot(mapping = aes(x=pilot24$nCount_Tg, y=pilot24$nCount_Hs))+
  geom_point(size=1, colour="grey30")+
  scale_x_log10(limits=c(10,100000), oob=squish, expand=c(0,0), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits=c(100,300000), oob=squish, expand=c(0,0), labels = trans_format("log10", math_format(10^.x)))+
  annotate(geom="rect", xmin=3000, xmax=50000, ymin=10000, ymax=200000, colour="slateblue", fill=NA, size=1.5)+
  theme_pubr(legend = "none")+
  labs(x="Tg UMI", y="Hs UMI", colour="sgRNA UMI")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=15, b = 5, l = 5))
ggplot(mapping=aes(x=pilot24$orig.ident, y=pilot24$percent.mt))+
  geom_jitter(size=1, colour="grey30")+
  scale_y_continuous(limits=c(0,100), expand=c(0,0))+
  theme_pubr(legend = "none")+
  labs(x="", y="% Mitochondrial UMIs")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=5, b = 5, l = 5))+
  annotate(geom="rect", xmin=0.5, xmax=1.5, ymin=0.5, ymax=10, colour="slateblue", fill=NA, size=1.5)
ggplot(mapping = aes(x=as.numeric(pilot24$nFeature_sgRNA)))+
  geom_histogram(binwidth = 1, colour="grey30", fill="grey50", bins = 8)+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 10, r=5, b = 5, l = 5))+
  labs(x="Unique sgRNAs per cell", y="Count")
pilot24 <- subset(pilot24, nCount_Hs>=10000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)

# Assign sgRNAs
pilot24 <- Assign_sgRNAs(object = pilot24, reference.file = "../inputs/seurat_sgrna_reference/PILOT_24_SGRNA_REFERENCE_221031.csv")
pilot24$guide.target.identity <- factor(pilot24$guide.target.identity, levels=gene.order)
Idents(pilot24) <- pilot24$guide.target.identity

#Plot number of cells per sgRNA
cells.per.gene.pilot24 <- CellsPerGene(pilot24)
cells.per.gene.pilot24$controls <- ifelse(cells.per.gene.pilot24$ident %in% c("MYR1", "HCE1", "IST", "GRA16", "GRA24", "GRA28", "ROP16", "GRA6", "GRA7", "ROP18"), yes = "Control", no = "Untested")

ggplot(cells.per.gene.pilot24, aes(x=ident, y=n, fill=controls))+
  geom_col()+
  scale_y_continuous(limits=c(0,90), expand=c(0,0))+
  theme_pubr(legend="none")+
  theme(text=element_text(size=15),
        axis.text.x = element_text(angle=45, hjust=1, size=15),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0))+
  labs(x="sgRNA", y="Number of cells")+
  scale_fill_manual(values=c("slateblue", "grey50"))
ggsave(filename = "../outputs/fig_s2/FIG_S2D_CELLS_PER_TARGET_GENE.png", width=8.2, height=4.3, dpi=300, scale=1)    

##### SCTransform, PCA, UMAP #####

# SCTransform, PCA, UMAP
pilot24 <- SCTransform(pilot24, assay="Hs")
pilot24 <- RunPCA(pilot24)
pilot24 <- RunUMAP(pilot24, dims=1:10)

# UMAP plot sg(UPRT)
ggplot()+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"]),
             size=2, colour="grey80")+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"][pilot24$guide.target.identity=="UPRT"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"][pilot24$guide.target.identity=="UPRT"]),
             size=2, colour="grey50")+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s2/FIG_S2E_UMAP_UPRT.png", width=3.7, height=3.7, dpi=300, scale=1)

# UMAP sg(MYR1)
ggplot()+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"]),
              size=1.5, colour="grey80")+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"][pilot24$guide.target.identity=="MYR1"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"][pilot24$guide.target.identity=="MYR1"]),
             size=2, colour="slateblue")+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s2/FIG_S2E_UMAP_MYR1.png", width=1.7, height=1.7, dpi=300, scale=1)

#sg(GRA16)
ggplot()+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"]),
             size=1.5, colour="grey80")+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"][pilot24$guide.target.identity=="GRA16"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"][pilot24$guide.target.identity=="GRA16"]),
             size=2, colour="slateblue")+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s2/FIG_S2E_UMAP_GRA16.png", width=1.7, height=1.7, dpi=300, scale=1)

#sg(GRA24)
ggplot()+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"]),
             size=1.5, colour="grey80")+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"][pilot24$guide.target.identity=="GRA24"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"][pilot24$guide.target.identity=="GRA24"]),
             size=2, colour="slateblue")+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s2/FIG_S2E_UMAP_GRA24.png", width=1.7, height=1.7, dpi=300, scale=1)

#sg(ROP16)
ggplot()+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"]),
             size=1.5, colour="grey80")+
  geom_point(mapping=aes(x=pilot24@reductions$umap@cell.embeddings[,"UMAP_1"][pilot24$guide.target.identity=="ROP16"], 
                         y=pilot24@reductions$umap@cell.embeddings[,"UMAP_2"][pilot24$guide.target.identity=="ROP16"]),
             size=2, colour="slateblue")+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s2/FIG_S2E_UMAP_ROP16.png", width=1.7, height=1.7, dpi=300, scale=1)

##### Find differentially expressed genes for each target #####

# Run FindMarkers for each sgRNA
markers <- data.frame(target=NA, gene=NA, avg_log2FC=NA, p_val=NA, p_val_bh=NA)
for (i in unique(Idents(pilot24))[unique(Idents(pilot24))!="UPRT"]) {
  markers.i <- FindMarkers(pilot24, ident.1 = i, ident.2 = "UPRT", assay = "SCT", logfc.threshold = 0.1, min.pct = 0)
  markers.i$p_val_bh <- p.adjust(markers.i$p_val, method="BH")
  markers.i <- rownames_to_column(markers.i, var="gene")
  markers.i <- data.frame(target=i, gene=markers.i$gene, avg_log2FC=markers.i$avg_log2FC, p_val=markers.i$p_val, p_val_bh=markers.i$p_val_bh)
  markers.i <- filter(markers.i, p_val_bh<=0.05)
  markers <- rbind(markers, markers.i)
  print(i)  
}
markers <- filter(markers, target!="NA")
write.csv(markers, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_2_24_GUIDE_DE_GENES.csv")
markers <- read.csv(file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_2_24_GUIDE_DE_GENES.csv")
# Plot number of DE genes for each sgRNA
markers.n <- markers %>%
  group_by(target) %>%
  dplyr::summarise(n=n())
markers.zero <- data.frame(target=gene.order[!gene.order%in%markers.n$target&gene.order!="UPRT"], n=0)
markers.n <- rbind(markers.n, markers.zero) %>%
  mutate(ident=factor(target, levels=gene.order)) %>%
  mutate(control=ifelse(target%in%c("MYR1", "HCE1", "IST", "GRA16", "GRA24", "GRA28", "ROP16"), yes = "Control", no = "Untested"))

ggplot(markers.n, aes(x=ident, y=n, fill=control))+
  geom_col(colour=NA)+
  theme_pubr(legend="none")+
  scale_y_continuous(trans=pseudo_log_trans(), breaks=c(1,10,100,1000, 10000), expand=c(0,0), limits=c(0,9000))+
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        text=element_text(size=15),
        plot.margin = margin(t = 0, b = 2, r = 0, l = 0))+
  labs(x="sgRNA", y="Number of DE genes")+
  scale_fill_manual(values=c("slateblue", "grey50"))
ggsave(filename = "../outputs/fig_1/FIG_1L_24_GUIDE_DE_GENES.png", width=8.2, height=4, dpi=300, scale=1)    

##### VISION Signature Scores #####

# Extract PCA embeddings and metadata for VISION
pca.embeddings <- as.data.frame(Embeddings(pilot24, reduction = "pca"))
pilot24$guide.target.identity <- factor(pilot24$guide.target.identity)
metadata <- pilot24[[]]

# Run Vision
pilot.vis.controls <- Vision(pilot24@assays$SCT@data[,], 
                         latentSpace=pca.embeddings, 
                         meta=metadata, 
                         projection_methods=NULL, 
                         signatures=c("../inputs/gene_sets/GRA16_MARKERS_BOUGDOUR_2013.txt",
                                      "../inputs/gene_sets/GRA24_MARKERS_BRAUN_2013.txt",
                                      "../inputs/gene_sets/GRA28_MARKERS_RUDZKI_2021.txt",
                                      "../inputs/gene_sets/HCE1_MARKERS_PANAS_2019.txt",
                                      "../inputs/gene_sets/IST_MARKERS_MATTA_2019.txt",
                                      "../inputs/gene_sets/MYR1_MARKERS_NAOR_2018.txt",
                                      "../inputs/gene_sets/ROP16_MARKERS_SAEIJ_2007.txt")) %>% calcSignatureScores()

# Test and plot MYR1 signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"MYR1_MARKERS_NAOR"][pilot.vis.controls@metaData$guide.target.identity=="MYR1"],
            y = pilot.vis.controls@SigScores[,"MYR1_MARKERS_NAOR"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "MYR1")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "MYR1")],
                   y=pilot.vis.controls@SigScores[,"MYR1_MARKERS_NAOR"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "MYR1")]))+
  geom_violin(adjust=1.5, scale = "width", draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"MYR1_MARKERS_NAOR"]), max(pilot.vis.controls@SigScores[,"MYR1_MARKERS_NAOR"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(MYR1)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_MYR1.png", width=2.2, height=4, dpi=300, scale=1)    

# Test and plot HCE1 signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"HCE1_MARKERS_PANAS"][pilot.vis.controls@metaData$guide.target.identity=="HCE1"],
            y = pilot.vis.controls@SigScores[,"HCE1_MARKERS_PANAS"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "HCE1")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "HCE1")],
                   y=pilot.vis.controls@SigScores[,"HCE1_MARKERS_PANAS"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "HCE1")]))+
  geom_violin(adjust=1.5, scale = "width", draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"HCE1_MARKERS_PANAS"]), max(pilot.vis.controls@SigScores[,"HCE1_MARKERS_PANAS"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(HCE1)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_HCE1.png", width=2.2, height=4, dpi=300, scale=1)    

# Test and plot IST signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"IST_MARKERS_MATTA"][pilot.vis.controls@metaData$guide.target.identity=="IST"],
            y = pilot.vis.controls@SigScores[,"IST_MARKERS_MATTA"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "IST")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "IST")],
                   y=pilot.vis.controls@SigScores[,"IST_MARKERS_MATTA"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "IST")]))+
  geom_violin(adjust=1.5, scale = "width", draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"IST_MARKERS_MATTA"]), max(pilot.vis.controls@SigScores[,"IST_MARKERS_MATTA"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(IST)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_IST.png", width=2.2, height=4, dpi=300, scale=1)    

# Test and plot GRA16 signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"GRA16_MARKERS_BOUGDOUR"][pilot.vis.controls@metaData$guide.target.identity=="GRA16"],
            y = pilot.vis.controls@SigScores[,"GRA16_MARKERS_BOUGDOUR"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA16")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA16")],
                   y=pilot.vis.controls@SigScores[,"GRA16_MARKERS_BOUGDOUR"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA16")]))+
  geom_violin(adjust=1.5, scale = "width", draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"GRA16_MARKERS_BOUGDOUR"]), max(pilot.vis.controls@SigScores[,"GRA16_MARKERS_BOUGDOUR"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(GRA16)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_GRA16.png", width=2.2, height=4, dpi=300, scale=1)    

# Test and plot GRA24 signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"GRA24_MARKERS_BRAUN"][pilot.vis.controls@metaData$guide.target.identity=="GRA24"],
           y = pilot.vis.controls@SigScores[,"GRA24_MARKERS_BRAUN"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA24")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA24")],
                   y=pilot.vis.controls@SigScores[,"GRA24_MARKERS_BRAUN"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA24")]))+
  geom_violin(adjust=1.5, scale = "width", draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"GRA24_MARKERS_BRAUN"]), max(pilot.vis.controls@SigScores[,"GRA24_MARKERS_BRAUN"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(GRA24)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_GRA24.png", width=2.2, height=4, dpi=300, scale=1)    

# Test and plot GRA28 signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"GRA28_MARKERS_RUDZKI"][pilot.vis.controls@metaData$guide.target.identity=="GRA28"],
            y = pilot.vis.controls@SigScores[,"GRA28_MARKERS_RUDZKI"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA28")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA28")],
                   y=pilot.vis.controls@SigScores[,"GRA28_MARKERS_RUDZKI"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "GRA28")]))+
  geom_violin(adjust=1.5, scale = "width",  draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"GRA28_MARKERS_RUDZKI"]), max(pilot.vis.controls@SigScores[,"GRA28_MARKERS_RUDZKI"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(GRA28)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_GRA28.png", width=2.2, height=4, dpi=300, scale=1)    

# Test and plot ROP16 signature scores
wilcox.test(x = pilot.vis.controls@SigScores[,"ROP16_MARKERS_SAEIJ"][pilot.vis.controls@metaData$guide.target.identity=="ROP16"],
            y = pilot.vis.controls@SigScores[,"ROP16_MARKERS_SAEIJ"][pilot.vis.controls@metaData$guide.target.identity=="UPRT"])
ggplot(mapping=aes(x=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "ROP16")],
                   fill=pilot.vis.controls@metaData$guide.target.identity[pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "ROP16")],
                   y=pilot.vis.controls@SigScores[,"ROP16_MARKERS_SAEIJ"][pilot.vis.controls@metaData$guide.target.identity%in%c("UPRT", "ROP16")]))+
  geom_violin(adjust=1.5, scale = "width", draw_quantiles = 0.5, size=0.5)+
  theme_pubr(legend = "none")+
  labs(x=NULL, y=NULL)+
  scale_fill_manual(values = c("grey50", "slateblue"))+
  theme(text = element_text(size=14),
        plot.margin = margin(t = 0, b = 5, r = 0, l = 5),
        axis.ticks.y = element_blank(), axis.text.y=element_blank(), axis.ticks.length.y = unit(0,"pt"))+
  scale_y_continuous(limits=c(min(pilot.vis.controls@SigScores[,"ROP16_MARKERS_SAEIJ"]), max(pilot.vis.controls@SigScores[,"ROP16_MARKERS_SAEIJ"])),
                     expand = expansion(mult = c(0, 0.2)))+
  scale_x_discrete(labels=c("sg(UPRT)", "sg(ROP16)"))
ggsave(filename = "../outputs/fig_1/FIG_1J_SIGNATURE_SCORE_ROP16.png", width=2.2, height=4, dpi=300, scale=1)    
