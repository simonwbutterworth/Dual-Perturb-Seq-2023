##### Further analysis of hit from screen: SOS1/TGGT1_222100 #####

##### Setup #####
rm(list=ls())
gc()
memory.limit(size = 80000)
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
controls <- c("MYR1", "MYR2", "MYR3", "MYR4", "ROP17", "GRA45", "GRA16", "GRA24", "GRA28", "IST", "HCE1", "ROP16", "PPM3C", "MAF1B")

##### Find differentially expressed genes for SOS1 #####

# Set idents for test
screen$test <- ifelse(screen$guide.target.identity%in%controls, yes="control", no="other")
screen$test <- ifelse(screen$guide.target.identity=="222100", yes="222100", no=screen$test)
Idents(screen) <- screen$test

# Find DE genes
markers.222100 <- FindMarkers(screen, ident.1 = "222100", ident.2 = "other", assay = "SCT", slot = "data", logfc.threshold = -Inf, min.pct = 0)
markers.222100$p_val_bh <- p.adjust(markers.222100$p_val, method="BH")
markers.222100 <- rownames_to_column(markers.222100, var="gene")

# Plot DE genes
ggplot(markers.222100, aes(x=avg_log2FC, y=-log10(p_val_bh), label=gene))+
  geom_point(size=3, colour="grey30")+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin=margin(c(0,0,0,0)))+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_continuous(limits=c(0,45), expand=c(0,0))+
  labs(x="Average L2FC", y="-Log10(p-adj)")
ggsave(filename = "../outputs/fig_4/FIG_4A_SOS1_DE_GENES.png", width = 3.9, height = 3.9, dpi=300, scale = 1)

##### Plot most significantly DE PID gene sets #####

# Import gene set data and filter for 222100
gene.sets.wilcox.222100 <- read.csv("../outputs/supplementary_data/SUPPLEMENTARY_DATA_7_PID_GENE_SETS.csv") %>% filter(ident=="222100")

# Plot DE gene sets
ggplot(gene.sets.wilcox.222100, aes(x=mean.score.z, y=-log10(p.value.bh.all), label=gene.set))+
  geom_point(size=3, colour="grey30")+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin=margin(c(0,0,0,0)))+
  scale_x_continuous(limits=c(-5,5))+
  scale_y_continuous(limits=c(0,24), expand=c(0,0))+
  labs(x="Signature Z-score", y="-Log10(p-adj)")
ggsave(filename = "../outputs/fig_4/FIG_4B_SOS1_DE_GENE_SETS.png", width = 3.9, height = 3.9, dpi=300, scale = 1)
