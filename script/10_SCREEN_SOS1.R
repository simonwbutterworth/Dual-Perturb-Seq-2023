##### Further analysis of hit from screen: SOS1/TGGT1_222100 #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")
#screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
controls <- c("MYR1", "MYR2", "MYR3", "MYR4", "ROP17", "GRA45", "GRA16", "GRA24", "GRA28", "IST", "HCE1", "ROP16", "PPM3C", "MAF1B")

##### Plot differentially expressed genes #####
markers.222100 <- read.csv("../outputs/supplementary_data/SUPPLEMENTARY_DATA_9_DE_GENES.csv") %>% filter(target=="222100")
ggplot(markers.222100, aes(x=avg_log2FC, y=-log10(p_val_bh), label=gene))+
  geom_point(size=2, colour="grey50")+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin=margin(c(0,0,0,0)))+
  scale_x_continuous(limits=c(-2.2,2.2))+
  scale_y_continuous(limits=c(0,45), expand=c(0,0))+
  labs(x="Average L2FC", y="-Log10(p-adj)")
ggsave(filename = "../outputs/fig_4/FIG_4A_SOS1_DE_GENES.png", width = 3.9, height = 3.9, dpi=300, scale = 1)

##### Plot differentially expressed PID gene sets #####
gene.sets.wilcox.222100 <- read.csv("../outputs/supplementary_data/SUPPLEMENTARY_DATA_8_PID_GENE_SETS.csv") %>% filter(target=="222100")
ggplot(gene.sets.wilcox.222100, aes(x=mean.score.z, y=-log10(p.value.bh), label=gene.set))+
  geom_point(size=2, colour="grey50")+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin=margin(c(0,0,0,0)))+
  scale_x_continuous(limits=c(-5,5))+
  scale_y_continuous(limits=c(0,22), expand=c(0,0))+
  labs(x="Signature Z-score", y="-Log10(p-adj)")
ggsave(filename = "../outputs/fig_4/FIG_4B_SOS1_DE_GENE_SETS.png", width = 3.9, height = 3.9, dpi=300, scale = 1)