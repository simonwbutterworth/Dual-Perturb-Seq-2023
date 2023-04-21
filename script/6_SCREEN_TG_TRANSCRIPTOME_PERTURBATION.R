##### Investigation whether any effector protein knockouts alter the Tg parasite transcriptome #####

##### Setup #####
rm(list=ls())
gc()
memory.limit(size = 80000)
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")

##### SCTransform, PCA, UMAP using Tg gene expression

# SCTransform, PCA, UMAP using Tg gene expression
screen <- SCTransform(screen, assay="Tg", vst.flavor="v2", vars.to.regress="orig.ident")
screen <- RunPCA(screen, assay="SCT")
screen <- RunUMAP(screen, dims=1:10)

# Plot unstimulated vs IFNG
data.frame(UMAP_1=screen@reductions$umap@cell.embeddings[,"UMAP_1"], 
           UMAP_2=screen@reductions$umap@cell.embeddings[,"UMAP_2"],
           ifng=screen$ifng,
           rand=runif(length(screen$ifng), 0, 1)) %>%
  arrange(rand) %>%
  ggplot()+
  geom_point(mapping=aes(x=UMAP_1, 
                         y=UMAP_2,
                         colour=ifng), size=0.1, alpha=1)+
  scale_colour_manual(values=c("royalblue4", "grey80"))+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s4/FIG_S4A_TG_UMAP_IFNG.png", width=5.1, height=5.1, dpi=300, scale=1)

##### Perturbation scores for Tg transcriptome #####

# Calculate perturbation scores
t2.results.tg <- t2.test(screen, npcs = 20)
t2.results.tg <- left_join(t2.results.tg, CellsPerGene(screen), by="ident")
t2.results.tg$significant <- ifelse(t2.results.tg$p.value.bh<=0.01, yes="Y", no="N")
write.csv(t2.results.tg, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_5_TG_PERTURBATION_SCORES.csv")

# Plot perturbation scores
ggplot(t2.results.tg, aes(x=t2.stat, y=pmin(-log10(p.value.bh),15), colour=significant))+
  geom_hline(yintercept = 2, colour="grey80", size=0.5, linetype="dashed")+
  geom_point(aes(size=pmax(pmin(log10(n),3),1)), alpha=0.75)+
  scale_x_log10(limits=c(8,1000), oob=squish)+
  scale_y_continuous(limits=c(0,16), breaks=c(0,5,10,15), labels=c(0, 5, 10, ">15"), expand=expansion(add=c(1,0)))+
  scale_colour_manual(values=c("grey70", "grey50"))+
  scale_radius(limits=c(1,3), breaks=c(1,2,3), labels = c(10,100,1000), range=c(2,10))+
  theme_pubr(border=T, legend="none")+
  labs(x="Perturbation score (t\u00b2 statistic)", y="-Log10(adjusted p-value)")+
  theme(text=element_text(size=15), 
        plot.margin = margin(t = 0,0,0,0), line=element_line(size = 1), panel.border = element_rect(size=1, fill = NA))
ggsave("../outputs/fig_s4/FIG_S4B_TG_PERTURBATION_SCORES.png", width = 5.4, height = 5.6, scale = 1, dpi = 300)

##### Find differentially expressed Tg genes in GRA3 knockout #####

# Find markers
markers.gra3 <- FindMarkers(screen, ident.1 = "GRA3", assay = "SCT", logfc.threshold = 0, min.pct = 0) %>% rownames_to_column(var="gene")
markers.gra3$gene <- gsub(x=markers.gra3$gene, pattern = "TGME49", rep="")
markers.gra3$gene <- gsub(x=markers.gra3$gene, pattern = "227280", rep="GRA3")
markers.gra3$p_val_bh <- p.adjust(markers.gra3$p_val, method="BH")
write.csv(markers.gra3, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_6_GRA3_DE_TG_GENES.csv")

ggplot(markers.gra3, aes(x=avg_log2FC, y=pmin(-log10(p_val_bh),100)))+
  geom_point(colour="grey50")+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_continuous(limits=c(0, 105), labels=c(0, 25, 75, 50, ">100"), expand=c(0,0))+
  theme_pubr()+
  theme(text=element_text(size=15), plot.margin = margin(t = 0, r = 0, b = 0, l = 0))+
  labs(x="Average Log2FC", y="-Log10(adjusted p-value)")
ggsave("../outputs/fig_s4/FIG_S4C_GRA_DE_TG_GENES.png", width = 5.4, height = 5.6, scale = 1, dpi = 300)
