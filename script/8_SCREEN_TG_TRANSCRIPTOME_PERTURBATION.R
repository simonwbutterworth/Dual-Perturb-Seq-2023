##### Investigation whether any effector protein knockouts alter the Tg parasite transcriptome #####

##### Setup #####
rm(list=ls())
gc()
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
t2.results.tg <- t2.test(screen, npcs = 20) %>%
  left_join(CellsPerGene(screen), by="ident") %>%
  rename(target=ident, n.cells=n) %>%
  mutate(significance=ifelse(p.value.bh<=0.01, yes="Significant", no="Non-significant"))
write.csv(t2.results.tg, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_5_TG_PERTURBATION_SCORES.csv", row.names = F)

# Plot perturbation scores
ggplot(t2.results.tg, aes(x=t2.stat, y=pmin(-log10(p.value.bh),15), colour=significance))+
  geom_hline(yintercept = 2, colour="grey80", linewidth=0.5, linetype="dashed")+
  geom_point(aes(size=pmax(pmin(log10(n.cells),3),1)), alpha=0.75)+
  scale_x_log10(limits=c(8,1000), oob=squish)+
  scale_y_continuous(limits=c(0,16), breaks=c(0,5,10,15), labels=c(0, 5, 10, ">15"), expand=expansion(add=c(1,0)))+
  scale_colour_manual(values=c("grey70", "grey50"))+
  scale_radius(limits=c(1,3), breaks=c(1,2,3), labels = c(10,100,1000), range=c(2,10))+
  theme_pubr(border=T, legend="none")+
  labs(x="Perturbation score (t\u00b2 statistic)", y="-Log10(adjusted p-value)")+
  theme(text=element_text(size=15), 
        plot.margin = margin(t = 0,0,0,0), line=element_line(size = 1), panel.border = element_rect(size=1, fill = NA), axis.title.y = element_text(margin=margin(0,11,0,0)))
ggsave("../outputs/fig_s4/FIG_S4B_TG_PERTURBATION_SCORES.png", width = 5.4, height = 5.6, scale = 1, dpi = 300)

##### Find differentially expressed Tg genes in GRA3 knockout #####

# Find markers
markers.gra3 <- FindMarkers(screen, ident.1 = "GRA3", assay = "SCT", logfc.threshold = -Inf, min.pct = 0) %>% 
  rownames_to_column(var="gene") %>%
  select(-p_val_adj) %>%
  mutate(gene=gsub(x=gene, pattern="-", rep="_")) %>%
  mutate(p_val_bh=p.adjust(p_val, method="BH")) %>%
  mutate(target="GRA3") %>%
  relocate(target, .before=gene) %>%
  relocate(p_val, .before=p_val_bh)
write.csv(markers.gra3, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_6_GRA3_DE_TG_GENES.csv", row.names = F)

ggplot(markers.gra3, aes(x=avg_log2FC, y=pmin(-log10(p_val_bh),100)))+
  geom_point(colour="grey50", size=2)+
  scale_x_continuous(limits=c(-3,3))+
  scale_y_continuous(limits=c(0, 105), labels=c(0, 25, 75, 50, ">100"), expand=c(0,0))+
  theme_pubr()+
  theme(text=element_text(size=15), plot.margin = margin(t = 0, r = 0, b = 0, l = 0))+
  labs(x="Average Log2FC", y="-Log10(adjusted p-value)")
ggsave("../outputs/fig_s4/FIG_S4C_GRA3_DE_TG_GENES.png", width = 5.4, height = 5.6, scale = 1, dpi = 300)

##### Test whether sgRNA target genes are downregulated #####

#Get list of targets represented by at least three cells and that are present in SCT data (genes that aren't in SCT data have all zero values)
Idents(screen) <- screen$guide.target.gene
cells.per.gene.3 <- CellsPerGene(screen) %>%filter(n>=3)
target.gene.ids <- data.frame(target.gene.id=screen$guide.target.gene, target.gene.name=screen$guide.target.identity) %>% 
  distinct() %>% 
  filter(target.gene.name%in%cells.per.gene.3$ident) %>% 
  filter(gsub(x=target.gene.id, pattern="_", replacement="-") %in% rownames(screen@assays$SCT@data))

#Run Seurat DE test for target genes in cells with cogate sgRNA
target.gene.foldchanges <- data.frame(gene=NA, p_val=NA, avg_log2FC=NA, pct.1=NA, pct.2=NA, p_val_adj=NA)
for (i in target.gene.ids$target.gene.id) {
  target.gene.foldchanges.i <- FindMarkers(screen, ident.1 = i, features = gsub(x = i, pattern = "_", replacement = "-"), assay = "SCT", logfc.threshold = -Inf, min.pct = 0)
  target.gene.foldchanges.i$gene <- i
  target.gene.foldchanges <- rbind(target.gene.foldchanges, target.gene.foldchanges.i)
}  
target.gene.foldchanges <- target.gene.foldchanges %>%
  filter(!is.na(gene)) %>%
  mutate(target=target.gene.ids$target.gene.name) %>%
  select(-p_val_adj) %>%
  mutate(p_val_bh=p.adjust(p_val, method="BH")) %>%
  relocate(target, .before=gene) %>%
  relocate(p_val, .before=p_val_bh)
write.csv(target.gene.foldchanges, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_7_DE_TG_TARGETED_GENES.csv", row.names = F)

#Plot DE genes
target.gene.foldchanges$p_val_bh <- ifelse(target.gene.foldchanges$p_val_bh<10^-100, yes = 10^-100, no=target.gene.foldchanges$p_val_bh)
p1 <- ggplot(filter(target.gene.foldchanges, pct.2>=0.25), aes(x=avg_log2FC))+
  geom_density(adjust=1, fill="grey50", colour="grey50", alpha=0.5, size=1.5)+
  scale_x_continuous(limits=c(-3.1,3.1), expand = c(0,0), oob=squish)+
  scale_y_continuous(limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.5))+
  theme_pubr()+
  labs(x=NULL, y="Density")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text=element_text(size=15), plot.margin = margin(5,0,10,13))
p2 <- ggplot(filter(target.gene.foldchanges, pct.2>=0.25), aes(x=avg_log2FC, y=-log10(p_val_bh), label=target))+
  geom_point(size=2, colour="grey50")+
  scale_x_continuous(limits=c(-3.1,3.1), expand = c(0,0), oob=squish)+
  scale_y_continuous(limits=c(0, 105), labels=c(0, 25, 75, 50, ">100"), expand=c(0,0))+
  geom_text_repel(data=filter(target.gene.foldchanges, pct.2>=0.25&p_val_bh<10^-25), nudge_x = 0.4, min.segment.length = 10, size=5, colour="grey50")+
  theme_pubr()+
  theme(text = element_text(size=15), plot.margin = margin(0,0,0,0))+
  labs(x="Average Log2FC", y="-Log10(adjusted p-value)")
ggarrange(plotlist=list(p1, p2), ncol=1, nrow=2, heights = c(1,4))
ggsave("../outputs/fig_s4/FIG_S4D_GRA3_DE_TG_GENES.png", width = 5.4, height = 5.6, scale = 1, dpi = 300)


