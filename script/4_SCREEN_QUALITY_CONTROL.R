##### Quality control analysis of dual perturb-seq screen data #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")

##### UMAP plots ##### 

# Create new ident for experimental replicate
screen$date.ident <- NA
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_1A", yes="2021_01_10_UNSTIM_A", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_1B", yes="2021_01_10_UNSTIM_B", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_2A", yes="2021_06_25_UNSTIM_C", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_2B", yes="2021_06_25_UNSTIM_D", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_3A", yes="2022_03_23_UNSTIM_E", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_3B", yes="2022_03_23_UNSTIM_F", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_3C", yes="2022_03_23_UNSTIM_G", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_UNSTIMULATED_3D", yes="2022_03_23_UNSTIM_H", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_1A", yes="2021_10_28_IFNG_A", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_1B", yes="2021_10_28_IFNG_B", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_1C", yes="2021_10_28_IFNG_C", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_1D", yes="2021_10_28_IFNG_D", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_2A", yes="2022_04_06_IFNG_E", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_2B", yes="2022_04_06_IFNG_F", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_2C", yes="2022_04_06_IFNG_G", no=screen$date.ident)
screen$date.ident <- ifelse(screen$orig.ident=="SCREEN_IFNG_2D", yes="2022_04_06_IFNG_H", no=screen$date.ident)             
screen$date.ident <- factor(screen$date.ident, levels=c("2021_01_10_UNSTIM_A",
                                                        "2021_01_10_UNSTIM_B",
                                                        "2021_06_25_UNSTIM_C",
                                                        "2021_06_25_UNSTIM_D",
                                                        "2022_03_23_UNSTIM_E",
                                                        "2022_03_23_UNSTIM_F",
                                                        "2022_03_23_UNSTIM_G",
                                                        "2022_03_23_UNSTIM_H",
                                                        "2021_10_28_IFNG_A",
                                                        "2021_10_28_IFNG_B",
                                                        "2021_10_28_IFNG_C",
                                                        "2021_10_28_IFNG_D",
                                                        "2022_04_06_IFNG_E",
                                                        "2022_04_06_IFNG_F",
                                                        "2022_04_06_IFNG_G",
                                                        "2022_04_06_IFNG_H"))

# UMAP plot coloured by replicate
data.frame(UMAP_1=screen@reductions$umap@cell.embeddings[,"UMAP_1"],
           UMAP_2=screen@reductions$umap@cell.embeddings[,"UMAP_2"],
           date=screen$date.ident,
           rand=runif(length(screen$date.ident), 0, 1)) %>%
  arrange(rand) %>%
  ggplot()+
  geom_point(mapping=aes(x=UMAP_1, y=UMAP_2, colour=date), size=0.1)+
  theme_pubr(border = T, legend="none")+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s3/FIG_S3A_SCREEN_UMAP_REPLICATE.png", width=3.7, height=3.7, dpi=300, scale=1)

# UMAP plot coloured by replicate split into panels
data.frame(UMAP_1=screen@reductions$umap@cell.embeddings[,"UMAP_1"],
           UMAP_2=screen@reductions$umap@cell.embeddings[,"UMAP_2"],
           date=screen$date.ident) %>%
  ggplot()+
  geom_point(mapping=aes(x=UMAP_1, y=UMAP_2, colour=date), size=0.1, alpha=1)+
  theme_pubr(border = T, legend="none")+
  labs(x=NULL, y=NULL, title=NULL)+
  scale_x_continuous(expand = expansion(mult=0.05))+
  scale_y_continuous(expand = expansion(mult=c(0.05,0.4)))+
  facet_wrap(~ date, dir = "v", nrow = 4)+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"),
        aspect=1,
        panel.spacing = unit(0, unit = "pt"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
ggsave(filename = "../outputs/fig_s3/FIG_S3B_SCREEN_UMAP_REPLICATE_PANELS.png", width=3.7, height=3.7, dpi=300, scale=1)

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
ggsave(filename = "../outputs/fig_s3/FIG_S3C_SCREEN_UMAP_UNSTIMULATED_IFNG.png", width=3.7, height=3.7, dpi=300, scale=1)

# Plot IRF1 expression
ggplot()+
  geom_point(mapping=aes(x=screen@reductions$umap@cell.embeddings[,"UMAP_1"], 
                         y=screen@reductions$umap@cell.embeddings[,"UMAP_2"],
                         colour=screen@assays$SCT@data["IRF1",]), size=0.1, alpha=1)+
  theme_pubr(border = T, legend="none")+
  scale_colour_gradient(low = "grey80", high = "royalblue4", oob=squish)+
  theme(axis.text = element_blank(), axis.ticks=element_blank(), text=element_text(size=15), axis.ticks.length = unit(0, "pt"), axis.title = element_blank(),
        plot.margin = margin(t = 0, r=0, b = 0, l = 0, unit = "pt"), aspect=1)+
  labs(x=NULL, y=NULL, title=NULL)
ggsave(filename = "../outputs/fig_s3/FIG_S3D_SCREEN_UMAP_IRF1.png", width=3.7, height=3.7, dpi=300, scale=1)

##### Guide metrics #####

# Histogram of cells per guide
data.frame(guide=screen$guide) %>%
  group_by(guide) %>%
  dplyr::summarise(cells.per.guide=n()) %>%
  ggplot(aes(x=cells.per.guide))+
  geom_histogram(binwidth = 10, boundary=0, fill="grey50", colour="grey50")+
  scale_x_continuous(limits=c(0,110), oob = squish, expand=c(0,0), breaks=seq(0,100,20), labels=c(0,20,40,60,80,"100+"))+
  scale_y_continuous(limits=c(0,550), oob = squish, expand=c(0,0), breaks=seq(0,1000,100))+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin = margin(t=0, b=5, l=5, r=0))+
  labs(y="Count", x="Cells per sgRNA")
ggsave(filename = "../outputs/fig_s3/FIG_S3E_HISTOGRAM_CELLS_PER_GUIDE.png", width=5.4, height=4, dpi=300, scale=1)

length(unique(screen$guide)) 
median(table(screen$guide))

# Histogram of cells per gene
data.frame(guide=screen$guide.target.identity) %>%
  group_by(guide) %>%
  dplyr::summarise(cells.per.gene=n()) %>%
  ggplot(aes(x=cells.per.gene))+
  geom_histogram(binwidth = 10, boundary=0, fill="grey50", colour="grey50")+
  scale_x_continuous(limits=c(0,510), oob = squish, expand=c(0,0), breaks=seq(0,500,100), labels=c(0,100,200,300,400,"500+"))+
  scale_y_continuous(limits=c(0,32), expand=c(0,0), breaks=seq(0,30,10))+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin = margin(t=0, b=5, l=5, r=10))+
  labs(y="Count", x="Cells per target gene (sgRNAs aggregated)")
ggsave(filename = "../outputs/fig_s3/FIG_S3F_HISTOGRAM_CELLS_PER_GENE.png", width=5.4, height=4, dpi=300, scale=1)

length(unique(screen$guide.target.identity))
median(table(screen$guide.target.identity))

# Histogram sgRNAs per target gene
data.frame(guide=screen$guide,
           gene=screen$guide.target.identity) %>%
  distinct() %>%
  group_by(gene) %>%
  dplyr::summarise(n=n()) %>%
  ggplot()+
  geom_histogram(mapping=aes(x=n), binwidth = 1, fill="grey50", colour="grey50")+
  scale_x_continuous(limits=c(0.5,5.5), oob=squish, breaks=c(1,2,3,4,5), labels=c(1,2,3,4,5))+
  scale_y_continuous(limits=c(0,160), expand=c(0,0), breaks=seq(0,150,50))+
  theme_pubr()+
  theme(text=element_text(size=15),
        plot.margin = margin(t=0, b=5, l=5, r=0))+
  labs(y="Count", x="sgRNAs detected per target gene")
ggsave(filename = "../outputs/fig_s3/FIG_S3G_HISTOGRAM_SGRNAS_PER_GENE.png", width=5.4, height=4, dpi=300, scale=1)

data.frame(guide=screen$guide,
           gene=screen$guide.target.identity) %>%
  distinct() %>%
  group_by(gene) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::summarise(median(n))

temp <- data.frame(guide=screen$guide,
           gene=screen$guide.target.identity) %>%
  distinct() %>%
  group_by(gene) %>%
  dplyr::summarise(n=n())
length(which(temp$n>=3))/length(temp$n)

##### Correlation of guide with bulk library counts #####

# Guide-level correlation
cells.per.guide <- data.frame(table(screen$guide)) %>% rename(guide.r=Var1, n.cells=Freq)
phenotype <- read.csv("../inputs/phenotype_scores/PHENOTYPE_SCORES_SIDIK_2016.csv") %>%
  group_by(TGME49.ID) %>%
  dplyr::summarise(score=mean(score)) %>%
  select(TGME49.ID, score) %>%
  rename(gene=TGME49.ID)
counts.guide <- read.csv("../inputs/seurat_sgrna_reference/SCREEN_SEURAT_REFERENCE_220428.csv") %>% 
  rename(ident=identity) %>%
  left_join(cells.per.guide, by="guide.r") %>%
  left_join(phenotype, by="gene")

length(unique(counts.guide$guide))
length(unique(counts.guide$gene))

ggplot(counts.guide, aes(x=plasmid.counts, y=n.cells, colour=score))+
  geom_point(size=2)+
  scale_x_continuous(trans=log_trans(), limits=c(10,3000000), breaks=c(10,100,1000,10000,100000,1000000), expand = c(0,0), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans=log_trans(), limits=c(1,500), breaks=c(1,10,100,1000,10000,100000,1000000), expand=c(0,0), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_colour_gradient(low="royalblue4", high="grey80", na.value = "grey80", limits=c(-5,2.5))+
  theme_pubr(legend="none")+
  labs(x="Plasmid library count", y="Captured cell count")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 0,b = 5,r = 0,l = 5))
ggsave(filename = "../outputs/fig_s3/FIG_S3H_CORRELATION_PLASMID_VS_CELLS_PER_GUIDE.png", width=8.2, height=3.7, dpi=300, scale=1)

cor.test(x=counts.guide$plasmid.counts, y=counts.guide$n.cells, method = "pearson")  

write.csv(counts.guide, file="../outputs/supplementary_data/SUPPLEMENTARY_DATA_3_PROTOSPACER_SEQUENCES.csv")


# Import plasmid counts for guides and sum counts by gene
cells.per.gene <- data.frame(table(screen$guide.target.gene)) %>% rename(gene=Var1, n.cells=Freq)
phenotype <- read.csv("../inputs/phenotype_scores/PHENOTYPE_SCORES_SIDIK_2016.csv") %>%
  group_by(TGME49.ID) %>%
  dplyr::summarise(score=mean(score)) %>%
  select(TGME49.ID, score) %>%
  rename(gene=TGME49.ID)
counts.gene <- read.csv("../inputs/seurat_sgrna_reference/SCREEN_SEURAT_REFERENCE_220428.csv") %>%
  group_by(gene) %>%
  mutate(sum.plasmid.counts=sum(plasmid.counts),
         sum.prefreezing.counts=sum(prefreezing.counts),
         sum.postfreezing.counts=sum(postfreezing.counts)) %>%
  select(gene, sum.plasmid.counts, sum.prefreezing.counts, sum.postfreezing.counts) %>%
  ungroup() %>%
  distinct() %>% 
  left_join(cells.per.gene, by="gene") %>%
  left_join(phenotype, by="gene")

ggplot(counts.gene, aes(x=sum.plasmid.counts, y=n.cells, colour=score))+
  geom_point(size=2)+
  scale_x_continuous(trans=log_trans(), limits=c(10000,1100000), breaks=c(10,100,1000,10000,100000,1000000), expand = c(0,0), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous(trans=log_trans(), limits=c(1,2000), breaks=c(1,10,100,1000,10000,100000,1000000), expand=c(0,0), labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_colour_gradient(low="royalblue4", high="grey80", na.value = "grey80", limits=c(-5,2.5))+
  theme_pubr(legend="none")+
  labs(x="Plasmid library count", y="Captured cell count")+
  theme(text=element_text(size=15),
        plot.margin = margin(t = 0,b = 5,r = 0,l = 5))
ggsave(filename = "../outputs/fig_s3/FIG_S3H_CORRELATION_PLASMID_VS_CELLS_PER_GENE.png", width=8.2, height=3.7, dpi=300, scale=1)

cor.test(x=counts.gene$sum.plasmid.counts, y=counts.gene$n.cells, method = "pearson")  
