##### Further analysis of hit from screen: GRA59 #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
controls <- c("MYR1", "MYR2", "MYR3", "MYR4", "ROP17", "GRA45", "GRA16", "GRA24", "GRA28", "IST", "HCE1", "ROP16", "PPM3C", "MAF1B")

##### Run VISION for MYR1-infected-like signatures score #####

# Get embeddings & metadata
pca.embeddings <- as.data.frame(Embeddings(screen))
screen$guide.target.identity <- factor(screen$guide.target.identity)
metadata <- screen[[]]

# Run Vision with MYR1 signature
screen.vision.myr1 <- Vision(data=screen@assays$SCT@data[,], 
                             latentSpace=pca.embeddings, 
                             meta=metadata, 
                             projection_methods=NULL, 
                             signatures="../inputs/gene_sets/MYR1_MARKERS_NAOR_2018.txt")
screen.vision.myr1 <- calcSignatureScores(object=screen.vision.myr1)


# Extract scores and set idents
screen.vision.myr1 <- data.frame(ident=screen.vision.myr1@metaData$guide.target.identity, signature="MYR1_MARKERS_NAOR_2018", score=screen.vision.myr1@SigScores[,1])
screen.vision.myr1$test <- ifelse(screen.vision.myr1$ident%in%controls, yes=NA,      no="OTHER")
screen.vision.myr1$test <- ifelse(screen.vision.myr1$ident=="MYR1",     yes="MYR1",  no=screen.vision.myr1$test)
screen.vision.myr1$test <- ifelse(screen.vision.myr1$ident=="GRA59",    yes="GRA59", no=screen.vision.myr1$test)
screen.vision.myr1$test <- factor(screen.vision.myr1$test, levels=c("OTHER", "GRA59", "MYR1"))

# Wilcox text
wilcox.test(x=screen.vision.myr1$score[screen.vision.myr1$test=="GRA59"],
            y=screen.vision.myr1$score[screen.vision.myr1$test=="OTHER"])$p.value*3
wilcox.test(x=screen.vision.myr1$score[screen.vision.myr1$test=="GRA59"],
            y=screen.vision.myr1$score[screen.vision.myr1$test=="MYR1"])$p.value*3
wilcox.test(x=screen.vision.myr1$score[screen.vision.myr1$test=="MYR1"],
            y=screen.vision.myr1$score[screen.vision.myr1$test=="OTHER"])$p.value*3

# Plot signature scores
ggplot(filter(screen.vision.myr1, !is.na(test)), aes(x=test, y=score, fill=test))+
  geom_violin(adjust=1.5, draw_quantiles = 0.5, size=0.5, colour="black")+
  scale_fill_manual(values = c("grey50", "goldenrod", "slateblue"))+
  theme_pubr(legend="none")+
  scale_x_discrete(labels=c("Other", "sg(GRA59)", "sg(MYR1)"))+
  labs(x=NULL, y="Signature score")+
  theme(text=element_text(size=15), plot.margin = margin(0,0,0,2))+
  scale_y_continuous(limits=c(-0.8, 0.65), breaks=seq(-1,1,0.5))+
  geom_signif(xmin=1, xmax=1.95, y_position = 0.45, tip_length = 0, size=0.5, textsize = 5, annotations = "< 0.001")+
  geom_signif(xmin=2.05, xmax=3, y_position = 0.45, tip_length = 0, size=0.5, textsize = 5, annotations = "< 0.001")+
  geom_signif(xmin=1, xmax=3, y_position = 0.6, tip_length = 0, size=0.5, textsize = 5, annotations = "< 0.001")
ggsave(filename = "../outputs/fig_3/FIG_3A_GRA59_MYR1_SIGNATURE_SCORE.png", width = 4.4, height = 2.9, dpi=600)

