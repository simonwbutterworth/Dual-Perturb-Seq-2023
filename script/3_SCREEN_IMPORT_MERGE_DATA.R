##### Import dual perturb-seq screen data, filter cell barcodes, merge replicates, normalise, run PCA/UMAP #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")

##### Import, filter, and merge data from unstimulated screen samples #####

# Import
screen.unstim.1a <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A21_BUT290A22_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A21_BUT290A22_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A21_BUT290A22_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_1A", sample.id.gex = "BUT290A21", sample.id.crispr = "BUT290A22")
screen.unstim.1b <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A23_BUT290A24_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A23_BUT290A24_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A23_BUT290A24_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_1B", sample.id.gex = "BUT290A23", sample.id.crispr = "BUT290A24")
screen.unstim.2a <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A25_BUT290A26_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A25_BUT290A26_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A25_BUT290A26_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_2A", sample.id.gex = "BUT290A25", sample.id.crispr = "BUT290A26")
screen.unstim.2b <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A27_BUT290A28_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A27_BUT290A28_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A27_BUT290A28_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_2B", sample.id.gex = "BUT290A27", sample.id.crispr = "BUT290A28")
screen.unstim.3a <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A37_BUT290A38_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A37_BUT290A38_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A37_BUT290A38_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_3A", sample.id.gex = "BUT290A37", sample.id.crispr = "BUT290A38")
screen.unstim.3b <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A39_BUT290A40_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A39_BUT290A40_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A39_BUT290A40_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_3B", sample.id.gex = "BUT290A39", sample.id.crispr = "BUT290A40")
screen.unstim.3c <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A41_BUT290A42_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A41_BUT290A42_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A41_BUT290A42_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_3C", sample.id.gex = "BUT290A41", sample.id.crispr = "BUT290A42")
screen.unstim.3d <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A43_BUT290A44_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A43_BUT290A44_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A43_BUT290A44_features.tsv.gz", sample.name = "SCREEN_UNSTIMULATED_3D", sample.id.gex = "BUT290A43", sample.id.crispr = "BUT290A44")

# Check nCount_Tg and nCount_Hs
nCountPlot(screen.unstim.1a)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 1000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.1b)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 1000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.2a)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 4000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.2b)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 1000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.3a)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 2000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.3b)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 2000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.3c)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 2000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.unstim.3d)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 2000, xmax = 1000000, fill=NA, colour="grey20")

# Check percent MT
VlnPlot(screen.unstim.1a, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.1b, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.2a, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.2b, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.3a, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.3b, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.3c, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.unstim.3d, features=c("percent.mt"))+geom_hline(yintercept=10)

# Filter cell barcodes
screen.unstim.1a <- subset(screen.unstim.1a, subset=nCount_Hs>=5000&nCount_Tg>=1000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.1b <- subset(screen.unstim.1b, subset=nCount_Hs>=5000&nCount_Tg>=1000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.2a <- subset(screen.unstim.2a, subset=nCount_Hs>=5000&nCount_Tg>=4000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.2b <- subset(screen.unstim.2b, subset=nCount_Hs>=5000&nCount_Tg>=1000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.3a <- subset(screen.unstim.3a, subset=nCount_Hs>=5000&nCount_Tg>=2000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.3b <- subset(screen.unstim.3b, subset=nCount_Hs>=5000&nCount_Tg>=2000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.3c <- subset(screen.unstim.3c, subset=nCount_Hs>=5000&nCount_Tg>=2000&percent.mt<=10&nFeature_sgRNA==1)
screen.unstim.3d <- subset(screen.unstim.3d, subset=nCount_Hs>=5000&nCount_Tg>=2000&percent.mt<=10&nFeature_sgRNA==1)

##### Import, filter, and merge data from IFNG-stimulated screen samples #####

# Import 
screen.ifng.1a <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A29_BUT290A30_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A29_BUT290A30_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A29_BUT290A30_features.tsv.gz", sample.name = "SCREEN_IFNG_1A", sample.id.gex = "BUT290A29", sample.id.crispr = "BUT290A30")
screen.ifng.1b <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A31_BUT290A32_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A31_BUT290A32_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A31_BUT290A32_features.tsv.gz", sample.name = "SCREEN_IFNG_1B", sample.id.gex = "BUT290A31", sample.id.crispr = "BUT290A32")
screen.ifng.1c <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A33_BUT290A34_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A33_BUT290A34_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A33_BUT290A34_features.tsv.gz", sample.name = "SCREEN_IFNG_1C", sample.id.gex = "BUT290A33", sample.id.crispr = "BUT290A34")
screen.ifng.1d <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A35_BUT290A36_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A35_BUT290A36_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A35_BUT290A36_features.tsv.gz", sample.name = "SCREEN_IFNG_1D", sample.id.gex = "BUT290A35", sample.id.crispr = "BUT290A36")
screen.ifng.2a <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A45_BUT290A46_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A45_BUT290A46_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A45_BUT290A46_features.tsv.gz", sample.name = "SCREEN_IFNG_2A", sample.id.gex = "BUT290A45", sample.id.crispr = "BUT290A46")
screen.ifng.2b <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A47_BUT290A48_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A47_BUT290A48_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A47_BUT290A48_features.tsv.gz", sample.name = "SCREEN_IFNG_2B", sample.id.gex = "BUT290A47", sample.id.crispr = "BUT290A48")
screen.ifng.2c <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A49_BUT290A50_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A49_BUT290A50_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A49_BUT290A50_features.tsv.gz", sample.name = "SCREEN_IFNG_2C", sample.id.gex = "BUT290A49", sample.id.crispr = "BUT290A50")
screen.ifng.2d <- ImportFeatureBC_Mtx(matrix = "../data/GSE229505_BUT290A51_BUT290A52_matrix.mtx.gz", cells = "../data/GSE229505_BUT290A51_BUT290A52_barcodes.tsv.gz", features = "../data/GSE229505_BUT290A51_BUT290A52_features.tsv.gz", sample.name = "SCREEN_IFNG_2D", sample.id.gex = "BUT290A51", sample.id.crispr = "BUT290A52")

# Check nCount_Tg and nCount_Hs
nCountPlot(screen.ifng.1a)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.1b)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.1c)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.1d)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.2a)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.2b)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.2c)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")
nCountPlot(screen.ifng.2d)+annotate(geom="rect", ymin = 5000, ymax = 1000000, xmin = 3000, xmax = 1000000, fill=NA, colour="grey20")

# Check percent MT
VlnPlot(screen.ifng.1a, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.1b, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.1c, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.1d, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.2a, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.2b, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.2c, features=c("percent.mt"))+geom_hline(yintercept=10)
VlnPlot(screen.ifng.2d, features=c("percent.mt"))+geom_hline(yintercept=10)

# Filter cell barcodes
screen.ifng.1a <- subset(screen.ifng.1a, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.1b <- subset(screen.ifng.1b, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.1c <- subset(screen.ifng.1c, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.1d <- subset(screen.ifng.1d, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.2a <- subset(screen.ifng.2a, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.2b <- subset(screen.ifng.2b, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.2c <- subset(screen.ifng.2c, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)
screen.ifng.2d <- subset(screen.ifng.2d, subset=nCount_Hs>=5000&nCount_Tg>=3000&percent.mt<=10&nFeature_sgRNA==1)

##### Merge, Assign sgRNAs, SCTransform, PCA, UMAP, saveRDS #####

# Merge replicates
screen <- merge(x =   screen.unstim.1a, 
                y = c(screen.unstim.1b,
                      screen.unstim.2a,
                      screen.unstim.2b,
                      screen.unstim.3a,
                      screen.unstim.3b, 
                      screen.unstim.3c,
                      screen.unstim.3d,
                      screen.ifng.1a,
                      screen.ifng.1b,
                      screen.ifng.1c,
                      screen.ifng.1d,
                      screen.ifng.2a,
                      screen.ifng.2b,
                      screen.ifng.2c,
                      screen.ifng.2d))

# Remove old objects
screen.unstim.1a <- NULL
screen.unstim.1b <- NULL
screen.unstim.2a <- NULL
screen.unstim.2b <- NULL
screen.unstim.3a <- NULL
screen.unstim.3b <- NULL
screen.unstim.3c <- NULL
screen.unstim.3d <- NULL
screen.ifng.1a <- NULL
screen.ifng.1b <- NULL
screen.ifng.1c <- NULL
screen.ifng.1d <- NULL
screen.ifng.2a <- NULL
screen.ifng.2b <- NULL
screen.ifng.2c <- NULL
screen.ifng.2d <- NULL

# Assign unstimalated/IFNG
screen$ifng <- ifelse(screen$orig.ident%in%c("SCREEN_UNSTIMULATED_1A",
                                             "SCREEN_UNSTIMULATED_1B",
                                             "SCREEN_UNSTIMULATED_2A",
                                             "SCREEN_UNSTIMULATED_2B",
                                             "SCREEN_UNSTIMULATED_3A",
                                             "SCREEN_UNSTIMULATED_3B",
                                             "SCREEN_UNSTIMULATED_3C",
                                             "SCREEN_UNSTIMULATED_3D"), yes="Unstimulated", no="IFNG")

# Assign sgRNAs
screen <- Assign_sgRNAs(screen, reference.file = "../inputs/seurat_sgrna_reference/SCREEN_SEURAT_REFERENCE_220428.csv")
screen$guide.target.identity <- gsub(screen$guide.target.identity, pattern="TGME49_", rep="")
Idents(screen) <- screen$guide.target.identity

# SCTransform, PCA, UMAP
screen <- SCTransform(screen, vst.flavor = "v2", assay = "Hs", vars.to.regress = "orig.ident")
screen <- RunPCA(screen, assay="SCT")
screen <- RunUMAP(screen, dims=1:10)

# Save RDS
saveRDS(screen, file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
