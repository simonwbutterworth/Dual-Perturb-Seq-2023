##### Identification of differentially expressed genes sets in dual perturb-seq screen #####

##### Setup #####
rm(list=ls())
gc()
source("../script/FUNCTIONS.R")
screen <- readRDS(file="../outputs/seurat_objects/SCREEN_MERGED_SEURAT.rds")
controls <- c("MYR1", "MYR2", "MYR3", "MYR4", "ROP17", "GRA45", "GRA16", "GRA24", "GRA28", "IST", "HCE1", "ROP16", "PPM3C", "MAF1B")

##### Run Seurat Wilcox test for all genes for each identified effector#####

#Select hits from screen
hits <- read.csv("../outputs/supplementary_data/SUPPLEMENTARY_DATA_4A_PERTURBATION_SCORES_COMBINED.csv") %>%
  filter(p.value.bh<=0.01) %>%
  pull(target)

#Create a temporary directory to hold DE results for each hit
dir.create("../temp/")

#Run Wilcox test for all genes for each identified effector
clock <- 1
for (i in hits) {
  print(paste("Starting", i, "..."))
  Idents(screen) <- ifelse(screen$guide.target.identity==i,
                           yes=i,
                           no=ifelse(screen$guide.target.identity%in%controls, yes="Control", no="Other"))
  markers.i <- FindMarkers(object=screen, ident.1=i, ident.2="Other", assay="SCT", slot="data", logfc.threshold=-Inf, min.pct=0, verbose = T)
  markers.i <- markers.i %>%
    rownames_to_column(var="gene") %>%
    mutate(target=i, .before=gene) %>%
    select(-p_val_adj) %>%
    mutate(p_val_bh=p.adjust(p_val, method="BH")) %>%
    relocate(p_val, .before=p_val_bh)
  write.csv(markers.i, file = paste("../temp/SCREEN_MARKERS_", clock, "_", i, ".csv", sep=""), row.names = F)
  print(paste(i, " complete (", clock, " of ", length(hits), ")", sep="" ))#
  clock=clock+1
}

#Read in files for each effector and rbind together
markers.all <- do.call(rbind,
                       lapply(list.files(path = "../temp/", pattern = "SCREEN_MARKERS_*", full.names = T), read.csv))
write.csv(markers.all, file = "../outputs/supplementary_data/SUPPLEMENTARY_DATA_9_DE_GENES.csv", row.names = F)

