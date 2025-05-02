# Script to extract count data & metadata from Seurat object

library(Seurat) 
library(cowplot)
library(dplyr)

CLUSTER <- 8
ALLOW_CONDS <- c("WT", "CRE", "PRE", "MLEU")

# Read Seurat Data
Bcells.combined <- readRDS(file ='seurat_obj.rds')

temp <- SubsetData(object=Bcells.combined, subset.name = "condition", accept.value = ALLOW_CONDS )

selected_data <- SubsetData(object=temp, subset.name = "seurat_clusters", accept.value = CLUSTER )

rm(Bcells.combined, temp)

# Extract data
write.csv(x=t(as.matrix(GetAssayData(object = selected_data, slot = "counts"))) , file = "counts.csv", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)
# or to select RNA
write.csv(x=t(as.matrix(GetAssayData(object = selected_data, assay = "RNA", slot = "data"))) , file = "rna.csv", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE)

write.csv(x=selected_data$seurat_clusters, file = "meta.csv", sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
write.csv(x=selected_datacondition, file = "condition.csv", sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
write.csv(x=selected_data$replicate, file = "replicate.csv", sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
write.csv(x=selected_data$sequencing, file = "sequencing.csv", sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)


