# Part 1 takes us from feature-barcode matrices to a unified Seurat object. 
# It also filters for standard QC metrics and creates a metadata table with CRISPR guide information. 

library(Seurat)
library(dplyr)
library(ggplot2)

file_directories <- c('L415_410' = 'L415_410/per_sample_outs/L415_410/count/sample_filtered_feature_bc_matrix', 
                      'L415_411' = 'L415_411/per_sample_outs/L415_411/count/sample_filtered_feature_bc_matrix', 
                      'L415_412' = 'L415_412/per_sample_outs/L415_412/count/sample_filtered_feature_bc_matrix', 
                      'L415_413' = 'L415_413/per_sample_outs/L415_413/count/sample_filtered_feature_bc_matrix', 
                      'L415_414' = 'L415_414/per_sample_outs/L415_414/count/sample_filtered_feature_bc_matrix', 
                      'L415_415' = 'L415_415/per_sample_outs/L415_415/count/sample_filtered_feature_bc_matrix', 
                      'L415_416' = 'L415_416/per_sample_outs/L415_416/count/sample_filtered_feature_bc_matrix')
Seurat_data <- Read10X(data.dir=file_directories)
Microglia_Perturb_Seq <- CreateSeuratObject(counts=Seurat_data$`Gene Expression`,
                                            min.features=100,
                                            project="Microglia_Perturb_Seq")
saveRDS(Microglia_Perturb_Seq,file="MPS_Seurat.RDS")

# Add number of genes per UMI for each cell to metadata
Microglia_Perturb_Seq$log10GenesPerUMI <- log10(Microglia_Perturb_Seq$nFeature_RNA) / log10(Microglia_Perturb_Seq$nCount_RNA)
# Compute percent mito ratio
Microglia_Perturb_Seq$mitoRatio <- PercentageFeatureSet(object = Microglia_Perturb_Seq, pattern = "^MT-")
Microglia_Perturb_Seq$mitoRatio <- Microglia_Perturb_Seq@meta.data$mitoRatio / 100
# Create metadata dataframe
metadata <- Microglia_Perturb_Seq@meta.data
# Fix metadata format 
metadata$orig.ident <- gsub('_................-1','', rownames(metadata))
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Add metadata back to Seurat object
Microglia_Perturb_Seq@meta.data <- metadata
# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = Microglia_Perturb_Seq, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
# Extract counts
counts <- LayerData(filtered_seurat, layer = "counts")
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

### Compare filtered object vs original
View(filtered_seurat@meta.data)
filtered_metadata <- filtered_seurat@meta.data
dim(metadata) # 105509 cells
dim(filtered_metadata) # 99601 cells

# Normalize data
MPS_filtered <- NormalizeData(MPS_filtered)
MPS_filtered <- FindVariableFeatures(MPS_filtered, 
                                     selection.method = "vst", 
                                     nfeatures=2000, 
                                     verbose=TRUE)
top100_MPS <- head(VariableFeatures(MPS_filtered), 100)
top100_MPS
plot1 <- VariableFeaturePlot(MPS_filtered)
LabelPoints(plot=plot1, points = top100_MPS, repel=TRUE, xnudge=0, ynudge=0)
MPS_filtered <- ScaleData(MPS_filtered)

#### CRISPR guide data
CRISPR_guide_data <- Seurat_data$`CRISPR Guide Capture`
## Subset/filter CRISPR data to match Seurat data
CRISPR_guide_data <- CRISPR_guide_data[,rownames(MPS_filtered@meta.data)]
## Correct columns in CRISPR data? 
all(colnames(CRISPR_guide_data) %in% rownames(MPS_filtered@meta.data)) # TRUE
all(colnames(CRISPR_guide_data) == rownames(MPS_filtered@meta.data)) # TRUE

# Generate new metadata for each guide RNA
append_to_metadata <- as.data.frame(matrix('Does_not_contain', 
                                           nrow = dim(CRISPR_guide_data)[1], 
                                           ncol = dim(CRISPR_guide_data)[2]))
rownames(append_to_metadata) <- rownames(CRISPR_guide_data)
colnames(append_to_metadata) <- colnames(CRISPR_guide_data)

# Get any cells with positive transduction
append_to_metadata[as.matrix(CRISPR_guide_data) >= 3] = 'Contains'
View(append_to_metadata)

# Add metadata
metadata_with_CRISPR <- MPS_filtered@meta.data

append_to_metadata <- t(append_to_metadata)
all(rownames(append_to_metadata) == rownames(metadata_with_CRISPR)) # TRUE
metadata_with_CRISPR <- cbind(metadata_with_CRISPR,append_to_metadata)
View(metadata_with_CRISPR)
MPS_filtered@meta.data <- metadata_with_CRISPR

## MOI - Calculate and plot transduction per cell
contains_count <- rowSums(metadata_with_CRISPR == "Contains")
rows_with_one_contains <- sum(contains_count == 1)
rows_with_one_contains # 24955 cells with only one guide
MOI_numbers <- c(24955,17857,10771,6011,6312)
MOI_labels <- c('1','2','3','4','5+')
barplot(MOI_numbers,names.arg=MOI_labels,xlab='Number of gRNAs',ylab='Number of Cells',
        ylim=c(0,30000),main='MOI of Microglia Perturb-seq')

## Cells with or without guides
rows_with_guides <- sum(contains_count >= 1)
rows_without_guides <- sum(contains_count == 0)
rows_with_guides # 65906 cells have at least one guide
rows_without_guides # 33695 cells have no guide assigned
