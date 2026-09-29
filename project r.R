install.packages("dplyr")
install.packages("ggplot2")
install.packages("Seurat")
install.packages("BiocManager")
install.packages("harmony")
BiocManager::install("AUCell")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scDblFinder")
BiocManager::install("bluster")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("miloR")
BiocManager::install("infercnv")
BiocManager::install("scater")
BiocManager::install("decoupleR")
BiocManager::install("dorothea")
BiocManager::install("tidyr")
BiocManager::install("tibble")
BiocManager::install("patchwork")
BiocManager::install("pheatmap")
BiocManager::install("CellChat") 
BiocManager::install("monocle3") 
devtools::install_github("jinworks/CellChat")
devtools::install_github('immunogenomics/presto')

library(dplyr)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(bluster)
library(SingleR)
library(celldex)
library(miloR)
library(MouseGastrulationData)
library(scater)
library(harmony)
library(AUCell)


#######QC & doublet detection###################################################
sc_data=Read10X(data.dir = "E:/single cell course/project")
project= CreateSeuratObject(
  counts = sc_data,
  project = "scRNAseq",
  min.cells=3 , min.features=200)
# Calculate Mitochondrial Percentage
project[["percent.mt"]] = PercentageFeatureSet(project, pattern = "^MT-")
# Visualization BEFORE Filtering
# Violin Plot
VlnPlot(project, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
# Scatter Plot
FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "percent.mt")
# Cell Filtering (QC)
n_before = ncol(project)
cat("Cells before filtering:", n_before, "\n")
project = subset(project, 
                 subset = nFeature_RNA > 200 & 
                   nFeature_RNA < 6000 & 
                   percent.mt < 10)
n_after = ncol(project)
cat("Cells after filtering:", n_after, "\n")
cat("Cells removed:", n_before - n_after, "\n")
cat("Percentage removed:", round((n_before - n_after) / n_before * 100, 2), "%\n")
# Visualization AFTER Filtering
VlnPlot(project, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

sce_d1=as.SingleCellExperiment(project)
set.seed(100)
sce_d1= scDblFinder(sce_d1)
project$doublet_score<-
  colData(sce_d1) $scDblFinder.score
project$doublet_class<-
  colData(sce_d1) $scDblFinder.class
table(project$doublet_class)
# Doublet Visualization
VlnPlot(project, features = "doublet_score",
        group.by = "doublet_class")
ggplot(project@meta.data, 
       aes(x = nCount_RNA, y = doublet_score, colour = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(x = "nCount_RNA", y = "Doublet Score") +
  theme_classic()
# Remove Doublets (Keep Singlets Only)
project<-subset(project,subset=doublet_class ==
                  "singlet")

#######normalization & scaling & feature selection##############################
project <- NormalizeData(
  project,
  normalization.method = "LogNormalize",
  scale.factor = 10000)
project <- FindVariableFeatures(
  project,
  selection.method = "vst",
  nfeatures = 2000)
top10 <- head(VariableFeatures(project), 10)
plot1 <- VariableFeaturePlot(project)
LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE)
VariableFeaturePlot(project) +
  LabelPoints(points = top10, repel = TRUE)
plot1 <- VariableFeaturePlot(project)
plot1 <- LabelPoints(
  plot   = plot1,
  points = top10,
  repel  = TRUE
)
project = ScaleData(project, features = VariableFeatures(project))

#######dimensionality reduction & clustering & UMAP#############################
#run PCA 
project <- RunPCA(
  project,
  features = VariableFeatures(object = project)
)
#visualization of PCA 
DimPlot(project, reduction = "pca")
#elbow plot
ElbowPlot(project)
#heatmap
DimHeatmap(project, dims = 1:10, cells = 500,
           balanced = TRUE)
# Find Neibhor Graph 
project <- FindNeighbors(project, dims = 1:10)
# Clustering the cells
project <- FindClusters(project, resolution = 0.2)
# RUn UMAP
project <- RunUMAP(project, dims = 1:10)


#########cell annotation######################################################## 
#1
canonical_markers <- c(
  "CD3D", "CD3E",     
  "MS4A1",            
  "NKG7",              
  "LYZ", "FCGR3A",    
  "CD1C",              
  "PPBP"               
)
FeaturePlot(project, features = canonical_markers)
VlnPlot(project, features = canonical_markers)
#2
project.markers = FindAllMarkers(project,
      only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(project.markers, file = "project.markers.cvs", row.names = FALSE)
top_5_project.markers = project.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  arrange(cluster, desc(avg_log2FC))

top5_genes <- top_5_project.markers$gene
top_markers <- project.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
View(top_markers)
project$cell_type <- dplyr::recode(
  as.character(project$seurat_clusters),
  "0"  = "Plasma cells (IgG lambda)",
  "1"  = "Cytotoxic T cells",
  "2"  = "Plasma cells (IgG kappa)",
  "3"  = "Fibroblasts",
  "4"  = "Gastric Epithelial cells",
  "5"  = "Macrophages",
  "6"  = "Plasma cells (IgG kappa)",
  "7"  = "Endothelial cells",
  "8"  = "B cells",
  "9"  = "Enteroendocrine cells",
  "10" = "Pericytes/Stromal cells",
  "11" = "Mast cells"
)
DimPlot(
  project,
  reduction = "umap",
  group.by  = "cell_type",
  label     = TRUE,
  repel     = TRUE
) +
  ggtitle("UMAP - Manual Marker-Based Annotation")

DoHeatmap(
  project,
  features = top5_genes,
  raster = TRUE,           
  disp.min = -2.5,
  disp.max = 2.5,
  group.bar = TRUE,        
  draw.lines = FALSE,
  size = 3                 
) + 
  scale_fill_gradientn(colors = c("magenta", "black", "yellow")) + 
  ggtitle("Top Differentially Expressed Marker Genes - Sample 33") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),      
    legend.text = element_text(size = 11),
    plot.margin = margin(15, 40, 15, 40)       
  )
ggsave("Heatmap_Top5_Markers.png", width = 20, height = 10, dpi = 300)

#Visualization
hpca_ref <- celldex::HumanPrimaryCellAtlasData()
expr_matrix <- GetAssayData(project, layer = "data")  
cluster_singler <- SingleR(
  test = expr_matrix,
  ref = hpca_ref,
  labels = hpca_ref$label.main,
  clusters = project$seurat_clusters
)
cluster_singler$labels   
cluster_map2 <- setNames(cluster_singler$labels, rownames(cluster_singler))
project$cluster_labels <- unname(cluster_map2[as.character(project$seurat_clusters)])


#two plots 
DimPlot(
  project,
  reduction = "umap",
  group.by  = "cluster_labels",
  label     = TRUE,
  repel     = TRUE
) +
  ggtitle("UMAP - Cell Type Names (Cluster-level)")

head(project@meta.data)

DimPlot(
  project,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE,
  repel     = TRUE
) +
  ggtitle("UMAP - Seurat Clusters")

#################downstream analysis############################################
################################################################################







# Macrophages
library(Seurat)

combined <- merge(project, y = s40, add.cell.ids = c("project", "s40"))
combined <- JoinLayers(combined)
combined <- NormalizeData(combined)

table(combined$cell_type)  

mac_cells <- subset(combined, subset = cell_type == "Macrophages")
DefaultAssay(mac_cells) <- "RNA"
mac_cells <- FindVariableFeatures(mac_cells, nfeatures = 2000)
mac_cells <- ScaleData(mac_cells)
mac_cells <- RunPCA(mac_cells, npcs = 20)
mac_cells <- FindNeighbors(mac_cells, reduction = "pca", dims = 1:20)
mac_cells <- FindClusters(mac_cells, resolution = 0.8)
mac_cells <- RunUMAP(mac_cells, reduction = "pca", dims = 1:20)

VlnPlot(mac_cells, features = c("SPP1", "C1QA", "C1QB", "C1QC"), group.by = "seurat_clusters")

# CD8+ T cells (combine both CD8-flavored clusters from your earlier annotation)
combined <- merge(project, y = s40, add.cell.ids = c("project", "s40"))
combined <- JoinLayers(combined)
combined <- NormalizeData(combined)

cd8_cells <- subset(combined, subset = cell_type == "Cytotoxic T cells")
DefaultAssay(cd8_cells) <- "RNA"
cd8_cells <- FindVariableFeatures(cd8_cells, nfeatures = 2000)
cd8_cells <- ScaleData(cd8_cells)
cd8_cells <- RunPCA(cd8_cells, npcs = 20)

table(cd8_cells$orig.ident)

cd8_cells <- RunHarmony(cd8_cells, group.by.vars = "orig.ident")
cd8_cells <- FindNeighbors(cd8_cells, reduction = "harmony", dims = 1:20)
cd8_cells <- FindClusters(cd8_cells, resolution = 0.8)
cd8_cells <- RunUMAP(cd8_cells, reduction = "harmony", dims = 1:20)

VlnPlot(cd8_cells, features = c("CCR7", "IL7R", "GZMK", "GZMB", "LAG3", "PDCD1", "HAVCR2", "TOX"), group.by = "seurat_clusters")

###############################################3

table(mac_cells$seurat_clusters)
mac_cells$mac_subtype <- dplyr::recode(as.character(mac_cells$seurat_clusters),
                                       "0" = "Macro_C1QC",
                                       "1" = "Macro_C1QC",
                                       "2" = "Macro_SPP1",
                                       "3" = "Macro_C1QC",
                                       "4" = "Macro_SPP1")

table(mac_cells$mac_subtype)   

# Merge fine labels back into the main object
combined$cell_type_fine <- as.character(combined$cell_type)     
combined$cell_type_fine[colnames(mac_cells)] <- mac_cells$mac_subtype
combined$cell_type_fine[colnames(cd8_cells)] <- cd8_cells$cd8_subtype


#############################################3

# BiocManager::install("AUCell")
library(AUCell)

gene_sets <- list(
  Exhaustion    = c("LAG3", "PDCD1", "HAVCR2", "CTLA4", "TOX", "TIGIT"),
  Cytotoxicity  = c("GZMB", "GZMA", "PRF1", "NKG7", "GNLY"),
  M2_macrophage = c("SPP1", "C1QA", "C1QB", "C1QC", "APOE", "MRC1")
)
# تأكدي الأول من الجينات المفقودة (لو فيه)
lapply(gene_sets, function(g) g[!g %in% rownames(combined)])
expr_matrix <- GetAssayData(combined, assay = "RNA", layer = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)

combined <- AddMetaData(combined, metadata = as.data.frame(t(getAUC(cells_AUC))))
VlnPlot(combined, features = c("Exhaustion", "Cytotoxicity", "M2_macrophage"), group.by = "cell_type_fine")
####################################################
library(CellChat)

cellchat <- createCellChat(object = combined, group.by = "cell_type_fine")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

netVisual_bubble(cellchat, sources.use = c("Macro_SPP1", "Macro_C1QC"),
                 targets.use = "CD8_Tex", 
                 signaling = c("MIF", "GALECTIN"), 
                 remove.isolate = FALSE)

##########################################
# BiocManager::install("monocle")
library(monocle3)

cd8_expr <- as.matrix(GetAssayData(cd8_cells, assay = "RNA", layer = "counts"))
cell_metadata <- cd8_cells@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(cd8_expr), row.names = rownames(cd8_expr))

cds <- new_cell_data_set(cd8_expr,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds, root_cells = colnames(cds)[cd8_cells$cd8_subtype == "CD8_Tnaive"][1])

plot_cells(cds, color_cells_by = "cd8_subtype", label_groups_by_cluster = FALSE)
plot_cells(cds, genes = c("LAG3", "PDCD1", "HAVCR2"), color_cells_by = "pseudotime")
