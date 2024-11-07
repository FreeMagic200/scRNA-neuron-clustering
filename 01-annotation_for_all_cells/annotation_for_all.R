# Load required libraries
library(Seurat)
library(ggplot2)
library(stringr)
library(future)
library(patchwork)
library(clustree)
library(dplyr)
library(scCustomize)
library(Nebulosa)

# Set random seed and future plan for parallel processing
set.seed(42)
plan("sequential", workers = 1)
options(future.globals.maxSize = Inf, future.seed = TRUE)

# Read preprocessed data
obj <- readRDS("./umap_tunned_103530_cells.rds")

# 1. Clustering Analysis ---------------------------------------------------
# Find neighbors and perform clustering at different resolutions
obj <- FindNeighbors(obj, dims = 1:12, reduction = "integrated.jointpca", 
                     k.param = 5, annoy.metric = "euclidean")
obj <- FindClusters(obj, method = "igraph", algorithm = 4,
                    resolution = seq(0, 2, 0.1), verbose = TRUE)

# Visualize clustering tree
k5_tree <- clustree(obj, prefix = "RNA_snn_res.")

# Analyze cluster stability
plot_data <- k5_tree$data %>%
  group_by(RNA_snn_res.) %>%
  summarise(
    avg_stability = mean(sc3_stability, na.rm = TRUE),
    cluster_count = n_distinct(cluster)
  ) %>%
  mutate(resolution = as.numeric(as.character(RNA_snn_res.)))

# Plot stability metrics
stability_plot <- ggplot(plot_data, aes(x = resolution, y = avg_stability)) +
  geom_line() +
  geom_point() +
  labs(x = "Resolution", y = "Average SC3 Stability",
       title = "Average SC3 Stability vs. Resolution") +
  theme_minimal() +
  scale_x_continuous(breaks = plot_data$resolution) +
  geom_text(aes(label = cluster_count), vjust = -0.5, size = 3)

# Save clustering analysis plots
ggsave("k5_clustree_indent_1.png", plot = k5_tree, width = 12, height = 14, dpi = 200)

# 2. Find Marker Genes ---------------------------------------------------
obj <- JoinLayers(obj)
scRNA.markers <- FindAllMarkers(obj,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.6,
                               slot = "data")

# Process markers
all.markers <- scRNA.markers %>%
  dplyr::select(gene, everything()) %>%
  dplyr::filter(avg_log2FC > 0, p_val < 0.05)

# Get top 50 markers per cluster
top50 <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

write.csv(top50, file = "./scRNA_top50_markers.csv", row.names = FALSE)

# 3. Cell Type Annotation ---------------------------------------------------
# Define cell type mapping
jointpca_cell_type_mapping <- list(
  GABA = c(1,2,7,8,9,13,21,23,25),
  GLU = c(3,4,5,6,14,16,18,24,33,35),
  IPC1 = c(12,20),
  IPC2 = c(28,29),
  RGC = c(11,15,26,27,31),
  AS = c(17,22),
  TC = c(30),
  Mig = 34,
  OPC = 19,
  low_quality1 = 10,
  low_quality2 = 32,
  Ependymal = 37,
  OD = 36
)

# Assign cell types
obj$jointpca_cell_type <- NA
for (i in seq_along(jointpca_cell_type_mapping)) {
  clusters <- jointpca_cell_type_mapping[[i]]
  cell_type <- names(jointpca_cell_type_mapping)[i]
  mask <- obj$jointpca.clust %in% clusters & is.na(obj$jointpca_cell_type)
  if (any(mask)) {
    obj$jointpca_cell_type[mask] <- cell_type
  }
}

# 4. Quality Control and Filtering ---------------------------------------------------
# Remove low quality cells
obj <- subset(obj, subset = jointpca_cell_type != "low_quality1" & 
                          jointpca_cell_type != "low_quality2")

# Add developmental stage information
obj$stage <- sub("_.*", "", obj$Batch)

# 5. Visualization ---------------------------------------------------
# UMAP plots
jointpca_cell_type_umap <- DimPlot(obj,
                                  reduction = "jointpca.umap",
                                  group.by = "jointpca_cell_type",
                                  raster = FALSE,
                                  label = TRUE) + 
                          ggtitle("jointpca celltype")

# Generate stage-specific plots
stage_plot_list <- list()
for (dev_stage in unique(obj$stage)) {
  plot <- DimPlot(obj,
                  reduction = "jointpca.umap",
                  group.by = "jointpca_cell_type",
                  raster = FALSE,
                  label = TRUE,
                  label.size = 5,
                  cells.highlight = WhichCells(obj, expression = stage == dev_stage),
                  cols.highlight = "red",
                  cols = "gray") + 
    ggtitle(paste("jointpca stage -", dev_stage)) +
    NoLegend()
  
  stage_plot_list[[dev_stage]] <- plot
}

# Combine and save stage plots
stage_combined_plot <- wrap_plots(stage_plot_list, ncol = 4)
stage_final_plot <- stage_combined_plot + 
  plot_annotation(title = "jointpca stages",
                 theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))

# Save final objects and plots
saveRDS(obj, "./filtered_jointpca_annotated_99252_cells.rds")
ggsave("jointpca_cell_type_umap.png", jointpca_cell_type_umap, 
       width = 18.15, height = 14.34, units = "cm")
ggsave("jointpca_stages.png", stage_final_plot, 
       width = 24, height = 12, limitsize = FALSE)