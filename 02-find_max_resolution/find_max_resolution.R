# Load required libraries
library(Seurat)  
library(dplyr)   
library(ggplot2) 
library(SeuratWrappers) 
library(future)
library(future.apply) # For applying functions in parallel
set.seed(42)    # Set seed for reproducibility

# Set up parallel processing plan
plan("multicore", workers = 15, gc = T) 

# Increase memory limits for future processing
options(
  future.globals.maxSize = 60128960000, 
  future.seed = TRUE, 
  future.stdout = TRUE,  
  future.rng.onMisuse = "ignore"  
)

# Load reticulate for Python integration
library(reticulate)

# Set Python environment for reticulate
reticulate::use_python("/opt/miniforge3/envs/rstudio_reticulate/bin/python", required = TRUE)

# Display Python configuration
py_config()

# Set output directory for saving results
OUTPUT_DIR = "/mnt/data/projects/sub_clustering/find_optimized_neu_res/"

# Read in the processed Seurat object
all_data <- readRDS(
    "/mnt/data/projects/scRNA-0722-Reviewed/data/processed/seurat/30d/jointpca_annotated_103530_cells.rds"
)
dim(all_data)  # Check dimensions of the data

# Subset data to include only neurons (GABA and GLU types)
neurons <- subset(all_data, subset = jointpca_cell_type %in% c("GABA", "GLU"))
table(neurons$jointpca_cell_type)  # Count the number of cells per type

# Join layers in the Seurat object
neurons <- JoinLayers(neurons)

# Garbage collection to free up memory
gc()

# Extract the expression matrix from the Seurat object
expr_matrix <- GetAssayData(neurons, slot = "counts")

# Get cell names and gene names
cell_names <- colnames(neurons)
gene_names <- rownames(neurons)

# Extract relevant metadata from the Seurat object
metadata <- neurons@meta.data[, c("nCount_RNA", "nFeature_RNA", "Batch", 
                                  "percent.mt", "percent.hb", "percent.ribo", "percent.Malat1",
                                  "S.Score", "G2M.Score", "Phase",
                                  "jointpca_cell_type", "jointpca.clust")]

# Ensure that row names of metadata match the column names of expression matrix
rownames(metadata) <- cell_names

# Check dimensions of extracted data
print(dim(expr_matrix))
print(head(metadata))

# Save the expression matrix to a file (using sparse matrix format)
library(Matrix)  # For sparse matrix functionality
sparse_expr_matrix <- as(expr_matrix, "dgCMatrix")
saveRDS(sparse_expr_matrix, file = paste0(OUTPUT_DIR, "neurons_expression_matrix.rds"))

# Save metadata to a CSV file
write.csv(metadata, file = paste0(OUTPUT_DIR, "neurons_metadata.csv"))

# Save cell names to a text file
writeLines(cell_names, paste0(OUTPUT_DIR, "neurons_cell_names.txt"))

# Save gene names to a text file
writeLines(gene_names, paste0(OUTPUT_DIR, "neurons_gene_names.txt"))

# Create a new Seurat object with the sparse expression matrix and metadata
neurons_trimmed <- CreateSeuratObject(counts = sparse_expr_matrix, 
                                      meta.data = metadata,
                                      project = "neurons")

# Check if the new Seurat object is valid
print(validObject(neurons_trimmed))

# Save the trimmed Seurat object to a file
saveRDS(neurons_trimmed, file = paste0(OUTPUT_DIR, "neurons_trimmed.rds"))

# Normalize the data using LogNormalize method
neurons_trimmed <- NormalizeData(neurons_trimmed, normalization.method = "LogNormalize")

# Identify highly variable features (genes)
neurons_trimmed <- FindVariableFeatures(neurons_trimmed, selection.method = "vst", nfeatures = 2000)

# Scale the data and regress out unwanted sources of variation
neurons_trimmed <- ScaleData(neurons_trimmed, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "Malat1", "S.Score", "G2M.Score"))

# Run PCA on the scaled data
neurons_trimmed <- RunPCA(neurons_trimmed, layer = "scale.data", verbose = F)

# Generate an elbow plot to determine the optimal number of principal components
ElbowPlot(neurons_trimmed)

# Save the object after regression
saveRDS(neurons_trimmed, file = paste0(OUTPUT_DIR, "neurons_trimmed_regressed.rds"))

# Split RNA data by Batch for integration
neurons_trimmed[["RNA"]] <- split(neurons_trimmed[["RNA"]], f = neurons_trimmed$Batch)

# Integrate layers using Joint PCA Integration
neurons_trimmed <- IntegrateLayers(
    neurons_trimmed,
    method = JointPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.jointpca",
    verbose = F
)

# Run UMAP for dimensionality reduction using the integrated data
neurons_trimmed <- RunUMAP(neurons_trimmed, 
                           reduction = "integrated.jointpca",
                           dims = 1:12,
                           n.neighbors = 60,
                           reduction.name = "integrated_regressed_umap",
                           verbose = F)

# Plot features on the UMAP
FeaturePlot(neurons_trimmed, reduction = "integrated_regressed_umap", features = c("Slc17a6", "Slc32a1"))

# Find neighbors for clustering
neurons_trimmed <- FindNeighbors(neurons_trimmed, dims = 1:12, reduction = "integrated.jointpca", verbose = F)

neurons_trimmed<-FindClusters(neurons_trimmed,resolution = c(0.1,seq(1,50,1)),algorithm = 4,method = "igraph",verbose = F)

saveRDS(neurons_trimmed, file = paste0(OUTPUT_DIR, "neurons_trimmed_max_clustered_interval_1.rds"))

# Define a function to process the data for each resolution
process_resolution <- function(seurat_obj, resolution) {
  # Get the clustering results for the specified resolution
  cluster_col <- paste0("RNA_snn_res.", resolution)
  
  # Create a cross-table using 'stage' instead of 'Batch'
  cross_table <- table(seurat_obj@meta.data[[cluster_col]], seurat_obj@meta.data$stage)
  
  # Calculate the number of stages in each cluster (i.e., non-zero values in each row)
  cluster_stage_counts <- apply(cross_table, 1, function(x) sum(x > 0))
  
  # Find the cluster with the minimum number of associated stages
  min_stage_count <- min(cluster_stage_counts)
  
  # Return a data frame with resolution and the minimum stage count
  return(data.frame(resolution = resolution, min_stage_count = min_stage_count))
}

# Define a sequence of resolutions to test
resolutions <- c(0.1, seq(1, 50, 1))  # Resolutions from 0.1 and then 1 through 50

# Apply the function to all specified resolutions
results <- lapply(resolutions, function(res) process_resolution(neurons_trimmed, res))

# Combine the results into a single data frame
all_results <- do.call(rbind, results)

# Create a line plot showing minimum stage count across resolutions
p <- ggplot(all_results, aes(x = resolution, y = min_stage_count, group = 1)) +
  geom_line() +  # Add a line for the resolution vs. min stage count
  geom_point() +  # Add points to highlight individual data points
  scale_x_continuous(breaks = resolutions) +  # Set x-axis breaks to the specified resolutions
  scale_y_continuous(breaks = seq(min(all_results$min_stage_count), max(all_results$min_stage_count), by = 1)) +  # Set y-axis breaks
  labs(title = "Minimum Stage Count Across Resolutions (0.1 + x, x = 1 to 50)",
       x = "Resolution",
       y = "Minimum Stage Count") +
  theme_minimal() +  # Use a minimal theme for the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Save the plot as a PNG image
ggsave(filename = paste0(OUTPUT_DIR, "min_stage_count_across_resolutions_0.1_plus_x.png"),
       plot = p,
       width = 12,  # Width of the image in inches
       height = 8,  # Height of the image in inches
       dpi = 300)  # Set the resolution of the image to 300 DPI    
                 
analyze_resolution <- function(resolution) {
  cat("\n分析分辨率:", resolution, "\n")
  
  # 设置分辨率
  Idents(neurons_trimmed) <- paste0("RNA_snn_res.", resolution)
  
  # 执行FindAllMarkers
  all_markers <- FindAllMarkers(neurons_trimmed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 保存所有markers（筛选前）
  write.csv(all_markers, file = paste0(OUTPUT_DIR, "res", resolution, "_all_markers_unfiltered.csv"), row.names = FALSE)
  
  # 筛选markers
  filtered_markers <- all_markers %>%
    filter(pct.1 > 0.25) %>%
    filter(p_val_adj < 0.05) %>%
    filter(avg_log2FC > 0.5)
  
  # 保存筛选后的markers
  write.csv(filtered_markers, file = paste0(OUTPUT_DIR, "res", resolution, "_all_markers_filtered.csv"), row.names = FALSE)
  
  # 打印一些统计信息
  cat("筛选前的marker基因总数：", nrow(all_markers), "\n")
  cat("筛选后的marker基因总数：", nrow(filtered_markers), "\n")
  
  # 打印每个聚类的marker基因数量
  cluster_marker_counts <- filtered_markers %>%
    group_by(cluster) %>%
    summarise(marker_count = n()) %>%
    arrange(desc(marker_count))
  
  print(cluster_marker_counts)
  
  # 找出marker基因数量最少的聚类
  min_marker_cluster <- cluster_marker_counts %>%
    filter(marker_count == min(marker_count))
  
  cat("marker基因数量最少的聚类：", min_marker_cluster$cluster, "\n")
  cat("最少的marker基因数量：", min_marker_cluster$marker_count, "\n")
  
  # 计算每个聚类的细胞数
  cluster_sizes <- table(Idents(neurons_trimmed))
  cluster_sizes_df <- as.data.frame(cluster_sizes)
  colnames(cluster_sizes_df) <- c("cluster", "cell_count")
  cluster_sizes_df <- cluster_sizes_df %>% arrange(desc(cell_count))
  
  print(cluster_sizes_df)
  
  # 找出细胞数量最少的聚类
  min_cell_cluster <- cluster_sizes_df %>%
    filter(cell_count == min(cell_count))
  
  cat("细胞数量最少的聚类：", min_cell_cluster$cluster, "\n")
  cat("最少的细胞数量：", min_cell_cluster$cell_count, "\n")
}

# 对所有指定的分辨率进行分析
resolutions <- c(39, 39.2, 39.4, 39.5, 39.6, 39.7, 39.8)

for (res in resolutions) {
  analyze_resolution(res)
}
                  
# 定义函数来处理每个分辨率的数据
process_resolution <- function(resolution) {
  # 读取筛选后的差异基因结果
  markers_filtered <- read.csv(paste0(OUTPUT_DIR, "res", resolution, "_all_markers_filtered.csv"))
  
  # 统计每个类群的DEGs数量
  cluster_deg_counts <- markers_filtered %>%
    group_by(cluster) %>%
    summarise(deg_count = n()) %>%
    arrange(deg_count)
  
  # 找出DEGs数量最少的类群
  min_deg_count <- min(cluster_deg_counts$deg_count)
  
  return(data.frame(resolution = resolution, min_deg_count = min_deg_count))
}

# 定义所有的分辨率
resolutions <- c(39, seq(39.1, 39.8, by = 0.1))

# 应用函数到所有分辨率
results <- lapply(resolutions, process_resolution)

# 合并结果
all_results <- do.call(rbind, results)

# 创建折线图
p <- ggplot(all_results, aes(x = resolution, y = min_deg_count, group = 1)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = resolutions) +
  labs(title = "Minimum DEG Count Across Resolutions (39 - 39.8)",
       x = "Resolution",
       y = "Minimum DEG Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存图片
ggsave(filename = paste0(OUTPUT_DIR, "min_deg_count_across_resolutions_39_to_39.8.png"),
       plot = p,
       width = 12,
       height = 8,
       dpi = 300)
                  
                  
# 获取整个数据
neu_jointpca_emb <- neurons_trimmed@reductions$integrated.jointpca@cell.embeddings

# 保存为CSV文件
write.csv(neu_jointpca_emb, "neu_jointpca_emb.csv", row.names = TRUE)
