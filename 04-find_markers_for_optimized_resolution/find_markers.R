# Load required libraries
library(Seurat)  
library(dplyr)   
library(ggplot2) 
set.seed(42)    # Set seed for reproducibility

# Set output directory for saving results
OUTPUT_DIR = "/mnt/data/projects/sub_clustering/find_optimized_neu_res/"

neurons_trimmed<-readRDS(paste0(OUTPUT_DIR, "neurons_trimmed_max_clustered_interval_1.rds"))

neurons_trimmed <- JoinLayers(neurons_trimmed)

# 设置分辨率
resolution <- 27

# 设置Seurat对象的标识为指定分辨率的聚类结果
Idents(neurons_trimmed) <- paste0("RNA_snn_res.", resolution)

# 使用默认参数进行差异基因分析
all_markers_27 <- FindAllMarkers(neurons_trimmed, only.pos = F)

# 保存默认参数的差异基因结果
write.csv(all_markers_27, file = paste0(OUTPUT_DIR, "res_27_all_markers.csv"), row.names = FALSE)

----

neu_markers<-read.csv("../find_optimized_neu_res/res_27_all_markers.csv")

# 先过滤掉符合特定模式的行
neu_markers_filtered <- neu_markers %>%
    filter(!grepl("^mt-", gene) & 
           !grepl("^Rp[sl]", gene) & 
           !grepl("^Hb[^(egf)|^(s1l)|^(p1)].+", gene) & 
           !grepl("\\dRik", gene) & 
           !grepl("^Gm\\d", gene))

# 然后筛选marker基因
neu_markers_filtered <- neu_markers_filtered %>%
    filter(pct.1 > 0.2 & p_val_adj < 0.05 & avg_log2FC > 0.6)

head(neu_markers_filtered)

neu_markers_filtered <- neu_markers_filtered %>%
  arrange(desc(avg_log2FC))

head(neu_markers_filtered)

# 按照cluster分组并筛选每个组的avg_log2FC排名前50基因
top50_genes <- neu_markers_filtered %>%
    group_by(cluster) %>%
    top_n(50, avg_log2FC) %>%
    ungroup()

# 查看每个cluster的top 50基因
print(top50_genes)

# 保存top50基因为CSV文件
write.csv(top50_genes, "clust27_top50_genes.csv", row.names = FALSE)
