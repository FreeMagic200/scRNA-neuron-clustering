import os
import dask.dataframe as dd
import cuml
from cuml.cluster import KMeans
from cuml.metrics.cluster import silhouette_score, adjusted_rand_score, homogeneity_score, completeness_score, mutual_info_score, v_measure_score
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from cuml.datasets import make_blobs

# 定义输入目录和文件名
INPUT_DIR = "/mnt/data/projects/sub_clustering/find_optimized_neu_res/"
file_path = os.path.join(INPUT_DIR, "neu_jointpca_emb.csv")

# 使用 Dask 读取 CSV 文件
cell_emb = dd.read_csv(file_path)

# 将第一列设置为行索引
cell_emb = cell_emb.set_index('Unnamed: 0')  # 假设第一列的列名为 'Unnamed: 0'

cell_emb.head()

# 读取 metadata
metadata_path = os.path.join(INPUT_DIR, "neurons_scaled_metadata.csv")
metadata = dd.read_csv(metadata_path, assume_missing=True)

# 将第一列设置为行索引
metadata = metadata.set_index('Unnamed: 0')  # 假设第一列的列名为 'Unnamed: 0'

metadata.head()

# 假设你的 DataFrame 是 metadata
na_counts = metadata.isnull().sum().compute()
print(na_counts)

# 初始化结果列表
results = []

# 计算每个分辨率的聚类指标
for resolution in range(1, 51):  # 假设分辨率范围是2到50
    # 获取当前分辨率的聚类标签
    labels = metadata[f'RNA_snn_res.{resolution}'].compute().to_numpy()
    
    # 计算聚类指标
    silhouette = silhouette_score(cell_emb, labels)
    
    # 将结果保存到列表
    results.append((resolution, silhouette))
    
    # 打印结果
    print(f"Resolution: {resolution}, Silhouette: {silhouette}")

# 可选：打印所有结果
print(results)

# 可视化结果
resolutions, silhouettes = zip(*results)  # 解压结果列表

plt.figure(figsize=(10, 6))
plt.plot(resolutions, silhouettes, marker='o')
plt.title('Silhouette Score vs Resolution')
plt.xlabel('Resolution')
plt.ylabel('Silhouette Score')
plt.xticks(resolutions)  # 设置 x 轴刻度
plt.grid()
plt.show()

# 初始化结果列表
x_results = []

# 计算每个分辨率的聚类指标
previous_labels = None  # 用于存储前一个分辨率的标签

for resolution in range(1, 51):  # 假设分辨率范围是2到50
    # 获取当前分辨率的聚类标签
    labels = metadata[f'RNA_snn_res.{resolution}'].compute().to_numpy()
    
    # 将 labels 转换为 int32 类型
    labels_int = labels.astype(np.int32)  # 确保 labels 是 int32 类型
    
    if previous_labels is not None:
        # 使用前一个分辨率的标签作为 y_true
        y_true = previous_labels.astype(np.int32)  # 确保 y_true 是 int32 类型
        
        # 计算聚类指标
        adjusted_rand = adjusted_rand_score(y_true, labels_int)
        homogeneity = homogeneity_score(y_true, labels_int)
        completeness = completeness_score(y_true, labels_int)
        mutual_info = mutual_info_score(y_true, labels_int)
        v_measure = v_measure_score(y_true, labels_int)
        
        # 将结果保存到 x_results 列表
        x_results.append({
            'resolution': resolution,
            'adjusted_rand': adjusted_rand,
            'homogeneity': homogeneity,
            'completeness': completeness,
            'mutual_info': mutual_info,
            'v_measure': v_measure
        })
        
        # 打印结果
        print(f"Resolution: {resolution}, Adjusted Rand Index: {adjusted_rand}, Homogeneity: {homogeneity}, Completeness: {completeness}, Mutual Information: {mutual_info}, V-measure: {v_measure}")
    
    # 更新前一个分辨率的标签
    previous_labels = labels_int
    
# 假设 x_results 已经计算并填充了数据
# 提取指标
resolutions = [result['resolution'] for result in x_results]
adjusted_rand_scores = [result['adjusted_rand'] for result in x_results]
homogeneity_scores = [result['homogeneity'] for result in x_results]
completeness_scores = [result['completeness'] for result in x_results]
v_measure_scores = [result['v_measure'] for result in x_results]

# 可视化结果
plt.figure(figsize=(12, 8))

plt.plot(resolutions, adjusted_rand_scores, marker='o', label='Adjusted Rand Index')
plt.plot(resolutions, homogeneity_scores, marker='o', label='Homogeneity')
plt.plot(resolutions, completeness_scores, marker='o', label='Completeness')
plt.plot(resolutions, v_measure_scores, marker='o', label='V-measure')

plt.title('Clustering Metrics vs Resolution (Excluding Mutual Information)')
plt.xlabel('Resolution')
plt.ylabel('Score')
plt.xticks(resolutions)  # 设置 x 轴刻度
plt.legend()
plt.grid()
plt.show()