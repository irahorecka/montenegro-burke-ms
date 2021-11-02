library('dplyr')
library('janitor')
library('pheatmap')


# Read and process CSV for plotting on cluster map
data <- read.csv('data/log2_nutrient_mean.csv')
# Ward Hierarchical Clustering
rownames(data) <- data$Sample.Group
data$Sample.Group <- NULL
data <- mutate_all(data, function(x) as.numeric(as.character(x)))

# Plot heatmap using ward.D2 clustering method employing Manhattan clustering by row
pheatmap(as.matrix(t(data)),
    cutree_rows=5,
    cluster_rows=T,
    cluster_cols=F,
    clustering_method="ward.D2",
    clustering_distance_rows="manhattan")
