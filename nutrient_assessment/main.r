library('dplyr')
library('janitor')
library('pheatmap')

data <- read.csv('data/log2_nutrient_mean.csv')

# Ward Hierarchical Clustering
rownames(data) <- data$Sample.Group
data$Sample.Group <- NULL
# data_t <- t(data)
# data_t %>%
#   col_to_names(col_number = 1)
data <- mutate_all(data, function(x) as.numeric(as.character(x)))
pheatmap(as.matrix(t(data)),
    cutree_rows=5,
    cluster_rows=T,
    cluster_cols=F,
    clustering_method="ward.D2",
    clustering_distance_rows="manhattan")
# d <- dist(t(data), method = "euclidean") # distance matrix
# print(d)
# fit <- hclust(d, method="ward")
# groups <- cutree(fit, k=5) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters
# rect.hclust(fit, k=5, border="red")