library('dplyr')
library('pheatmap')


# Read and process CSV for plotting on cluster map
data <- read.csv('data/log2_nutrient_mean.csv')
# Index Sample.Group as row names and drop Sample.Group column.
rownames(data) <- data$Sample.Group
data$Sample.Group <- NULL
data <- mutate_all(data, function(x) as.numeric(as.character(x)))

# Plot heatmap using ward.D2 clustering method employing Manhattan clustering by row
pheatmap(as.matrix(t(data)),
    cutree_rows=8,
    cluster_rows=T,
    cluster_cols=F,
    clustering_method="ward.D2",
    clustering_distance_rows="manhattan",
    display_numbers=T,
    main="Pertubations in the yeast metabolic profile\nas a result of different nutrient compositions"
    )
