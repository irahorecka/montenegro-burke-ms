library('optparse')
library('dplyr')
library('pheatmap')


# Get external arguments passed by the command line
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
                help="dataset file name", metavar="character")); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check valid file path is provided
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}


# Read and process CSV for plotting on cluster map
data <- read.csv(opt$file)
# Index Sample.Group as row names and drop Sample.Group column.
rownames(data) <- data$Sample.Group
data$Sample.Group <- NULL
data <- mutate_all(data, function(x) as.numeric(as.character(x)))

# Plot heatmap using ward.D2 clustering method employing Manhattan clustering by row
pheatmap(as.matrix(t(data)),
    cutree_rows=5,
    cluster_rows=T,
    cluster_cols=F,
    clustering_method="ward.D2",
    clustering_distance_rows="manhattan",
    display_numbers=T,
    main="Upregulated metabolites in yeast\nas a result of varying nutrient conditions"
    )
