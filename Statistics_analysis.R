# Genomic Analysis in golub data
# Author: Shalini 

# Load required packages
library(multtest)   # For loading the dataset
library(dplyr)      # For data manipulation
library(ggplot2)    # For data visualization
library(cluster)    # For clustering analysis
library(corrplot)   # For correlation plots

### Data Loading and Preprocessing ###

# Load the Golub dataset
data(golub, package = "multtest")
golub_gnames <- golub.gnames[, 2]

### Problem 1: Clustering Analysis ###

# Select specific genes for clustering analysis
CCND3 <- grep("CCND3 Cyclin D3", golub_gnames)
CCND3_data <- golub[CCND3, ]

# Perform hierarchical clustering
hc_single <- hclust(dist(CCND3_data))
hc_ward <- hclust(dist(CCND3_data), method = "ward.D2")

# Plot dendrograms
plot(hc_single, main = "Single Linkage Hierarchical Clustering", xlab = "Patients", ylab = "Distance")
plot(hc_ward, main = "Ward Linkage Hierarchical Clustering", xlab = "Patients", ylab = "Distance")

# Perform k-means clustering
k <- 2  # Number of clusters
kmeans_clusters <- kmeans(CCND3_data, centers = k)

# Visualize k-means clusters
clusplot(CCND3_data, kmeans_clusters$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)

### Problem 2: Correlation Analysis ###

# Correlation matrix of all genes in the dataset
cor_matrix <- cor(golub)

# Plot correlation matrix using corrplot
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", addrect = 2)

### Problem 3: Differential Expression Analysis ###

# Perform differential expression analysis
group <- golub.cl
gene_expression <- golub
de_results <- apply(gene_expression, 1, function(x) t.test(x ~ group)$p.value)

# Adjust p-values for multiple testing using FDR
de_results_adj <- p.adjust(de_results, method = "fdr")

# Identify significant genes
significant_genes <- which(de_results_adj < 0.05)
significant_genes_names <- golub_gnames[significant_genes]

# Output significant genes
cat("Significant genes:", significant_genes_names, "\n")

### Visualization 

# Boxplot of gene expression across patient groups
 ggplot(data = golub, aes(x = golub.cl, y = golub[1, ])) +
   geom_boxplot() +
   labs(title = "Gene Expression by Patient Group", x = "Group", y = "Expression")

# Save plots
 ggsave("gene_expression_boxplot.png")