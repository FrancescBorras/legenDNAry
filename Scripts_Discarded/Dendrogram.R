# PR Tree species dendrogram

# Welcome to my attempt to visualize how the different tree species are related using their P6-loop sequenced data.
# I would expect species within the same family/order to be grouped together. Will this be the case? LET'S FIND OUT


# Load required libraries (for now I don't need it tho)
library(ggplot2)


# Read the data
data <- read.csv("PRreflib_haplotype_data2.csv", stringsAsFactors = FALSE, sep=";")


# Extract the species names and their most abundant OTUs and remove the extension from their name
data$Column1 <- gsub("_F_filt.fastq.gz", "", data$Column1)
species <- data$Column1
most_abundant_otu <- apply(data[, -1], 1, function(x) which.max(x))


# Create a matrix of the most abundant OTUs
otu_matrix <- t(sapply(data[, -1], function(x) {
  max_otu <- which.max(x)
  x / x[max_otu]
}))


# Compute distance matrix based on Euclidean distance
dist_matrix <- dist(otu_matrix, method = "euclidean")


# Hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "complete") #There's several method options, type hclust in help for more info


# Plot the dendrogram
plot(hclust_result, hang = -1, labels = F, main = "Species Dendrogram")


# Quick check on what the clustering is informing us about 
str(hclust_result)











          # WORKING ON IT Create a circular dendrogram
library(circlize)
par(mar = c(1, 1, 1, 1))  # Adjust margin for better plotting
dend <- as.dendrogram(hclust_result)
circos.dendrogram(dend, labels = species, main = "Circular Dendrogram")


install.packages("dendextend")
library(dendextend)
library(circlize)



# And here is the second attempt at making this 

# Load required libraries
library(ggplot2)

# Read the data
data <- read.csv("PRreflib_haplotype_data2.csv", stringsAsFactors = FALSE, sep=";")

# Extract the species names and their most abundant OTUs and remove the extension from their name
data$Column1 <- gsub("_F_filt.fastq.gz", "", data$Column1)
species <- data$Column1
most_abundant_otu <- apply(data[, -1], 1, function(x) which.max(x))

# Create a matrix of the most abundant OTUs
otu_matrix <- t(sapply(data[, -1], function(x) {
  max_otu <- which.max(x)
  x / x[max_otu]
}))

# Compute distance matrix based on Euclidean distance
dist_matrix <- dist(otu_matrix, method = "euclidean")



# Hierarchical clustering
require(stats)
hclust_result <- hclust(dist_matrix, method = "complete")

# Plot the dendrogram (DOESN'T WORK PROPERLY)
plot(x = hclust_result)
require(factoextra)
require(vctrs)
fviz_dend(x = hclust_result, cex = 0.5, lwd = 0.7)


plot(hclust_result, hang = -1, labels = F, main = "Species Dendrogram")

dend <- as.dendrogram(hclust_result(gdist1, method = "complete"))




