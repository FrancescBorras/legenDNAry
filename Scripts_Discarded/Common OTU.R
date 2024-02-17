# PR Tree Species common OTU
# Some species share their most common OTU, with the following script I'll attempt to visualize it


# Read data
data <- read.csv("PRreflib_haplotype_data2.csv", stringsAsFactors = FALSE, sep=";")


# Extract the species names (without their extension) and their most abundant OTUs 
data$Column1 <- gsub("_F_filt.fastq.gz", "", data$Column1)
species <- data$Column1
most_abundant_otu <- apply(data[, -1], 1, function(x) which.max(x))


# Identify species that share the same most abundant OTU
shared_otu_counts <- table(most_abundant_otu)
shared_otu_species <- names(shared_otu_counts)[shared_otu_counts > 1]


# Visualize the results
for (otu in shared_otu_species) {
  species <- names(which(most_abundant_otu == otu))
  cat(paste("Species sharing most abundant OTU", otu, ": ", species, "\n"))
}




# SECOND ATTEMPT (Not working atm...)

# Replace 'your_data.csv' with the actual file name and path
data <- read.csv("PRreflib_haplotype_data2.csv", header = TRUE, sep = ";")

# Assuming 'data' is your data frame with columns 'Species' and OTUs

# Identify Most Abundant OTU for Each Species
most_abundant_otu_index <- apply(data[, -1, drop = FALSE], 1, which.max)

# Create a new data frame with results
result_df <- data.frame(
  Species = data$Species,
  MostAbundantOTU = colnames(data)[-1][most_abundant_otu_index],
  stringsAsFactors = FALSE  # Ensure character columns are not converted to factors
)

# Remove rows with NA values in MostAbundantOTU
result_df <- result_df[complete.cases(result_df), ]

# View Results
print(result_df)



