library(vegan)

libfilts<- readRDS("Processed_data/Biggrid_libsize_filter.RDA")
library(vegan)

# Define a function to create rarefaction plots
create_rarefaction_plot <- function(phylo_obj, col, min_reads, xlim) {
  # Convert the OTU table to a matrix
  ttab <- otu_table(phylo_obj)
  class(ttab) <- "matrix"
  
  # Create the rarefaction plot
  rarecurve(ttab, step = 10, col = col, label = FALSE, xlim = xlim)
  
  # Add a vertical line for the minimum reads
  abline(v = min_reads, col = "red", lty = 2)
}

# Set up a 2x2 panel layout
par(mfrow = c(2, 2))

# Define thresholds for the titles
thresholds <- c(10, 40, 70, 200)

# Iterate over the libfilts list
for (i in seq_along(libfilts)) {
  # Extract the phyloseq object
  phylo_obj <- libfilts[[i]]
  
  # Calculate the minimum reads and number of samples
  min_reads <- min(sample_sums(phylo_obj))
  num_samples <- nsamples(phylo_obj)
  
  # Create the rarefaction plot
  create_rarefaction_plot(phylo_obj, col = "blue", min_reads = min_reads, xlim = c(0, 200000))
  
  # Add a title with the threshold and number of samples
  title(paste(">= ", thresholds[i], "k Reads (", num_samples, " samples)", sep = ""))
}
