##read in data for plots and summaries

prepost<- readRDS("Processed_data/Prereflibcut_Postreflibcut.RDA")

precutsums <- rowSums(otu_table(prepost$precut))
postcutsums <- rowSums(otu_table(prepost$postcut))

##plots of sequences remaining after filtering for only reference library OTUs
## total sequences per sample before V after
plot(precutsums, postcutsums, xlab ="All Plants OTUs - library size", ylab ="Only Ref. Lib. OTUs - library size ")
abline(0,1, col="red")
## proportion remaining for each sample sequence count
boxplot(postcutsums/precutsums, ylab = "Prop. soil reads matching reference library taxa")

summary(postcutsums/precutsums)
sd(postcutsums/precutsums)

## total OTUs  per sample
precutotusums <- apply(otu_table(prepost$precut), MARGIN=1, FUN=function(x) {length( x[x > 0] )} ) 
postcutotusums <- apply(otu_table(prepost$postcut), MARGIN=1, FUN=function(x) {length( x[x > 0] )} ) 


plot(precutotusums, postcutotusums, xlab ="All Plants OTUs count", ylab ="Only Ref. Lib. OTUs count ")
abline(0,1, col="red")

## proportion OTUs remaining in each sample after filtering for what is in reference library
boxplot(postcutotusums/precutotusums, ylab = "Prop. soil OTUs matching reference library taxa")

plot(postcutsums/precutsums, postcutotusums/precutotusums, xlab ="Prop. reads remaining after ref. lib. filtering" , ylab ="Prop. OTUs remaining after ref. lib. filtering")
## reads
summary(postcutsums/precutsums)
## OTUs
summary(postcutotusums/precutotusums)
# so in general 91-99% (range of mean & median) of reads are from the reference library taxa 
# equating to 65 - 66% of OTUs (from the reference library taxa)
size <- rowSums(otu_table(prepost$precut))
finsum <- as.data.frame(cbind(postcutsums/precutsums, postcutotusums/precutotusums, size))
colnames(finsum) <- c("libprop", "otuprop", "libsize")

### making a final summary plot of the proportion of soil reads and OTUs that reference library taxa
# Load necessary libraries
library(ggplot2)
library(patchwork)

# Prepare data for the combined boxplot
combined_data <- data.frame(
  value = c(postcutsums / precutsums, postcutotusums / precutotusums),
  type = rep(c("Reads", "OTUs"), each = length(postcutsums))
)

# Combined boxplot for Reads and OTUs (no legend)
plot_a <- ggplot(combined_data, aes(x = type, y = value, fill = type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Reads" = "blue", "OTUs" = "orange")) + # Custom box colors
  labs(
    y = "Prop. matching reference library taxa",
    x = NULL,
    fill = NULL # Removes the legend title
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title.y = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = "none", # Remove legend
    plot.margin = margin(5, 5, 2, 5) # Reduce top margin
  ) +
  annotate("text", x = 0.5, y = max(combined_data$value) + 0.05, label = "A", size = 5, fontface = "bold", hjust = 0)

# Scatter plot with default vertical legend
plot_b <- ggplot(finsum, aes(x = libprop, y = otuprop, color = libsize)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(5)) + # Custom color gradient
  labs(
    x = "Prop. soil reads matching reference library",
    y = "Prop. soil OTUs matching reference library",
    colour = "Seq. Depth" # Legend label
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "right", # Default vertical legend on the right
    plot.margin = margin(5, 5, 2, 5) # Reduce top margin
  ) +
  annotate("text", x = min(finsum$libprop), y = max(finsum$otuprop) + 0.05, label = "B", size = 5, fontface = "bold", hjust = 0)

# Combine the combined boxplot and scatter plot into a single row
combined_plot <- plot_a + plot_b + plot_layout(ncol = 2)
combined_plot 
# Save the combined plot as a PDF
pdf("Figures/FigureS8.pdf", width = 17 / 2.54, height = 8 / 2.54) # Convert cm to inches
print(combined_plot)
dev.off()

## So most reads and about 60% of OTUs remain when only using the reference library matches
## and there is no clear relationship here with individual sample sequencing depth (the colours)
