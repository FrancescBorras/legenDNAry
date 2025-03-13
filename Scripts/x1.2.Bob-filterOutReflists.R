library(phyloseq)
library(microViz)
load("Raw_data/CescPRdata_v2.RData")

## The object with the reflib / Genbank ID indexing is a list called all.taxs.source.otu
## it is larger than the OTU tables oin your phlyoseq objects because it includes all OTUs
## even ones that had no classification from Genbank that were thrown out immediately
## The last column of each data frame in the list is $idsource - here you can see where the OTU was identified from

# View(all.taxs.source.otu)
# View(all.taxs.source.otu$dada.nopool.nochim)
## You can use the OTU code to index and join the data frames.. 

R1Ref.lib.ids <- lapply(R1all.taxs.source.otu, function(x) subset(x, idsource == "RefLib"))
R2Ref.lib.ids <- lapply(R2all.taxs.source.otu, function(x) subset(x, idsource == "RefLib"))

str(R1Ref.lib.ids)
View(R1Ref.lib.ids$dada.pspool.nc.lulu)
View(R2Ref.lib.ids$dada.pspool.nc.lulu)

## Note that OTU34 needs removal - it was from an unidentified "vine" in the plot (i.e. paraphyletic)
## and thus is not in the census data as well as having no identification use at all
## Ref.lib.ids <- lapply(Ref.lib.ids, function(x) x[!(row.names(x) %in% "OTU34"), ])

# View(Ref.lib.ids)
# View(Ref.lib.ids$dada.pspool.nc.lulu)

## To filter your phyloseq objects 
R1reflib.otus <- lapply(R1Ref.lib.ids, function(x) rownames(x))
R2reflib.otus <- lapply(R2Ref.lib.ids, function(x) rownames(x))

## OK so now the names your Taxa found in the soil that was identified in the reference libraries 
## is in the object reflib.otus - a list of 6 charcater vectors for: dada.no.pool, dada.pooled, dada.psuedo pooled, lulu.no.pool, lulu.pooled, lulu.psuedo.pooled

## So you make your new lists with (for example, on the R1phylo.plants.tax object:
## bob I know you would write a loop because you are a wild and crazy guy...
library(phyloseq)
cutlist1 <- prune_taxa(R1reflib.otus[[1]], R1phylo.plants.tax[[1]])
cutlist2 <- prune_taxa(R1reflib.otus[[2]], R1phylo.plants.tax[[2]])
cutlist3 <- prune_taxa(R1reflib.otus[[3]], R1phylo.plants.tax[[3]])
cutlist4 <- prune_taxa(R1reflib.otus[[4]], R1phylo.plants.tax[[4]])
cutlist5 <- prune_taxa(R1reflib.otus[[5]], R1phylo.plants.tax[[5]])
cutlist6 <- prune_taxa(R1reflib.otus[[6]], R1phylo.plants.tax[[6]])

rep2cutlist1 <- prune_taxa(R2reflib.otus[[1]], R2phylo.plants.tax[[1]])
rep2cutlist2 <- prune_taxa(R2reflib.otus[[2]], R2phylo.plants.tax[[2]])
rep2cutlist3 <- prune_taxa(R2reflib.otus[[3]], R2phylo.plants.tax[[3]])
rep2cutlist4 <- prune_taxa(R2reflib.otus[[4]], R2phylo.plants.tax[[4]])
rep2cutlist5 <- prune_taxa(R2reflib.otus[[5]], R2phylo.plants.tax[[5]])
rep2cutlist6 <- prune_taxa(R2reflib.otus[[6]], R2phylo.plants.tax[[6]])

R1ref.lib.list <- list(cutlist1, cutlist2, cutlist3, cutlist4, cutlist5, cutlist6 )
names(R1ref.lib.list) <- c("R1.dada.nopool.reflib", "R1.dada.pool.reflib" , "R1.dada.pspool.reflib", 
                           "R1.lulu.nopool.reflib", "R1.lulu.pool.reflib", "R1.lulu.pspool.reflib" )

R2ref.lib.list <- list(rep2cutlist1, rep2cutlist2, rep2cutlist3, rep2cutlist4, rep2cutlist5, rep2cutlist6 )
names(R2ref.lib.list) <- c("R2.dada.nopool.reflib", "R2.dada.pool.reflib" , "R2.dada.pspool.reflib", 
                           "R2.lulu.nopool.reflib", "R2.lulu.pool.reflib", "R2.lulu.pspool.reflib" )

## Getting how many sequences are cut out when only including OTUs - R1.lulu.pspool.reflib
## For paper, excluding all samples except big grid (40 points across LFDP) and small grid (2x2m plot)
library(microViz)
precut <- R1phylo.plants.tax[[6]]
precut <- prune_samples(grepl("bigplot", precut@sam_data$experiment), precut)
precut <- subset_samples(precut, experiment != "bigplot/samptreat")
precut <- subset_samples(precut, DNAtreat == "clean")
precut@sam_data$origreflibsize <- rowSums(otu_table(precut))
precut@sam_data$origreflibotus <- rowSums(otu_table(precut) != 0)
precut <-phyloseq_validate(precut, remove_undetected = TRUE)

postcut <- R1ref.lib.list[[6]]
postcut <- prune_samples(grepl("bigplot", postcut@sam_data$experiment), postcut)
postcut <- subset_samples(postcut, experiment != "bigplot/samptreat")
postcut <- subset_samples(postcut, DNAtreat == "clean")
postcut@sam_data$origreflibsize <- rowSums(otu_table(precut))
postcut@sam_data$origreflibotus <- rowSums(otu_table(precut) != 0)
postcut@sam_data$afterlibsize <- rowSums(otu_table(postcut))
postcut@sam_data$afterlibotus <- rowSums(otu_table(postcut) != 0)
postcut <-phyloseq_validate(postcut, remove_undetected = TRUE)


precutsums <- rowSums(otu_table(precut))
postcutsums <- rowSums(otu_table(postcut))

##plots of sequences remaining after filtering for only reference library OTUs
## total sequences per sample before V after
plot(precutsums, postcutsums, xlab ="All Plants OTUs - library size", ylab ="Only Ref. Lib. OTUs - library size ")
abline(0,1, col="red")
## proportion remaining for each sample sequence count
boxplot(postcutsums/precutsums, ylab = "Prop. soil reads matching reference library taxa")

summary(postcutsums/precutsums)
sd(postcutsums/precutsums)

## total OTUs  per sample
precutotusums <- apply(otu_table(precut), MARGIN=1, FUN=function(x) {length( x[x > 0] )} ) 
postcutotusums <- apply(otu_table(postcut), MARGIN=1, FUN=function(x) {length( x[x > 0] )} ) 


plot(precutotusums, postcutotusums, xlab ="All Plants OTUs count", ylab ="Only Ref. Lib. OTUs count ")
abline(0,1, col="red")

## proportion OTUs remaining in each sample after filtering for what is in reference library
boxplot(postcutotusums/precutotusums, ylab = "Prop. soil OTUs matching reference library taxa")

plot(postcutsums/precutsums, postcutotusums/precutotusums, xlab ="Prop. reads remaining after ref. lib. filtering" , ylab ="Prop. OTUs remaining after ref. lib. filtering")
## reads
summary(postcutsums/precutsums)
## OTUs
summary(postcutotusums/precutotusums)
# so in general 91-98% (range of mean & median) of reads are from the reference library taxa 
# equating to 60 - 61.5% of OTUs (from the reference library taxa)
size <- rowSums(otu_table(precut))
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

# Save the combined plot as a PDF
pdf("FigureS8.pdf", width = 17 / 2.54, height = 8 / 2.54) # Convert cm to inches
print(combined_plot)
dev.off()

## So most reads and about 60% of OTUs remain when only using the reference library matches
## and there is no clear relationship here with individual sample sequencing depth (the colours)

<<<<<<< Updated upstream
=======
### Now separating just the BIG GRID and SMALL GRID samples separately to summarize

biggridprecut <- subset_samples(precut, experiment != "bigplot/smallss")
biggridprecut <- prune_samples(!grepl("single", biggridprecut@sam_data$discretehomo), biggridprecut)
biggridprecut@sam_data$origreflibsize <- rowSums(otu_table(biggridprecut))
biggridprecut@sam_data$origreflibotus <- rowSums(otu_table(biggridprecut) != 0)
biggridpostcut <- subset_samples(postcut, experiment != "bigplot/smallss")
biggridpostcut <- prune_samples(!grepl("single", biggridpostcut@sam_data$discretehomo), biggridpostcut)
biggridpostcut@sam_data$origreflibsize <- rowSums(otu_table(biggridpostcut))
biggridpostcut@sam_data$origreflibotus <- rowSums(otu_table(biggridpostcut) != 0)
biggridprecut <-phyloseq_validate(biggridprecut, remove_undetected = TRUE)
biggridpostcut <- phyloseq_validate(biggridpostcut, remove_undetected = TRUE)
biggridprecut
biggridpostcut
smallgridprecut <- subset_samples(precut, experiment == "bigplot/smallss")
smallgridprecut <-phyloseq_validate(smallgridprecut, remove_undetected = TRUE)
smallgridpostcut <- subset_samples(postcut, experiment == "bigplot/smallss")
smallgridpostcut <-phyloseq_validate(smallgridpostcut, remove_undetected = TRUE)

###### Checking library size effects with postcut samples
tab <- biggridpostcut
ttab <- otu_table(tab)
class(ttab) <- "matrix" 
library(vegan)
raremin <- min(rowSums(ttab))
raremin
rarecurve(ttab, step = 10, col = "blue", label = FALSE, xlim=c(0,200000))
rarecurve(ttab, step = 10, sample = 12302, col = "blue", label = FALSE, xlim=c(0,200000))

## Now to link the OTU names to the POTU table that cesc has, just access the object:

## to get the link between otu names and POTU (cesc's code for the reference libraries, load this list)
# otu.potu.link <- readRDS("/Users/glennd/Documents/GitHub/legenDNAry/Raw_data/Reference_library_filtering/OTU-to-RefIDs-List.rds")
otu.potu.link <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_v1.rds")
unique(otu.potu.link$dada.pspool.nc.lulu$OTU)
>>>>>>> Stashed changes
### Now to get some filtered variants for analyses representing lenient and stringent filtering
### Note that first sub-setting data so that only samples from the biggrid (no incubation experiments etc)
### are included
library(microViz)
lenient <- R1ref.lib.list[[6]]
lenient <- prune_samples(grepl("bigplot", lenient@sam_data$experiment), lenient)
lenient <- subset_samples(lenient, experiment != "bigplot/samptreat")
lenient <- subset_samples(lenient, DNAtreat == "clean")
lenient@sam_data$origreflibsize <- rowSums(otu_table(lenient))
lenient@sam_data$origreflibotus <- rowSums(otu_table(lenient) != 0)
lenient <-phyloseq_validate(lenient, remove_undetected = TRUE)

tabjoin <- merge(sample_data(biggridpostcut), otu.potu.link$dada.pspool.nc.lulu$OTU
                 , by.x =)
colnames(otu_table(biggridpostcut))
## Getting stats on library sizes
biggridsums <- rowSums(otu_table(lenient))

## Removing samples with >10,000 sequences in the sample
n <- 10000

# Check the sample sums
sample_counts <- sample_sums(lenient)

# Filter samples that meet the threshold
lenient <- prune_samples(sample_counts >= n, lenient)

# Check the number of samples remaining
nsamples(lenient)

### getting average and SD of library sizes
summary(sample_sums(lenient))
##se of lenient
sd(sample_sums(lenient))/sqrt(length(sample_sums(lenient)))

repfiltered <- R2ref.lib.list[[6]] ## same PS object as lenient, except only census OTUs occuring at at least 2/3 PCRs
repfiltered <- prune_samples(grepl("bigplot", repfiltered@sam_data$experiment), repfiltered)
repfiltered <- subset_samples(repfiltered, experiment != "bigplot/samptreat")
repfiltered <- subset_samples(repfiltered, DNAtreat == "clean")
repfiltered@sam_data$origreflibsize <- rowSums(otu_table(repfiltered))
repfiltered@sam_data$origreflibotus <- rowSums(otu_table(repfiltered) != 0)
repfiltered <-phyloseq_validate(repfiltered, remove_undetected = TRUE)

library(metagMisc)
##filtering lenient removeing OTUs per sample with < 1 sequence per 10,000 reads
lenientf1 <- phyloseq_filter_sample_wise_abund_trim(lenient, minabund = 0.0001, relabund = TRUE)
repfilteredf1 <- phyloseq_filter_sample_wise_abund_trim(repfiltered, minabund = 0.0001, relabund = TRUE)
##filtering lenient removeing OTUs per sample with < 1 sequence per 1,000 reads

## now we have 4 filters in the data from lenient to stringent
#1 lenient (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering) = 75 OTUs
#2 lenientf1 (lulu curated, dada psuedopooled option with no replicate bu rel. abund. threshold: OTU must have  >1 read per 10000 sequences) = 69 OTUs
#3 repfiltered (lulu curated, dada psuedopooled option with no relative abund. filtering, but OTU needs to be in at least 2/3 PCR replicates) = 58 OTUs)
#4 repfilteredf1 (lulu curated, dada psuedopooled option with rel. abund. threshold: OTU must have  >1 read per 10000 sequences and be in at least 2/3 PCR replicates) = 58 OTUs

## now applying rarefaction normalization of libraries to 12,302 sequences (the next largest library size over a 10,000 cut off threshold)
rarfun  <- function(x) {
  rfy <- 12302
  sar.rats.rarefy <- rarefy_even_depth(x, sample.size=rfy, replace=FALSE, rngseed = 1)
  return(sar.rats.rarefy)
}

rare_lenient <- rarfun(lenient)
rare_lenientf1 <- rarfun(lenientf1)
rare_repfiltered <- rarfun(repfiltered)
rare_repfilteredf1 <- rarfun(repfilteredf1)

## SO here are filtering varients - purely for the combined BIG GRID and NESTED samples:
## BIG GRID = (the main 40 samples with a few poitns dropped due to mall library size)
## NESTED = the 2x2m plot

#1 lenient:  (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering) = 62 OTUs
#1b rare_lenient: lenient (#1) but all libraries normalized by rarefying to an even depth = 56 OTUs
#2 lenientf1 (lulu curated, dada psuedopooled option with no replicate bu rel. abund. threshold: before rarefying OTU must have  >1 read per 10000 sequences) = 55 OTUs
#2b rare_lenientf1: lenientf1 (#2) but all libraries normalized by rarefying to an even depth = 54 OTUs
#3 repfiltered (lulu curated, dada psuedopooled option with no relative abund. filtering, but OTU needs to be in at least 2/3 PCR replicates) = 56 OTUs)
#3b rare_repfiltered: repfiltered (#3) but all libraries normalized by rarefying to an even depth = 55 OTUs
#4 repfilteredf1 (lulu curated, dada psuedopooled option with rel. abund. threshold: before rarefying OTU must have  >1 read per 10000 sequences and be in at least 2/3 PCR replicates) = 55 OTUs
#4b rare_repfilteredf1: repfilteredf1 (#4) but all libraries normalized by rarefying to an even depth = 55 OTUs


data <- list(lenient=lenient,
             rare_lenient=rare_lenient,
             lenientf1=lenientf1,
             rare_lenientf1=rare_lenientf1,
             repfiltered=repfiltered,
             rare_repfiltered=rare_repfiltered,
             repfilteredf1=repfilteredf1,
             rare_repfilteredf1=rare_repfilteredf1)

saveRDS(data, "Processed_data/PR_eDNA-for-analysis_2024-10-29.RData")
