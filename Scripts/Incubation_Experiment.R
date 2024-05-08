# If loading from here al libraries used above need reloading...
library(dada2)
library(textshape)
library(tibble)
library(purrr)
library(Biostrings)
library(phyloseq)
library(microViz)
library(speedyseq)
library(phylosmith)
library(phyloseq)
library(vegan)
library(metagMisc)
library(ggplot2)
library(ggExtra)
library(ggplotify)
library(vegan)
library(iNEXT)
library(dplyr)
## Loading all data
tab <- readRDS("Processed_data/dilu_exp.rds")
## tab <- R1phylo.plants.tax[[5]] - note tab is THIS - R1phylo.plants.tax[[5]]

tabinc <- subset_samples(tab, factorlevel == "ambient" |  factorlevel == "ice" | factorlevel == "silica")
tabinc <- subset_samples(tabinc, DNAtreat == "clean")
tabinc <- phyloseq_validate(tabinc, remove_undetected = TRUE)
tabinc
View(sample_data(tabinc))
## reordering samples to make plots more informative
remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)
samorder <-sample.data.frame(tabinc)
samorder <- arrange(samorder, i.e.site, factorlevel)

n <- rownames(samorder)
n
##  Visualising paired-dilution relative compositions (top 20 taxa by total seq count)
a<- tabinc %>%
  comp_barplot(
    tax_level = "speciesOTU",
    sample_order = n ,
    #label = "realsample", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    #taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other OTUs", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) +
  coord_flip()
a

## rarefying data for diversity etc comparisons between samples
rfy <- min(rowSums(otu_table(tabinc)))
tabinc_scaled <- rarefy_even_depth(tabinc, sample.size=rfy, replace=FALSE, rngseed = 1) 
tabinc_scaled  %>%
  tax_transform("rclr", rank = "speciesOTU") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "factorlevel", shape = "i.e.site", size = 4) +
  #stat_ellipse(aes(colour = "i.e.site")) +
  scale_colour_brewer(palette = "Dark2")

rfy <- min(rowSums(otu_table(tabinc)))
tabinc_scaled <- rarefy_even_depth(tabinc, sample.size=rfy, replace=FALSE, rngseed = 1) 
grid <- get_variable(tabinc_scaled, "homo")
dists <- vegdist(otu_table(tabinc_scaled), binary=FALSE, method="robust.aitchison") 
# Adonis test - Permanova on rarefied data for if treatment had an effect (including site and site*treatment interaction to 
# look for a treatment effect while accounting for differences betwen sites )
sampledf <- data.frame(sample_data(tabinc_scaled))

adonis2(dists ~ i.e.site*factorlevel, data = sampledf)


## Quick look at OTU totals across all samples
otucounts <- phyloseq_ntaxa_by_tax(tabinc_scaled, TaxRank = "speciesOTU")
otucounts1 <- phyloseq_ntaxa_by_tax(tabinc_scaled, TaxRank = "phylum")
otusummary <- otucounts %>%
  group_by(name) %>%
  dplyr::summarise(
    OTUs = sum(N.OTU),
    .groups = "drop"
  )
ie.richness <- merge(otucounts1, otusummary, by = "name", all.x = TRUE)
View(ie.richness)
