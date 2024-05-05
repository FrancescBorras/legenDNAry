# DNA extraction treatment experiment 
# load("CescPRdata_v1.RData") # - remove first # to load environment.. 

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
## Loading all data
tab <- readRDS("Processed_data/dilu_exp.rds")
#########################################
## pulling out the DILUTION EXPERIMENT
## tab <- R1phylo.plants.tax[[5]] - note tab is THIS - R1phylo.plants.tax[[5]]
#tab <- k2
tabinc <- subset_samples(tab, dilution == "yes")
tabinc <- phyloseq_validate(tabinc, remove_undetected = TRUE)
tabinc
## reordering samples to make plots more informative
remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)
samorder <-sample.data.frame(tabinc)
library(dplyr)
samorder <- arrange(samorder, realsample, paircombo)

n <- rownames(samorder)

##  Visualising paired-dilution relative compositions
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
## Seems a Consistent difference in th incubation experiment samples only...
## lets see if the number of taxa detected differs between treatments

otucounts <- phyloseq_ntaxa_by_tax(tabinc, TaxRank = "speciesOTU")
otucounts1 <- phyloseq_ntaxa_by_tax(tabinc, TaxRank = "phylum")
otusummary <- otucounts %>%
  group_by(paircombo) %>%
  dplyr::summarise(
    OTUs = sum(N.OTU),
    .groups = "drop"
  )
dilu.richness <- merge(otucounts1, otusummary, by = "paircombo", all.x = TRUE)
View(dilu.richness)
## in order to plot paired samples, removing second repeat of leS2A3 (ddiltue) - it has an equal number of OTUs (22) anyway
dilu.richness <- dilu.richness[-c(31), ]
library(ggpubr)
library(dplyr)
library(tidyr)
df <- dilu.richness
df1 <- df %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = OTUs)
fab <- otu_table(tabinc)
fab <- rowSums(fab)
fab <- sort(fab)
mdilu <- merge(dilu.richness, data.frame(fab, newname2 = names(fab)), by = "newname2", all.x = TRUE)
df <- mdilu
df2 <- df %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = fab)

ggpaired(df1, cond1 = "clean", cond2 = "ddilute",
         fill = "condition", palette = "npg", line.color = "gray")
colnames(df2) <- c("realsample", "clean.lib", "ddilute.lib")
## important to remember that this is all confounded by library size -- we can add that to this paired table.. and normalize OTUs according to it. 
fab <- otu_table(tabinc)
fab <- rowSums(fab)
fab <- sort(fab)
mdilu <- merge(dilu.richness, data.frame(fab, newname2 = names(fab)), by = "newname2", all.x = TRUE)
str(mdilu)
richlib <- merge(df1, df2, by = "realsample", all.x = TRUE)
## normalizing OTUs by librarysize by OTUs/librarysize - OTUs per 10,000 reads
richlib$cleanNO <- 10000*(richlib$clean/richlib$clean.lib)
richlib$diluteNO <- 10000*(richlib$ddilute/richlib$ddilute.lib)

## boxplot on the normalized OTUs
ggpaired(richlib, cond1 = "cleanNO", cond2 = "diluteNO",
         fill = "condition", palette = "npg", line.color = "gray", legend.title="", tickslab= FALSE) + 
  theme(axis.text.x=element_blank())+
  xlab("Extraction Protocol") +
  ylab("OTU's per 10000 reads") +
  scale_fill_hue(labels = c("Cleaned", "Phenol.Dilu")) 

## seeing if we can T-test this difference -- normality in groups first 
shapiro.test(richlib$cleanNO)
shapiro.test(richlib$diluteNO)
## looks like not normal ( P < 0.05 = reject the null hypothesis that the data are normally distributed)
## we can use non-parametric tests
wilcox.test(richlib$cleanNO, richlib$diluteNO, paired = TRUE, alternative = "two.sided")

## We can also check this with a diversity index which uses the distribution of OTUs across libraries and sample sizes within libraies to calculate it
## (The shannon index) - we use shannon here because chao is really sensitive to the lack of singletons (a problem with sequence data)
plot_richness(tabinc, x="DNAtreat",  measures=c( "Shannon"))
plot_richness(tabinc, x="DNAtreat", measures="Shannon", color = "DNAtreat")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

## Looking at diversity difference using shannon and inverse simpson measures
alpha.diversity1 <- estimate_richness(tabinc, measures=c("Shannon"))
alpha.diversity2 <- estimate_richness(tabinc, measures=c("InvSimpson"))

chaodilu1 <- merge(dilu.richness, data.frame(alpha.diversity1, newname2 = rownames(alpha.diversity1)), by = "newname2", all.x = TRUE)
chaodilu2 <- merge(dilu.richness, data.frame(alpha.diversity2, newname2 = rownames(alpha.diversity2)), by = "newname2", all.x = TRUE)

df31 <- chaodilu1 %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = Shannon)
df32 <- chaodilu2 %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = InvSimpson)
colnames(df31) <- c("realsample", "clean.shan", "dilute.shan")
colnames(df32) <- c("realsample", "clean.InvSimpson", "dilute.InvSimpson")
richlib1 <- merge(richlib, df31, by = "realsample", all.x = TRUE)
richlib2 <- merge(richlib, df32, by = "realsample", all.x = TRUE)
richlib1
richlib2

shapiro.test(richlib1$clean.shan)
shapiro.test(richlib1$dilute.shan)

shapiro.test(richlib2$clean.InvSimpson)
shapiro.test(richlib2$dilute.InvSimpson)
## data are for shan, not for invsimpson normal so we can use a paired sample t test for shan and wilcox for Invsimp

t.test(richlib1$clean.shan, richlib1$dilute.shan, paired = TRUE, alternative = "two.sided")

wilcox.test(richlib2$clean.InvSimpson,  richlib2$dilute.InvSimpson, paired = TRUE, alternative = "two.sided")

## plot of differences
p1 <- ggpaired(richlib1, cond1 = "clean.shan", cond2 = "dilute.shan",
               fill = "condition", palette = "npg", line.color = "gray", legend.title="", tickslab= FALSE) + 
  theme(axis.text.x=element_blank())+
  xlab("Extraction Protocol") +
  ylab("Shannon Diversity Index") +
  scale_fill_hue(labels = c("Cleaned", "Chloro.Dilu")) 

p2 <- ggpaired(richlib2, cond1 = "clean.InvSimpson", cond2 = "dilute.InvSimpson",
               fill = "condition", palette = "npg", line.color = "gray", legend.title="", tickslab= FALSE) + 
  theme(axis.text.x=element_blank())+
  xlab("Extraction Protocol") +
  ylab("Inverse Simpson Diversity Index") +
  scale_fill_hue(labels = c("Cleaned", "Chloro.Dilu"))
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)

## ok so with the raw normalised data we have no evidence of a difference in alpha diversity (meaning cleaning the DNA made no difference)..but 
## tests with shannon index indicate significantly greater A diversity with the unclean DNA. FML...
## so buying those expensive kits and running through an entire column extract protocol was a waste of time... in fact counter productive
## excellent ....