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
R1Ref.lib.ids <- lapply(R1Ref.lib.ids, function(x) x[!(row.names(x) %in% "OTU34"), ])
R2Ref.lib.ids <- lapply(R2Ref.lib.ids, function(x) x[!(row.names(x) %in% "OTU34"), ])

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

postcut_r2 <- R2ref.lib.list[[6]]
postcut_r2 <- prune_samples(grepl("bigplot", postcut_r2@sam_data$experiment), postcut_r2)
postcut_r2 <- subset_samples(postcut_r2, experiment != "bigplot/samptreat")
postcut_r2 <- subset_samples(postcut_r2, DNAtreat == "clean")
postcut_r2@sam_data$origreflibsize <- rowSums(otu_table(precut))
postcut_r2@sam_data$origreflibotus <- rowSums(otu_table(precut) != 0)
postcut_r2@sam_data$afterlibsize <- rowSums(otu_table(postcut_r2))
postcut_r2@sam_data$afterlibotus <- rowSums(otu_table(postcut_r2) != 0)
postcut_r2 <-phyloseq_validate(postcut_r2, remove_undetected = TRUE)
postcut_r2
## For downstream graphs summarizing the proportions of reads and OTUs discarded in each sample
## because they did not appear in the reference library
saveRDS(list(precut=precut, postcut=postcut), 
         "Processed_data/Prereflibcut_Postreflibcut.RDA")

### Now separating just the BIG GRID and SMALL GRID samples separately to summarize

biggridprecut <- subset_samples(precut, experiment != "bigplot/smallss")
biggridprecut <- prune_samples(!grepl("single", biggridprecut@sam_data$discretehomo), biggridprecut)
biggridprecut@sam_data$origreflibsize <- rowSums(otu_table(biggridprecut))
biggridprecut@sam_data$origreflibotus <- rowSums(otu_table(biggridprecut) != 0)
## post reference library filtering
biggridpostcut <- subset_samples(postcut, experiment != "bigplot/smallss")
biggridpostcut <- prune_samples(!grepl("single", biggridpostcut@sam_data$discretehomo), biggridpostcut)
biggridpostcut@sam_data$origreflibsize <- rowSums(otu_table(biggridpostcut))
biggridpostcut@sam_data$origreflibotus <- rowSums(otu_table(biggridpostcut) != 0)
biggridprecut <-phyloseq_validate(biggridprecut, remove_undetected = TRUE)
biggridpostcut <- phyloseq_validate(biggridpostcut, remove_undetected = TRUE)
## post reference library filtering for R2 (OTUs in 2/3 replicates)
biggridpostcut_r2 <- subset_samples(postcut_r2, experiment != "bigplot/smallss")
biggridpostcut_r2 <- prune_samples(!grepl("single", biggridpostcut_r2@sam_data$discretehomo), biggridpostcut_r2)
biggridpostcut_r2@sam_data$origreflibsize <- rowSums(otu_table(biggridpostcut_r2))
biggridpostcut_r2@sam_data$origreflibotus <- rowSums(otu_table(biggridpostcut_r2) != 0)
biggridpostcut_r2 <- phyloseq_validate(biggridpostcut_r2, remove_undetected = TRUE)

biggridprecut
biggridpostcut
biggridpostcut_r2
## checking library sizes after excluding non-reference library sequences
## and then selecting three cutoffs for further analysis
rowSums(otu_table(biggridpostcut))
sum(rowSums(otu_table(biggridpostcut)) > 10000)
sum(rowSums(otu_table(biggridpostcut)) > 20000)
sum(rowSums(otu_table(biggridpostcut)) > 30000)
sum(rowSums(otu_table(biggridpostcut)) > 40000)
sum(rowSums(otu_table(biggridpostcut)) > 50000)
sum(rowSums(otu_table(biggridpostcut)) > 60000)
sum(rowSums(otu_table(biggridpostcut)) > 70000)
sum(rowSums(otu_table(biggridpostcut)) > 80000)
sum(rowSums(otu_table(biggridpostcut)) > 90000)
sum(rowSums(otu_table(biggridpostcut)) > 200000)

## So, cutoffs can be >10,000 = retain 95% of samples, >40,000 retain 70% of samples, 
## >70,000 retain over 50% of samples, >200,000 retain 25% of samples

biggrid_10k <- prune_samples(sample_sums(biggridpostcut) >= 10000, biggridpostcut)
biggrid_10k <- phyloseq_validate(biggrid_10k, remove_undetected = TRUE)
biggrid_10k_r2 <- prune_samples(sample_sums(biggridpostcut_r2) >= 10000, biggridpostcut_r2)
biggrid_10k_r2 <- phyloseq_validate(biggrid_10k_r2, remove_undetected = TRUE)
biggrid_40k <- prune_samples(sample_sums(biggridpostcut) >= 40000, biggridpostcut)
biggrid_40k <- phyloseq_validate(biggrid_40k, remove_undetected = TRUE)
biggrid_70k <- prune_samples(sample_sums(biggridpostcut) >= 70000, biggridpostcut)
biggrid_70k <- phyloseq_validate(biggrid_70k, remove_undetected = TRUE)
biggrid_200k <- prune_samples(sample_sums(biggridpostcut) >= 200000, biggridpostcut)
biggrid_200k <- phyloseq_validate(biggrid_200k, remove_undetected = TRUE)

saveRDS(list(biggrid_10k=biggrid_10k, biggrid_40k=biggrid_40k,
             biggrid_70k=biggrid_70k, biggrid_200k=biggrid_200k), 
        "Processed_data/Biggrid_libsize_filter.RDA")

smallgridprecut <- subset_samples(precut, experiment == "bigplot/smallss")
smallgridprecut <-phyloseq_validate(smallgridprecut, remove_undetected = TRUE)
smallgridpostcut <- subset_samples(postcut, experiment == "bigplot/smallss")
smallgridpostcut <-phyloseq_validate(smallgridpostcut, remove_undetected = TRUE)

###### Go to xxx script to make plot of this biggrid minimum library size filtering

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
lenient <- biggrid_10k
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
biggridsums <- rowSums(otu_table(biggrid_10k))

### getting together data for analysis and plots of sub-samples and small grid samples
ssdat <- R1ref.lib.list[[6]]
ssdat <- prune_samples(grepl("bigplot", ssdat@sam_data$experiment), ssdat)
ssdat <- subset_samples(ssdat, experiment != "bigplot/samptreat")
ssdat <- subset_samples(ssdat, DNAtreat == "clean")
ssdat@sam_data$origreflibsize <- rowSums(otu_table(ssdat))
ssdat@sam_data$origreflibotus <- rowSums(otu_table(ssdat) != 0)
ssdat <-phyloseq_validate(ssdat, remove_undetected = TRUE)
rarfun  <- function(x) {
  rfy <- min(rowSums(otu_table(x)))
  sar.rats.rarefy <- rarefy_even_depth(x, sample.size=rfy, replace=FALSE, rngseed = 1)
  return(sar.rats.rarefy)
}
repfil_ssdat <- R2ref.lib.list[[6]]
repfil_ssdat <- prune_samples(grepl("bigplot", repfil_ssdat@sam_data$experiment), repfil_ssdat)
repfil_ssdat <- subset_samples(repfil_ssdat, experiment != "bigplot/samptreat")
repfil_ssdat <- subset_samples(repfil_ssdat, DNAtreat == "clean")
repfil_ssdat@sam_data$origreflibsize <- rowSums(otu_table(repfil_ssdat))
repfil_ssdat@sam_data$origreflibotus <- rowSums(otu_table(repfil_ssdat) != 0)
repfil_ssdat <-phyloseq_validate(repfil_ssdat, remove_undetected = TRUE)

rare_ssdat <- rarfun(ssdat)
repfil_ssdattf1 <- phyloseq_filter_sample_wise_abund_trim(repfil_ssdat, minabund = 0.0001, relabund = TRUE)

dat1 <- list(ssdat = ssdat,
             rare_ssdat = rare_ssdat,
             repfil_ssdat = repfil_ssdat,
             repfil_ssdattf1 = repfil_ssdattf1)

saveRDS(dat1, "Processed_data/subsample-smallplot.RData")

### getting average and SD of library sizes
summary(sample_sums(lenient))
##se of lenient
sd(sample_sums(lenient))/sqrt(length(sample_sums(lenient)))

### Getting the replicate filtered data (OTU in at least 2/3 PCR replicates)

repfiltered <- biggrid_10k_r2 ## same PS object as lenient, except only census OTUs occuring at at least 2/3 PCRs
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
#1 lenient (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering) = 61 ASVs
#2 lenientf1 (lulu curated, dada psuedopooled option with no replicate bu rel. abund. threshold: OTU must have  >1 read per 10000 sequences) = 53 OTUs
#3 repfiltered (lulu curated, dada psuedopooled option with no relative abund. filtering, but OTU needs to be in at least 2/3 PCR replicates) = 50 OTUs)
#4 repfilteredf1 (lulu curated, dada psuedopooled option with rel. abund. threshold: OTU must have  >1 read per 10000 sequences and be in at least 2/3 PCR replicates) = 49 OTUs

## now applying rarefaction normalization of libraries to the minimum sequence count of each dataset
##  (the next largest library size over a 10,000 cut off threshold)
## Since rarefying wihtout replcaement fails if there are too few non-zero OTUs in a sample, thes
rarfun  <- function(x) {
  rfy <- min(rowSums(otu_table(x)))
  sar.rats.rarefy <- rarefy_even_depth(x, sample.size=rfy, replace=FALSE, rngseed = 1)
  return(sar.rats.rarefy)
}

rare_lenient <- rarfun(lenient)
rare_lenientf1 <- rarfun(lenientf1) ## Rarefaction without replacement does not function - ignore
rare_repfiltered <- rarfun(repfiltered)
rare_repfilteredf1 <- rarfun(repfilteredf1)  ## Rarefaction without replacement does not function - ignore
rare_lenient_40k <- rarfun(biggrid_40k)
rare_lenient_70k <- rarfun(biggrid_70k)
rare_lenient_200k <- rarfun(biggrid_200k)


## SO here are filtering varients - purely for the combined BIG GRID samples:
## BIG GRID = (the main 40 samples with a some poitns dropped due to library size filtering)

#1: lenient:  (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering, min library size = 10k) = 61 OTUs
#1b: rare_lenient: lenient (#1) but all libraries normalized by rarefying to an even depth = 54 OTUs
#2: lenientf1 (lulu curated, dada psuedopooled option with no replicate bu rel. abund. threshold, min library size = 10k: before rarefying OTU must have  >1 read per 10000 sequences) = 53 OTUs
#3: repfiltered (lulu curated, dada psuedopooled option with no relative abund. filtering, min library size = 10k, but OTU needs to be in at least 2/3 PCR replicates) = 50 OTUs)
#3b: rare_repfiltered: repfiltered (#3) but all libraries normalized by rarefying to an even depth = 49 OTUs
#4: repfilteredf1 (lulu curated, dada psuedopooled option with rel. abund. threshold: , min library size = 10k, before rarefying OTU must have  >1 read per 10000 sequences and be in at least 2/3 PCR replicates) = 49 OTUs
#5:lenient_40k (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering, min library size = 40k) = 55 OTUs
#6:rare_lenient_40k: lenient_40k (#5) but all libraries normalized by rarefying to an even depth = 54 OTUs
#7:lenient_70k (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering, min library size = 70k) = 50 OTUs
#8:rare_lenient_70k lenient_70k (#7) but all libraries normalized by rarefying to an even depth = 50 OTUs
#9:lenient_200k  (lulu curated, dada psuedopooled option with no replicate & relative abund. filtering, min library size = 200k) = 45 OTUs
#10:rare_lenient_200k lenient_200k (#9) but all libraries normalized by rarefying to an even depth = 50 OTUs
data <- list(lenient_10k=lenient,
             rare_lenient_10k=rare_lenient,
             lenientf1=lenientf1,
             repfiltered=repfiltered,
             rare_repfiltered=rare_repfiltered,
             repfilteredf1=repfilteredf1, 
             lenient_40k=biggrid_40k, 
             rare_lenient_40k=rare_lenient_40k,
             lenient_70k=biggrid_70k, 
             rare_lenient_70k=rare_lenient_70k,
             lenient_200k=biggrid_200k, 
             rare_lenient_200k=rare_lenient_200k)

saveRDS(data, "Processed_data/PR_eDNA-for-analysis_2024-12-04.RData")
