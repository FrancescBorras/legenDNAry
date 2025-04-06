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
R1Ref.lib.ids$dada.pspool.nc.lulu
R2Ref.lib.ids$dada.pspool.nc.lulu

## Note that OTU34 in R1 , and OTU 33 in R2 needs removal - it was from an unidentified "vine" in the plot (i.e. paraphyletic)
## and thus is not in the census data as well as having no identification use at all
R1Ref.lib.ids <- lapply(R1Ref.lib.ids, function(x) x[!(row.names(x) %in% "OTU34"), ])
R2Ref.lib.ids <- lapply(R2Ref.lib.ids, function(x) x[!(row.names(x) %in% "OTU33"), ])

str(R1Ref.lib.ids$dada.pspool.nc.lulu)
str(R2Ref.lib.ids$dada.pspool.nc.lulu)

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
postcut
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
#otu.potu.link <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_v1.rds")
#unique(otu.potu.link$dada.pspool.nc.lulu$OTU)
otu.potu.linkR1 <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R1.rds")
otu.potu.linkR2 <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R2.rds")
str(otu.potu.linkR1)
unique(otu.potu.linkR1$dada.pspool.nc.lulu$OTU)
unique(otu.potu.linkR2$dada.pspool.nc.lulu$OTU)

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
install.packages("remotes")
remotes::install_github("vmikk/metagMisc")
library(metagMisc)
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
#rare_lenientf1 <- rarfun(lenientf1) ## Rarefaction without replacement does not function - ignore
rare_repfiltered <- rarfun(repfiltered)
#rare_repfilteredf1 <- rarfun(repfilteredf1)  ## Rarefaction without replacement does not function - ignore
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

saveRDS(data, "Processed_data/PR_eDNA-for-analysis_2025-04-01.RData")

## also need an object with the above 12 sampling iterations for the intensive plot for confusion matrix analysis
inten_lenient_10k <- R1ref.lib.list[[6]]
inten_lenient_10k <- subset_samples(inten_lenient_10k, experiment == "bigplot/smallss")
inten_lenient_10k <-phyloseq_validate(inten_lenient_10k, remove_undetected = TRUE)
inten_rare_lenient_10k <- rarfun(inten_lenient_10k)
inten_lenientf1_10k <- phyloseq_filter_sample_wise_abund_trim(inten_lenient_10k, minabund = 0.0001, relabund = TRUE)
inten_repfil_10k <- R2ref.lib.list[[6]]
inten_repfil_10k <- subset_samples(inten_repfil_10k, experiment == "bigplot/smallss")
inten_repfil_10k <-phyloseq_validate(inten_repfil_10k, remove_undetected = TRUE)
inten_rare_repfiltered_10k <- rarfun(inten_repfil_10k)
inten_repfilteredf1_10k <- phyloseq_filter_sample_wise_abund_trim(inten_repfil_10k, minabund = 0.0001, relabund = TRUE)
inten_lenient_40k <- prune_samples(sample_sums(smallgridpostcut) >= 40000, smallgridpostcut)
inten_rare_lenient_40k <- rarfun(inten_lenient_40k)
inten_lenient_70k <- prune_samples(sample_sums(smallgridpostcut) >= 70000, smallgridpostcut)
inten_rare_lenient_70k <- rarfun(inten_lenient_70k)
inten_lenient_200k <- prune_samples(sample_sums(smallgridpostcut) >= 200000, smallgridpostcut)
inten_rare_lenient_200k <- rarfun(inten_lenient_70k)

datainten <- list(inten_lenient_10k=inten_lenient_10k,
                  inten_rare_lenient_10k=inten_rare_lenient_10k,
                  inten_lenientf1_10k=inten_lenientf1_10k,
                  inten_repfil_10k=inten_repfil_10k,
                  inten_rare_repfiltered_10k =inten_rare_repfiltered_10k ,
                  inten_repfilteredf1_10k=inten_repfilteredf1_10k, 
                  inten_lenient_40k=inten_lenient_40k, 
                  inten_rare_lenient_40k=inten_rare_lenient_40k,
                  inten_lenient_70k =inten_lenient_70k , 
                  inten_rare_lenient_70k=inten_rare_lenient_70k,
                  inten_lenient_200k=inten_lenient_200k, 
                  inten_rare_lenient_200k=inten_rare_lenient_200k)

saveRDS(datainten, "Processed_data/PR_eDNA-for-analysis_intenplot_2025-04-01.RData")

#### Looking at matching OTU - POTU - gOTU for intensive grid supporting information plots
## Matching

## Link between OTU and pOTU codes - old object incorrect for R2 iterations
otupotu <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_v1.rds")
otupotu <- otupotu[[6]]
str(otupotu)
## Link between OTU and pOTU codes
otupotuR1 <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R1.rds")
otupotuR2 <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R2.rds")
otupotuR1 <- otupotuR1[[6]] # only using pseudopooled, lulu filtered data (6th object on list)
otupotuR2 <- otupotuR2[[6]]
str(otupotuR1)
## sanity check are old OTU-to-RefIDs-List_v1.rds and OTU-to-RefIDs-List_R1.rds identical as they should be?
str(otupotu)
str(otupotuR1)
identical(otupotu, otupotuR1)

## need to remove OTU34 (unnecessary identified plot taxon) - OTU34 for R1 and OTU33 for R2
otupotuR1 <- otupotuR1[otupotuR1$OTU != "OTU34", ]
otupotuR2 <- otupotuR2[otupotuR2$OTU != "OTU33", ]
View(otupotuR1)

## Now getting pOTU to gOTU
codes <- readxl::read_xlsx("Raw_data/LFDP-SPcodes-Wseq.xlsx")
str(codes)
codes <- subset(codes, LFDP2023==1) # making sure only 2023 data is used
##object that has all three fields, OTU, POTU and gOTU
otu_gotu_mapR1 <- merge(otupotuR1, codes, by.x = "Genus", by.y = "POTU")[, c("OTU", "gOTU")]
otu_gotu_mapR2 <- merge(otupotuR2, codes, by.x = "Genus", by.y = "POTU")[, c("OTU", "gOTU")]

#otu_gotu_mapR1 <- merge(otupotuR1, codes, by.x = "sequence", by.y = "sequence")[, c("OTU", "gOTU")]
#otu_gotu_mapR2 <- merge(otupotuR2, codes, by.x = "sequence", by.y = "sequence")[, c("OTU", "gOTU")]

## Now looping through all PS objects in alldat (the different bioinformatic filterings of the 
## 38 biggrid points) and making new PS objects with just GOTU codes

alldat <- readRDS("Processed_data/PR_eDNA-for-analysis_2025-04-01.RData")
alldat
# Create named vectors for mapping (R1 and R2)
# Create named vectors for mapping (R1 and R2)
otu_to_gotu_R1 <- setNames(otu_gotu_mapR1$gOTU, otu_gotu_mapR1$OTU)
otu_to_gotu_R2 <- setNames(otu_gotu_mapR2$gOTU, otu_gotu_mapR2$OTU)

# Initialize list for updated phyloseq objects
alldatgotu <- list()

# Loop through each phyloseq object in the list
for (name in names(alldat)) {
  message("🔧 Processing: ", name)
  
  ps_obj <- alldat[[name]]
  
  # Choose appropriate OTU-to-gOTU mapping
  if (grepl("repfiltered", name, ignore.case = TRUE)) {
    otu_to_gotu <- otu_to_gotu_R2
    message("📌 Using R2 mapping")
  } else {
    otu_to_gotu <- otu_to_gotu_R1
    message("📌 Using R1 mapping")
  }
  
  # Extract components
  otu_tab <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(otu_table(ps_obj))) {
    otu_tab <- t(otu_tab)
  }
  sample_dat <- sample_data(ps_obj)
  tax_tab <- tax_table(ps_obj)
  
  # Match and rename OTUs
  matched_otus <- intersect(rownames(otu_tab), names(otu_to_gotu))
  renamed_otus <- otu_to_gotu[matched_otus]
  rownames(otu_tab)[match(matched_otus, rownames(otu_tab))] <- renamed_otus
  
  # Collapse by gOTU
  otu_tab_collapsed <- as.matrix(rowsum(otu_tab, group = rownames(otu_tab)))
  
  # Keep only valid gOTUs
  valid_gOTUs <- intersect(rownames(otu_tab_collapsed), unique(otu_to_gotu))
  otu_tab_collapsed <- otu_tab_collapsed[valid_gOTUs, , drop = FALSE]
  
  # Update taxonomy table to match collapsed gOTUs
  if (!is.null(tax_tab)) {
    rownames(tax_tab)[rownames(tax_tab) %in% matched_otus] <- renamed_otus
    tax_tab <- tax_tab[!duplicated(rownames(tax_tab)), , drop = FALSE]
    tax_tab <- tax_tab[valid_gOTUs, , drop = FALSE]
  }
  
  # Rebuild the phyloseq object
  new_ps <- phyloseq(
    otu_table(otu_tab_collapsed, taxa_are_rows = TRUE),
    sample_dat,
    if (!is.null(tax_tab)) tax_table(tax_tab)
  )
  
  # Store in output list
  alldatgotu[[name]] <- new_ps
  message("✅ Finished: ", name)
}


alldatgotu

#################################


####### Now working with gOTUs for just the big plot and smallplot samples for combined
####### figure 2 plots (rarefaction / hill numbers) etc 
subsample_etc<- readRDS("Processed_data/subsample-smallplot.RData")
subsample_etc
# Initialize list for updated phyloseq objects
# Create named vectors for mapping (R1 and R2)
otu_to_gotu_R1 <- setNames(otu_gotu_mapR1$gOTU, otu_gotu_mapR1$OTU)
otu_to_gotu_R2 <- setNames(otu_gotu_mapR2$gOTU, otu_gotu_mapR2$OTU)

# Initialize list for updated phyloseq objects
subsample_etc_gotus <- list()

# Loop through all subsample_etc phyloseq objects
for (name in names(subsample_etc)) {
  ps_obj <- subsample_etc[[name]]
  
  # Select correct OTU-to-gOTU mapping
  if (grepl("repfil", name)) {
    otu_to_gotu <- otu_to_gotu_R2
  } else {
    otu_to_gotu <- otu_to_gotu_R1
  }
  
  # Extract and process components
  otu_tab <- as(otu_table(ps_obj), "matrix")
  sample_dat <- sample_data(ps_obj)
  tax_tab <- tax_table(ps_obj)
  
  # Transpose OTU table if needed
  if (!taxa_are_rows(otu_table(ps_obj))) {
    otu_tab <- t(otu_tab)
  }
  
  # Match and rename OTUs
  matched_otus <- intersect(rownames(otu_tab), names(otu_to_gotu))
  renamed_otus <- otu_to_gotu[matched_otus]
  rownames(otu_tab)[match(matched_otus, rownames(otu_tab))] <- renamed_otus
  
  # Collapse duplicate gOTUs
  otu_tab_collapsed <- as.matrix(rowsum(otu_tab, group = rownames(otu_tab)))
  
  # Remove OTUs that weren't mapped to gOTUs
  valid_gOTUs <- intersect(rownames(otu_tab_collapsed), unique(otu_to_gotu))
  otu_tab_collapsed <- otu_tab_collapsed[valid_gOTUs, , drop = FALSE]
  
  # Update taxonomy table
  if (!is.null(tax_tab)) {
    rownames(tax_tab)[rownames(tax_tab) %in% matched_otus] <- renamed_otus
    tax_tab <- tax_tab[!duplicated(rownames(tax_tab)), , drop = FALSE]
    tax_tab <- tax_tab[valid_gOTUs, , drop = FALSE]
  }
  
  # Create updated phyloseq object
  new_ps <- phyloseq(
    otu_table(otu_tab_collapsed, taxa_are_rows = TRUE),
    sample_dat,
    if (!is.null(tax_tab)) tax_table(tax_tab)
  )
  
  # Store result
  subsample_etc_gotus[[name]] <- new_ps
  message("✅ Processed: ", name)
}




saveRDS(subsample_etc_gotus, "Processed_data/subsample-smallplot-gotus.RData")

## Getting gotus from 12 different bioinformatic filtering options

## lets extract the intensive plot and get some summary stats - total LFDP OTUs and mean / SD per sample
tab <- subsample_etc_gotus$ssdat
tabinc1 <- subset_samples(tab, experiment == "bigplot/smallss")
tabinc1 <- phyloseq_validate(tabinc1, remove_undetected = TRUE)
tabinc1 ## total LFDP OTUs

###Checking if all taxa in the unpooled and intensive plot samples are in the big grid samples
###########################################

alldatgotu$lenient_10k # used big grid data
tabinc1 # intentisive sample data
tab <- subsample_etc_gotus$ssdat
c12 <- subset_samples(tab, realsample1 == "C12")
c12 <- subset_samples(c12, DNAtreat == "clean")
c12<- phyloseq_validate(c12, remove_undetected = TRUE) # c12 pooled and unpooled sample data

tab <- subsample_etc_gotus$ssdat
c18 <- subset_samples(tab, realsample1 == "C18")
c18 <- subset_samples(c18, DNAtreat == "clean")
c18<- phyloseq_validate(c18, remove_undetected = TRUE) # c18 pooled and unpooled sample data

tab <- subsample_etc_gotus$ssdat
c8 <- subset_samples(tab, realsample1 == "C8")
c8 <- subset_samples(c8, DNAtreat == "clean")
c8<- phyloseq_validate(c8, remove_undetected = TRUE) # c8 pooled and unpooled sample data

# Function to check for missing taxa
check_missing_taxa <- function(small_phy, big_phy, label = "Dataset") {
  small_taxa <- taxa_names(small_phy)
  big_taxa <- taxa_names(big_phy)
  
  missing <- setdiff(small_taxa, big_taxa)
  
  if (length(missing) == 0) {
    cat(paste0("✅ All taxa in ", label, " are present in alldatgotu$lenient_10k.\n"))
  } else {
    cat(paste0("❌ ", label, " is missing ", length(missing), " taxa from alldatgotu$lenient_10k:\n"))
    print(missing)
  }
}

# Run the checks
check_missing_taxa(tabinc1, alldatgotu$lenient_10k, "tabinc1")
check_missing_taxa(c12, alldatgotu$lenient_10k, "c12")
tax_table(c12)["gOTU27", ]
check_missing_taxa(c18, alldatgotu$lenient_10k, "c18")
check_missing_taxa(c8, alldatgotu$lenient_10k, "c8")


