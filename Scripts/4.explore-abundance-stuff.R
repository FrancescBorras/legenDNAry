########################################################
########################################################
##################
### Explore false presences and abundance
##################
########################################################
########################################################

### Load libraries (not sure how many if any of these are needed below...)
library(vegan)
library(caret)
library(compositions)
library(viridis)
library(ecodist)
library(sf)
library(spdep)

datafile <- paste0("Processed_data/stem-soil-40pt-data-lenient_10k-20241205.RDA")
  
  ### Load data
  dnamat <- readRDS(datafile)[[1]]
  stem.23 <- readRDS(datafile)[[2]]
  stem.16 <- readRDS(datafile)[[4]]
  traits <- readRDS(datafile)[[3]]
  
  ### Count number of valid eDNA sample points (out of 40)
  (npts <- nrow(stem.23$abund[[1]]) - length(grep('random', rownames(stem.23$abund[[1]]))))
  
  ### gOTU full plot summaries
  lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
  
  ### Make presence absence matrices
  dnamat.pa <- 1*(dnamat>0)
  stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
  stem.pa.list16 <- lapply(stem.16$abund, function(x) 1*(x>0))
  
  ### Color palette for plotting
  cp <- rev(viridis::viridis(20))
  
  # gOTUs in the DNA data
  gotu_in <- colnames(dnamat.pa)[colSums(dnamat.pa)>0]
  
  sum(lfdp23$total_abund[rownames(lfdp23) %in% gotu_in]) / sum(lfdp23$total_abund)
  sum(lfdp23$total_ba[rownames(lfdp23) %in% gotu_in]) / sum(lfdp23$total_ba)
  length(gotu_in)/nrow(lfdp23)
  