# 0. LOAD THE LFDP CENSUS DATA AS 'census'
head(census)
census$Mnemonic <- as.factor(census$Mnemonic)

# 1. READ SAMPLE POINT LOCATIONS
sample_xy <- read.csv("LFDP-eDNA-xy-V2.csv")

# 2. COMPUTE BASAL AREA OF EACH STEM
census$BA <- pi * (census$DBH/2000)^2

# 3. GET FUNCTION TO FILTER STEM DATA TO TREES IN A CERTAIN RADIUS AROUND SAMPLE PLOT
tree_dist <- function(sample_xy, census, radius=10) {
  d <- sqrt((sample_xy$PX - census$PX)^2 + (sample_xy$PY - census$PY)^2)
  foctrees <- data.frame(census[d <= radius & census$Status=="alive",])
  return(foctrees)
}

# 4. EXTRACT NEIGHBORHOODS AT A SET OF DIFFERENT RADII
radii <- seq(5, 100, by=5)
abund <- list()
ba <- list()

for(r in seq_along(radii)){
  neighborhoods <- list()
  for(i in 1:nrow(sample_xy)){
    neighborhoods[[i]] <- tree_dist(sample_xy[i,], census, radius=radii[r])
  }
  
  abund[[r]] <- do.call(rbind, lapply(neighborhoods, function(x) tapply(x$Mnemonic, x$Mnemonic, length)))
  abund[[r]][is.na(abund[[r]])] <- 0
  rownames(abund[[r]]) <- sample_xy$siteID
  
  ba[[r]] <- do.call(rbind, lapply(neighborhoods, function(x) tapply(x$BA, x$Mnemonic, sum, na.rm=T)))
  ba[[r]][is.na(ba[[r]])] <- 0
  rownames(ba[[r]]) <- sample_xy$siteID
}

names(abund) <- radii
names(ba) <- radii

# 5. GET NEAREST NEIGHBOR TREE TO EACH POINT
nearestsp <- function(xy, census) {
  d <- sqrt((xy$PX - census$PX)^2 + (xy$PY - census$PY)^2)
  nearest <- tapply(d, census$Mnemonic, min)
  return(nearest)
}

nntree_list <- list()
for(pt in 1:nrow(sample_xy)){
  xy <- sample_xy[pt, ]
  nntree_list[[pt]] <- nearestsp(xy, census)
}

nntrees <- do.call(rbind, nntree_list)

# 6. COMBINE AND SAVE RESULTS
results <- list(abund=abund, 
                ba=ba,
                nearest_sp=nntrees)

saveRDS(results, "Raw_data/LFDP2023-extract-20240510.RDA")

# 6. PLEASE SEND THE OUTPUT FILE to robert.muscarella@ebc.uu.se


