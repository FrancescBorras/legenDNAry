### FAVA for compositional similarity between samples (eDNA and stem) at different scales

library("FAVA")

data_selector <- 1

outfiles <- ("stem-soil-40pt-data-lenient_10k-20250312.RDA")

### Load selected data file
datafile <- paste0("Processed_data/", outfiles[data_selector])

### Load data
dnamat <- readRDS(datafile)[[1]]
stem.23 <- readRDS(datafile)[[2]]

?fava


# Normalize to relative abundance (each row sums to 1)
dnamat_ra <- dnamat / rowSums(dnamat)


obs <- matrix(nrow=nrow(dnamat), ncol=20)
ses <- matrix(nrow=nrow(dnamat), ncol=20)

for(r in 1:20){
  
  print(r)

  npts <- nrow(dnamat_ra)
  
  stem_ra_mat <- stem.23$abund[[r]] / rowSums(stem.23$abund[[r]])
  
  for(site in 1:npts){
    
    dna_ra <- dnamat_ra[site,]
    
    obs[site, r] <- fava(rbind(dna_ra, stem_ra_mat[site,]))
   
    rand <- matrix(nrow=nrow(stem.23$abund[[r]]) - (npts+1))

    for(i in 1:99){ # nrow(rand)){
      
      rand[i,] <- fava(rbind(dna_ra, stem_ra_mat[npts+i,]))
    
    }
    
    ses[site, r] <- (obs[site, r] - mean(rand))/sd(rand)
    
  }
  
}


### Color palette for plotting
cp <- rev(viridis::viridis(20))

par(mfrow=c(2,2), mar=c(4,4,1,1))

boxplot(output, col=cp, names=seq(5,100,5),
        ylab="FAVA (OBS)")

boxplot(ses, col=cp, names=seq(5,100,5),
        ylab="FAVA (SES)")
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(ses, col=cp, names=seq(5,100,5), add=T)





######## FACTORED VERSION ############

# obs <- matrix(nrow = nrow(dnamat), ncol = 20)
# ses <- matrix(nrow = nrow(dnamat), ncol = 20)


npts <- nrow(dnamat_ra)

for (r in 1:20) {
  
  cat("Distance class:", r, "\n")  # Faster printing in loops
  flush.console()
  
  stem_abund_r <- stem.23$abund[[r]]
  stem_ra_mat <- stem_abund_r / rowSums(stem_abund_r)
  
  for (site in 1:npts) {
    
    dna_ra <- dnamat_ra[site, ]
    
    obs[site, r] <- fava(rbind(dna_ra, stem_ra_mat[site, ]))
    
    # Precompute indices for random samples
    rand_indices <- (npts + 1):(npts + 100) # nrow(stem_abund_r)
    
    # Compute all random values in a single operation (vectorized)
    rand <- apply(stem_ra_mat[rand_indices, , drop = FALSE], 1, function(x) fava(rbind(dna_ra, x)))
    
    ses[site, r] <- (obs[site, r] - mean(rand)) / sd(rand)
  }
}



######## PARALLEL VERSION ###############

library(parallel)

# Detect the number of available cores (adjust if needed)
n_cores <- detectCores() - 1  # Leave 1 core free for system stability

# obs <- matrix(nrow = nrow(dnamat), ncol = 20)
# ses <- matrix(nrow = nrow(dnamat), ncol = 20)

npts <- nrow(dnamat_ra)

# Parallel loop over 'r'
results <- mclapply(1:20, function(r) {
  
  stem_abund_r <- stem.23$abund[[r]]
  stem_ra_mat <- stem_abund_r / rowSums(stem_abund_r)
  
  obs_r <- numeric(npts)
  ses_r <- numeric(npts)
  
  for (site in 1:npts) {
    
    dna_ra <- dnamat_ra[site, ]
    
    obs_r[site] <- fava(rbind(dna_ra, stem_ra_mat[site, ]))
    
    rand_indices <- (npts + 1):nrow(stem_abund_r)
    
    # Compute random values in parallel (vectorized)
    rand <- apply(stem_ra_mat[rand_indices, , drop = FALSE], 1, function(x) fava(rbind(dna_ra, x)))
    
    ses_r[site] <- (obs_r[site] - mean(rand)) / sd(rand)
  }
  
  list(obs_r, ses_r)  # Return results for this 'r'
}, mc.cores = n_cores)

# Combine results
for (r in 1:20) {
  obs[, r] <- results[[r]][[1]]
  ses[, r] <- results[[r]][[2]]
}


