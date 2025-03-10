##########################################
##########################################
##########################################
### COMPUTE SPECIES-SPECIFIC CONFUSION MATRIX STATS
##########################################
##########################################
##########################################


#####################
### Load data
#####################

samp.dna <- readRDS("Processed_data/stem-soil-40pt-data-lenient_10k-20241205.RDA")[[1]]
samp.dna.pa <- 1 * (samp.dna>0)

stem.23 <- readRDS("Processed_data/stem-soil-40pt-data-lenient_10k-20241205.RDA")[[2]]

# Get number of sample points to use as "observed", leaving 999 "random" points
npts <- nrow(stem.23$abund$`5`) - sum(grepl("random", rownames(stem.23$abund$`5`)))
nrand <- nrow(stem.23$abund$`5`) - npts
nsp <- ncol(samp.dna.pa)
  
# Convert stem data to presence absence
samp.abun.pa <- lapply(stem.23$abund, function(x){ 1 * (x[1:npts,] > 0)})

### gOTU full plot summaries
lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")

#####################
### Helper function
#####################
# Function to compute standardized effect size
ses <- function(obs, rand){
  return((obs - colMeans(rand))/apply(rand, 2, sd))
}

#####################
### Run the calculations
#####################

# Initialize results holders
sp_conf_stats_obs_list <- list()
sp_conf_stats_ses_list <- list()

for (r in 1) { #20) {
  message(paste('radius', r))
  # Subset the stem and dna communities to the focal radius
  tmp.samp.abun.pa <- samp.abun.pa[[r]]
  
  # Initialize matrices to store observation and SES results
  conf_stats_obs <- matrix(nrow = nsp, ncol = 12)
  conf_stats_ses <- matrix(nrow = nsp, ncol = 12)
  
  # Loop through sites
  for (sp in 1:nsp) {
    message(paste('= species', sp))
    
    # Initialize a results holder
    if(sp == 1) { confus_rand <- confus_obs <- list() }
    
    confus_obs[[sp]] <- caret::confusionMatrix(table(factor(+(!samp.dna.pa[,sp]),
                                                            levels=c(0,1)),
                                                       factor(+(!tmp.samp.abun.pa[,sp]), 
                                                              levels=c(0,1))))
    
    # Loop through randomizations
    for (rand in 1:999){
      
      # Shuffle species across sites
      confus_rand[[rand]] <- caret::confusionMatrix(table(factor(+(!samp.dna.pa[,sample(1:nsp, 1)]), 
                                                                 levels=c(0,1)),
                                                          factor(+(!tmp.samp.abun.pa[,sp]), 
                                                                 levels=c(0,1))))

      # Within species shuffle
    #   confus_rand[[rand]] <- caret::confusionMatrix(table(factor(sample(+(!samp.dna.pa[,sp])),
    #                                                              levels=c(0,1)),
    #                                                       factor(+(!tmp.samp.abun.pa[,sp]), 
    #                                                              levels=c(0,1))))
    }

    # Calculate SES
    conf_stats_ses[sp,1:11] <- ses(confus_obs[[sp]]$byClass,
                                     do.call(rbind, lapply(confus_rand, function(x) x$byClass)))
    
    obs_mcc <- mltools::mcc(confusionM = matrix(confus_obs[[sp]]$table, nrow=2))
    
    rand_mmc <- do.call(c, lapply(confus_rand, function(x)
                            mltools::mcc(confusionM = matrix(x$table, nrow=2))))

    conf_stats_ses[sp,12] <- (obs_mcc - mean(rand_mmc))/(sd(rand_mmc))
  }
  
  conf_stats_obs[,1:11] <- do.call(rbind, lapply(confus_obs, function(x) x$byClass))
  
  conf_stats_obs[,12] <- do.call(rbind, lapply(confus_obs, function(x)
                            mltools::mcc(confusionM = matrix(x$table, nrow=2))
                            ))
  
  colnames(conf_stats_ses) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_ses <- as.data.frame(conf_stats_ses)

  colnames(conf_stats_obs) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_obs <- as.data.frame(conf_stats_obs)
  
  sp_conf_stats_obs_list[[r]] <- conf_stats_obs
  sp_conf_stats_ses_list[[r]] <- conf_stats_ses
}


# Make colors to be relative to log basal area (x-axis is log abundance)
cp <- (viridis::viridis_pal(option = "A")(78))

# Get a simple vector of total abundance in same order as 

lfdp23 <- lfdp23[match(colnames(samp.dna.pa), rownames(lfdp23)),]

cpsp <- cp[round(rank(lfdp23$total_abund))]


par(mfcol=c(3,2), mar=c(4,4,2,1))

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$Sensitivity))
plot(lfdp23$total_abund, x[,1], log='x', 
     pch=21, bg=cpsp, cex=2,
     xlab="Total abundance", ylab="Sensitivity")
mtext("A", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$Specificity))
plot(lfdp23$total_abund, x[,1], log='x', 
     pch=21, bg=cpsp, cex=2,
     xlab="Total abundance", ylab="Specificity")
mtext("B", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$MCC))
plot(lfdp23$total_abund, x[,1], log='x', 
     pch=21, bg=cpsp, cex=2, 
     xlab="Total abundance", ylab="MCC")
mtext("C", adj=0, line=0.5)


x <- do.call(cbind, lapply(sp_conf_stats_ses_list, function(x) x$Sensitivity))

plot(lfdp23$total_abund, x[,1], log='x', pch=NA, bg=cp, 
     xlab="Total abundance", ylab="SES Sensitivity")
polygon(x=c(1e-6,5e5,5e5,1e-6), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(lfdp23$total_abund, x[,1], 
       pch=21, bg=cpsp, cex=2)
mtext("D", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_ses_list, function(x) x$Specificity))

plot(a, x[,1], log='x', pch=NA, bg=cp, 
     xlab="Total abundance", ylab="SES Specificity")
polygon(x=c(1e-6,5e5,5e5,1e-6), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(lfdp23$total_abund, x[,1], log='x', 
       pch=21, bg=cpsp, cex=2)
mtext("E", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_ses_list, function(x) x$MCC))

plot(a, x[,1], log='x', pch=NA, bg=cp, 
     xlab="Total abundance", ylab="SES MCC")
polygon(x=c(1e-6,5e5,5e5,1e-6), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(lfdp23$total_abund, x[,1], log='x', 
       pch=21, bg=cpsp, cex=2)
mtext("F", adj=0, line=0.5)









