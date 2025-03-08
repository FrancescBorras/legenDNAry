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
lfdp23 <- lfdp23[order(lfdp23$total_abund, decreasing=F),]
lfdp23$cumprop <- round(cumsum(lfdp23$total_abund)/sum(lfdp23$total_abund), 3)
lfdp23$cumprop_ba <- round(cumsum(lfdp23$total_ba)/sum(lfdp23$total_ba), 3)

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

for (r in 1:20) {
  message(paste('radius', r))
  # Subset the stem and dna communities to the focal radius
  tmp.samp.abun.pa <- samp.abun.pa[[r]]
  
  # Initialize matrices to store observation and SES results
  conf_stats_obs <- matrix(nrow = nsp, ncol = 12)
  conf_stats_ses <- matrix(nrow = nsp, ncol = 12)
  
  # Loop through sites
  for (sp in 1:nsp) {
    message(paste('-species', sp))
    # Initialize a results holder
    if(sp == 1) { confus_rand <- confus_obs <- list() }
    
    confus_obs[[sp]] <- caret::confusionMatrix(table(factor(+(!samp.dna.pa[,sp]),
                                                            levels=c(0,1)),
                                                       factor(+(!tmp.samp.abun.pa[,sp]), 
                                                              levels=c(0,1))))
    
    # Loop through randomizations
    # for (rand in 1:nrand){
    #   confus_rand[[rand]] <- caret::confusionMatrix(table(+(!samp.dna.pa[,sp]),
    #                                                       +(!stem.23$abund[[r]][npts + rand,sp])))
    # }
    # 
    
    # Calculate SES
    # conf_stats_ses[site,1:11] <- ses(confus_obs[[site]]$byClass,
    #                                  do.call(rbind, lapply(confus_rand, function(x) x$byClass)))
    
    obs_mcc <- mltools::mcc(confusionM = matrix(confus_obs[[sp]]$table, nrow=2))
    
    # rand_mmc <- do.call(c, lapply(confus_rand, function(x)
    #   mltools::mcc(confusionM = matrix(x$table, nrow=2)
    #   )))
    
    # conf_stats_ses[site,12] <- (obs_mcc - mean(rand_mmc))/(sd(rand_mmc))
    
  }
  
  conf_stats_obs[,1:11] <- do.call(rbind, lapply(confus_obs, function(x) x$byClass))
  
  conf_stats_obs[,12] <- do.call(rbind, lapply(confus_obs, function(x){
    mltools::mcc(confusionM = matrix(x$table, nrow=2))
  }))
  
  
  # colnames(conf_stats_ses) <- c(names(confus_obs[[1]]$byClass), "MCC")
  # conf_stats_ses <- as.data.frame(conf_stats_ses)
  
  colnames(conf_stats_obs) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_obs <- as.data.frame(conf_stats_obs)
  
  sp_conf_stats_obs_list[[r]] <- conf_stats_obs
  # sp_conf_stats_ses_list[[r]] <- conf_stats_ses
  
}



cp <- (viridis::viridis_pal(option = "A")(78))

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$`Balanced Accuracy`))

plot(seq(5,100,5), x[1,], type='l', ylim=c(0,1), col=cp[sp],
     xlab="Scale (m)", ylab="Sensitivity")
for(sp in 2:78) {
  lines(seq(5,100,5), x[sp,], type='l', col=cp[sp])
}


hist(x[,1], breaks=20)
abline(v=0.5, col='blue', lwd=3)


a <- lfdp23$total_abund[match(colnames(samp.dna.pa), rownames(lfdp23))]
cp <- cp[match(colnames(samp.dna.pa), rownames(lfdp23))]


par(mfrow=c(2,2), mar=c(4,4,1,1))
x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$Sensitivity))
plot(a, x[,1], log='x', pch=21, bg=cp, 
     xlab="Total abundance", ylab="Sensitivity")
x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$Specificity))
plot(a, x[,1], log='x', pch=21, bg=cp,
     xlab="Total abundance", ylab="Specificity")
x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$`Balanced Accuracy`))
plot(a, x[,1], log='x', pch=21, bg=cp,
     xlab="Total abundance", ylab="Balanced Accuracy")













