##########################################
##########################################
##########################################
### SIMULATION TO SEE HOW GOOD EDNA NEEDS TO BE IN ORDER TO GET 'GOOD' CONFUSION MATRIX STATS
##########################################
##########################################
##########################################


#####################
### Load data
###### Bioinformatic method doesn't matter since we simulate eDNA
###### Note <40 sample points are represented because start data has been pruned earlier...
#####################

stem.23 <- readRDS("Processed_data/stem-soil-40pt-data-lenient_10k-20241205.RDA")[[2]]

# Get number of sample points to use as "observed", leaving 999 "random" points
nrand <- 999
npts <- nrow(stem.23$abund$`5`) - nrand

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
### Find gOTUS to drop
#####################

### Let's drop all gOTUs that make up < 10% of total cumulative abundance
# drops <- rownames(lfdp23)[lfdp23$cumprop < 0.1]

### Let's drop all gOTUs that make up < 25% of total cumulative abundance (69 species)
# drops <- rownames(lfdp23)[lfdp23$cumprop < 0.25]

### Let's drop the 5 most common gOTUs
# drops <- rownames(lfdp23)[order(lfdp23$total_abund, decreasing=T)][1:5]

### Let's drop all gOTUs that make up the lower 25% of total cumulative basal area (69 species)
# lfdp23 <- lfdp23[order(lfdp23$total_ba, decreasing=F),]
# lfdp23$cumprop_ba <- round(cumsum(lfdp23$total_ba)/sum(lfdp23$total_ba), 3)
# drops <- rownames(lfdp23)[lfdp23$cumprop_ba < 0.25]

### Let's drop the 5 most common gOTUs in terms of basal area
# drops <- rownames(lfdp23)[order(lfdp23$total_ba, decreasing=T)][1:5]

### Let's not drop anything!
drops <- NULL

#####################
### Find gOTUS to add
#####################

### Let's add the 5 most common gOTUs in terms of basal area to all samples
# adds <- rownames(lfdp23)[order(lfdp23$total_ba, decreasing=T)][1:5]

### Let's not add anything!
adds <- NULL


#####################
### Setup eDNA community
#####################

# Choose the DNA community (100% match with the 5m census)
samp.dna.pa <- samp.abun.pa$`30`

# Impose undetected species (drops from above)
samp.dna.pa[,drops] <- samp.dna.pa[,drops] * 0

# Impose the falsely detected species (adds from above)
samp.dna.pa[,adds] <- 1

#####################
### Run the simulation
#####################

# Initialize results holders
conf_stats_obs_list <- list()
conf_stats_ses_list <- list()

# Loop through radii
# seq(5,100,5)[10]

for (r in 1:20) {
  message(paste('radius', r))
  # Subset the stem and dna communities to the focal radius
  tmp.samp.abun.pa <- samp.abun.pa[[r]]
  
  # Initialize matrices to store observation and SES results
  conf_stats_obs <- matrix(nrow = npts, ncol = 12)
  conf_stats_ses <- matrix(nrow = npts, ncol = 12)
  
  # Loop through sites
  for (site in 1:npts) {
    message(paste('-site', site))
    # Initialize a results holder
    if(site == 1) { confus_rand <- confus_obs <- list() }

    confus_obs[[site]] <- caret::confusionMatrix(table(+(!samp.dna.pa[site,]),
                                                 +(!tmp.samp.abun.pa[site,])))

    # Loop through randomizations
    for (rand in 1:nrand){
      confus_rand[[rand]] <- caret::confusionMatrix(table(+(!samp.dna.pa[site,]),
                                                          +(!stem.23$abund[[r]][npts + rand,])))
      }
    
    
    # Calculate SES
    conf_stats_ses[site,1:11] <- ses(confus_obs[[site]]$byClass,
                                     do.call(rbind, lapply(confus_rand, function(x) x$byClass)))
    
    obs_mcc <- mltools::mcc(confusionM = matrix(confus_obs[[site]]$table, nrow=2))
    
    rand_mmc <- do.call(c, lapply(confus_rand, function(x)
      mltools::mcc(confusionM = matrix(x$table, nrow=2)
      )))
    
    conf_stats_ses[site,12] <- (obs_mcc - mean(rand_mmc))/(sd(rand_mmc))

  }

  conf_stats_obs[,1:11] <- do.call(rbind, lapply(confus_obs, function(x) x$byClass))
  
  conf_stats_obs[,12] <- do.call(rbind, lapply(confus_obs, function(x){
    mltools::mcc(confusionM = matrix(x$table, nrow=2))
  }))
  
  
  colnames(conf_stats_ses) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_ses <- as.data.frame(conf_stats_ses)

  colnames(conf_stats_obs) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_obs <- as.data.frame(conf_stats_obs)
  
  conf_stats_obs_list[[r]] <- conf_stats_obs
  conf_stats_ses_list[[r]] <- conf_stats_ses
  
}


#####################
### Save output
#####################

saveRDS(list(conf_stats_obs_list=conf_stats_obs_list, 
             conf_stats_ses_list=conf_stats_ses_list), 
        file="Processed_data/Sim_output_30meDNA-5to100m-20250307.RDA")


###### Names of previous saved simulation output files ######

## NO DROPS NOR ADDS (perfect sampling simulation)
"Processed_data/Sim_output_5meDNA-5to100m-20250307.RDA" # Match to 5 m community
"Processed_data/Sim_output_30meDNA-5to100m-20250307.RDA" # Match to 30 m community

## DROPS simulations
"Processed_data/Sim_output_25m-nodrops-20250305.RDA"
"Processed_data/Sim_output_5m-nodrops-20250305.RDA"
"Processed_data/Sim_output_5m-droplower10cumsum-5to30m-20250305.RDA"
"Processed_data/Sim_output_5m-droplower25cumsum-5to30m-20250305.RDA"
"Processed_data/Sim_output_5m-drop5topgOTUs-5to30m-20250305.RDA"
"Processed_data/Sim_output_5m-droplower25cumsumBA-5to30m-20250305.RDA"
"Processed_data/Sim_output_5m-drop5topgOTUs-BA-5to30m-20250305.RDA"

## ADDS simulations
"Processed_data/Sim_output_5m-add5topgOTUs-BA-5to30m-20250305.RDA"

## ADD + DROP simulations
"Processed_data/Sim_output_5m-add5topgOTUsBA-droplower10pctBA-5to30m-20250305.RDA"
"Processed_data/Sim_output_5m-add5topgOTUsBA-droplower25pctBA-5to30m-20250305.RDA"
"Processed_data/Sim_output_5m-add5topgOTUsBA-droplower25pctBA-5to50m-20250305.RDA"
"Processed_data/Sim_output_5m-add5topgOTUsBA-droplower25pctBA-5to100m-20250305.RDA"

conf_stats_obs_list <- readRDS("Processed_data/Sim_output_5meDNA-5to100m-20250307.RDA")[[1]]
conf_stats_ses_list <- readRDS("Processed_data/Sim_output_5meDNA-5to100m-20250307.RDA")[[2]]

#####################
### Plot it!
#####################

par(mfcol=c(3,2), mar=c(4,4,1,1))

# cols <- rev(viridis::viridis(20))
cols <- rev(viridis::viridis(length(conf_stats_obs_list)))

# "Observed" (raw simulated) values
boxplot(lapply(conf_stats_obs_list, function(x) x$Sensitivity),
        ylab="Sensitivity", 
        xlab="Radius (m)", axes=F, col=cols)
axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
mtext("A", adj=0, line=0.5)
abline(v=1, col='red', lwd=2, lty=2)

boxplot(lapply(conf_stats_obs_list, function(x) x$Specificity), 
        ylab="Specificity", 
        xlab="Radius (m)", axes=F, col=cols)
axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
mtext("B", adj=0, line=0.5)
abline(v=1, col='red', lwd=2, lty=2)

boxplot(lapply(conf_stats_obs_list, function(x) x$MCC),
               ylab="MCC", 
               xlab="Radius (m)", axes=F, col=cols)
axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
mtext("C", adj=0, line=0.5)
abline(v=1, col='red', lwd=2, lty=2)

# SES values
boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
        ylab="Sensitivity (SES)", 
        xlab="Radius (m)", axes=F, col=cols)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
        axes=F, col=cols, add=T)
axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
mtext("D", adj=0, line=0.5)
abline(v=1, col='red', lwd=2, lty=2)

boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
        ylab="Specificity (SES)", 
        xlab="Radius (m)", axes=F, col=cols)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
        axes=F, col=cols, add=T)
axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
mtext("E", adj=0, line=0.5)
abline(v=1, col='red', lwd=2, lty=2)

boxplot(lapply(conf_stats_ses_list, function(x) x$MCC), 
        ylab="MCC (SES)", 
        xlab="Radius (m)", axes=F, col=cols)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$MCC), 
        axes=F, col=cols, add=T)
axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
mtext("F", adj=0, line=0.5)
abline(v=1, col='red', lwd=2, lty=2)
