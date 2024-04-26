########################################################
########################################################
##################
### Load packages and summary data
##################
########################################################
########################################################

### Load libraries
library(vegan)
library(caret)
library(compositions)
library(viridis)
library(ecodist)

### Load data
dnamat <- readRDS("Processed_data/stem-soil-39pt-data-20240411.RDA")[[1]]
stem.23 <- readRDS("Processed_data/stem-soil-39pt-data-20240411.RDA")[[2]]
stem.16 <- readRDS("Processed_data/stem-soil-39pt-data-20240423.RDA")[[4]]
traits <- readRDS("Processed_data/stem-soil-39pt-data-20240411.RDA")[[3]]
sampxy <- read.csv("Raw_data/LFDP-sample39-coordinates.csv", row.names = 1)

### Make presence absence matrices
dnamat.pa <- 1*(dnamat>0)
stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
stem.pa.list16 <- lapply(stem.16$abund, function(x) 1*(x>0))

### Color palette for plotting
cp <- rev(viridis(20))


########################################################
########################################################
##################
### SUMMARY DATA AT WHOLE PLOT LEVEL
##################
########################################################
########################################################

# Total number of gOTUs in the DNA data
(no_gotu_dna <- sum(colSums(dnamat.pa)>0))

# Total gOTUs in the stem data:
(no_gotu_stem <- sum(colSums(stem.23$abund[[20]])>0))

### Percent of gOTUs from stem data that are detected in DNA data
(100 * (no_gotu_dna/no_gotu_stem))

# Total abundance of stems per gOTU (as log count in 100 m radius around 39 points) vs.number of sites where gOTU was detected in the DNA
plot(colSums(stem.23$abund[[20]]), colSums(dnamat.pa), log='x')

# Significant positive correlation (i.e., more abundant gOTUs occur in more DNA samples)
cor.test(colSums(dnamat.pa), log10(colSums(stem.23$abund[[20]])))

# Total abundance of stems per gOTU (as count in 100 m radius around 39 points) vs.
plot(colSums(stem.23$abund[[20]]), 1*(colSums(dnamat.pa)>0), log='x')

# Significant logistic regression (i.e., more abundant gOTUs are more likely to have been detected overall in DNA)
a <- log10(colSums(stem.23$abund[[20]]))
m1 <- glm(1*(colSums(dnamat.pa)>0) ~ a, family="binomial")
summary(m1)

# Rank total abundance in stem data (based on 100 m radii) vs DNA read rank abundance
plot(rank(colSums(stem.23$abund[[20]])), rank(colSums(dnamat.pa)))

# Rank stem abundance is significantly positively correlated with rank DNA read abundance
cor.test(rank(colSums(stem.23$abund[[20]])), rank(colSums(dnamat.pa)))



### Figure 1 - there are a variety of plots below; need to decide which to include

pdf("Figures/Fig1.rank-abundance-plot.pdf")

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(rank(colSums(stem.23$abund[[20]])), rank(colSums(dnamat.pa)), 
     ylab="Rank abundance of DNA reads",
     xlab="Rank abundance of stems", 
     pch=21, bg='grey')
# abline(lm(rank(colSums(dnamat.pa)) ~ rank(colSums(stem.23$abund[[20]]))),
#        col='blue', lwd=2)

plot(colSums(stem.23$ba[[20]]),
     colSums(dnamat.pa), log='x',
     ylab="Sites detected in DNA",
     xlab="Total basal area (m^2) [log10]",
     pch=21, bg='grey')

# plot(colSums(stem.23$ba[[20]][1:39,]>0),
#      colSums(dnamat.pa),
#      ylab="Sites detected in DNA",
#      xlab="Sites detected in stem data",
#      pch=21, bg='grey')

# corres <- vector()
# for(i in 1:20){
#   corres[i] <- cor(colSums(stem.23$ba[[i]][1:39,]>0),
#                    colSums(dnamat.pa))
# }
# abline(lm(colSums(dnamat.pa) ~ colSums(stem.23$ba[[20]])),
#        col='blue', lwd=2)

# plot(colSums(stem.23$abund[[20]]),
#      jitter(1*(colSums(dnamat.pa)>0), 0.1), log='x',
#      ylab="Detected in DNA",
#      xlab="Total number of stems",
#      pch=21, bg='grey')

# Add logistic prediction
# nd <- data.frame(a=seq(min(log10(colSums(stem.23$abund[[20]]))),
#                        max(log10(colSums(stem.23$abund[[20]]))),
#                        length.out=1000))
# ypred <- predict(m1, newdata=nd, type="response")
# lines(10^(nd$a), ypred, col='blue', lwd=2)

# Species accumulation in stem data with radius around the sample points
b <- boxplot(sapply(stem.23$ba, function(x) rowSums(x[1:39,]>0)), col=cp,
             xlab="Radius (m)", 
             ylab="Taxon richness", 
             xlim=c(-0.5,20), ylim=c(0,80))
boxplot(rowSums(dnamat.pa), add=T, at=-0.5, width=2, col=2)
abline(v=0.25)
abline(h=79, lty=2)
axis(1, at=-0.5, labels="DNA")

# Species accumulation in DNA samples
# plot(specaccum(dnamat.pa),
#      xlab="Number of samples",
#      ylab="Taxon richness", 
#      ylim=c(0,80))
# plot(specaccum(dnamat.pa)$richness, pch=21, bg='grey',
#      xlab="Number of samples",
#      ylab="Taxon richness", 
#      ylim=c(0,80))
# lines(specaccum(dnamat.pa)$richness+specaccum(dnamat.pa)$sd)
# lines(specaccum(dnamat.pa)$richness-specaccum(dnamat.pa)$sd)
# abline(h=79, lty=2)

# Species accumulation in spatial aggregations of DNA samples
df <- data.frame(sapply(out, "length<-", max(lengths(out))))
colnames(df) <- rads
boxplot(df, col=rev(viridis(ncol(df))), 
        xlab="Aggregating Distance (m)", 
        ylab="Taxon richness",
        ylim=c(0,80))
abline(h=79, lty=2)

dev.off()



########################################################
########################################################
##################
### 39 Point Analysis
##################
########################################################
########################################################

########################################################
### 1. Correlation of species richness between DNA and stem data
########################################################

# Correlation of DNA and stem species richness across spatial scales
# Using both raw and rarefied data
corrs <- corrs_rare <- list()
for(r in seq_along(stem.23$abund)){
  corrs[[r]] <- cor.test(rowSums(dnamat>0), rowSums(stem.23$abund[[r]][1:39,]>0))
  
  corrs_rare[[r]] <- cor.test(rarefy(stem.23$abund[[r]][1:39,], 
                                        min(rowSums(stem.23$abund[[r]][1:39,]))),
                                 rarefy(dnamat, min(rowSums(dnamat))))
}

### Figure 2
pdf("Figures/Fig2.Richness-correlations.pdf", width = 9, height = 4)

par(mfrow=c(1,2))

plot(seq_along(stem.23$abund), sapply(corrs, function(x) x$estimate), 
     ylim=c(-1,1), axes=F,
     pch=21, bg='grey', xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness (raw)")
segments(seq_along(stem.23$abund), y0=sapply(corrs, function(x) x$conf.int[1]),
         y1=sapply(corrs, function(x) x$conf.int[2]))
abline(h=0, lty=2)
points(seq_along(stem.23$abund), sapply(corrs, function(x) x$estimate), 
       pch=21, bg='grey')
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)

plot(seq_along(stem.23$abund), sapply(corrs_rare, function(x) x$estimate), 
     ylim=c(-1,1), axes=F,
     pch=21, bg='grey', xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness (rarefied)")
segments(seq_along(stem.23$abund), y0=sapply(corrs_rare, function(x) x$conf.int[1]),
         y1=sapply(corrs_rare, function(x) x$conf.int[2]))
abline(h=0, lty=2)
points(seq_along(stem.23$abund), sapply(corrs_rare, function(x) x$estimate), 
       pch=21, bg='grey')
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)

dev.off()


########################################################
### 2. Procrusties test of stem ordination and dna ordinatino
########################################################

# Similar to Mantel test, lots to consider in terms of transformation, dist metric, ordination...

prores <- list()
for(r in 1:20){
  dnaord <- metaMDS(dnamat.pa, distance="jaccard")
  stemord <- metaMDS(1*(stem.23$abund[[r]][1:39,]>0), distance="jaccard")
  prores[[r]] <- protest(dnamat.pa, stemord, symmetric=T)
}

# A smaller Procrustes SS indicates a better fit or alignment between the two sets of points. It essentially tells you how well one set of points can be adjusted to match another set, considering only rotation, scaling, and translation as transformation operations.
# The correlation (x$t0) tells you the goodness-of-fit between the two configurations of multivariate data.
# A 'significant' p value (<0.05) indicates that the two matrices are similar

vals <- sapply(prores, function(x) summary(permustats(x))$z)
sig <- sapply(prores, function(x) x$signif)

### Figure 3
pdf("Figures/Fig3.Procrustes_correlations.pdf", width = 8, height = 8)

par(mar=c(4,4,1,1))
plot(vals, 
     pch=21, bg=cp, cex=2, axes=F,
     xlab="Radius (m)",
     ylab="Procrustes Correlation SES")
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(vals, pch=21, bg=cp, cex=2)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)
graphics::box()

dev.off()



########################################################
### 3. Confusion matrix
########################################################

### COMPUTE BALANCED ACCURACY 
### (INCLUDING STANDARDIZED EFFECT SIZE OF BALANCED ACCURACY)
### **Note that this takes while to run**

library(parallel)

# Set the number of cores to use for parallel processing
num_cores <- parallel::detectCores()

# Create clusters for parallel processing
cl <- makeCluster(num_cores)

# Initialize lists to store results
conf_stats_obs_list <- vector("list", length = 20)
conf_stats_ses_list <- vector("list", length = 20)

nruns <- 1000

# Export the caret package to all nodes
clusterEvalQ(cl, library(caret))

# Export necessary objects to all nodes
clusterExport(cl, c("dnamat.pa", "stem.23", "nruns", "ses", "nruns"))

# Parallel loop over radius
results <- parLapply(cl, 1:20, function(r) {
  message(paste("radius =", r))
  
  # Initialize matrices to store observation and SES results
  conf_stats_obs <- matrix(nrow = 39, ncol = 11)
  conf_stats_ses <- matrix(nrow = 39, ncol = 11)
  
  # Perform computations for each site
  for (site in 1:39) {
    # Compute confusion matrix for observation
    confus_obs <- caret::confusionMatrix(table(dnamat.pa[site,], 
                                               1 * (stem.23$abund[[r]][site,] > 0)))
    conf_stats_obs[site,] <- confus_obs$byClass
    
    # Generate randomizations
    randomizations <- 1:nruns + 39
    confus_rands <- lapply(randomizations, function(s) {
      confusionMatrix(table(dnamat.pa[site,],
                            1 * (stem.23$abund[[r]][s,] > 0)))
    })
    
    # Calculate SES for each class
    conf_stats_ses[site,] <- ses(conf_stats_obs[site,], 
                                 do.call(rbind, lapply(confus_rands, function(x) x$byClass)))
  }
  
  # Convert matrices to data frames
  conf_stats_obs <- as.data.frame(conf_stats_obs)
  conf_stats_ses <- as.data.frame(conf_stats_ses)
  names(conf_stats_ses) <- names(conf_stats_obs) <- names(confus_obs$byClass)
  
  # Return results as a list
  list(conf_stats_obs = conf_stats_obs, conf_stats_ses = conf_stats_ses)
})

# Close clusters
stopCluster(cl)

# Unpack results
for (i in 1:20) {
  conf_stats_obs_list[[i]] <- results[[i]]$conf_stats_obs
  conf_stats_ses_list[[i]] <- results[[i]]$conf_stats_ses
}


### There is a significant difference of observed balanced accuracy with increasing scale.
### Generally, the values decline as you get bigger scales
ba_obs <- do.call(cbind, lapply(conf_stats_obs_list, function(x) x$`Balanced Accuracy`))
ba_ses <- do.call(cbind, lapply(conf_stats_ses_list, function(x) x$`Balanced Accuracy`))

m2_obs <- aov(as.numeric(ba_obs) ~ as.factor(rep(1:20, each=39)))
anova(m2_obs)

### However, there is no difference in SES balanced accuracy across the scales
### A slight hump appers in ~ 30-40 m range but these are not sig different from 
m2_ses <- aov(as.numeric(ba_ses) ~ as.factor(rep(1:20, each=39)))
anova(m2_ses)
TukeyHSD(m2_ses)

boxplot(ba_ses[,1], ba_ses[,7])
t.test(ba_ses[,1], ba_ses[,7], paired=T)


### Figure 4

pdf("Figures/Fig4.Confusion-matrix-stats.pdf", width = 8, height = 10)

par(mfrow=c(3,2), mar=c(4,4,1,1))


# Sensitivity (true positive rate)
boxplot(lapply(conf_stats_obs_list, function(x) x$Sensitivity), 
        ylab="Sensitivity (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()

boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
        ylab="Sensitivity (SES)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
# abline(h=c(-1.96, 1.96), lwd=2, lty=2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
        axes=F, col=cp, add=T)
graphics::box()


# Specificity (true negative rate)
boxplot(lapply(conf_stats_obs_list, function(x) x$Specificity), 
        ylab="Specificity (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()

boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
        ylab="Specificity (SES)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
        axes=F, col=cp, add=T)
graphics::box()

# Balanced Accuracy
boxplot(lapply(conf_stats_obs_list, function(x) x$`Balanced Accuracy`), 
        ylab="Balanced Accuracy (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()

boxplot(lapply(conf_stats_ses_list, function(x) x$`Balanced Accuracy`), 
        ylab="Balanced Accuracy (SES)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
# abline(h=c(-1.96, 1.96), lwd=2, lty=2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$`Balanced Accuracy`), 
        axes=F, col=cp, add=T)
graphics::box()

dev.off()