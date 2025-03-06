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
library(sf)
library(spdep)

### Make labels for saving plots
label <- c("lenient_10k",
           "lenient_200k",
           "lenient_40k",
           "lenient_70k",
           "lenientf1",
           "rare_lenient_10k",
           "rare_lenient_200k",
           "rare_lenient_40k",
           "rare_lenient_70k",
           "rare_repfiltered",
           "repfiltered",
           "repfilteredf1")

### Name data files
outfiles <- c( "stem-soil-40pt-data-lenient_10k-20241205.RDA",
               "stem-soil-40pt-data-lenient_200k-20241205.RDA",
               "stem-soil-40pt-data-lenient_40k-20241205.RDA",
               "stem-soil-40pt-data-lenient_70k-20241205.RDA",
               "stem-soil-40pt-data-lenientf1-20241205.RDA",
               "stem-soil-40pt-data-rare_lenient_10k-20241205.RDA",
               "stem-soil-40pt-data-rare_lenient_200k-20241205.RDA",
               "stem-soil-40pt-data-rare_lenient_40k-20241205.RDA",
               "stem-soil-40pt-data-rare_lenient_70k-20241205.RDA",
               "stem-soil-40pt-data-rare_repfiltered-20241205.RDA",
               "stem-soil-40pt-data-repfiltered-20241205.RDA",
               "stem-soil-40pt-data-repfilteredf1-20241205.RDA")

### Select data file to load (for testing)
# data_selector <- 1

### Select data file to load (for full loop)
for(data_selector in seq_along(outfiles)){
  
  message(paste("Working on dataset", data_selector, ":", label[data_selector]))
  
### Load selected data file
datafile <- paste0("Processed_data/", outfiles[data_selector])

### Load data
dnamat <- readRDS(datafile)[[1]]
stem.23 <- readRDS(datafile)[[2]]
stem.16 <- readRDS(datafile)[[4]]
traits <- readRDS(datafile)[[3]]
# sampxy <- read.csv("Raw_data/LFDP-sample39-coordinates.csv", row.names = 1)

### Count number of valid eDNA sample points (out of 40)
(npts <- nrow(stem.23$abund[[1]]) - length(grep('random', rownames(stem.23$abund[[1]]))))

### Species code full plot summaries
lfdp <- readRDS("Raw_data/LFDP2023-extract-v2-20240427.RDA")
df <- readRDS("Raw_data/LFDP2016-extract-v2-20240427.RDA")

### gOTU full plot summaries
lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
lfdp16 <- readRDS("Raw_data/LFDP2016-extract-v2-20240427-gOTUs.RDA")

### Make presence absence matrices
dnamat.pa <- 1*(dnamat>0)
stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
stem.pa.list16 <- lapply(stem.16$abund, function(x) 1*(x>0))

### Color palette for plotting
cp <- rev(viridis::viridis(20))

### TRY DELETING THE EDGE SAMPLES...
# focsites <- rownames(sampxy)[sampxy$X<300 & sampxy$Y<450]
# dnamat.pa <- dnamat.pa[rownames(dnamat.pa) %in% focsites,]
# stem.23$abund <- lapply(stem.23$abund, function(x) x[rownames(x) %in% focsites,])
# stem.23$ba <- lapply(stem.23$ba, function(x) x[rownames(x) %in% focsites,])
# stem.23$nn <- stem.23$nn[rownames(stem.23$nn) %in% focsites,]

# gOTUs in the DNA data
gotu_in <- colnames(dnamat.pa)[colSums(dnamat.pa)>0]

sum(lfdp23$total_abund[rownames(lfdp23) %in% gotu_in]) / sum(lfdp23$total_abund)
sum(lfdp23$total_ba[rownames(lfdp23) %in% gotu_in]) / sum(lfdp23$total_ba)
length(gotu_in)/nrow(lfdp23)

########################################################
########################################################
##################
### SUMMARY DATA AT WHOLE PLOT LEVEL
##################
########################################################
########################################################

message(paste("Doing plot-level summary, Fig. 1"))

### Initiate a table to hold plot-level summary stats
if(data_selector==1){
plot_level_summary_table <- data.frame(dataset=NA,
                                       npts=NA,
                                       no_gotu_dna=NA,
                                       no_gotu_stem=NA,
                                       mdna=NA, sddna=NA, mindna=NA, maxdna=NA,
                                       m5=NA, sd5=NA, min5=NA, max5=NA,
                                       m100=NA, sd100=NA, min100=NA, max100=NA,
                                       rank_corr_est=NA, rank_corr_p=NA,
                                       ba_corr_est=NA, ba_corr_p=NA,
                                       log_p=NA)
}


mean(plot_level_summary_table$mdna)
range(plot_level_summary_table$mdna)

# Total number of gOTUs in the DNA data
(no_gotu_dna <- sum(colSums(dnamat.pa)>0))

# Total gOTUs in the stem data:
(no_gotu_stem <- sum(lfdp23$total_abund>0))

### Percent of gOTUs from stem data that are detected in DNA data
round(100 * (no_gotu_dna/no_gotu_stem), 1)

### Mean and SD of taxon richness in eDNA data
(mdna <- mean(rowSums(dnamat.pa)))
(sddna <- sd(rowSums(dnamat.pa)))
(mindna <- min(rowSums(dnamat.pa)))
(maxdna <- max(rowSums(dnamat.pa)))

### Mean and SD of taxon richness in stem data at 5 m
(m5 <- mean(rowSums(stem.23$abund[[1]]>0)[1:npts]))
(sd5 <- sd(rowSums(stem.23$abund[[1]]>0)[1:npts]))
(min5 <- min(rowSums(stem.23$abund[[1]]>0)[1:npts]))
(max5 <- max(rowSums(stem.23$abund[[1]]>0)[1:npts]))

### Mean and SD of taxon richness in stem data at 100 m
(m100 <- mean(rowSums(stem.23$abund[[20]]>0)[1:npts]))
(sd100 <- sd(rowSums(stem.23$abund[[20]]>0)[1:npts]))
(min100 <- min(rowSums(stem.23$abund[[20]]>0)[1:npts]))
(max100 <- max(rowSums(stem.23$abund[[20]]>0)[1:npts]))


##################
### Figure 1 - there are a variety of plots below; need to decide which to include

# pdf("Figures/Fig1.rank-abundance-plot.pdf")

pdf(paste0("Figures/compare_filtering/Fig1.rank-abundance-plot_", 
           label[data_selector],
           ".pdf"), width = 7, height = 6)


par(mfrow=c(2,2), mar=c(4,4,1,1))

# TAXON ACCUMULATION IN STEM DATA WITH INCREASING RADII AROUND SAMPLE POINTS
b <- boxplot(sapply(stem.23$abund, function(x) rowSums(x>0)[1:npts]), 
             col=cp,
             xlab="Radius (m)", 
             ylab="Taxon richness", 
             xlim=c(-0.5,20), ylim=c(0,no_gotu_stem))
b2 <- boxplot(rowSums(dnamat.pa), add=T, at=-0.5, width=2, col=2)
axis(1, at=-0.5, labels="DNA")
abline(v=0.25)
abline(h=no_gotu_stem, lty=2)
mtext("A", adj=0.1, line=-2)

# Full plot abundance
plot(rank(lfdp23$total_abund), rank(colSums(dnamat.pa)), 
     ylab="Rank abundance of DNA reads",
     xlab="Rank abundance of stems", 
     pch=21, bg='grey')
mtext("B", adj=0.05, line=-1.5)

# Rank stem abundance is significantly positively correlated with rank DNA read abundance
(rankcor <- cor.test(rank(lfdp23$total_abund), rank(colSums(dnamat.pa))))
rank_corr_est <- rankcor$estimate
rank_corr_p <- rankcor$p.value


plot(lfdp23$total_ba * 100,
     (colSums(dnamat.pa)/npts), log='x',
     ylab="Prop. of sites detected in DNA",
     xlab="Total basal area (m^2) [log10]",
     pch=21, bg='grey', ylim=c(0,1))
mtext("C", adj=0.05, line=-1.5)

# Total BA is significantly positively correlated with prop. of DNA samples where detected
(bacor <- cor.test(lfdp23$total_ba * 100, colSums(dnamat.pa)/npts))
ba_corr_est <- bacor$estimate
ba_corr_p <- bacor$p.value


# plot(colSums(stem.23$abund[[20]][1:39,]>0),
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


### Logistic regression with distance to nearest tree and detection
y <- as.vector(dnamat.pa[1:npts,])
x <- log10(as.vector(unlist(stem.23$nn[1:npts,])))

plot(x, jitter(y, 0.1), pch=16, cex=0.75, col=rgb(0,0,0,0.05),
     ylab="Presence in DNA",
     xlab="Dist. to nearest tree (m)", #[log10]", 
     xlim=c(0,3), cex.lab=1.25, 
     axes=F)
axis(1, c(0,1,2,3), labels=c(1, 10, 100,1000))
axis(2)
box()

mod <- glm(as.vector(dnamat.pa) ~ x, family="binomial")
nd <- data.frame(x=seq(0, max(x), length.out=100))
ypred <- predict(mod, nd, type='response')
# 
# # Overall fit
lines(nd$x, ypred, col="blue", lwd=3)
# 
# gOTU-specific fits
coeffs <- vector()
for(sp in 1:ncol(stem.23$nn)){
  xx <- log(unlist(stem.23$nn[1:npts,sp]))
  yy <- unlist(as.data.frame(1*(dnamat[,sp]>1)))
  mod <- glm(yy ~ xx, family=binomial)
  coeffs[sp] <- coef(mod)[2]
  nd <- data.frame(xx=seq(0, max(x), length.out=100))
  ypred <- predict(mod, nd, type='response')
  lines(nd$xx, ypred, col=rgb(0,0,0,0.2), lwd=0.5)
  #col=ifelse(coeffs[sp]>0, rgb(1,0,0,0.4), rgb(0,0,1,0.4)))
}
mtext("D", adj=0.05, line=-1.5)

log_p <- summary(mod)$coefficients["xx","Pr(>|z|)"]

# plot(lfdp23$total_abund,
#      jitter(1*(colSums(dnamat.pa)>0), 0.1), log='x',
#      xlim=c(1,10000),
#      ylab="Detected in DNA",
#      xlab="Total number of stems",
#      pch=21, bg='grey')
# mtext("C", adj=0.05, line=-1.5)
# 
# # Add logistic prediction
# nd <- data.frame(a=seq(log10(min(lfdp23$total_abund)),
#                        log10(max(lfdp23$total_abund)),
#                        length.out=1000))
# ypred <- predict(m1, newdata=nd, type="response")
# lines(10^(nd$a), ypred, col='blue', lwd=2)
# 
# lfdp23$eDNA <- 1*(colSums(dnamat.pa)>0)
# lfdp23[order(lfdp23$eDNA, lfdp23$total_abund),]



### TAXON ACCUMULATION CURVE WHEN YOU AGGREGATE DNA SAMPLES RANDOMLY
# plot(specaccum(dnamat.pa),
#      xlab="Number of samples",
#      ylab="Taxon richness",
#      ylim=c(0,79))
# abline(h=no_gotu_stem, lty=2)
# mtext("D", adj=0.05, line=-2)

### TAXON ACCUMULATION CURVE WHEN YOU AGGREGATE DNA SAMPLES SPATIALLY
# df <- data.frame(sapply(out, "length<-", max(lengths(out))))
# colnames(df) <- seq(5,100,5)
# boxplot(df, col=rev(viridis(ncol(df))), 
#         xlab="Aggregating Distance (m)", 
#         ylab="Taxon richness",
#         ylim=c(0, no_gotu_stem))
# abline(h=no_gotu_stem, lty=2)

dev.off()

plot_level_summary_table[data_selector,] <- c(label[data_selector],
                                              npts,
                                              no_gotu_dna,
                                              no_gotu_stem,
                                              mdna, sddna, mindna, maxdna,
                                              m5, sd5, min5, max5,
                                              m100, sd100, min100, max100,
                                              rank_corr_est, rank_corr_p,
                                              ba_corr_est, ba_corr_p,
                                              log_p)



plot_level_summary_table[,-1] <- round(apply(plot_level_summary_table[,-1], 2, as.numeric), 3)

# plot_level_summary_table

# x <- plot_level_summary_table
# cs <- RColorBrewer::brewer.pal(12, "Paired")
# plot(x$npts, x$no_gotu_dna, pch=1:12, col=cs, cex=2, lwd=3,
#      xlab="Number of sample points retained",
#      ylab="Number of gOTUs detected in eDNA")
# legend("topleft", legend=x$dataset, cex=0.75, pch=1:12, col=cs, bty='n', pt.cex=1, pt.lwd = 2)

write.csv(plot_level_summary_table, "Results/plot_level_summary_table.csv", row.names = F)


########################################################
########################################################
##################
### 40 Point Analysis
##################
########################################################
########################################################

########################################################
### 1. Correlation of species richness between DNA and stem data
########################################################

# Correlation of DNA and stem species richness across spatial scales
# Using both raw and rarefied data
corrs <- corrs_rare <- list()
corrs16 <- corrs_rare16 <- list()

## This block was from previous data - maybe able to delete now...
# For some reason, with the 7th file, the rarefy function is throwing error
# It seems to be some kind of strange formatting of the `dnamat` numbers but I cannot solve it
# If I save it and read back in then it seems to convert properly and works as expected.
# if(data_selector==7){
#   write.csv(dnamat, file="Processed_data/temp_dnamat.csv")
#   dnamat <- read.csv("Processed_data/temp_dnamat.csv", row.names = 1)
# }

for(r in seq_along(stem.23$abund)){
  
  corrs[[r]] <- cor.test(rowSums(dnamat > 0), rowSums(stem.23$abund[[r]][1:npts,]>0))
  corrs16[[r]] <- cor.test(rowSums(dnamat > 0), rowSums(stem.16$abund[[r]][1:npts,]>0))
  
  corrs_rare[[r]] <- cor.test(rarefy(stem.23$abund[[r]][1:npts,],
                                     min(rowSums(stem.23$abund[[r]][1:npts,]))),
                              rarefy(dnamat, min(rowSums(dnamat))))
  corrs_rare16[[r]] <- cor.test(rarefy(stem.16$abund[[r]][1:npts,],
                                     min(rowSums(stem.16$abund[[r]][1:npts,]))),
                              rarefy(dnamat, min(rowSums(dnamat))))
  
}



# corrs <- list()
# for(r in seq_along(stem.23$abund)){
#   corrs[[r]] <- cor.test(rowSums(dnamat>0), 
#                          renyi(stem.23$abund[[r]][1:39,], hill = T, scales=2))
# }

### Figure 2
# pdf("Figures/Fig2.Richness-correlations.pdf", width = 9, height = 4)

pdf(paste0("Figures/compare_filtering/Fig2.Richness-correlations_", 
    label[data_selector],
    ".pdf"), width = 9, height = 4)

par(mfrow=c(1,2))

plot(seq_along(stem.23$abund), sapply(corrs, function(x) x$estimate), 
     ylim=c(-1,1), axes=F,
     pch=21, bg=cp, xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness (raw)")
segments(seq_along(stem.23$abund), y0=sapply(corrs, function(x) x$conf.int[1]),
         y1=sapply(corrs, function(x) x$conf.int[2]))
abline(h=0, lty=2)
points(seq_along(stem.23$abund), sapply(corrs, function(x) x$estimate), 
       pch=21, bg=cp)
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
       pch=21, bg=cp)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)

dev.off()

########################################################
### 2. Procrusties test of stem ordination and dna ordination
########################################################

### Similar to Mantel test, lots to consider in terms of transformation, dist metric, ordination...
# 
# prores <- list()
# for(r in 1:20){
#   dnaord <- metaMDS(dnamat.pa, distance="jaccard")
#   
#   dnaord <- prcomp(dnamat)
#   
#   stemord <- metaMDS(1*(stem.23$abund[[r]][1:npts,]>0), distance="jaccard")
#   prores[[r]] <- protest(dnaord, stemord, symmetric=T)
# }

# A smaller Procrustes SS indicates a better fit or alignment between the two sets of points. It essentially tells you how well one set of points can be adjusted to match another set, considering only rotation, scaling, and translation as transformation operations.
# The correlation (x$t0) tells you the goodness-of-fit between the two configurations of multivariate data.
# A 'significant' p value (<0.05) indicates that the two matrices are similar

# corrvals <- sapply(prores, function(x) x$t0)

# vals <- sapply(prores, function(x) summary(permustats(x))$z)

# vals <- sapply(prores, function(x) x$t0)
# 
# quants <- sapply(prores, function(x) {
#   quantile(x$t0 - x$t, probs=c(0.025, 0.975)) + median(x$t)}
#   )

# sig <- sapply(prores, function(x) x$signif)

# ### Figure 3
# pdf("Figures/Fig3.Procrustes_correlations.pdf", width = 8, height = 8)
# 
# par(mfrow=c(1,1))
# plot(vals, 
#      pch=21, bg=cp, cex=2, axes=F,
#      xlab="Radius (m)",
#      ylab="Procrustes Correlation", ylim=c(-0.1, 0.7))
# segments(1:20, quants[1,], 1:20, quants[2,], lwd=2)
# # polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
# points(vals, pch=21, bg=cp, cex=2)
# axis(1, labels=seq(5, 100, 5), at=1:20)
# axis(2)
# abline(h=0, lty=2)
# graphics::box()
# 
# dev.off()



### COMPARE WITH RANDOM POINTS
# prores_rand <- list()
# for(r in 1:20){
#   dnaord <- metaMDS(dnamat.pa, distance="jaccard")
#   stemord <- metaMDS(1*(stem.23$abund[[r]][151:(150+npts),]>0), distance="jaccard")
#   prores_rand[[r]] <- protest(dnaord, stemord, symmetric=T)
# }
# 
# vals <- sapply(prores_rand, function(x) summary(permustats(x))$z)
# par(mar=c(4,4,1,1))
# plot(vals, 
#      pch=21, bg=cp, cex=2, axes=F,
#      xlab="Radius (m)",
#      ylab="Procrustes Correlation SES")
# polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
# points(vals, pch=21, bg=cp, cex=2)
# axis(1, labels=seq(5, 100, 5), at=1:20)
# axis(2)
# graphics::box()


########################################################
### 3. Confusion matrix
########################################################

### COMPUTE BALANCED ACCURACY 
### (INCLUDING STANDARDIZED EFFECT SIZE OF BALANCED ACCURACY)
### **Note that this takes while to run**

message(paste("working on confusion matrix for", datafile, "; dataset", data_selector))


library(parallel)

# Set the number of cores to use for parallel processing
num_cores <- parallel::detectCores()

# Create clusters for parallel processing
cl <- makeCluster(num_cores)

# Initialize lists to store results
conf_stats_obs_list <- vector("list", length = 20)
conf_stats_ses_list <- vector("list", length = 20)

nruns <- 999

# Function to compute standardized effect size
ses <- function(obs, rand){
  return((obs - colMeans(rand))/apply(rand, 2, sd))
}

# Export the caret package to all nodes
clusterEvalQ(cl, library(caret))
clusterEvalQ(cl, library(mltools))

# Export necessary objects to all nodes
clusterExport(cl, c("dnamat.pa", "stem.23", "nruns", "ses", "nruns", "npts"))

# Parallel loop over radius
results <- parLapply(cl, 1:20, function(r) {
  message(paste("radius =", r))
  
  # Initialize matrices to store observation and SES results
  conf_stats_obs <- matrix(nrow = npts, ncol = 12)
  conf_stats_ses <- matrix(nrow = npts, ncol = 12)
  
  # Perform computations for each site
  for (site in 1:npts) {
    # Compute confusion matrix for observation
    # confus_obs <- caret::confusionMatrix(table(dnamat.pa[site,], 
    #                                            1*(stem.23$abund[[r]][site,]>0)))

        confus_obs <- caret::confusionMatrix(table(+(!dnamat.pa[site,]), 
                                               +(!(1*(stem.23$abund[[r]][site,]>0)))
                                               ))
    
    conf_stats_obs[site,1:11] <- confus_obs$byClass
    
    # Compute MCC
    conf_stats_obs[site,12] <- mltools::mcc(TP=confus_obs$table[2,2],
                                            FP=confus_obs$table[2,1],
                                            TN=confus_obs$table[1,1],
                                            FN=confus_obs$table[1,2])
    
    # Generate randomizations
    randomizations <- 1:nruns + npts
    
    confus_rands <- lapply(randomizations, function(s) {
      confusionMatrix(table(+(!dnamat.pa[site,]),
                            +(!(1 * (stem.23$abund[[r]][s,] > 0)))))
    })
    
    # Calculate SES for each class
    conf_stats_ses[site,1:11] <- ses(conf_stats_obs[site,1:11], 
                                     do.call(rbind, lapply(confus_rands, 
                                                           function(x) x$byClass)))
    
    # Compute MCC of randomized points
    rand_mcc <- unlist(lapply(confus_rands, function(tab) {
      mltools::mcc(confusionM = matrix(tab$table, nrow = 2))
    }))
    
    conf_stats_ses[site,12] <- (conf_stats_obs[site,12] - mean(rand_mcc))/sd(rand_mcc)
    
  }
  
  # Convert matrices to data frames
  conf_stats_obs <- as.data.frame(conf_stats_obs)
  conf_stats_ses <- as.data.frame(conf_stats_ses)
  names(conf_stats_ses) <- names(conf_stats_obs) <- c(names(confus_obs$byClass), "MCC")
  
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


m2_obs <- aov(as.numeric(ba_obs) ~ as.factor(rep(1:20, each=npts)))
anova(m2_obs)


### However, there is no difference in SES balanced accuracy across the scales
### A slight hump appers in ~ 30-40 m range but these are not sig different from 
m2_ses <- aov(as.numeric(ba_ses) ~ as.factor(rep(1:20, each=npts)))
anova(m2_ses)
TukeyHSD(m2_ses)

### Median Sensitivity, Specificity & Balanced Accuracy at 5 m
median(lapply(conf_stats_obs_list, function(x) x$Sensitivity)[[1]])
median(lapply(conf_stats_obs_list, function(x) x$Specificity)[[1]])
median(lapply(conf_stats_obs_list, function(x) x$`Balanced Accuracy`)[[1]])

### Median Sensitivity, Specificity & Balanced Accuracy at 100 m
median(lapply(conf_stats_obs_list, function(x) x$Sensitivity)[[20]])
median(lapply(conf_stats_obs_list, function(x) x$Specificity)[[20]])
median(lapply(conf_stats_obs_list, function(x) x$`Balanced Accuracy`)[[5]])
median(lapply(conf_stats_obs_list, function(x) x$`MCC`)[[20]])



### Figure 4

# pdf("Figures/Fig4.Confusion-matrix-stats.pdf", width = 8, height = 10)

pdf(paste0("Figures/compare_filtering/Fig4.Confusion-matrix-stats_", 
           label[data_selector],
           ".pdf"), width = 8, height = 10)


par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

## Sensitivity (true positive rate)
boxplot(lapply(conf_stats_obs_list, function(x) x$Sensitivity), 
        ylab="Sensitivity (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()
mtext("A", adj=0, line=0.5)

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
mtext("D", adj=0, line=0.5)

## Specificity (true negative rate)
boxplot(lapply(conf_stats_obs_list, function(x) x$Specificity), 
        ylab="Specificity (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()
mtext("B", adj=0, line=0.5)

boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
        ylab="Specificity (SES)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
        axes=F, col=cp, add=T)
graphics::box()
mtext("E", adj=0, line=0.5)

## Balanced Accuracy
# boxplot(lapply(conf_stats_obs_list, function(x) x$`Balanced Accuracy`), 
#         ylab="Balanced Accuracy (observed)", 
#         xlab="Radius (m)", axes=F, col=cp)
# axis(1, labels=names(stem.23$abund), at=1:20)
# axis(2)
# graphics::box()
# mtext("C", adj=0, line=0.5)
# 
# boxplot(lapply(conf_stats_ses_list, function(x) x$`Balanced Accuracy`), 
#         ylab="Balanced Accuracy (SES)", 
#         xlab="Radius (m)", axes=F, col=cp)
# axis(1, labels=names(stem.23$abund), at=1:20)
# axis(2)
# polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
# boxplot(lapply(conf_stats_ses_list, function(x) x$`Balanced Accuracy`), 
#         axes=F, col=cp, add=T)
# graphics::box()
# mtext("F", adj=0, line=0.5)

## MCC
boxplot(lapply(conf_stats_obs_list, function(x) x$`MCC`), 
        ylab="MCC (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()
mtext("C", adj=0, line=0.5)

boxplot(lapply(conf_stats_ses_list, function(x) x$`MCC`),
        ylab="MCC (SES)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$`MCC`), 
        axes=F, col=cp, add=T)
graphics::box()
mtext("F", adj=0, line=0.5)


dev.off()

}

plot(sapply(conf_stats_ses_list, function(x) median(x$`Balanced Accuracy`)))


plot(do.call(c, lapply(conf_stats_ses_list, function(x) x$`Balanced Accuracy`)),
     do.call(c, lapply(conf_stats_ses_list, function(x) x$`MCC`)))

### Explore the plots where SES Balanced Accuracy is higher than expected

stem.23$abund[[1]][which(conf_stats_ses_list[[1]]$`Balanced Accuracy`>1.96),]



plot(jitter(rowSums(stem.23$abund[[1]]>0)[1:npts]),
     conf_stats_ses_list[[1]]$`Balanced Accuracy`,
     pch=21, bg=ifelse(conf_stats_ses_list[[1]]$`Balanced Accuracy`>1.96, 'red', 1),
     xlab='Taxon richness',
     ylab="SES Balanced accuracy")
abline(h=1.96)

t.test(conf_stats_ses_list[[1]]$`Balanced Accuracy`,
       rep(1.96, times=length(conf_stats_ses_list[[1]]$`Balanced Accuracy`)))















########################################################
########################################################
##################
### LAND USE ZONE ANALYSIS
##################
########################################################
########################################################

### NOT TOTALLY SURE WHAT TO DO WITH THIS PART.... LET'S REVISIT AFTER CESC'S THESIS

# Identify LU history for each sample site
low <- rownames(sampxy)[sampxy$CoverClass < 4]
high <- rownames(sampxy)[sampxy$CoverClass == 4]

# DNA species richness in different land-use categories
sum(1 * (colSums(dnamat[sampxy$CoverClass < 4,])>0))
sum(1 * (colSums(dnamat[sampxy$CoverClass == 4,])>0))

# Stem species richness in different land-use categories
sum(lfdp23$low_LU_abund>0)
sum(lfdp23$high_LU_abund>0)


rbind(table(rownames(lfdp23), lfdp23$low_LU_abund>0)[,2],
      table(rownames(lfdp23), lfdp23$high_LU_abund>0)[,2])






########################################################
########################################################
##################
### OLD STUFF FOR FULL PLOT
##################
########################################################
########################################################


# Total abundance of stems per gOTU (as log count in 100 m radius around 39 points) vs.number of sites where gOTU was detected in the DNA
# plot(lfdp23$total_abund, colSums(dnamat.pa), log='x')
# plot(lfdp16$total_abund, colSums(dnamat.pa), log='x')
# text(lfdp23$total_abund, colSums(dnamat.pa), labels=colnames(dnamat.pa))

# Significant positive correlation (i.e., more abundant gOTUs occur in more DNA samples)
(cor.test(colSums(dnamat.pa), log10(lfdp23$total_abund)))

# Total abundance of stems per gOTU (as count in 100 m radius around 39 points) vs.
plot(lfdp23$total_abund, 1*(colSums(dnamat.pa)>0), log='x')

# Significant logistic regression (i.e., more abundant gOTUs are more likely to have been detected overall in DNA)
a <- log10(lfdp23$total_abund)
m1 <- glm(1*(colSums(dnamat.pa)>0) ~ a, family="binomial")
summary(m1)

# Rank total abundance in stem data (based on 100 m radii) vs DNA read rank abundance
plot(rank(lfdp23$total_abund), rank(colSums(dnamat.pa)))

# Rank stem abundance is significantly positively correlated with rank DNA read abundance
(rankcor <- cor.test(rank(lfdp23$total_abund), rank(colSums(dnamat.pa))))
rank_corr_est <- rankcor$estimate
rank_corr_p <- rankcor$p.value

### Species accumulation in spatial aggregations of DNA samples
# dxy <- dist(sampxy, 'euclidean')
# dxy <- st_as_sf(sampxy, coords = c("X", "Y"), remove=FALSE)

# Find neighbors in different distance radii
# out <- list()
# rads <- seq(50,600,25)
# for(r in seq_along(rads)){
#   dxynn <- dnearneigh(dxy, 0, rads[r])
#   grps <- unique(lapply(include.self(dxynn), sort))
#   tmp <- vector()
#   for(g in seq_along(grps)){
#     if(length(grps[[g]])>1){
#       tmp[g] <- sum(colSums(dnamat.pa[grps[[g]],])>0)
#     } else {
#       tmp[g] <- sum(dnamat.pa[grps[[g]],]>0)
#     }
#   }
#   out[[r]] <- tmp
# }


# dnaaccum <- specaccum(dnamat.pa)
# 
# fitsp1 <- fitspecaccum(dnaaccum, model="michaelis-menten")
# fitsp2 <- fitspecaccum(dnaaccum, model="gleason")
# fitsp3 <- fitspecaccum(dnaaccum, model="asymp")
# fitsp4 <- fitspecaccum(dnaaccum, model="gompertz")
# fitsp5 <- fitspecaccum(dnaaccum, model="weibull")
# fitsp6 <- fitspecaccum(dnaaccum, model="michaelis-menten")
# fitsp7 <- fitspecaccum(dnaaccum, model="gitay")
# 
# plot(predict(fitsp1, newdata = 1:1000), type='l', ylim=c(0, 80))
# abline(h=78, lty=2)
# lines(predict(fitsp2, newdata = 1:1000), col=2)
# lines(predict(fitsp3, newdata = 1:1000), col=3)
# lines(predict(fitsp4, newdata = 1:1000), col=4)
# lines(predict(fitsp5, newdata = 1:1000), col=5)
# lines(predict(fitsp6, newdata = 1:1000), col=6)
# lines(predict(fitsp7, newdata = 1:1000), col=7)


# library(fossil)
# fossil::chao2(dnamat.pa, taxa.row = F)

