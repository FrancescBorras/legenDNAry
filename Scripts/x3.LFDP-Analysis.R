### Load libraries
library(vegan)
library(caret)
library(compositions)
library(viridis)
library(ecodist)


### Load summary data
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

sampxy$quad <- sample_data$largeGridlocation[match(rownames(dnamat.pa), rownames(sample_data))]


##################
### NOTES
##################

# Some 'first pass' exploratory analyses are below, separated by code blocks
# The order is meaningless

# SOME IDEAS OF THINGS TO DO

# - rank abundance of stem vs. dna data

# - Is there a land-use signature? 
#     - Pool sites by land use and look at summary stats
#     - PCA on DNA data and then look at land-use, soil type


########################################################
### SAC of DNA data
########################################################

library(spdep)
library(sf)
library(dplyr)
library(ggplot2)

plot(sampxy)

dxy <- dist(sampxy, 'euclidean')

# convert dataframe to sf object,
dxy <- st_as_sf(sampxy, coords = c("X", "Y"), remove=FALSE)

# find distance-based neigbours, upper distance bound is set ti 0.2km,
out <- list()
rads <- seq(50,600,25)
for(r in seq_along(rads)){
  
  dxynn <- dnearneigh(dxy, 0, rads[r])
  
  grps <- unique(lapply(include.self(dxynn), sort))
  
  tmp <- vector()
  for(g in seq_along(grps)){
    if(length(grps[[g]])>1){
      tmp[g] <- sum(colSums(dnamat.pa[grps[[g]],])>0)
    } else {
      tmp[g] <- sum(dnamat.pa[grps[[g]],]>0)
    }
  }
  out[[r]] <- tmp
}


par(mfrow=c(2,1), mar=c(4,4,1,1))

df <- data.frame(sapply(out, "length<-", max(lengths(out))))
colnames(df) <- rads
boxplot(df, col=rev(viridis(ncol(df))), 
        xlab="Aggregating Distance (m)", 
        ylab="Species richness",
        main="DNA sample gOTU accumulation curve")


b <- boxplot(sapply(stem.23$ba, function(x) rowSums(x[1:39,]>0)), col=cp,
             xlab="Radius (m)", ylab="Species richness (2023)", 
             xlim=c(-0.5,20), ylim=c(0,80))
boxplot(rowSums(dnamat.pa), add=T, at=-0.5, width=2, col=2)
abline(v=0.25)
abline(h=79, lty=2)
axis(1, at=-0.5, labels="DNA")

b <- boxplot(sapply(stem.16$ba, function(x) rowSums(x[1:39,]>0)), col=cp,
             xlab="Radius (m)", ylab="Species richness (2016)", 
             xlim=c(-0.5,20), ylim=c(0,80))
boxplot(rowSums(dnamat.pa), add=T, at=-0.5, width=2, col=2)
abline(v=0.25)
abline(h=79, lty=2)
axis(1, at=-0.5, labels="DNA")






########################################################
### Ordination of DNA matrix (testing)
########################################################

dna.ord <- metaMDS(dnamat.pa, distance = 'raup')

plot(dna.ord, display="sites", 
     xlim=range(dna.ord$species, na.rm=T),
     ylim=range(dna.ord$species, na.rm=T))

par(mfrow=c(1,2), mar=c(4,4,1,1))
plot(sampxy[,1:2], pch=15, 
     cex=4, col=viridis(10)[cut(dna.ord$points[,1], 10)])
plot(sampxy[,1:2], pch=15, 
     cex=4, col=viridis(10)[cut(dna.ord$points[,2], 10)])

lu <- factor(c(rep('a', 20), rep('b', 19)))
stem.ord <- cca(1*(stem.23$abund[[3]][1:39,]>0) ~ lu)
plot(stem.ord, display=c("wa", 'cn'))

dna.ord <- metaMDS(dnamat.pa)
stem.ord <- metaMDS(1*(stem.23$abund[[3]][1:39,]>0))


# Remove columns for species that were not observed
dnamat.pa2 <- dnamat.pa[,colSums(dnamat.pa)!=0]

# Create a matrix of pairwise dissimiliarity between quadrats
sp.dist <- vegdist(1*(stem.23$abund[[3]][1:39,]>0),
                   method = "bray") # distance metric

sp.dist <- vegdist(dnamat.pa2,
                   method = "bray") # distance metric



dnamat.pa


### GET LAND-USE HISTORY
physenvir <- read.table("/Users/au529793/Projects/Puerto Rico/Spatial Aggregation/DATA/LFDPEnvironment20.txt", header=T, row.names = 1)
physenvir$history <- as.factor(ifelse(physenvir$CoverClass<4, "H", "L"))

sampxy$history <- physenvir$history[match(sampxy$quad, physenvir$Quadrat)]
# sampxy$history <- as.numeric(sampxy$history)

### SIGNIFICANTLY HIGHER gOTU RICHNESS DETECTED IN 

boxplot(rowSums(dnamat.pa[LU=="L",]),
        rowSums(dnamat.pa[LU=="H",]), col=c(4,2), names=c("low","high"))

t.test(rowSums(dnamat.pa[LU=="L",]),
       rowSums(dnamat.pa[LU=="H",]))




### NO SIGNIFICANT DIFFERENCE BETWEEN SPECIES RICHNESS
x <- do.call(rbind, lapply(stem.pa.list, function(x) rowSums(x[1:39,][LU=="L",])))
y <- do.call(rbind, lapply(stem.pa.list, function(x) rowSums(x[1:39,][LU=="H",])))

for(i in 1:20){
  if(i==1){
    b <- boxplot(x[1,],y[1,], xlim=c(1, 40), ylim=c(0,80), col=c(4,2), axes=F,
                 xlab="Radius (m)", ylab="Species richness", main="2023")
  } else {
    ats <- seq(1, 44, 2)
    boxplot(x[i,],y[i,], add=T, at=c(ats[i], ats[i]+1), col=c(4,2), axes=F)
  }
}
axis(1, labels = seq(5, 100, 5), at=seq(1.5, 39.5, 2))
axis(2)

legend("topleft", legend=c("Low LU", "High LU"), bty='n', 
       pch=22, pt.bg=c(4,2), pt.cex=2)


### NO SIGNIFICANT DIFFERENCE BETWEEN SPECIES RICHNESS
x <- do.call(rbind, lapply(stem.pa.list16, function(x) rowSums(x[1:39,][LU=="L",])))
y <- do.call(rbind, lapply(stem.pa.list16, function(x) rowSums(x[1:39,][LU=="H",])))

for(i in 1:20){
  if(i==1){
    b <- boxplot(x[1,],y[1,], xlim=c(1, 40), ylim=c(0,80), col=c(4,2), axes=F,
                 xlab="Radius (m)", ylab="Species richness", main="2016")
  } else {
    ats <- seq(1, 44, 2)
    boxplot(x[i,],y[i,], add=T, at=c(ats[i], ats[i]+1), col=c(4,2), axes=F)
  }
}
axis(1, labels = seq(5, 100, 5), at=seq(1.5, 39.5, 2))
axis(2)

legend("topleft", legend=c("Low LU", "High LU"), bty='n', 
       pch=22, pt.bg=c(4,2), pt.cex=2)






### 2023
par(mar=c(5,4,1.5,1.5))
x <- do.call(rbind, lapply(stem.pa.list, function(x) rowSums(x[1:39,][LU=="L",])))
y <- do.call(rbind, lapply(stem.pa.list, function(x) rowSums(x[1:39,][LU=="H",])))
boxplot(rowSums(dnamat.pa[LU=="L",]),
        rowSums(dnamat.pa[LU=="H",]), col=c(4,2),
        xlim=c(1.5, 42.5), ylim=c(0,70), axes=F,
        xlab="Radius (m)", ylab="Species richness")
for(i in 1:20){
    ats <- seq(4, 46, 2)
    boxplot(x[i,],y[i,], add=T, at=c(ats[i], ats[i]+1), col=c(4,2), axes=F)
  }
axis(1, labels = seq(5, 100, 5), at=seq(4.5, 42.5, 2))
axis(1, labels = "DNA", at=1.5)
axis(2)
graphics::box()
abline(v=3)
legend("bottomright", legend=c("Low LU", "High LU"), bty='n', 
       pch=22, pt.bg=c(4,2), pt.cex=2)

### 2016
boxplot(rowSums(dnamat.pa[LU=="L",]),
        rowSums(dnamat.pa[LU=="H",]), col=c(4,2),
        xlim=c(1.5, 42.5), ylim=c(0,70), axes=F,
        xlab="Radius (m)", ylab="Species richness")
x <- do.call(rbind, lapply(stem.pa.list16, function(x) rowSums(x[1:39,][LU=="L",])))
y <- do.call(rbind, lapply(stem.pa.list16, function(x) rowSums(x[1:39,][LU=="H",])))
for(i in 1:20){
  ats <- seq(4, 46, 2)
  boxplot(x[i,],y[i,], add=T, at=c(ats[i], ats[i]+1), col=c(4,2), axes=F)
}
axis(1, labels = seq(5, 100, 5), at=seq(4.5, 42.5, 2))
axis(1, labels = "DNA", at=1.5)
axis(2)
graphics::box()
abline(v=3)
legend("bottomright", legend=c("Low LU", "High LU"), bty='n', 
       pch=22, pt.bg=c(4,2), pt.cex=2)







dna.ord <- metaMDS(dnamat.pa)
stem.ord <- metaMDS(1*(stem.23$abund[[1]][1:39,]>0))

dna.efit <- envfit(dna.ord ~ sampxy$history)
stem.efit <- envfit(stem.ord ~ sampxy$history)

plot(dna.ord)
plot(dna.efit, add=T)

plot(stem.ord)
plot(stem.efit, add=T)

LU <- sampxy$history
dna.ord <- cca(dnamat.pa ~ LU)
plot(dna.ord, display=c("cn"))

stem.ord <- cca(1*(stem.23$abund[[3]][1:39,]>0) ~ LU)
plot(stem.ord, display=c("wa","cn"))




########################################################
### Mantel test of stem distance matrix and dna distance matrix
########################################################

# Note: Lots to consider about how to transform the data and which distance metric to use...

# Note: If ecodist::mantel() function gives errors about dimensions, then you have to restart R and load ecodist package again. Another package ('proxies'?) seems to change the handling of distance objects.

# dnamat.dist <- vegdist(clr(dnamat), method='euclidean')
dnamat.dist <- vegdist(dnamat.pa, method='raup')
# dnamat.dist <- vegdist(dnamat.pa[,colSums(dnamat.pa)>0], method='raup')

manres <- matrix(nrow=20, ncol=6)
for(r in 1:20){
  stem.dist <- vegdist(stem.23$ba[[r]][1:39,]>0, method='raup')
  # stem.dist <- vegdist(decostand(stem.23$ba[[r]][1:39,], "total", 1), method='jaccard')
  manres[r,] <- round(ecodist::mantel(dnamat.dist ~ stem.dist), 3)
}

colnames(manres) <- names(ecodist::mantel(dnamat.dist ~ stem.dist))
rownames(manres) <- seq(5, 100, 5)
manres <- as.data.frame(manres)

plot(rownames(manres), manres$mantelr, ylim=range(manres[,5:6]), 
     pch=ifelse(manres$pval3<0.05, 21, 16), 
     xlab="Radius (m)", ylab="Mantel correlation",
     main="eDNA p/a dissimilarity vs. Stem basal area dissimilarity")
segments(seq(5, 100, 5), manres$`llim.2.5%`, seq(5, 100, 5), manres$`ulim.97.5%`)
points(rownames(manres), manres$mantelr, ylim=range(manres[,5:6]), 
     pch=ifelse(manres$pval3<0.05, 21, 16), bg=ifelse(manres$pval3<0.05, 'white', 1),
     xlab="Radius (m)", ylab="Mantel correlation")
abline(h=0, lty=2)

# Essentially, it seems that sites that are more dissimilar from one another in terms of DNA presence/absence data are also more dissimilar in terms of stem presence/absence for spatial radii of 10,15,25,30 m.  These dissimilarities are not significantly correlated for other spatial scales.


########################################################
### Procrusties test of stem ordination and dna ordinatino
########################################################

# Similar to Mantel test, lots to consider in terms of transformation, dist metric, ordination...

prores <- list()
prores16 <- list()

for(r in 1:20){
  
  print(r)
  dnaord <- rda(dnamat.pa, scale=T)
  # dnaord <- metaMDS(clr(dnamat), dist='euclidean')
  
  stemord <- rda(1*(stem.23$abund[[r]][1:39,]>0), scale=T)
  stemord16 <- rda(1*(stem.16$abund[[r]][1:39,]>0), scale=T)
  
  # stemord <- metaMDS(1*(stem.23$abund[[r]][1:39,]>0), distance="raup")
  # stemord <- metaMDS(stem.23$ba[[r]][1:39,])

  prores[[r]] <- protest(dnaord, stemord, permutations = 999)
  prores16[[r]] <- protest(dnaord, stemord16, permutations = 999)
  # prores[[r]] <- protest(dnamat.pa, 1*(stem.23$abund[[r]][1:39,]>0))
  # prores[[r]] <- protest(clr(dnamat), stem.23$abund[[r]][1:39,])
  
}

# A 'significant' p value (<0.05) indicates that the two matrices are similar
# A 'nonsignificant' p value (>0.05) indicates that the two matrices are different

# A smaller Procrustes SS indicates a better fit or alignment between the two sets of points. It essentially tells you how well one set of points can be adjusted to match another set, considering only rotation, scaling, and translation as transformation operations.

# vals <- sapply(prores, function(x) x$t0)
vals <- sapply(prores, function(x) summary(permustats(x))$z)
# low <- sapply(prores, function(x) quantile(summary(permustats(x))$permutations, 0.025))
# hi <- sapply(prores, function(x) quantile(summary(permustats(x))$permutations, 0.975))
sig <- sapply(prores, function(x) x$signif)

# vals16 <- sapply(prores16, function(x) x$t0)
vals16 <- sapply(prores16, function(x) summary(permustats(x))$z)
low16 <- sapply(prores16, function(x) quantile(summary(permustats(x))$permutations, 0.025))
hi16 <- sapply(prores16, function(x) quantile(summary(permustats(x))$permutations, 0.975))
sig16 <- sapply(prores16, function(x) x$signif)

plot(vals, vals16, xlab="Procrusties correlation 2023",
     ylab="Procrusties correlation 2016", pch=21, bg=scales::alpha(cp, 0.8),
     cex=2)
segments(vals, low16, vals, hi16)
segments(low, vals16, hi, vals16)
abline(0,1, lty=2)

abline(-1.96, 1)
abline(1.96, 1)



legend('topleft', legend=c("5m", "50m", "100m"), 
       pt.bg=cp[c(1,10,20)], pch=21, bty='n')

summary(permustats(prores[[1]]))


# With the presence-absence data matrices (both dna and stem), p<0.05 for all spatial scales
# With the metaMDS ordinations on presence-absence data, p>0.05 for 5m
# Basically most evidence points to the two matrices being similar across scales with p/a data.
# With abundance data (clr(dna) and raw stem abund), all p vals are >0.05 (matrices are diff)


### DEMO OF HOW THE PROTEST WORKS...
# set.seed(123) 
# self <- matrix(rnorm(2000), ncol = 200, nrow = 10) 
# other <- matrix(rnorm(2000), ncol = 200, nrow = 10) 
# protest(self, self, permutations=100)
# protest(self, other, permutations=100)




########################################################
### Whole plot rank abundance vs rank abundance of DNA data
########################################################

# "Whole plot" stem rank abundance and rank DNA read abundance is similar (~0.5) across spatial scales
# This *should* be based on the sum total stems for the full plot 
# but we don't have that for 2023 census so this uses sum of stems in our summary data

est <- matrix(nrow=20, ncol=3)

for(r in 1:20){
  focdat <- stem.23$ba[[r]]
  res <- cor.test(rank(colSums(focdat)), rank(colSums(dnamat)))
  est[r,] <- c(res$estimate, res$conf.int)
}

plot(est[,1], ylim=c(-1,1))
segments(1:20, est[,2], 1:20, est[,3])
abline(h=0, lty=2)




########################################################
### COMPUTE BALANCED ACCURACY 
### (INCLUDING STANDARDIZED EFFECT SIZE OF BALANCED ACCURACY)
### **Note that this takes ~30 min to run on Bob's computer!!!**
########################################################

# balaccuracy.obs <- matrix(nrow=39, ncol=20)
# accuracy.obs <- matrix(nrow=39, ncol=20)
# balaccuracy.ses <- matrix(nrow=39, ncol=20)

conf.stats.obs.list <- list()
conf.stats.ses.list <- list()

### Radius loop
for(r in 1:20){
  message(paste("radius =", r))
  conf.stats.obs <- matrix(nrow=39, ncol=11)
  conf.stats.ses <- matrix(nrow=39, ncol=11)

  ### Site loop  
  for(site in 1:39){
    confus.obs <- confusionMatrix(table(dnamat.pa[site,], 
                                    1*(stem.23$abund[[r]][site,] > 0)))
    # balaccuracy.obs[site,r] <- (confus.obs$byClass[1] + confus.obs$byClass[2])/2
    # accuracy.obs[site,r] <- confus.obs$overall[1]
    conf.stats.obs[site,] <- confus.obs$byClass
    
    conf.stats.rand <- matrix(nrow=1000, ncol=ncol(conf.stats.obs))
    
    ### Randomization loop
    for(s in 1:50){
      confus.rand <- confusionMatrix(table(dnamat.pa[site,],
                                           1*(stem.23$abund[[r]][s+39,] > 0)))
      conf.stats.rand[s,] <- confus.rand$byClass
    }
    
    conf.stats.ses[site,] <- (conf.stats.obs[site,] - apply(conf.stats.rand, 2, mean, na.rm=T)) / apply(conf.stats.rand, 2, sd, na.rm=T)
    
    # balaccuracy.ses[site,r] <- (balaccuracy.obs[site,r] - mean(balaccuracy.rand)) / sd(balaccuracy.rand)
  }
  
  conf.stats.obs <- as.data.frame(conf.stats.obs)
  names(conf.stats.obs) <- names(confus.obs$byClass)
  conf.stats.ses <- as.data.frame(conf.stats.ses)
  names(conf.stats.ses) <- names(confus.rand$byClass)
  
  conf.stats.obs.list[[r]] <- conf.stats.obs
  conf.stats.ses.list[[r]] <- conf.stats.ses
  
}


### PLOT SES OF A FEW METRICS
par(mfrow=c(2,2))
boxplot(lapply(conf.stats.obs.list, function(x) x$Sensitivity), 
        main="Obs Sensitivity", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)

boxplot(lapply(conf.stats.obs.list, function(x) x$Specificity), 
        main="Obs Specificity", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)

boxplot(lapply(conf.stats.obs.list, function(x) x$`Balanced Accuracy`), 
        main="Obs Balanced Accuracy", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)

boxplot(lapply(conf.stats.ses.list, function(x) x$Sensitivity), 
        main="SES Sensitivity", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
abline(h=c(-1.96, 1.96), lwd=2, lty=2)

boxplot(lapply(conf.stats.ses.list, function(x) x$Specificity), 
        main="SES Specificity", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
abline(h=c(-1.96, 1.96), lwd=2, lty=2)

boxplot(lapply(conf.stats.ses.list, function(x) x$`Balanced Accuracy`), 
        main="SES Balanced Accuracy", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
# abline(h=c(-1.96, 1.96), lwd=2, lty=2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf.stats.ses.list, function(x) x$`Balanced Accuracy`), 
        axes=F, col=cp, add=T)








# Plot the results
boxplot(balaccuracy.ses, main="SES of Balanced Accuracy", 
        xlab="Radius (m)", axes=F, col=cp)
polygon(c(-1,200,200,-1), c(-1.96,-1.96, 1.96, 1.96), col='grey', lty=2)
boxplot(balaccuracy.ses, col=cp, add=T, axes=F)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
box()

par(mfrow=c(2,3))

# Look at raw sensitivity results
boxplot(lapply(conf.stats.obs.list, function(x) x$`Sensitivity`), 
        main="Observed Sensitivity", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
box()

# Look at raw specificity results
boxplot(lapply(conf.stats.obs.list, function(x) x$`Specificity`), 
        main="Observed Specificity", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
box()

# Look at obs. balanced accuracy results
names(conf.stats.obs.list[[1]])
boxplot(lapply(conf.stats.obs.list, function(x) x$`Balanced Accuracy`), 
        main="Observed Balanced Accuracy", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
box()








########################################################
## COMMUNITY DISSIMILIARITY ANALYSIS
########################################################

## Gives warnings with 'raup' because there are 'empty' species
res <- matrix(nrow=39, ncol=20)

# Loop through spatial radii
for(r in 1:20){
  # Remove taxa with all zeros in both soil and stem data for this radius
  tmp <- rbind(dnamat, stem.23$abund[[r]][1:39,])
  tmp <- tmp[,which(colSums(tmp)!=0)]
  # Compute the index
  tmpres <- as.matrix(vegdist(tmp>0, method="raup", binary = TRUE))
  tmpres <- tmpres[1:39, 40:ncol(tmpres)]
  res[,r] <- diag(tmpres)
}

boxplot(res, ylab="Raup Dissimilarity", xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)


### Dissimilarity using Chase et al. 2011 null model and Raup 
### NOTE: Different null models (r1, r0) give different responses
res_rcnull <- matrix(nrow=39, ncol=20)
for(r in 1:20){
  message(paste("Computing r=", r))
  # Remove taxa with all zeros in both soil and stem data for this radius
  tmp <- rbind(dnamat, stem.23$abund[[r]][1:39,])
  tmp <- tmp[,which(colSums(tmp)!=0)]
  tmp <- as.matrix(raupcrick(tmp>0, null = "r0"))
  tmp <- as.matrix(tmp)[1:39, 40:ncol(tmp)]
  res_rcnull[,r] <- diag(tmp)
}

boxplot(res_rcnull, ylab="Raup-Crick Null Dissimilarity", xlab="Radius (m)", axes=F)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)


########################################################
### LOOK AT DETECTION PROBABILITY VS BASAL AREA IN A NEIGHBORHOOD
########################################################

## We need to think about how to 

# Choose a certain spatial scale for convenience
focdat <- stem.23$ba[[20]]

plot(range(focdat), c(0,1), pch=NA,
     ylab="Probability of eDNA detection",
     xlab="Basal area (units?)", xlim=c(0,0.1))

for(sp in 1:ncol(focdat)){
  y <- (1*(dnamat>0)[,sp])
  x <- focdat[1:39,sp]
  mod <- glm(y ~ x, family="binomial")
  coeffs[sp] <- coef(mod)[2]
  nd <- data.frame(x=seq(0, max(focdat), length.out=1000))
  ypred <- predict(mod, nd, type='response')
  lines(nd$x, ypred, col=ifelse(coeffs[sp]>0, "blue", "red"))
}


#####
### RELATE DETECTION TO  IS DETECTION SLOPE RELATED TO LEAF TRAITS?
#####

# One way to do this would be to say, e.g., at what basal area does each species 
# have a 50% predicted prob of being detected?
# Then we can correlate those values with traits.

# Another way would be to include traits directly in the regression as interaction term

# We need more thought on this...

library(lme4)

y <- as.vector((1*(dnamat>0)))
x <- scale(as.vector(as.matrix(focdat[1:39,])))
tz <- scale(rep(traits$WD, each=39))
tz[is.na(tz)] <- 0

summary(glm(y ~ x + tz + x*tz, family="binomial"))





########################################################
### RAW ABUNDANCE OF gOTU sequences vs. stems
########################################################
focdat <- stem.23$abund[[5]][1:39,]

focdat.ra <- decostand(focdat, 1, method = "total")
dnamat.ra <- decostand(dnamat, 1, method = "total")

x <- data.frame(a=unlist(focdat.ra),
                b=unlist(dnamat.ra))
par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(x, 
     xlab="Stem abundance (25 m radius)",
     ylab="Number of reads",
     pch=21, bg=rgb(0,0,0,0.5),
     main="untransformed"
)
# abline(lm((x2$b) ~ (x2$a)), col='blue', lwd=3)

plot(x, 
     xlab="Stem abundance (25 m radius)",
     ylab="Number of reads",
     pch=21, bg=rgb(0,0,0,0.5), log='xy',
     main="log-log"
)
# abline(lm(log10(x2$b) ~ log10(x2$a)), col='blue', lwd=3)


# Look across spatial scales
out <- matrix(nrow=20, ncol=5)
for(i in seq_along(stem.23$abund)){
  focdat <- stem.23$abund[[i]][201:239,]
  
  focdat.ra <- decostand(focdat, 1, method = "total")
  dnamat.ra <- decostand(dnamat, 1, method = "total")
  
  x <- data.frame(a=unlist(focdat.ra), 
                  b=unlist(dnamat.ra))
  
  # x2 <- x[rowSums(apply(x, 2, function(x)x==0))==0,]
  # s <- summary(lm(log10(x2$b) ~ log10(x2$a)))
  
  out[i,1] <- round(s$adj.r.squared, 2)
  out[i,2] <- round(s$coefficients[2,1], 2)
  out[i,3] <- round(cor.test(x$b, x$a, use='p')$estimate, 2)
  out[i,4:5] <- round(cor.test(x$b, x$a)$conf.int, 2)
}

out <- as.data.frame(out)
names(out) <- c("r2","coef", "pearsons")

par(mfrow=c(2,1))

plot(seq(5,100,5), out$pearsons, 
     xlab="Radius (m)", ylab="Pearsons correlation", 
     ylim=c(-0.01, 0.5), 
     main="Sequence vs. Stem relative abundance", col='red')
abline(h=0, lty=2)
segments(seq(5,100,5), out[,4], 
         seq(5,100,5), out[,5], col='red')



# Look across spatial scales
out <- matrix(nrow=20, ncol=5)
for(i in seq_along(stem.23$ba)){
  focdat <- stem.23$ba[[i]][201:239,]
  
  focdat.ra <- decostand(focdat, 1, method = "total")
  dnamat.ra <- decostand(dnamat, 1, method = "total")
  
  x <- data.frame(a=unlist(focdat.ra), 
                  b=unlist(dnamat.ra))
  
  out[i,1] <- round(s$adj.r.squared, 2)
  out[i,2] <- round(s$coefficients[2,1], 2)
  out[i,3] <- round(cor.test(x$b, x$a, use='p')$estimate, 2)
  out[i,4:5] <- round(cor.test(x$b, x$a)$conf.int, 2)
}

out <- as.data.frame(out)
names(out) <- c("r2","coef", "pearsons")

plot(seq(5,100,5), out$pearsons, 
     xlab="Radius (m)", ylab="Pearsons correlation", 
     ylim=c(-0.01, 0.5), 
     main="Sequence vs. Stem relative basal area", col='red')
abline(h=0, lty=2)
segments(seq(5,100,5), out[,4], 
         seq(5,100,5), out[,5], col='red')






########################################################
### SPECIES RICHNESS OF DNA AND STEM DATA
########################################################


# Quick plot for a given spatial scale (no rarefaction)
plot(rowSums(stem.23$abund[[20]][1:39,]>0),
     rowSums(dnamat.pa))

# Quick plot for a given spatial scale after rarefaction
plot(rarefy(stem.23$abund[[20]][1:39,], min(rowSums(stem.23$abund[[20]][1:39,]))),
     rarefy(dnamat, min(rowSums(dnamat))))

# Look at correlation coefficients across spatial scales for 2023 data
corrs <- list()
corrs_rare <- list()

for(dist in seq_along(stem.23$abund)){
  corrs[[dist]] <- cor.test(rowSums(dnamat>0), rowSums(stem.23$abund[[dist]][1:39,]>0))
  
  corrs_rare[[dist]] <- cor.test(rarefy(stem.23$abund[[dist]][1:39,], 
                                        min(rowSums(stem.23$abund[[dist]][1:39,]))),
                                 rarefy(dnamat, min(rowSums(dnamat))))
}


# Look at correlation coefficients across spatial scales for 2016 data
corrs16 <- list()
corrs_rare16 <- list()

for(dist in seq_along(stem.16$abund)){
  corrs16[[dist]] <- cor.test(rowSums(dnamat>0), rowSums(stem.16$abund[[dist]][1:39,]>0))
  
  corrs_rare16[[dist]] <- cor.test(rarefy(stem.16$abund[[dist]][1:39,], 
                                        min(rowSums(stem.16$abund[[dist]][1:39,]))),
                                 rarefy(dnamat, min(rowSums(dnamat))))
}


par(mfrow=c(1,2))

plot(seq_along(stem.23$abund), sapply(corrs, function(x) x$estimate), 
     ylim=c(-1,1), axes=F,
     pch=16, cex=1, xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness [raw]")
segments(seq_along(stem.23$abund), y0=sapply(corrs, function(x) x$conf.int[1]),
         y1=sapply(corrs, function(x) x$conf.int[2]))

points(seq_along(stem.16$abund)+0.5, sapply(corrs16, function(x) x$estimate), 
       pch=16, col=2)
segments(seq_along(stem.16$abund)+0.5, y0=sapply(corrs16, function(x) x$conf.int[1]),
         y1=sapply(corrs16, function(x) x$conf.int[2]), col=2)


abline(h=0, lty=2)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)
legend("top", legend=c("2023", "2016"), pch=16, col=1:2, bty='n')


plot(seq_along(stem.23$abund), sapply(corrs_rare, function(x) x$estimate), ylim=c(-1,1), axes=F,
     pch=16, xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness [rarefied]")
segments(seq_along(stem.23$abund), y0=sapply(corrs_rare, function(x) x$conf.int[1]),
         y1=sapply(corrs_rare, function(x) x$conf.int[2]))

points(seq_along(stem.16$abund)+0.5, sapply(corrs_rare16, function(x) x$estimate), 
       pch=16, col=2)
segments(seq_along(stem.16$abund)+0.5, 
         y0=sapply(corrs_rare16, function(x) x$conf.int[1]),
         y1=sapply(corrs_rare16, function(x) x$conf.int[2]), col=2)

abline(h=0, lty=2)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)

legend("top", legend=c("2023", "2016"), pch=16, col=1:2, bty='n')




########################################################
### PRESENCE IN DNA AS A FUNCTION OF DISTANCE TO NEAREST (KNOWN) STEM
########################################################

x <- unlist(log10(stem.23$nn[1:39,]))
y <- as.vector(dnamat.pa)

x <- unlist(log10(stem.23$nn[1:39,]))
y <- as.vector(dnamat.pa)


plot(x, jitter(y, 0.1), pch=16, col=rgb(0,0,0,0.2), 
     ylab="Presence in DNA",
     xlab="Distance to nearest individual (log10 m)")
mod <- glm(y ~ x, family="binomial")
nd <- data.frame(x=seq(0, max(x), length.out=100))
ypred <- predict(mod, nd, type='response')

# Overall fit
lines(nd$x, ypred, col="blue", lwd=3)

# gOTU-specific fits
coeffs <- vector()
for(sp in 1:ncol(stem.23$nn)){
  xx <- log(unlist(stem.23$nn[1:39,sp]))
  yy <- unlist(as.data.frame(1*(dnamat[,sp]>1)))
  mod <- glm(yy ~ xx, family=binomial)
  coeffs[sp] <- coef(mod)[2]
  nd <- data.frame(xx=seq(0, max(x), length.out=100))
  ypred <- predict(mod, nd, type='response')
  lines(nd$xx, ypred, col=ifelse(coeffs[sp]>0, rgb(1,0,0,0.4), rgb(0,0,1,0.4)))
}



