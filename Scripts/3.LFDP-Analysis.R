### Load libraries
library(vegan)
library(caret)
library(compositions)
library(viridis)
library(ecodist)

### Load summary data
dnamat <- readRDS("Processed_data/stem-soil-39pt-data-20240411.RDA")[[1]]
stem.otu <- readRDS("Processed_data/stem-soil-39pt-data-20240411.RDA")[[2]]
traits <- readRDS("Processed_data/stem-soil-39pt-data-20240411.RDA")[[3]]

### Make a DNA presence absence matrix
dnamat.pa <- 1*(dnamat>0)

### Color palette for plotting
cp <- rev(viridis(20))



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
### Mantel test of stem distance matrix and dna distance matrix
########################################################

# Lots to consider about how to transform the data and which distance metric to use...

# dnamat.dist <- vegdist(clr(dnamat), method='euclidean')
dnamat.dist <- vegdist(dnamat.pa, method='raup')
# dnamat.dist <- vegdist(dnamat.pa[,colSums(dnamat.pa)>0], method='raup')

manres <- matrix(nrow=20, ncol=6)
for(r in 1:20){
  stem.dist <- vegdist(stem.otu$ba[[r]][1:39,], method='raup')
  # stem.dist <- vegdist(decostand(stem.otu$ba[[r]][1:39,], "total", 1), method='jaccard')
  manres[r,] <- round(ecodist::mantel(dnamat.dist ~ stem.dist), 3)
}

colnames(manres) <- names(ecodist::mantel(dnamat.dist ~ stem.dist))
rownames(manres) <- seq(5, 100, 5)
manres <- as.data.frame(manres)

plot(rownames(manres), manres$mantelr, ylim=range(manres[,5:6]), pch=16, 
     xlab="Radius (m)", ylab="Mantel correlation")
segments(seq(5, 100, 5), manres$`llim.2.5%`, seq(5, 100, 5), manres$`ulim.97.5%`)
abline(h=0, lty=2)


########################################################
### Procrusties test of stem ordination and dna ordinatino
########################################################

# Similar to Mantel test, lots to consider in terms of transformation, dist metric, ordination...

prores <- list()

for(r in 1:20){
  
  dnaord <- metaMDS(dnamat.pa, dist="raup")
  # dnaord <- metaMDS(clr(dnamat), dist='euclidean')
  
  stemord <- metaMDS(vegdist(stem.otu$ba[[r]][1:39,], method='raup'))
  # stemord <- metaMDS(stem.otu$ba[[r]][1:39,])
                       
  prores[[r]] <- protest(stemord, dnaord)
}

plot(sapply(prores, function(x) x$signif), ylim=c(0,1))
abline(h=c(0.05, 0.95))


########################################################
### Whole plot rank abundance vs rank abundance of DNA data
########################################################

# "Whole plot" stem rank abundance and rank DNA read abundance is similar (~0.5) across spatial scales
# This *should* be based on the sum total stems for the full plot 
# but we don't have that for 2023 census so this uses sum of stems in our summary data

est <- matrix(nrow=20, ncol=3)

for(r in 1:20){
  focdat <- stem.otu$ba[[r]]
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

balaccuracy.obs <- matrix(nrow=1000, ncol=20)
balaccuracy.ses <- matrix(nrow=39, ncol=20)

for(r in 1:20){
  message(paste("radius =", r))
  
  for(site in 1:39){
    
    confus.obs <- confusionMatrix(table(dnamat.pa[site,], 
                                    1*(stem.otu$abund[[r]][site,] > 0)))
    balaccuracy.obs[site,r] <- (confus.obs$byClass[1] + confus.obs$byClass[2])/2
    
    balaccuracy.rand <- vector()

    for(s in 1:1000){
      confus.rand <- confusionMatrix(table(dnamat.pa[site,],
                                           1*(stem.otu$abund[[r]][-c(1:39),] > 0)[s,]))

      balaccuracy.rand[s] <- (confus.rand$byClass[1] + confus.rand$byClass[2])/2
    }

    balaccuracy.ses[site,r] <- (balaccuracy.obs[site,r] - mean(balaccuracy.rand)) / sd(balaccuracy.rand)
  }
}

# Plot the results
boxplot(balaccuracy.ses, main="SES of Balanced Accuracy", 
        xlab="Radius (m)", axes=F, col=cp)
polygon(c(-1,200,200,-1), c(-1.96,-1.96, 1.96, 1.96), col='grey', lty=2)
boxplot(balaccuracy.ses, col=cp, add=T, axes=F)
axis(1, labels=names(stem.otu$abund), at=1:20)
axis(2)
box()

# Look at raw balanced accuracy results
boxplot(balaccuracy.obs, main="Observed Balanced Accuracy", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.otu$abund), at=1:20)
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
  tmp <- rbind(dnamat, stem.otu$abund[[r]][1:39,])
  tmp <- tmp[,which(colSums(tmp)!=0)]
  # Compute the index
  tmpres <- as.matrix(vegdist(tmp>0, method="raup", binary = TRUE))
  tmpres <- tmpres[1:39, 40:ncol(tmpres)]
  res[,r] <- diag(tmpres)
}

boxplot(res, ylab="Raup Dissimilarity", xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.otu$abund), at=1:20)
axis(2)


### Dissimilarity using Chase et al. 2011 null model and Raup 
### NOTE: Different null models (r1, r0) give different responses
res_rcnull <- matrix(nrow=39, ncol=20)
for(r in 1:20){
  message(paste("Computing r=", r))
  # Remove taxa with all zeros in both soil and stem data for this radius
  tmp <- rbind(dnamat, stem.otu$abund[[r]][1:39,])
  tmp <- tmp[,which(colSums(tmp)!=0)]
  tmp <- as.matrix(raupcrick(tmp>0, null = "r0"))
  tmp <- as.matrix(tmp)[1:39, 40:ncol(tmp)]
  res_rcnull[,r] <- diag(tmp)
}

boxplot(res_rcnull, ylab="Raup-Crick Null Dissimilarity", xlab="Radius (m)", axes=F)
axis(1, labels=names(stem.otu$abund), at=1:20)
axis(2)


########################################################
### LOOK AT DETECTION PROBABILITY VS BASAL AREA IN A NEIGHBORHOOD
########################################################

## We need to think about how to 

# Choose a certain spatial scale for convenience
focdat <- stem.otu$ba[[20]]

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
focdat <- stem.otu$abund[[5]][1:39,]

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
for(i in seq_along(stem.otu$abund)){
  focdat <- stem.otu$abund[[i]][201:239,]
  
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
for(i in seq_along(stem.otu$ba)){
  focdat <- stem.otu$ba[[i]][201:239,]
  
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
plot(rowSums(stem.otu$abund[[20]][1:39,]>0),
     rowSums(dnamat.pa))

# Quick plot for a given spatial scale after rarefaction
plot(rarefy(stem.otu$abund[[20]][1:39,], min(rowSums(stem.otu$abund[[20]][1:39,]))),
     rarefy(dnamat, min(rowSums(dnamat))))

# Look at correlation coefficients across spatial scales
corrs <- list()
corrs_rare <- list()

for(dist in seq_along(stem.otu$abund)){
  corrs[[dist]] <- cor.test(rowSums(dnamat>0), rowSums(stem.otu$abund[[dist]][1:39,]>0))
  
  corrs_rare[[dist]] <- cor.test(rarefy(stem.otu$abund[[dist]][1:39,], 
                                        min(rowSums(stem.otu$abund[[dist]][1:39,]))),
                                 rarefy(dnamat, min(rowSums(dnamat))))
}

par(mfrow=c(1,2))

plot(seq_along(stem.otu$abund), sapply(corrs, function(x) x$estimate), ylim=c(-1,1), axes=F,
     pch=16, cex=1, xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness [raw]")
segments(seq_along(stem.otu$abund), y0=sapply(corrs, function(x) x$conf.int[1]),
         y1=sapply(corrs, function(x) x$conf.int[2]))
abline(h=0, lty=2)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)
box()


plot(seq_along(stem.otu$abund), sapply(corrs_rare, function(x) x$estimate), ylim=c(-1,1), axes=F,
     pch=16, xlab="Radius (m)",
     ylab="Pearson correlation",
     main="Stem vs. DNA richness [rarefied]")
segments(seq_along(stem.otu$abund), y0=sapply(corrs_rare, function(x) x$conf.int[1]),
         y1=sapply(corrs_rare, function(x) x$conf.int[2]))
abline(h=0, lty=2)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)
box()





########################################################
### PRESENCE IN DNA AS A FUNCTION OF DISTANCE TO NEAREST (KNOWN) STEM
########################################################

x <- unlist(log10(stem.otu$nn[1:39,]))
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
for(sp in 1:39){
  xx <- log(unlist(stem.otu$nn[1:39,sp]))
  yy <- unlist(as.data.frame(1*(dnamat[,sp]>1)))
  mod <- glm(yy ~ xx, family=binomial)
  coeffs[sp] <- coef(mod)[2]
  nd <- data.frame(xx=seq(0, max(x), length.out=100))
  ypred <- predict(mod, nd, type='response')
  lines(nd$xx, ypred, col=ifelse(coeffs[sp]>0, rgb(1,0,0,0.4), rgb(0,0,1,0.4)))
}



