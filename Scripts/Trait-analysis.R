### Trait analysis

### Load libraries
library(vegan)
library(caret)
library(compositions)
library(viridis)
library(ecodist)
library(sf)
library(spdep)

### Make labels for saving plots
label <- c("lenient_10k")

### Name data files
outfiles <- c( "stem-soil-40pt-data-lenient_10k-20250312.RDA")

### Select data file to load (for testing)
data_selector <- 1

### Load selected data file
datafile <- paste0("Processed_data/", outfiles[data_selector])
  
### Load data
dnamat <- readRDS(datafile)[[1]]
stem.23 <- readRDS(datafile)[[2]]
stem.16 <- readRDS(datafile)[[4]]
traits <- readRDS(datafile)[[3]]

### gOTU full plot summaries
lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
lfdp16 <- readRDS("Raw_data/LFDP2016-extract-v2-20240427-gOTUs.RDA")

### Make presence absence matrices
dnamat.pa <- 1*(dnamat>0)
stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
stem.pa.list16 <- lapply(stem.16$abund, function(x) 1*(x>0))

### Color palette for plotting
cp <- rev(viridis::viridis(20))


# Normalize to relative abundance (each row sums to 1)
dnamat_ra <- dnamat / rowSums(dnamat)



# DNA presence
y <- as.vector(dnamat.pa)

# Abundance of gOTUs in e.g. 10 m
x <- unlist(stem.23$ba$`10`[1:nrow(dnamat),])
x <- x + 0.0001

# Traits of gOTUs
t <- 1/as.numeric(scale(rep(traits$SLA.wp, times=nrow(dnamat.pa))))
  
m1 <- glm(y ~ x + t, family = "binomial")

summary(m1)

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(x+0.00001, jitter(y, 0.25), pch=16, col=rgb(0,0,0,0.1), 
     xlab="Stem abundance in 10 m",
     ylab="DNA detected", log='x')

# ===
nd <- data.frame(x=seq(0.0001, max(x), length.out=1000), t=mean(t, na.rm=T))
ypred <- predict(m1, nd, type='response')

# Overall fit
lines(nd$x, ypred, col="blue", lwd=3)




plot(t, jitter(y, 0.25), pch=16, col=rgb(0,0,0,0.1), 
     xlab="LMA (g / m2)",
     ylab="DNA detected", log='')

nd <- data.frame(x=mean(x, na.rm=T), t=seq(min(t, na.rm=T), 
                                           max(t, na.rm=T), 
                                           length.out=1000))
ypred <- predict(m1, nd, type='response')
# Overall fit
lines(nd$t, ypred, col="blue", lwd=3)


# === HOW DO LEAF TRAITS AFFECT FALSE PRESENCES?

x <- unlist(stem.23$ba$`10`[1:nrow(dnamat),])

fp <- 1 * (y==1 & x==0)

m2 <- glm(fp ~ t, family = "binomial")

summary(m2)

plot(t, jitter(fp, 0.25), pch=16, col=rgb(0,0,0,0.1), 
     ylab="False presences",
     xlab="LMA (g / m2)") 

# ===
nd <- data.frame(t=seq(min(t, na.rm=T), 
                         max(t, na.rm=T), 
                         length.out=1000))
ypred <- predict(m2, nd, type='response')

# Overall fit
lines(nd$t, ypred, col="blue", lwd=3)


# === HOW DO LEAF TRAITS AFFECT FALSE ABSENCES?

x <- unlist(stem.23$ba$`10`[1:nrow(dnamat),])

fa <- 1 * (y==0 & x > 0)

m3 <- glm(fa ~ t, family = "binomial")

summary(m3)

plot(t, jitter(fa, 0.25), pch=16, col=rgb(0,0,0,0.1), 
     ylab="False absences",
     xlab="LMA (g / m2)") 

# ===
nd <- data.frame(t=seq(min(t, na.rm=T), 
                       max(t, na.rm=T), 
                       length.out=1000))
ypred <- predict(m3, nd, type='response')

# Overall fit
lines(nd$t, ypred, col="blue", lwd=3)












dnamat
traits$THK


