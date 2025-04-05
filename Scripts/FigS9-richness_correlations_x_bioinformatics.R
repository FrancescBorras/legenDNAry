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
           "rare_lenient_10k",
           "lenientf1_10k",
           "repfiltered_10k",
           "rare_repfiltered_10k",
           "rare_repfilteredf1_10k",
           "lenient_40k",
           "rare_lenient_40k",
           "lenient_70k",
           "rare_lenient_70k",
           "lenient_200k",
           "rare_lenient_200k")

### Name data files
outfiles <- c( "stem-soil-40pt-data-lenient_10k-20250312.RDA",
               "stem-soil-40pt-data-rare_lenient_10k-20250312.RDA",
               "stem-soil-40pt-data-lenientf1-20250312.RDA",
               "stem-soil-40pt-data-repfiltered-20250312.RDA",
               "stem-soil-40pt-data-rare_repfiltered-20250312.RDA",
               "stem-soil-40pt-data-repfilteredf1-20250312.RDA",
               "stem-soil-40pt-data-lenient_40k-20250312.RDA",
               "stem-soil-40pt-data-rare_lenient_40k-20250312.RDA",
               "stem-soil-40pt-data-lenient_70k-20250312.RDA",
               "stem-soil-40pt-data-rare_lenient_70k-20250312.RDA",
               "stem-soil-40pt-data-lenient_200k-20250312.RDA",
               "stem-soil-40pt-data-rare_lenient_200k-20250312.RDA")

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
  
  ### Count number of valid eDNA sample points (out of 40)
  (npts <- nrow(stem.23$abund[[1]]) - length(grep('random', rownames(stem.23$abund[[1]]))))
  
  ### Species code full plot summaries
  lfdp <- readRDS("Raw_data/LFDP2023-extract-v2-20240427.RDA")
  
  ### gOTU full plot summaries
  lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
  
  ### Make presence absence matrices
  dnamat.pa <- 1*(dnamat>0)
  stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
  
########################################################
### 1. Correlation of species richness between DNA and stem data
########################################################

# Correlation of DNA and stem species richness across spatial scales
corrs <- corrs_rare <- list()

for(r in seq_along(stem.23$abund)){
  corrs_rare[[r]] <- cor.test(rarefy(stem.23$abund[[r]][1:npts,],
                                     min(rowSums(stem.23$abund[[r]][1:npts,]))),
                              rarefy(dnamat, min(rowSums(dnamat))))

  corrs[[r]] <- cor.test(rowSums(stem.23$abund[[r]][1:npts,]>0), rowSums(dnamat.pa))
}

cp <- rev(viridis::viridis(20))

if(data_selector==1){
  pdf(paste0("Figures/compare_filtering/FigSx.Alpha_div_correlations_x_bioinformatics.pdf"),
      width=10, height=8)
  par(mfrow=c(3,4), mar=c(4,4,2,1))
} 

plot(seq_along(stem.23$abund), sapply(corrs_rare, function(x) x$estimate), 
     ylim=c(-1,1), axes=F,
     pch=21, bg='grey', xlab="Radius (m)",
     ylab="Pearson correlation",
     )
segments(seq_along(stem.23$abund), 
         y0=sapply(corrs_rare, function(x) x$conf.int[1]),
         y1=sapply(corrs_rare, function(x) x$conf.int[2]))
abline(h=0, lty=2)
points(seq_along(stem.23$abund), 
       sapply(corrs_rare, function(x) x$estimate), 
       pch=21, bg=cp, cex=1.5)
axis(1, labels=seq(5, 100, 5), at=1:20)
axis(2)
box()

mtext(paste0(LETTERS[data_selector], ") ", label[data_selector]), adj=0, line=0.1)

}

dev.off()

