########################################################
########################################################
##################
### Explore false presences and abundance
##################
########################################################
########################################################

### Load data
dnamat <- readRDS("Processed_data/stem-soil-40pt-data-lenient_10k-20241205.RDA")[[1]]
stem.23 <- readRDS("Processed_data/stem-soil-40pt-data-lenient_10k-20241205.RDA")[[2]]

### Count number of valid eDNA sample points (out of 40)
(npts <- nrow(stem.23$abund[[1]]) - length(grep('random', rownames(stem.23$abund[[1]]))))

### gOTU full plot summaries
lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")

### Make presence absence matrices
dnamat.pa <- 1*(dnamat>0)
stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))

### Color palette for plotting
cp <- rev(viridis::viridis(20))

### Make a matrix of false positives at 5 m
fp <- 1 * (dnamat.pa==1 & stem.pa.list$`5`[1:npts,]==0)

### Make a matrix of false absences at 5 m
fa <- 1 * (dnamat.pa==0 & stem.pa.list$`5`[1:npts,]==1)

### Make a matrix of true presences at 5 m
tp <- 1 * (dnamat.pa==1 & stem.pa.list$`5`[1:npts,]==1)

### Make a matrix of false absences at 5 m
ta <- 1 * (dnamat.pa==0 & stem.pa.list$`5`[1:npts,]==0)

### Get proportion of sites where a gOTU was a false positive at 5 m
fp_prop <- colSums(fp)/npts
fa_prop <- colSums(fa)/npts
tp_prop <- colSums(tp)/npts
ta_prop <- colSums(ta)/npts


### Plot!

par(mfcol=c(2,3), mar=c(4,4,1,1))

cp <- viridis_pal(option = "A")(20)[cut(log10(lfdp23$total_ba+0.0001), 20)]

plot(lfdp23$total_abund, tp_prop, log='x',
     xlab="Total Abundance", ylab="Prop sites where TP @5m",
     pch=21, bg=cp, lwd=0.5, ylim=c(0,1))
# If you want to add text labels to know which OTU is which:
 # text(lfdp23$total_abund, tp_prop, rownames(lfdp23), cex=0.75)
mtext("A", adj=0, line=0.5)

plot(lfdp23$total_abund, fa_prop, log='x',
     xlab="Total Abundance", ylab="Prop sites where FA @5m",
     pch=21, bg=cp, lwd=0.5, ylim=c(0,1))
mtext("B", adj=0, line=0.5)

plot(lfdp23$total_abund, fp_prop, log='x',
     xlab="Total Abundance", ylab="Prop sites where FP @ 5m",
     pch=21, bg=cp, lwd=0.5, ylim=c(0,1))
mtext("C", adj=0, line=0.5)

plot(lfdp23$total_abund, ta_prop, log='x',
     xlab="Total Abundance", ylab="Prop sites where TA @5m",
     pch=21, bg=cp, lwd=0.5, ylim=c(0,1))
mtext("D", adj=0, line=0.5)

### LOOK AT TRUE PRESENCE IN 5M AND HOW THAT RELATES TO TP / FA ETC.
pres <- colSums(stem.pa.list$`5`[1:npts,])

plot(pres/npts, colSums(tp)/npts, cex=sqrt(lfdp23$total_abund)/30,
     pch=21, bg=cp, lwd=0.5, xlim=c(0,1), ylim=c(0,1),
     ylab="Prop sites present @5m", xlab="Prop sites where TP @5m")
abline(0,1, lty=2)
mtext("E", adj=0, line=0.5)

plot(pres/npts, colSums(fa)/npts, cex=sqrt(lfdp23$total_abund)/30,
     pch=21, bg=cp, lwd=0.5, xlim=c(0,1), ylim=c(0,1),
     ylab="Prop sites present @5m", xlab="Prop sites where FA @5m")
abline(0,1, lty=2)
mtext("F", adj=0, line=0.5)

mean(rowSums(stem.pa.list$`25`[1:npts,]))




### SAME BUT FOR TOTAL BASAL AREA
par(mfrow=c(3,2), mar=c(4,4,1,1))

plot(lfdp23$total_ba, tp_prop, log='x',
     xlab="Total Basal Area (m2)", ylab="Prop sites where FP @5m",
     pch=16)
# If you want to add text labels to know which OTU is which:
# text(lfdp23$total_abund, fp_prop, rownames(lfdp23), cex=0.75)

plot(lfdp23$total_ba, fa_prop, log='x',
     xlab="Total Basal Area (m2)", ylab="Prop sites where FA @5m",
     pch=16)

plot(lfdp23$total_ba, fa_prop, log='x',
     xlab="Total Basal Area (m2)", ylab="Prop sites where TP @ 5m",
     pch=16)

plot(lfdp23$total_ba, ta_prop, log='x',
     xlab="Total Basal Area (m2)", ylab="Prop sites where TA @5m",
     pch=16)

