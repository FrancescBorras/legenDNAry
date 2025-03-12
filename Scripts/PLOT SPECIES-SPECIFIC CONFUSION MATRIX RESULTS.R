### PLOT SPECIES-SPECIFIC CONFUSION MATRIX RESULTS

# Make colors to be relative to log basal area (x-axis is log abundance)
cp <- (viridis::viridis_pal(option = "A")(78))

# Get a simple vector of total abundance in same order as 

lfdp23 <- lfdp23[match(colnames(samp.dna.pa), rownames(lfdp23)),]

cpsp <- cp[round(rank(lfdp23$total_ba))]


par(mfcol=c(3,2), mar=c(4,4,2,1))

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$Sensitivity))
plot(lfdp23$total_ba, x[,1], log='x', 
     pch=21, bg=cpsp, cex=2,
     xlab="Total basal area (m2)", ylab="Sensitivity")
mtext("A", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$Specificity))
plot(lfdp23$total_ba, x[,1], log='x', 
     pch=21, bg=cpsp, cex=2,
     xlab="Total basal area (m2)", ylab="Specificity")
mtext("B", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_obs_list, function(x) x$MCC))
plot(lfdp23$total_ba, x[,1], log='x', 
     pch=21, bg=cpsp, cex=2, 
     xlab="Total basal area (m2)", ylab="MCC")
mtext("C", adj=0, line=0.5)


x <- do.call(cbind, lapply(sp_conf_stats_ses_list, function(x) x$Sensitivity))

plot(lfdp23$total_ba, x[,1], log='x', pch=NA, bg=cp, 
     xlab="Total basal area (m2)", ylab="SES Sensitivity")
polygon(x=c(1e-6,5e5,5e5,1e-6), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(lfdp23$total_ba, x[,1], 
       pch=21, bg=cpsp, cex=2)
mtext("D", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_ses_list, function(x) x$Specificity))

plot(lfdp23$total_ba, x[,1], log='x', pch=NA, bg=cp, 
     xlab="Total basal area (m2)", ylab="SES Specificity")
polygon(x=c(1e-6,5e5,5e5,1e-6), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(lfdp23$total_ba, x[,1],
       pch=21, bg=cpsp, cex=2)
mtext("E", adj=0, line=0.5)

x <- do.call(cbind, lapply(sp_conf_stats_ses_list, function(x) x$MCC))

plot(lfdp23$total_ba, x[,1], log='x', pch=NA, bg=cp, 
     xlab="Total basal area (m2)", ylab="SES MCC")
polygon(x=c(1e-6,5e5,5e5,1e-6), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
points(lfdp23$total_ba, x[,1],
       pch=21, bg=cpsp, cex=2)
mtext("F", adj=0, line=0.5)









