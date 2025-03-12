##########################################
##########################################
##########################################
### Plot confusion matrix boxplot figure with points added for small grid
##########################################
##########################################
##########################################

### Color palette for plotting
cp <- rev(viridis::viridis(20))

g_conf_stats_obs <- readRDS("Processed_data/Conf_matrix_output-Gplots-20250310.RDA")[[1]]
g_conf_stats_ses <- readRDS("Processed_data/Conf_matrix_output-Gplots-20250310.RDA")[[2]]

conf_stats_obs_list <- readRDS("Processed_data/Conf_matrix_output-lenient_10k-20250307.RDA")[[1]]
conf_stats_ses_list <- readRDS("Processed_data/Conf_matrix_output-lenient_10k-20250307.RDA")[[2]]


pdf(paste0("Figures/compare_filtering/Fig4.Confusion-matrix-stats_",
           label[data_selector],
           ".pdf"), width = 8, height = 10)


par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

## Sensitivity (true positive rate)
boxplot(lapply(conf_stats_obs_list, function(x) x$Sensitivity), 
        ylab="Sensitivity (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
points(1:20, g_conf_stats_obs$Sensitivity, pch=23, bg='red', lwd=0.5)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
graphics::box()
mtext("A", adj=0, line=0.5)

legend("topright", legend=c("Big grid samples", "Small grid samples"), 
       pch=c(22,23), pt.cex=c(2.5,1), 
       pt.bg=c(cols[1], 'red'), bty='n', lty=c(1,NA))

boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
        ylab="Sensitivity (SES)", 
        xlab="Radius (m)", axes=F, col=cp)
axis(1, labels=names(stem.23$abund), at=1:20)
axis(2)
# abline(h=c(-1.96, 1.96), lwd=2, lty=2)
polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
        axes=F, col=cp, add=T)
points(1:20, g_conf_stats_ses$Sensitivity, pch=23, bg='red', lwd=0.5)
graphics::box()
mtext("D", adj=0, line=0.5)

## Specificity (true negative rate)
boxplot(lapply(conf_stats_obs_list, function(x) x$Specificity), 
        ylab="Specificity (observed)", 
        xlab="Radius (m)", axes=F, col=cp, ylim=c(0.7,1))
points(1:20, g_conf_stats_obs$Specificity, pch=23, bg='red', lwd=0.5)
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
points(1:20, g_conf_stats_ses$Specificity, pch=23, bg='red', lwd=0.5)
graphics::box()
mtext("E", adj=0, line=0.5)

## MCC
boxplot(lapply(conf_stats_obs_list, function(x) x$`MCC`), 
        ylab="MCC (observed)", 
        xlab="Radius (m)", axes=F, col=cp)
points(1:20, g_conf_stats_obs$MCC, pch=23, bg='red', lwd=0.5)
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
points(1:20, g_conf_stats_ses$MCC, pch=23, bg='red', lwd=0.5)
graphics::box()
mtext("F", adj=0, line=0.5)


dev.off()
