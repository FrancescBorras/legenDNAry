### Bar plot of the publication counts through time

x <- read.csv("Raw_data/PubMed_Timeline_Results_by_Year.csv", skip=1)

x <- x[order(x$Year),]

x <- x[x$Year %in% 2000:2023,]

sum(x$Count)

pdf("Figures/Publication_count.pdf", height=5)

b <- barplot(x$Count, col=rev(viridis::mako(40)), ylim=c(0,1400), axes=F)

# axis(1, labels=2000:2023, at=b, las=2)
axis(1, labels=seq(2000, 2020, 5), at=b[seq(1, 23, 5)])
axis(1, labels=NA, at=b)

axis(2, las=2)

box()

dev.off()