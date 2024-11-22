
pdf("Figures/FigX.Maps-with-samp-points-20240923.pdf")

par(mar=rep(1, 4))

pts <- st_as_sf(census[census$Status=="alive", c("PX", "PY","DBH", "Mnemonic")],
                coords = c(1,2))

b <- st_as_sfc(st_bbox(pts))

set.seed(191)
cols <- sample(colors(distinct=T), length(unique(pts$Mnemonic)))
cols <- grey.colors(length(unique(pts$Mnemonic)))

plot(pts$geometry,
     cex=scales::rescale((pts$DBH), to=c(0.01,3)),
     lwd=scales::rescale((pts$DBH), to=c(0.5,1.5)),
     col=cols[as.factor(pts$Mnemonic)])
plot(b, add=T)
axis(1, at=seq(0,300,100))
axis(2, at=seq(0,500,100))

points(rbind(sampxy[,1:2], data.frame(X=220,Y=220)), pch=3, lwd=2, col='blue', cex=1)

dev.off()