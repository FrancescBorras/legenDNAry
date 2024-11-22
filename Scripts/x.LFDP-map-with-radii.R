library(viridis)
library(sf)
library(EDIutils)

sampxy <- read.csv("Raw_data/LFDP-sample39-coordinates.csv", row.names = 1)
sampxy <- rbind(sampxy, c(220,220))

### Load 2016 census data
raw6 <- read_data_entity(packageId="knb-lter-luq.119.1545979", 
                         entityId="325c43057e0dd4e1cd6a13fa5125a76d")
census <- readr::read_csv(file = raw6)
census <- census[census$Status %in% "alive",]


pts <- st_as_sf(census[census$Status=="alive", c("PX", "PY","DBH", "Mnemonic")],
                coords = c(1,2))

pdf("Figures/FigX.Map-40m-radius.pdf")


plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=scales::rescale(log(pts$DBH), to=c(0.25,0.75)),
     col=cols[as.factor(pts$Mnemonic)])
plot(st_as_sfc(st_bbox(pts)), add=T)
axis(1, at=seq(0,300,100))
axis(2, at=seq(0,500,100))

points(sampxy[,1:2], pch=3, lwd=3)

plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=0.5,
     xlim=c(100,200), ylim=c(100,200),
     col='grey')
box(); axis(1); axis(2)

s <- st_as_sf(data.frame(x=150, y=150),
              coords = c(1,2))
for(r in seq_along(seq(5,40,5))){
  plot(sf::st_buffer(s, seq(5,40,5)[r]), add=T,
       border=rev(viridis(8))[r], lwd=2)
}

points(150, 150, pch=3, lwd=3)

dev.off()