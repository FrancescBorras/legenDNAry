library(viridis)
library(sf)
library(EDIutils)


### Load data
dnamat <- readRDS("Processed_data/stem-soil-39pt-data-20240510.RDA")[[1]]
stem.23 <- readRDS("Processed_data/stem-soil-39pt-data-20240510.RDA")[[2]]
stem.16 <- readRDS("Processed_data/stem-soil-39pt-data-20240510.RDA")[[4]]
traits <- readRDS("Processed_data/stem-soil-39pt-data-20240510.RDA")[[3]]
sampxy <- read.csv("Raw_data/LFDP-sample39-coordinates.csv", row.names = 1)

### Make presence absence matrices
dnamat.pa <- 1*(dnamat>0)
stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
stem.pa.list16 <- lapply(stem.16$abund, function(x) 1*(x>0))

### Color palette for plotting
cp <- rev(viridis(20))

### Load 2016 census data
raw6 <- read_data_entity(packageId="knb-lter-luq.119.1545979", 
                         entityId="325c43057e0dd4e1cd6a13fa5125a76d")
census <- readr::read_csv(file = raw6)
census <- census[census$Status %in% "alive",]


### MAP PLOTS

pdf("Figures/FigX.Maps-without-samp-points.pdf")

layout(matrix(data=c(1,1,2,3), nrow=2, ncol=2))

par(mar=rep(1, 4))

pts <- st_as_sf(census[census$Status=="alive", c("PX", "PY","DBH", "Mnemonic")],
         coords = c(1,2))

b <- st_as_sfc(st_bbox(pts))

set.seed(191)
cols <- sample(colors(distinct=T), length(unique(pts$Mnemonic)))

plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=scales::rescale(log(pts$DBH), to=c(0.25,0.75)),
     col=cols[as.factor(pts$Mnemonic)])
plot(b, add=T)
axis(1, at=seq(0,300,100))
axis(2, at=seq(0,500,100))
# points(sampxy[,1:2], pch=5, lwd=5)


plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=0.5,
     xlim=c(50,250), ylim=c(50,250),
     col=cols[as.factor(pts$Mnemonic)])
box(); axis(1); axis(2)


plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=0.5,
     xlim=c(50,250), ylim=c(50,250),
     col='grey')
box(); axis(1); axis(2)

s <- st_as_sf(data.frame(x=150, y=150),
              coords = c(1,2))

for(r in seq_along(seq(5,100,5))){
  plot(sf::st_buffer(s, seq(5,100,5)[r]), add=T, 
       border=rev(viridis(20))[r], lwd=2)
}

points(150, 150, pch=4, col='red', cex=3, lwd=8)

dev.off()




###
# 
# pts <- st_as_sf(census[census$Status=="alive", c("PX", "PY","DBH", "Mnemonic")],
#                 coords = c(1,2))
# b <- st_as_sfc(st_bbox(pts))
# 
# set.seed(191)
# cols <- sample(colors(distinct=T), length(unique(pts$Mnemonic)))
# 
# plot(pts$geometry,
#      cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
#      lwd=scales::rescale(log(pts$DBH), to=c(0.25,0.75)),
#      col=cols[as.factor(pts$Mnemonic)])
# plot(b, add=T)
# axis(1, at=seq(0,300,100))
# axis(2, at=seq(0,500,100))
# points(sampxy[,1:2], pch=3, cex=1, lwd=3)
# 
# 
# s <- st_as_sf(data.frame(x=150, y=150),
#               coords = c(1,2))
# 
# for(r in seq_along(seq(5,100,5))){
#   plot(sf::st_buffer(s, seq(5,100,5)[r]), add=T, 
#        border=rev(viridis(20))[r], lwd=2)
# }
# 
# points(150, 150, pch=4, col='red', cex=3, lwd=8)








### Schematic diagram

pdf("Figures/FigX.map1-for-schematic.pdf")

pts <- st_as_sf(census[census$Status=="alive", c("PX", "PY","DBH", "Mnemonic")],
                coords = c(1,2))

b <- st_as_sfc(st_bbox(pts))

191
set.seed(192)
cols <- sample(colors(distinct=T), length(unique(pts$Mnemonic)))

plot(pts$geometry,
     cex=scales::rescale((pts$DBH), to=c(0.01, 2)),
     lwd=scales::rescale(log(pts$DBH), to=c(0.25,1)),
     col=scales::alpha(cols[as.factor(pts$Mnemonic)], 0.5))

plot(b, add=T)
axis(1, at=seq(0,300,100))
axis(2, at=seq(0,500,100))
points(sampxy[,1:2], pch=4, lwd=2, cex=0.75, col=1)
# points(sampxy[,1:2], pch=3, col='white', cex=0.5)

dev.off()


plot(pts$geometry,
     xlim=c(150,170), ylim=c(90,110),
     cex=scales::rescale((pts$DBH), to=c(0.5, 18)),
     lwd=scales::rescale(log(pts$DBH), to=c(1,5)),
     col=scales::alpha(cols[as.factor(pts$Mnemonic)], 0.75))
box()
axis(1, at=seq(0,300,10))
axis(2, at=seq(0,500,10))
points(sampxy[,1:2], pch=4, lwd=2, cex=0.75, col=1)



plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=0.5,
     xlim=c(50,250), ylim=c(50,250),
     col=cols[as.factor(pts$Mnemonic)])
box(); axis(1); axis(2)


plot(pts$geometry,
     cex=scales::rescale(log(pts$DBH), to=c(0.01,0.8)),
     lwd=0.5,
     xlim=c(50,250), ylim=c(50,250),
     col='grey')
box(); axis(1); axis(2)

s <- st_as_sf(data.frame(x=150, y=150),
              coords = c(1,2))

for(r in seq_along(seq(5,100,5))){
  plot(sf::st_buffer(s, seq(5,100,5)[r]), add=T, 
       border=rev(viridis(20))[r], lwd=2)
}
