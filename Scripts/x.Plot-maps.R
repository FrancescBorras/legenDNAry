### MAP PLOTS

par(mfrow=c(2,2))
hist(stem.otu$nn[1:39,][dnamat.pa==1], xlim=c(0,1000), breaks=50)
hist(stem.otu$nn[1:39,][dnamat.pa==0], xlim=c(0,1000), breaks=50)



plot(unlist(stem.otu$nn[1:39,]), jitter(dnamat.pa))


m0 <- glm(as.vector(dnamat.pa) ~ unlist(stem.otu$nn[1:39,]))
summary(m0)


library(sf)

pts <- st_as_sf(census[census$Status=="alive", c("PX", "PY","DBH","Mnemonic")],
         coords = c(1,2))

b <- st_as_sfc(st_bbox(pts))

plot(pts$geometry,
     cex=sqrt(pts$DBH)/30, 
     lwd=0.4,
     col=rgb(0,0,0,0.5))
plot(b, add=T)
axis(1, at=seq(0,300,100))
axis(2, at=seq(0,500,100))




plot(pts$geometry,
     cex=sqrt(pts$DBH)/5, 
     lwd=0.5,
     xlim=c(175,225), ylim=c(175,225),
     pch=21, bg=cols[as.numeric(as.factor(pts$Mnemonic))])
abline(v=seq(0,500,10), h=seq(0,500,10), col='grey')
box(); axis(1); axis(2)
points(200,200, pch=4, col='red', cex=3, lwd=8)

s <- st_as_sf(data.frame(x=200, y=200),
                coords = c(1,2))

for(r in seq_along(seq(5,100,5))){
  plot(sf::st_buffer(s, seq(5,100,5)[r]), add=T, 
       border=rev(viridis(20))[r], lwd=5)
}



plot(pts$geometry,
     cex=sqrt(pts$DBH)/10, 
     lwd=0.5,
     xlim=c(100,300), ylim=c(100,300),
     col='grey'
     # ,pch=21, bg=cols[as.numeric(as.factor(pts$Mnemonic))]
     )
# abline(v=seq(0,500,10), h=seq(0,500,10), col='grey')
box(); axis(1); axis(2)

s <- st_as_sf(data.frame(x=200, y=200),
              coords = c(1,2))

for(r in seq_along(seq(5,100,5))){
  plot(sf::st_buffer(s, seq(5,100,5)[r]), add=T, 
       border=rev(viridis(20))[r], lwd=2)
}

points(200, 200, pch=4, col='red', cex=3, lwd=8)

