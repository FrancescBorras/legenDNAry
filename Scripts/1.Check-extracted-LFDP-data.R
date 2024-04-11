### Work with extracted stem data

data <- readRDS("Data/LFDP2023-extract-20240327.RDA")

str(data)


r <- do.call(cbind, lapply(data$abund, function(x){
  rowSums(x>0)
}))



radii <- as.numeric(names(data$abund))

plot(range(radii), range(r), pch=NA, 
     ylab="Species richness", xlab="Radius (m)")
for(i in 1:nrow(r)){
  lines(radii, r[i,], col=rgb(0,0,1,0.05), lwd=2)
}






