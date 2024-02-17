# Summarize_stem_data.R


# READ THE CENSUS DATA (2016)
library(EDIutils)
raw6 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "325c43057e0dd4e1cd6a13fa5125a76d")
census <- readr::read_csv(file = raw6)


# GENERATE FAKE SAMPLE POINT COORDINATES FOR TESTING
sample_xy <- data.frame(matrix(c(50, 50, 100, 100, 200, 200), nrow=3, ncol=2))
colnames(sample_xy) <- c("PX", "PY")


### FILTER STEM DATA TO TREES IN A CERTAIN RADIUS AROUND SAMPLE PLOT
tree_dist <- function(sample_xy, census, radius=10) {
  d <- sqrt((sample_xy$PX - census$PX)^2 + (sample_xy$PY - census$PY)^2)
  foctrees <- data.frame(census[d <= radius & census$Status=="alive",])
  return(foctrees)
}


### SUMMARIZE THE NEIGHBORHOOD TREES IN SOME SAY (E.G., RICHNESS)
neighborhoods <- list()
for(i in 1:nrow(sample_xy)){
  neighborhoods[[i]] <- tree_dist(sample_xy[i,], census, radius=5)
}
neighborhoods

richness <- unlist(lapply(neighborhoods, function(x) length(unique(x$Mnemonic))))





