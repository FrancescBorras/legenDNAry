sp.conf.obs.list <- list()

### Radius loop
for(r in 1:20){
  message(paste("radius =", r))
  sp.conf.obs <- matrix(nrow=ncol(dnamat.pa), ncol=11)
  # sp.conf.ses <- matrix(nrow=ncol(dnamat.pa), ncol=11)
  
  ### Taxon loop
  for(g in 1:ncol(dnamat.pa)){

    pred <- abs(dnamat.pa[,g]-1)
    ref <- abs(1*(stem.otu$abund[[r]][1:39,g] > 0)-1)
    u <- c(0,1)
    
    t <- table(factor(pred, u), factor(ref, u))
    
    confus.obs <- confusionMatrix(t)
    sp.conf.obs[g,] <- confus.obs$byClass
       
  }
  sp.conf.obs <- as.data.frame(sp.conf.obs)
  names(sp.conf.obs) <- names(confus.obs$byClass)
  sp.conf.obs.list[[r]] <- sp.conf.obs
}
  
corest <- vector()
for(i in 1:20){
  sensitiv <- 1-sp.conf.obs.list[[i]]$Recall
  logabund <- log10(colSums(stem.otu$ba[[i]]))
  logabund[is.infinite(logabund)] <- NA
  plot(logabund, sensitiv)
  abline(lm(sensitiv~logabund))
  tmp <- cor.test(logabund, sensitiv, use='p')
  corest[i] <- tmp$estimate
}
plot(corest, xlab="Spatial scale (m)", ylab="Correlation (Abundance x Recall)",
     pch=16, type='l', axes=F)
axis(1, labels=names(stem.otu$abund), at=1:20)
axis(2)



rowSums(stem.otu$ba[[1]]>0)

b <- boxplot(sapply(stem.otu$ba, function(x) rowSums(x[1:39,]>0)), col=cp,
        xlab="Radius (m)", ylab="Species richness", xlim=c(-0.5,20), ylim=c(0,80))


boxplot(rowSums(dnamat.pa), add=T, at=-0.5, width=2, col=2)
abline(v=0.25)
abline(h=79, lty=2)
axis(1, at=-0.5, labels="DNA")




plot(specaccum(dnamat.pa))


# Precision is the proportion of times a taxon is detected in DNA among sites where it was present.
# More abundant species have lower Precision.
# 

At which spatial scale do we capture species in the dna data best when we know they are present in the stem data?
  i.e., what is the true presence proprtion?


That means the proportion of correctly predicted positive cases among all actual positive cases is lower for abundant species.  That makes sense; if they are common everywhere, then the have a greater chance to be missed.


### AVERAGE E.G. RECALL ACROSS TAXA AT DIFFERENT SPATIAL SCALES
boxplot(do.call(cbind, lapply(sp.conf.obs.list, function(x) x$F1)),
        col=cp)








