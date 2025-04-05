
"Raw_data/LFDP2016-extract-v2-20240427.RDA"

### PROCESS ROUND 2 OF LFDP 2023 DATA EXTRACT

res <- matrix(nrow=length(codes.collapse.list), ncol=6) # from 2.Matching...

for(i in 1:length(codes.collapse.list)){ 
  res[i,1:6] <- colSums(lfdp[lfdp$spcode %in% codes.collapse.list[[i]],2:7], na.rm=T)
}

rownames(res) <- names(codes.collapse.list)
res <- as.data.frame(res)
names(res) <- names(lfdp)[-1]

# saveRDS(res, "Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
# saveRDS(res, "Raw_data/LFDP2023-extract-v2-20250404-gOTUs.RDA")
# saveRDS(res, "Raw_data/LFDP2016-extract-v2-20250404-gOTUs.RDA")
# saveRDS(res, "Raw_data/LFDP2023-extract-v2-drop_noSeq_gOTUs-20240427-gOTUs.RDA")

# 'df' is a data.frame build by counting up 2016 stem / basal area data

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(df$total_abund[match(lfdp$spcode, df$spcode)], lfdp$total_abund, 
     log='xy', pch=21, bg="grey", 
     xlab="2016 Abundance", ylab="2023 Abundance")
abline(0, 1, col='blue')

plot(df$total_ba[match(lfdp$spcode, df$spcode)], lfdp$total_ba, 
     log='xy', pch=21, bg="grey", 
     xlab="2016 Basal area", ylab="2023 Basal area")
abline(0, 1, col='blue')

hist(log10(census$DBH[census$Status=='alive']), main=NA,
     axes=F, xlab="DBH mm (log10)")
axis(1, at=1:3, labels=c(10,100,1000))
axis(2)

plot(df$total_ba[match(lfdp$spcode, df$spcode)], 100*lfdp$total_ba, 
     log='xy', pch=21, bg="lightpink", 
     xlab="2016 Basal area", ylab="100 * 2023 Basal area")
abline(0, 1, col='blue')




