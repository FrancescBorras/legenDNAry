########################################################
### COMPUTE CONFUSION MATRIX INDICES WITH THE "G" (INTENSIVE SUB-SAMPLE) PLOTS
########################################################

library(phyloseq)
library(googlesheets4)
library(EDIutils)

  ### Load the LFDP 2023 census extract data
  tree <- readRDS("Raw_data/LFDP2023-extract-20240510.RDA")
  codes <- readxl::read_xlsx("Raw_data/LFDP-SPcodes.xlsx")
  
  ### Load the sample point coordinates
  sample_xy <- read.csv("Raw_data/LFDP-eDNA-xy-V2.csv")
  
  ### Get G sample points phyloseq object
  # g <- readRDS("Processed_data/G_sample_phyloseq_20250305.RDA")
  
  # Get coordinates from sample points included
  xy <- data.frame(X=140, Y=260)
  
  # ### Make the rows of the stem data match the rows of the eDNA otu-table data
  # 
  # ### Create stem data for 2023
  # stem_abund <- lapply(tree$abund, function(x) {
  #   rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
  #         x[88:1087,])})
  # for(i in seq_along(stem_abund)){
  #   rownames(stem_abund[[i]]) <- c("g", paste0('random', 1:1000))
  # }
  # 
  # stem_ba <- lapply(tree$ba, function(x) {
  #   rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
  #         x[88:1087,])})
  # for(i in seq_along(stem_ba)){
  #   rownames(stem_ba[[i]]) <- c("g", paste0('random', 1:1000))
  # }
  # 
  # stem_nn <- rbind(tree$nearest_sp[match(paste(xy$X, xy$Y),
  #                                        paste(sample_xy$PX, sample_xy$PY)),], 
  #                  tree$nearest_sp[88:1087,])
  # rownames(stem_nn) <- c("g", paste0('random', 1:1000))
  # 
  # stem <- list(abund=stem_abund, ba=stem_ba, nearest_sp=stem_nn)
  # 
  # ##### THIS BLOCK COLLAPSES THE DNA DATA TO THE gOTU CLUSTERS
  # ### Identify all POTUs/OTUs we need to keep because they are assigned to a gOTU
  # potus.keep <- codes$POTU[!is.na(codes$gOTU)]
  # otus.keep <- otu.potu.link[[6]]$OTU[match(potus.keep, otu.potu.link[[6]]$abundance)]
  # gotus.keep <- codes$gOTU[!is.na(codes$gOTU)]
  # codes.keep <- codes$`SPECIES CODE`[!is.na(codes$gOTU)]
  # 
  # # POTUs that we keep and to be collapsed
  # codes.collapse.list <- split(codes.keep, gotus.keep)
  # otus.collapse.list <- split(otus.keep, gotus.keep)
  # 
  # # loop to collapse taxa in phyloseq objects
  # for(i in seq_along(otus.collapse.list)){
  #   if(length(unique(otus.collapse.list[[i]][!is.na(otus.collapse.list[[i]])])) > 1){
  #     taxcollapse <- unique(otus.collapse.list[[i]][!is.na(otus.collapse.list[[i]])])
  #     taxcollapse <- taxcollapse[taxcollapse %in% rownames(g@tax_table)]
  #     g <- merge_taxa(g, taxcollapse, 1)
  #   }
  # }
  # 
  # ##### THIS BLOCK COLLAPSES THE STEM DATA TO THE gOTU CLUSTERS
  # stem.otu <- list()
  # for(r in seq_along(stem$abund)){
  #   
  #   # Initialize a matrix to hold collapsed results for a given radius
  #   tmat_abund <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
  #   tmat_ba <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
  #   tmat_nn <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
  #   
  #   # loop to collapse taxa in phyloseq objects
  #   for(i in seq_along(codes.collapse.list)){
  #     tmp_abund <- stem$abund[[r]][,which(colnames(stem$abund[[1]]) %in% codes.collapse.list[[i]])]
  #     tmp_ba <- stem$ba[[r]][,which(colnames(stem$ba[[1]]) %in% codes.collapse.list[[i]])]
  #     tmp_nn <- stem$nearest_sp[,which(colnames(stem$nearest_sp) %in% codes.collapse.list[[i]])]
  #     
  #     if(is.matrix(tmp_abund)){
  #       tmat_abund[,i] <- rowSums(tmp_abund)
  #       tmat_ba[,i] <- rowSums(tmp_ba)
  #       tmat_nn[,i] <- rowSums(tmp_nn)
  #     } else {
  #       tmat_abund[,i] <- tmp_abund
  #       tmat_ba[,i] <- tmp_ba
  #       tmat_nn[,i] <- tmp_nn
  #     }
  #   }
  #   
  #   tmat_abund <- as.data.frame(tmat_abund)
  #   tmat_ba <- as.data.frame(tmat_ba)
  #   tmat_nn <- as.data.frame(tmat_nn)
  #   
  #   names(tmat_nn) <- names(tmat_ba) <- names(tmat_abund) <- names(codes.collapse.list)
  #   rownames(tmat_nn) <- rownames(tmat_ba) <- rownames(tmat_abund) <- rownames(stem$abund[[1]])
  #   
  #   stem.otu$abund[[r]] <- tmat_abund
  #   stem.otu$ba[[r]] <- tmat_ba
  #   stem.otu$nn <- tmat_nn
  #   
  #   names(stem.otu$ba)[r] <- names(stem.otu$abund)[r] <- names(stem$abund)[r]
  # }
  # 
  # 
  # ### To make the columns match, and add zero columns to OTU data
  # 
  # # The new phyloseq object should have these OTUs as colnames
  # collapsed.otus <- unlist(lapply(otus.collapse.list, function(x) x[!is.na(x)][1]))
  # 
  # # The stem data has these names exactly
  # if(!all(names(collapsed.otus) == names(stem.otu$abund[[1]]))){
  #   message("Warning: all(names(collapsed.otus) != names(stem.otu$abund[[1]]))")
  # }
  # 
  # # Some gOTU names are in the soil data but not the collapsed list because they aren't trees
  # qs <- colnames(g@otu_table)[!colnames(g@otu_table) %in% collapsed.otus]
  # otu.potu.link[[6]][otu.potu.link[[6]]$OTU %in% qs,]
  # 
  # # Drop these from the OTU Table
  # g <- prune_taxa(colnames(g@otu_table)[colnames(g@otu_table) %in% collapsed.otus], g)
  # 
  # # Some gOTU names are not in the soil data (not detected).  For these we add a zero column.
  # (absent.gOTUs <- collapsed.otus[!collapsed.otus %in% colnames(g@otu_table)])
  # 
  # # Make the tables comparable
  # dnamat <- g@otu_table
  # colnames(dnamat) <- names(collapsed.otus)[match(colnames(dnamat), collapsed.otus)]
  # dnamat <- dnamat[,match(colnames(stem.otu$abund[[1]]), colnames(dnamat))]
  # colnames(dnamat) <- colnames(stem.otu$abund[[1]])
  # dnamat <- as.data.frame(dnamat)
  # dnamat[is.na(dnamat)] <- 0
  # 
  # 
  # #################################
  # ### SAVE DATA FOR DOWNSTREAM ANALYSES
  # #################################
  # 
  # dnamat.pa <- 1 * (dnamat > 0)
  # dnamat.pa.G <- 1 * (colSums(dnamat.pa) > 0)
  # 
  # saveRDS(list(dnamat.pa.G=dnamat.pa.G, 
  #              dnamat.pa=dnamat.pa, 
  #              dnamat=dnamat, 
  #              stem.otu=stem.otu), 
  #         "Processed_data/stem-soil-Gsubsample-data-20250305.RDA")
  # 
  # 
  # #################################
  # ### Read data
  # #################################
  # 
  g <- readRDS("Processed_data/stem-soil-Gsubsample-data-20250305.RDA")

  # Function to compute standardized effect size
  ses <- function(obs, rand){
    return((obs - colMeans(rand))/apply(rand, 2, sd))
  }
  
  # Initialize holders for results
  conf_stats_obs <- matrix(nrow = 20, ncol = 13)
  conf_stats_ses <- matrix(nrow = 20, ncol = 13)
  
  for(r in 1:20){
    message(paste("working on radius", r))
    confus_obs <- caret::confusionMatrix(table(+(!g$dnamat.pa.G),
                                               +(!(1 * (g$stem.otu$abund[[r]][1,] > 0)))))
   
    conf_stats_obs[r, 1:12] <- c(seq(5,100,5)[r], confus_obs$byClass)
    
    conf_stats_obs[r, 13] <- mltools::mcc(TP=confus_obs$table[2,2],
                                          FP=confus_obs$table[2,1],
                                          TN=confus_obs$table[1,1],
                                          FN=confus_obs$table[1,2])
    
    # Randomization loop
    for(rand in 1:999){
      
      if(rand==1) {confus_rand <- list()}
      
      confus_rand[[rand]] <- caret::confusionMatrix(table(+(!g$dnamat.pa.G),
                                                 +(!(1 * (g$stem.otu$abund[[r]][rand,] > 0)))))
    }
      # Calculate SES for each class
      conf_stats_ses[r, 2:12] <- ses(conf_stats_obs[r, 2:12], 
                                    do.call(rbind, lapply(confus_rand, function(x) x$byClass)))
      
      
      mcc_all <- unlist(lapply(confus_rand, function(x) {
        mltools::mcc(TP=x$table[2,2],
                     FP=x$table[2,1],
                     TN=x$table[1,1],
                     FN=x$table[1,2])}))
      
      conf_stats_ses[r, 13] <- (conf_stats_obs[r, 13] - mean(mcc_all))/sd(mcc_all)

    }

  
  conf_stats_ses[, 1] <- seq(5,100,5)
  
  colnames(conf_stats_obs) <- c("r", names(confus_obs$byClass), "MCC")
  colnames(conf_stats_ses) <- c("r", names(confus_obs$byClass), "MCC")
  
  conf_stats_obs <- data.frame(conf_stats_obs)
  conf_stats_ses <- data.frame(conf_stats_ses)
  
  
  
  saveRDS(list(conf_stats_obs=conf_stats_obs, 
               conf_stats_ses=conf_stats_ses), 
          file=paste0("Processed_data/Conf_matrix_output-Gplots-20250310.RDA"))
  
  
  par(mfrow=c(1,2))
  plot(conf_stats_obs$r, conf_stats_obs$Sensitivity, ylim=c(0,1), pch=16, type='b',
       ylab="Confusion matrix statisitic",
       xlab="Neighborhood radius (m)")
  points(conf_stats_obs$r, conf_stats_obs$Specificity, pch=16, col='red', type='b')
  points(conf_stats_obs$r, conf_stats_obs$Balanced.Accuracy, pch=16, col='blue', type='b')
  legend('bottom', legend=c("Sensitivity","Specificity","Balanced Accuracy"), pch=16, 
         col=c('black', 'red', 'blue'), bty='n', lty=1)
  
  plot(conf_stats_ses$r, conf_stats_ses$Sensitivity, ylim=c(-3,3), pch=16, type='b',
       ylab="SES of confusion matrix statisitic",
       xlab="Neighborhood radius (m)")
  polygon(x=c(0,120,120,0), y=c(-1.96,-1.96,1.96,1.96), col='grey', lty=0)
  points(conf_stats_ses$r, conf_stats_ses$Sensitivity, ylim=c(-3,3), pch=16, type='b')
  points(conf_stats_ses$r, conf_stats_ses$Specificity, pch=16, col='red', type='b')
  points(conf_stats_ses$r, conf_stats_ses$Balanced.Accuracy, pch=16, col='blue', type='b')
  legend('bottom', legend=c("SES Sensitivity","SES Specificity","SES Balanced Accuracy"), 
         pch=16, col=c('black', 'red', 'blue'), bty='n', lty=1)
  abline(h=0, lty=2)
  

  
  
  
  
  
  
  
  