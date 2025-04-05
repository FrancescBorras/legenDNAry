library(phyloseq)
library(googlesheets4)
library(EDIutils)

### Load data from Glenn
# List of 12 phyloseq objects (lenient, etc.)
data <- readRDS("Processed_data/PR_eDNA-for-analysis_2025-04-01.RData")

### Load link between OTUs and RefIDs
# otu.potu.link <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_v1.rds")[[6]]
otupotuR1 <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R1.rds")[[6]]
otupotuR2 <- readRDS("Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R2.rds")[[6]]

for(data_selector in seq_along(data)){
  d <- data[[data_selector]]
  
  ### Determine which otu.potu mapping to use
  if(!grepl("repfiltered", names(data)[data_selector])) {
    otu.potu.link <- otupotuR1
  } else {
    otu.potu.link <- otupotuR2
  }
  
  ### Prune phyloseq object to the OTUs that appear in the soil data AND correspond to the reference library - this should be redundant with previous processing by Glenn
  d <- prune_taxa(otu.potu.link$OTU, d)

  ### Load the LFDP 2023 census extract data
  tree <- readRDS("Raw_data/LFDP2023-extract-20240510.RDA")
  tree16 <- readRDS("Raw_data/LFDP2016-extract-20240510.RDA")
  
  # list of species codes observed in tree data
  # writexl::write_xlsx(data.frame(code=colnames(tree$abund$`100`)), "Raw_data/LFDP2023-codes.xlsx")
  
  # Access raw data from 2016 census
  # census6 <-  readr::read_csv(file=EDIutils::read_data_entity(packageId = "knb-lter-luq.119.1545979", 
  #                                    entityId = "325c43057e0dd4e1cd6a13fa5125a76d"))
  # Species 
  # sort(unique(census6$Mnemonic[!census6$Mnemonic %in% colnames(tree$abund$`100`)]))
  
  ### Load LFDP species codes
  # codes_raw <- EDIutils::read_data_entity(packageId = "knb-lter-luq.118.215602", 
  #                                    entityId = "cffeaeca64642878797cbb1621fcf07e")
  # codes <- readr::read_csv(file = codes_raw)
  
  # # Read in the leaf OTUs googlesheet
  # leaf.otus <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1Ihe1dXkENK4zZ0lBsN1_1tDnaIvVD3u9oiNHjtMRHzo/edit?usp=sharing", sheet="Leaf OTUs")
  # 
  # leaf.otus$code <- substring(leaf.otus$mnemonic, 1, 6)
  # 
  # codes$POTU <- leaf.otus$'FASTA#'[match(codes$`SPECIES CODE`, leaf.otus$code)]
  # 
  # writexl::write_xlsx(codes, "Raw_data/LFDP-SPcodes.xlsx")
  
  codes <- readxl::read_xlsx("Raw_data/LFDP-SPcodes.xlsx")
  
  ### Keep only species codes that are present in the 2023 tree census
  codes <- codes[codes$LFDP2023==1,]
  
  ### Load the sample point coordinates
  sample_xy <- read.csv("Raw_data/LFDP-eDNA-xy-V2.csv")
  
  ### Sample "C24_P2" (a "normal" sampling location) is not missing from soil data.
  ### Glen confirmed that this was dropped due to 'bad' eDNA data (small library size of <10000 sequences)
  
  ### Prune phyloseq object to the sample points from the '40 point analysis'
  d <- prune_samples(grepl("normal", d@sam_data$factorlevel), d)
  
  # Get coordinates from sample points included
  xy <- d@sam_data[grepl("normal", d@sam_data$factorlevel),c("X","Y")]
  
  # Correct coordinates (they were off by 20 m in both directions...)
  xy <- xy-20
  
  # write.csv(xy, file="Raw_data/LFDP-sample39-coordinates.csv")
  
  ### Make the rows of the stem data match the rows of the eDNA otu-table data
  
  ### Create stem data for 2023
  stem_abund <- lapply(tree$abund, function(x) {
    rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
          x[88:1087,])})
  for(i in seq_along(stem_abund)){
    rownames(stem_abund[[i]]) <- c(rownames(d@otu_table), paste0('random', 1:1000))
  }
  
  stem_ba <- lapply(tree$ba, function(x) {
    rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
          x[88:1087,])})
  for(i in seq_along(stem_ba)){
    rownames(stem_ba[[i]]) <- c(rownames(d@otu_table), paste0('random', 1:1000))
  }
  
  stem_nn <- rbind(tree$nearest_sp[match(paste(xy$X, xy$Y),
                                         paste(sample_xy$PX, sample_xy$PY)),], 
                   tree$nearest_sp[88:1087,])
  rownames(stem_nn) <- c(rownames(d@otu_table), paste0('random', 1:1000))
  
  stem <- list(abund=stem_abund, ba=stem_ba, nearest_sp=stem_nn)
  
  
  #### Create stem data for 2016
  stem_abund16 <- lapply(tree16$abund, function(x) {
    rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
          x[88:1087,])})
  
  for(i in seq_along(stem_abund16)){
    rownames(stem_abund16[[i]]) <- c(rownames(d@otu_table), paste0('random',1:1000))
  }
  
  stem_ba16 <- lapply(tree16$ba, function(x) {
    rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
          x[88:1087,])})
  
  for(i in seq_along(stem_ba16)){
    rownames(stem_ba16[[i]]) <- c(rownames(d@otu_table), paste0('random',1:1000))
  }
  
  stem_nn16 <- rbind(tree16$nearest_sp[match(paste(xy$X, xy$Y),
                                             paste(sample_xy$PX, sample_xy$PY)),], tree16$nearest_sp[88:1087,])
  rownames(stem_nn16) <- c(rownames(d@otu_table), paste0('random',1:1000))
  
  stem16 <- list(abund=stem_abund16, ba=stem_ba16, nearest_sp=stem_nn16)
  
  
  ##### THIS BLOCK COLLAPSES THE DNA DATA TO THE gOTU CLUSTERS
  ### Identify all POTUs/OTUs we need to keep because they are assigned to a gOTU
  potus.keep <- codes$POTU[!is.na(codes$gOTU)]
  otus.keep <- otu.potu.link$OTU[match(potus.keep, otu.potu.link$abundance)]
  gotus.keep <- codes$gOTU[!is.na(codes$gOTU)]
  codes.keep <- codes$`SPECIES CODE`[!is.na(codes$gOTU)]
  
  # POTUs that we keep and to be collapsed
  codes.collapse.list <- split(codes.keep, gotus.keep)
  otus.collapse.list <- split(otus.keep, gotus.keep)
  
  # loop to collapse taxa in phyloseq objects
  for(i in seq_along(otus.collapse.list)){
    if(length(unique(otus.collapse.list[[i]][!is.na(otus.collapse.list[[i]])])) > 1){
      taxcollapse <- unique(otus.collapse.list[[i]][!is.na(otus.collapse.list[[i]])])
      taxcollapse <- taxcollapse[taxcollapse %in% rownames(d@tax_table)]
      d <- merge_taxa(d, taxcollapse, 1)
    }
  }
  
  ##### THIS BLOCK COLLAPSES THE STEM DATA TO THE gOTU CLUSTERS
  stem.otu <- list()
  for(r in seq_along(stem$abund)){
    
    # Initialize a matrix to hold collapsed results for a given radius
    tmat_abund <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    tmat_ba <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    tmat_nn <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    
    # loop to collapse taxa in phyloseq objects
    for(i in seq_along(codes.collapse.list)){
      tmp_abund <- stem$abund[[r]][,which(colnames(stem$abund[[1]]) %in% codes.collapse.list[[i]])]
      tmp_ba <- stem$ba[[r]][,which(colnames(stem$ba[[1]]) %in% codes.collapse.list[[i]])]
      tmp_nn <- stem$nearest_sp[,which(colnames(stem$nearest_sp) %in% codes.collapse.list[[i]])]
      
      if(is.matrix(tmp_abund)){
        tmat_abund[,i] <- rowSums(tmp_abund)
        tmat_ba[,i] <- rowSums(tmp_ba)
        tmat_nn[,i] <- rowSums(tmp_nn)
      } else {
        tmat_abund[,i] <- tmp_abund
        tmat_ba[,i] <- tmp_ba
        tmat_nn[,i] <- tmp_nn
      }
    }
    
    tmat_abund <- as.data.frame(tmat_abund)
    tmat_ba <- as.data.frame(tmat_ba)
    tmat_nn <- as.data.frame(tmat_nn)
    
    names(tmat_nn) <- names(tmat_ba) <- names(tmat_abund) <- names(codes.collapse.list)
    rownames(tmat_nn) <- rownames(tmat_ba) <- rownames(tmat_abund) <- rownames(stem$abund[[1]])
    
    stem.otu$abund[[r]] <- tmat_abund
    stem.otu$ba[[r]] <- tmat_ba
    stem.otu$nn <- tmat_nn
    
    names(stem.otu$ba)[r] <- names(stem.otu$abund)[r] <- names(stem$abund)[r]
  }
  
  
  ##### THIS BLOCK COLLAPSES THE 2016 STEM DATA TO THE gOTU CLUSTERS
  stem.otu16 <- list()
  for(r in seq_along(stem16$abund)){
    
    # Initialize a matrix to hold collapsed results for a given radius
    tmat_abund <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    tmat_ba <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    tmat_nn <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    
    # loop to collapse taxa in phyloseq objects
    for(i in seq_along(codes.collapse.list)){
      tmp_abund <- stem16$abund[[r]][,which(colnames(stem16$abund[[1]]) %in% codes.collapse.list[[i]])]
      tmp_ba <- stem16$ba[[r]][,which(colnames(stem16$ba[[1]]) %in% codes.collapse.list[[i]])]
      tmp_nn <- stem16$nearest_sp[,which(colnames(stem16$nearest_sp) %in% codes.collapse.list[[i]])]
      
      if(is.matrix(tmp_abund)){
        tmat_abund[,i] <- rowSums(tmp_abund)
        tmat_ba[,i] <- rowSums(tmp_ba)
        tmat_nn[,i] <- rowSums(tmp_nn)
      } else {
        tmat_abund[,i] <- tmp_abund
        tmat_ba[,i] <- tmp_ba
        tmat_nn[,i] <- tmp_nn
      }
    }
    
    tmat_abund <- as.data.frame(tmat_abund)
    tmat_ba <- as.data.frame(tmat_ba)
    tmat_nn <- as.data.frame(tmat_nn)
    
    names(tmat_nn) <- names(tmat_ba) <- names(tmat_abund) <- names(codes.collapse.list)
    rownames(tmat_nn) <- rownames(tmat_ba) <- rownames(tmat_abund) <- rownames(stem16$abund[[1]])
    
    stem.otu16$abund[[r]] <- tmat_abund
    stem.otu16$ba[[r]] <- tmat_ba
    stem.otu16$nn <- tmat_nn
    
    names(stem.otu16$ba)[r] <- names(stem.otu16$abund)[r] <- names(stem16$abund)[r]
  }
  
  
  ### To make the columns match, and add zero columns to OTU data
  
  # The new phyloseq object should have these OTUs as colnames
  collapsed.otus <- unlist(lapply(otus.collapse.list, function(x) x[!is.na(x)][1]))
  
  # The stem data has these names exactly
  if(!all(names(collapsed.otus) == names(stem.otu$abund[[1]]))){
    message("Warning: all(names(collapsed.otus) != names(stem.otu$abund[[1]]))")
  }
  if(!all(names(collapsed.otus) == names(stem.otu16$abund[[1]]))){
    message("Warning: all(names(collapsed.otus) != names(stem.otu16$abund[[1]]))")
  }
  
  # Some gOTU names are in the soil data but not the collapsed list because they aren't trees
  qs <- colnames(d@otu_table)[!colnames(d@otu_table) %in% collapsed.otus]
  otu.potu.link[otu.potu.link$OTU %in% qs,]
  
  # Drop these from the OTU Table
  d <- prune_taxa(colnames(d@otu_table)[colnames(d@otu_table) %in% collapsed.otus], d)
  
  # Some gOTU names are not in the soil data (not detected).  For these we add a zero column.
  absent.gOTUs <- collapsed.otus[!collapsed.otus %in% colnames(d@otu_table)]
  
  # Make the tables comparable
  dnamat <- d@otu_table
  colnames(dnamat) <- names(collapsed.otus)[match(colnames(dnamat), collapsed.otus)]
  dnamat <- dnamat[,match(colnames(stem.otu$abund[[1]]), colnames(dnamat))]
  colnames(dnamat) <- colnames(stem.otu$abund[[1]])
  dnamat <- as.data.frame(dnamat)
  dnamat[is.na(dnamat)] <- 0
  
  
  #################################
  ### ADD TRAIT DATA
  #################################
  trait_data <- read.csv("Raw_data/Fixed_traits.csv", row.names=1)
  traits <- matrix(nrow=length(codes.collapse.list), ncol=7)
  colnames(traits) <- colnames(trait_data)
  rownames(traits) <- names(codes.collapse.list)
  for(i in seq_along(codes.collapse.list)){
    tmp <- trait_data[rownames(trait_data) %in% codes.collapse.list[[i]],]
    if(nrow(tmp) > 0){ traits[i,] <- colMeans(tmp) } 
  }
  traits <- as.data.frame(traits)
  
  #################################
  ### SAVE DATA FOR DOWNSTREAM ANALYSES
  #################################
  
  # saveRDS(list(dnamat=dnamat, stem.otu=stem.otu, traits=traits, stem.otu16=stem.otu16), 
  #         "Processed_data/stem-soil-39pt-data-20240510.RDA")
  
  outfile <- paste0("Processed_data/stem-soil-40pt-data-",
                    names(data)[data_selector], 
                    "-20250405.RDA")
  
  saveRDS(list(dnamat=dnamat, stem.otu=stem.otu, traits=traits, stem.otu16=stem.otu16), outfile)
}


