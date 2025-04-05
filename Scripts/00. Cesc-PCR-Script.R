# Make sure you have access to a unix system
# install conda (or miniconda etc)
# create an environment in conda, call it a name and install:


#### THIS IS ALL IN THE TERMINAL - PLATE

# cutadapt
conda install -c bioconda cutadapt
# sabre
conda install -c bioconda sabre
# dos2unix
conda install -c conda-forge dos2unix
# R
conda install -c r r

# Open the terminal and navigate to your sequencing raw_data folder. For me it is here (this will be different for you)
cd /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData

##  Make some directories to process the samples in the plate format that they were PCR'ed in (Plate 1 and PLate 2)
mkdir plate1
mkdir plate2

# Move samples into appropriate plate directory

mv /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/*P1 /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1

mv /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/*P2 /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2

## go into plate1 directory
cd plate1
## extract all sequence files from sub-directories into main directory
mv **/*.fq.gz /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1
## exit that directory 
cd ..
## etc
cd plate2
mv **/*.fq.gz /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2
cd ..

## Lets demultiplex our replicates - plate 1

cd plate1

## take a look at sabre, what it does and how you run it on github
## based on the forward barcode read for plate1 that is designated in another text file in the same directory (plate1_barcode_data.txt)
## make shell script to run sabre 
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/barcode_data.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
mv rep1f ${bn}_rep1f.fq;
mv rep1r ${bn}_rep1r.fq;
mv rep2f ${bn}_rep2f.fq;
mv rep2r ${bn}_rep2r.fq;
mv rep3f ${bn}_rep3f.fq;
mv rep3r ${bn}_rep3r.fq;
mv rep4f ${bn}_rep4f.fq;
mv rep4r ${bn}_rep4r.fq;  done' > sabreplate1.sh

## Using cut adapt to remove primers - Doing replicate by replicate and make reverse primer an exact match to the corresponding reverse barcode that was not removed by sabre
## obviously you will need to adapt these sequences for the specific index used for each replicate of each plate..  rc = reverse compliment

#### Trnl base primers (i.e. tag IDs etc) are:
# Trnlg: GGGCAATCCTGAGCCAA
# Trnlh: CCATTGAGTCTCTGCACCTATC

#Plate 1 rep 1 reverse barcode is AGGAA
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGTTCCT
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): AGGAACCATTGAGTCTCTGCACCTATC

echo 'for i in *rep1f.fq; do bn=${i/rep1f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTTCCT -A AGGAACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep1out1.fq.gz --untrimmed-paired-output ${bn}.rep1out2.fq.gz -o ${bn}.rep1.trim1.fq.gz -p ${bn}.rep1.trim2.fq.gz ${bn}rep1f.fq ${bn}rep1r.fq; done' > plate1rep1.sh
bash plate1rep1.sh

#Plate 1 rep 2 reverse barcode is GAGTGG
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGCCACTC
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): GAGTGGCCATTGAGTCTCTGCACCTATC

echo 'for i in *rep2f.fq; do bn=${i/rep2f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -A GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep2out1.fq.gz --untrimmed-paired-output ${bn}.rep2out2.fq.gz -o ${bn}.rep2.trim1.fq.gz -p ${bn}.rep2.trim2.fq.gz ${bn}rep2f.fq ${bn}rep2r.fq; done' > plate1rep2.sh
bash plate1rep2.sh

#Plate 1 rep 3 reverse barcode is CCACGTC
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGGACGTGG
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): CCACGTCCCATTGAGTCTCTGCACCTATC

echo 'for i in *rep3f.fq; do bn=${i/rep3f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -A CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep3out1.fq.gz --untrimmed-paired-output ${bn}.rep3out2.fq.gz -o ${bn}.rep3.trim1.fq.gz -p ${bn}.rep3.trim2.fq.gz ${bn}rep3f.fq ${bn}rep3r.fq; done' > plate1rep3.sh
bash plate1rep3.sh

#Plate 1 rep 3 reverse barcode is TTCTCAGC
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGGCTGAGAA
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): TTCTCAGCCCATTGAGTCTCTGCACCTATC

echo 'for i in *rep4f.fq; do bn=${i/rep4f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -A TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep4out1.fq.gz --untrimmed-paired-output ${bn}.rep4out2.fq.gz -o ${bn}.rep4.trim1.fq.gz -p ${bn}.rep4.trim2.fq.gz ${bn}rep4f.fq ${bn}rep4r.fq; done' > plate1rep4.sh
bash plate1rep4.sh

## moving trimmed files to their own directory for dada2 analysis
## make a new folder called p1trimmed
mkdir p1trimmed
## move all your demultiplexed, oligo-trimmed files there.
mv *trim1.fq.gz /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/p1trimmed
mv *trim2.fq.gz /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/p1trimmed


### SAME PROCESS FOR PLATE 2

## Lets demultiplex our replicates - plate 1

cd ..
cd plate2

## take a look at sabre, what it does and how you run it on github
## based on the forward barcode read for plate1 that is designated in another text file in the same directory (plate1_barcode_data.txt)
## make shell script to run sabre 
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/barcode_data.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
mv rep1f ${bn}_rep1f.fq;
mv rep1r ${bn}_rep1r.fq;
mv rep2f ${bn}_rep2f.fq;
mv rep2r ${bn}_rep2r.fq;
mv rep3f ${bn}_rep3f.fq;
mv rep3r ${bn}_rep3r.fq;
mv rep4f ${bn}_rep4f.fq;
mv rep4r ${bn}_rep4r.fq;  done' > sabreplate2.sh

## Using cut adapt to remove primers - Doing replicate by replicate and make reverse primer an exact match to the corresponding reverse barcode that was not removed by sabre
## obviously you will need to adapt these sequences for the specific index used for each replicate of each plate..  rc = reverse compliment

#### Trnl base primers (i.e. tag IDs etc) are:
# Trnlg: GGGCAATCCTGAGCCAA
# Trnlh: CCATTGAGTCTCTGCACCTATC

#Plate 1 rep 1 reverse barcode is AGGAA
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGTTCCT
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): AGGAACCATTGAGTCTCTGCACCTATC

echo 'for i in *rep1f.fq; do bn=${i/rep1f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTTCCT -A AGGAACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep1out1.fq.gz --untrimmed-paired-output ${bn}.rep1out2.fq.gz -o ${bn}.rep1.trim1.fq.gz -p ${bn}.rep1.trim2.fq.gz ${bn}rep1f.fq ${bn}rep1r.fq; done' > plate1rep1.sh
bash plate1rep1.sh

#Plate 1 rep 2 reverse barcode is GAGTGG
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGCCACTC
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): GAGTGGCCATTGAGTCTCTGCACCTATC

echo 'for i in *rep2f.fq; do bn=${i/rep2f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -A GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep2out1.fq.gz --untrimmed-paired-output ${bn}.rep2out2.fq.gz -o ${bn}.rep2.trim1.fq.gz -p ${bn}.rep2.trim2.fq.gz ${bn}rep2f.fq ${bn}rep2r.fq; done' > plate1rep2.sh
bash plate1rep2.sh

#Plate 1 rep 3 reverse barcode is CCACGTC
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGGACGTGG
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): CCACGTCCCATTGAGTCTCTGCACCTATC

echo 'for i in *rep3f.fq; do bn=${i/rep3f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -A CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep3out1.fq.gz --untrimmed-paired-output ${bn}.rep3out2.fq.gz -o ${bn}.rep3.trim1.fq.gz -p ${bn}.rep3.trim2.fq.gz ${bn}rep3f.fq ${bn}rep3r.fq; done' > plate1rep3.sh
bash plate1rep3.sh

#Plate 1 rep 3 reverse barcode is TTCTCAGC
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGGCTGAGAA
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): TTCTCAGCCCATTGAGTCTCTGCACCTATC

echo 'for i in *rep4f.fq; do bn=${i/rep4f.fq}; cutadapt -a GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -A TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output ${bn}.rep4out1.fq.gz --untrimmed-paired-output ${bn}.rep4out2.fq.gz -o ${bn}.rep4.trim1.fq.gz -p ${bn}.rep4.trim2.fq.gz ${bn}rep4f.fq ${bn}rep4r.fq; done' > plate1rep4.sh
bash plate1rep4.sh

## moving trimmed files to their own directory for dada2 analysis
## make a new folder called p1trimmed
mkdir p2trimmed
## move all your demultiplexed, oligo-trimmed files there.
mv *trim1.fq.gz /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/p2trimmed
mv *trim2.fq.gz /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/p2trimmed

##### SWITCHING TO THE R TERMINAL NOW

library(dada2)
## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/p1trimmed")
path <- "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/p1trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))

###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[290:300]) # 
plotQualityProfile(fnRs[290:300])  #



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

## Warning that some samples did not pass the filter - need to find them

df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)
## Get vector of files that dont exist after filtering from original files
drops <- as.numeric(rownames(subset(df.fe, theref == "FALSE")))
drops
##redefining 'out' matrix so these files are not included
out <- out[file.exists(filtFs),]

### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
### to tell which numbers to add, loo at the output from the above subset commands.
filtFs <- filtFs[-drops]
filtRs <- filtRs[-drops]

## Learning error rates on these data

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

## getting some weird dog-leg fits typical of 

## There is a problem with these error functions with NOVASEQ data. Thus we need to hack the function fitter
## Ok we have an issue here with the error models due to NOvaSeq's binning of quality scores. Implementing the solution of 
## JacobRPrice, Here: https://github.com/benjjneb/dada2/issues/1307
## by modifying the error model, his trial 1: alter loess arguments (weights and span) & enforce monotonicity
library(magrittr)
library(dplyr)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

errF <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

## Having a look at new error plots with altered learn errors function
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

## dereplicating our sequence data set to remove redundency before applying dada denoising algorithm
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names - use sample.names1 if you needed to remove samples after filtering
names(derepFs) <- sample.names
##   Update/Change the sample.names object with the removed samples
sample.names <- sample.names[-c(34,35,46,47,51,62, 110, 314, 343, 344, 348, 354, 355, 356, 358, 360, 372, 382)]

names(derepFs) <- sample.names
### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference

dadaFPPs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)

dadaFpsPPs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
P1seqtab <- makeSequenceTable(mergers)
P1seqtabPP <- makeSequenceTable(mergersPP)
P1seqtabpsPP <- makeSequenceTable(mergers_psPP)

dim(P1seqtab)
dim(P1seqtabPP)
dim(P1seqtabpsPP)

table(nchar(getSequences(P1seqtab)))
table(nchar(getSequences(P1seqtabPP)))
table(nchar(getSequences(P1seqtabpsPP)))

## Removing Chimeras
P1seqtab.nochim <- removeBimeraDenovo(P1seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(P1seqtab.nochim)
P1seqtabPP.nochim <- removeBimeraDenovo(P1seqtabPP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(P1seqtabPP.nochim)
P1seqtabpsPP.nochim <- removeBimeraDenovo(P1seqtabpsPP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(P1seqtabpsPP.nochim)

sum(P1seqtab.nochim)/sum(P1seqtab)
sum(P1seqtabPP.nochim)/sum(P1seqtabPP)
sum(P1seqtabpsPP.nochim)/sum(P1seqtabpsPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
p1track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(P1seqtab.nochim))
colnames(p1track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "Chi-removed")
head(p1track[c(100:150),])

#### Tracking read loss through the pipeline - total pooled
p1trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN), rowSums(P1seqtabPP.nochim))
colnames(p1trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "Chi-removed")
head(p1trackPP)

#### Tracking read loss through the pipeline - psuedo pooled
p1trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN), rowSums(P1seqtabpsPP.nochim))
colnames(p1trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "Chi-removed")
head(p1trackpsPP)

#######################################
#################################
## Curating datasets with LULU
library(devtools)
install_github("tobiasgf/lulu")  
library(lulu)

## Make a directory in the plate1 dir for this in the terminal (seqfiles)
mkdir seqfiles

## getting sequences
uniquesToFasta(P1seqtab.nochim, fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/seqfiles/P1seqtab.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(P1seqtab.nochim)))))
uniquesToFasta(P1seqtabPP.nochim, fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/seqfiles/P1seqtabPP.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(P1seqtabPP.nochim)))))
uniquesToFasta(P1seqtabpsPP.nochim, fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/seqfiles/P1seqtabpsPP.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(P1seqtabpsPP.nochim)))))

## Make LULU OTU tables (OTUs: rows, samples: columns)
npool.lulu <- P1seqtab.nochim
colnames(npool.lulu) <- paste0("OTU", seq(length(getSequences(P1seqtab.nochim))))
npool.lulu <- t(npool.lulu)

pool.lulu <- P1seqtabPP.nochim
colnames(pool.lulu) <- paste0("OTU", seq(length(getSequences(P1seqtabPP.nochim))))
pool.lulu <- t(pool.lulu)

pspool.lulu <- P1seqtabpsPP.nochim
colnames(pspool.lulu) <- paste0("OTU", seq(length(getSequences(P1seqtabpsPP.nochim))))
pspool.lulu <- t(pspool.lulu)

########### THIS NEXT PART IN THE BASH TERMINAL WITH BLAST INSTALLED IN PATH (or conda environment)
################### 

cd /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/seqfiles
  #First produce a blast databases with the OTUs
makeblastdb -in P1seqtab.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in P1seqtabPP.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in P1seqtabpsPP.nochim.fasta -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
blastn -db P1seqtab.nochim.fasta -outfmt '6 qseqid sseqid pident' -out NoPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query P1seqtab.nochim.fasta
blastn -db P1seqtabPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out Pool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query P1seqtabPP.nochim.fasta
blastn -db P1seqtabpsPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out psPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query P1seqtabpsPP.nochim.fasta

### Now running LULU algorithm IN R
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate1/seqfiles")
NoPool_match_list.txt <- read.table("NoPool_match_list.txt")
str(NoPool_match_list.txt)
str(npool.lulu)
Pool_match_list.txt <- read.table("Pool_match_list.txt")
psPool_match_list.txt <- read.table("psPool_match_list.txt")
nopool.nochim.curated_result <- lulu(as.data.frame(npool.lulu), NoPool_match_list.txt)
pool.nochim.curated_result <- lulu(as.data.frame(pool.lulu), Pool_match_list.txt)
pspool.nochim.curated_result <- lulu(as.data.frame(pspool.lulu), psPool_match_list.txt)

## Check out how many OTUs were collapsed:
print(paste0("Not Pooled: ", "OTUs after Lulu: ", nopool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(nopool.nochim.curated_result$original_table)))
print(paste0("Pooled: ", "OTUs after Lulu: ", pool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pool.nochim.curated_result$original_table)))
print(paste0("PsudoPooled: ", "OTUs after Lulu: ", pspool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pspool.nochim.curated_result$original_table)))

## Making sequence tables compatible with summary routines below...
## Making vector of row numbers of OTUs kept after curation - correspond to column numbers from the precurated data
rownames(nopool.nochim.curated_result$curated_table)

nopool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(nopool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pspool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pspool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))

p1nopool.lulu <- t(nopool.nochim.curated_result$curated_table)
colnames(p1nopool.lulu) <- colnames(P1seqtab.nochim[, nopool.kept.otus])
p1pool.lulu <- t(pool.nochim.curated_result$curated_table)
colnames(p1pool.lulu) <- colnames(P1seqtabPP.nochim[, pool.kept.otus])
p1pspool.lulu <- t(pspool.nochim.curated_result$curated_table)
colnames(p1pspool.lulu) <- colnames(P1seqtabpsPP.nochim[, pspool.kept.otus])


##Going through the same process with plate 2 samples DADA2 and LULU


## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/p2trimmed")
path <- "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/p2trimmed"
fnFs <- sort(list.files(path, pattern=".trim1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern=".trim2.fq.gz", full.names = TRUE))

###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
reps <- rep(c("r1", "r2", "r3", "r4"))

sample.names <-  paste0(sample.names, reps)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[290:300]) # 
plotQualityProfile(fnRs[290:300])  #

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, # truncLen=c(45,45),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

## Warning that some samples did not pass the filter - need to find them

df.fe <-  data.frame(theref = file.exists(filtFs), therer = file.exists(filtRs), filef = filtRs, filer = filtRs)

## Get vector of files that dont exist after filtering from original files
drops <- as.numeric(rownames(subset(df.fe, theref == "FALSE")))
drops
##redefining 'out' matrix so these files are not included
out <- out[file.exists(filtFs),]

### dropping samples that are empty - these numbers will change plate by plate - BE CAREFUL THESE ARE FOR PLATE 1 ONLY
### to tell which numbers to add, loo at the output from the above subset commands.
filtFs <- filtFs[-drops]
filtRs <- filtRs[-drops]

## Learning error rates on these data

errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST=20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST=20)

## Having a look at error plots
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

## getting some weird dog-leg fits typical of novaseq data

## There is a problem with these error functions with NOVASEQ data. Thus we need to hack the function fitter
## Ok we have an issue here with the error models due to NOvaSeq's binning of quality scores. Implementing the solution of 
## JacobRPrice, Here: https://github.com/benjjneb/dada2/issues/1307
## by modifying the error model, his trial 1: alter loess arguments (weights and span) & enforce monotonicity
library(magrittr)
library(dplyr)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

errF <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

## Having a look at new error plots with altered learn errors function
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

## dereplicating our sequence data set to remove redundency before applying dada denoising algorithm
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names - use sample.names1 if you needed to remove samples after filtering
names(derepFs) <- sample.names
##   Update/Change the sample.names object with the removed samples
sample.names <- sample.names[-drops]

names(derepFs) <- sample.names
### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

###  also applying core inference algorith - but using pooling (all samples pooled together for sample inference

dadaFPPs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRPPs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

###  also applying core inference algorith - but using pseudo pooling (samples processed independent after sharing info
### between samples)

dadaFpsPPs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRpsPPs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Merging the paired-ends - true pooling
mergersPP <- mergePairs(dadaFPPs, filtFs, dadaRPPs, filtRs, verbose=TRUE)

### Merging the paired-ends - pseudo-pooling
mergers_psPP <- mergePairs(dadaFpsPPs, filtFs, dadaRpsPPs, filtRs, verbose=TRUE)

## Making sequence table
P2seqtab <- makeSequenceTable(mergers)
P2seqtabPP <- makeSequenceTable(mergersPP)
P2seqtabpsPP <- makeSequenceTable(mergers_psPP)

dim(P2seqtab)
dim(P2seqtabPP)
dim(P2seqtabpsPP)

table(nchar(getSequences(P2seqtab)))
table(nchar(getSequences(P2seqtabPP)))
table(nchar(getSequences(P2seqtabpsPP)))

## Removing Chimeras
P2seqtab.nochim <- removeBimeraDenovo(P2seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(P2seqtab.nochim)
P2seqtabPP.nochim <- removeBimeraDenovo(P2seqtabPP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(P2seqtabPP.nochim)
P2seqtabpsPP.nochim <- removeBimeraDenovo(P2seqtabpsPP, method="consensus", multithread=TRUE, verbose=TRUE)
dim(P2seqtabpsPP.nochim)

sum(P2seqtab.nochim)/sum(P2seqtab)
sum(P2seqtabPP.nochim)/sum(P2seqtabPP)
sum(P2seqtabpsPP.nochim)/sum(P2seqtabpsPP)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
P2track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(P2seqtab.nochim))
colnames(P2track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "Chi-removed")
head(P2track[c(100:150),])

#### Tracking read loss through the pipeline - total pooled
P2trackPP <- cbind(out, sapply(dadaFPPs, getN), sapply(dadaRPPs, getN), sapply(mergersPP, getN), rowSums(P2seqtabPP.nochim))
colnames(P2trackPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "Chi-removed")
head(P2trackPP[c(110:120),], 20)

#### Tracking read loss through the pipeline - psuedo pooled
P2trackpsPP <- cbind(out, sapply(dadaFpsPPs, getN), sapply(dadaRpsPPs, getN), sapply(mergers_psPP, getN), rowSums(P2seqtabpsPP.nochim))
colnames(P2trackpsPP) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "Chi-removed")
head(P2trackpsPP)

#######################################
#################################
## Curating datasets with LULU
library(devtools)
install_github("tobiasgf/lulu")  
library(lulu)

## Make a directory in the plate2 dir for this in the terminal (seqfiles)
mkdir seqfiles

## getting sequences
uniquesToFasta(P2seqtab.nochim, fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/seqfiles/P2seqtab.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(P2seqtab.nochim)))))
uniquesToFasta(P2seqtabPP.nochim, fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/seqfiles/P2seqtabPP.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(P2seqtabPP.nochim)))))
uniquesToFasta(P2seqtabpsPP.nochim, fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/seqfiles/P2seqtabpsPP.nochim.fasta", ids=paste0("OTU",  seq(length(getSequences(P2seqtabpsPP.nochim)))))

## Make LULU OTU tables (OTUs: rows, samples: columns)
npool.lulu <- P2seqtab.nochim
colnames(npool.lulu) <- paste0("OTU", seq(length(getSequences(P2seqtab.nochim))))
npool.lulu <- t(npool.lulu)

pool.lulu <- P2seqtabPP.nochim
colnames(pool.lulu) <- paste0("OTU", seq(length(getSequences(P2seqtabPP.nochim))))
pool.lulu <- t(pool.lulu)

pspool.lulu <- P2seqtabpsPP.nochim
colnames(pspool.lulu) <- paste0("OTU", seq(length(getSequences(P2seqtabpsPP.nochim))))
pspool.lulu <- t(pspool.lulu)

########### THIS NEXT PART IN THE BASH TERMINAL WITH BLAST INSTALLED IN PATH (or conda environment)
################### 

cd /Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/seqfiles
#First produce a blast databases with the OTUs
makeblastdb -in P2seqtab.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in P2seqtabPP.nochim.fasta -parse_seqids -dbtype nucl
makeblastdb -in P2seqtabpsPP.nochim.fasta -parse_seqids -dbtype nucl

# Then blast the OTUs against the database
blastn -db P2seqtab.nochim.fasta -outfmt '6 qseqid sseqid pident' -out NoPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query P2seqtab.nochim.fasta
blastn -db P2seqtabPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out Pool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query P2seqtabPP.nochim.fasta
blastn -db P2seqtabpsPP.nochim.fasta -outfmt '6 qseqid sseqid pident' -out psPool_match_list.txt -qcov_hsp_perc 80 -perc_identity 84 -query P2seqtabpsPP.nochim.fasta

### Now running LULU algorithm
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/plate2/seqfiles")
NoPool_match_list.txt <- read.table("NoPool_match_list.txt")
str(NoPool_match_list.txt)
str(npool.lulu)
Pool_match_list.txt <- read.table("Pool_match_list.txt")
psPool_match_list.txt <- read.table("psPool_match_list.txt")
nopool.nochim.curated_result <- lulu(as.data.frame(npool.lulu), NoPool_match_list.txt)
pool.nochim.curated_result <- lulu(as.data.frame(pool.lulu), Pool_match_list.txt)
pspool.nochim.curated_result <- lulu(as.data.frame(pspool.lulu), psPool_match_list.txt)

## Check out how many OTUs were collapsed:
print(paste0("Not Pooled: ", "OTUs after Lulu: ", nopool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(nopool.nochim.curated_result$original_table)))
print(paste0("Pooled: ", "OTUs after Lulu: ", pool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pool.nochim.curated_result$original_table)))
print(paste0("PsudoPooled: ", "OTUs after Lulu: ", pspool.nochim.curated_result$curated_count, " --- ", "OTUs before Lulu: ", nrow(pspool.nochim.curated_result$original_table)))

## Making sequence tables compatible with summary routines below...
## Making vector of row numbers of OTUs kept after curation - correspond to column numbers from the precurated data
rownames(nopool.nochim.curated_result$curated_table)

nopool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(nopool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))
pspool.kept.otus <-  as.numeric(gsub("OTU", "", rownames(pspool.nochim.curated_result$curated_table), ignore.case =FALSE, perl = TRUE))

p2nopool.lulu <- t(nopool.nochim.curated_result$curated_table)
colnames(p2nopool.lulu) <- colnames(P2seqtab.nochim[, nopool.kept.otus])
p2pool.lulu <- t(pool.nochim.curated_result$curated_table)
colnames(p2pool.lulu) <- colnames(P2seqtabPP.nochim[, pool.kept.otus])
p2pspool.lulu <- t(pspool.nochim.curated_result$curated_table)
colnames(p2pspool.lulu) <- colnames(P2seqtabpsPP.nochim[, pspool.kept.otus])

## now we have LULU curated dataframes. We will subtract the max # sequences in No Template Controls for each OTU from sample OTUs
# for plate one we have the matrices: "P1seqtab.nochim, P1seqtabPP.nochim, P1eqtabpsPP.nochim" (non-pooled, pooled and psuedopooled, chimera cleaned)
# and p1nopool.lulu, p1pool.lulu, p1pspool.lulu (non-pooled, pooled and psuedopooled, chimera cleaned and lulu curated)



## Generate summary info for the data to use later - can be modified as needed according to naming convention
library(stringr)
library(magrittr)

# a little function to summarize each matrix:

index.info <- function(x){
  y <- data.frame(matrix(NA,    
                         nrow = nrow(x),
                         ncol = 0))
  j <- as.data.frame(x) %>%          
    mutate(sampletype=case_when(startsWith(rownames(x), "E_") ~ "ExtractBlank", #Search for 000 followed by a digit from 1-6 followed by 0
                                startsWith(rownames(x), "N_") ~ "NTC", 
                                startsWith(rownames(x), "S_") ~ "Sample",
                                startsWith(rownames(x), "P_") ~ "Positive",
                                startsWith(rownames(x), "T_") ~ "Tagcatch",
                                startsWith(rownames(x), "T1_") ~ "Blankwell")) 
  y$rep <- str_sub(rownames(x), start= -1)
  samples <- rownames(x) %<>%
    gsub("E_", "", .) %>%
    gsub("N_", "", .) %>%
    gsub("S_", "", .) %>%
    gsub("P_", "", .) %>%
    gsub("T_", "", .) %>%
    gsub("T1_", "", .)
  y$plate <- substr(samples, nchar(samples)-3, nchar(samples)-2)
  y$sample <- str_sub(samples,1,nchar(samples)-5)
  rownames(x) <-  gsub("__.rep1.trim1.fq.gz", "", rownames(x), ignore.case =FALSE, perl = TRUE)
  y$sampletype <- j$sampletype
  y$full <- rownames(x)
  y$totseq <- rowSums(x)
  y$OTUs <- rowSums(x>0)
  rownames(y) <- y$full
  return(y)
}


## getting summary info for raw data processed and lulu OTU matrices
##dada processed P1 (not pooled, pooled, Psudopooled)
p1seqtab.nochim.index <- index.info(P1seqtab.nochim)
p1seqtabPP.nochim.index <- index.info(P1seqtabPP.nochim)
p1seqtabpsPP.nochim.index <- index.info(P1seqtabpsPP.nochim)

##lulu P1 (not pooled, pooled, Psudopooled)
p1nopool.lulu.index <- index.info(p1nopool.lulu)
p1pool.lulu.index <- index.info(p1pool.lulu)
p1pspool.lulu.index <- index.info(p1pspool.lulu)

##dada processed P2 (not pooled, pooled, Psudopooled)
p2seqtab.nochim.index <- index.info(P2seqtab.nochim)
p2seqtabPP.nochim.index <- index.info(P2seqtabPP.nochim)
p2seqtabpsPP.nochim.index <- index.info(P2seqtabpsPP.nochim)

##lulu P2 (not pooled, pooled, Psudopooled)
p2nopool.lulu.index <- index.info(p2nopool.lulu)
p2pool.lulu.index <- index.info(p2pool.lulu)
p2pspool.lulu.index <- index.info(p2pspool.lulu)

## lets look at some general diagnostics for each sequencing library (the different replicate PCRs and control samples)
## using callaghan decontam package to control for contamination - inspecting individual PCR library sizes
## only need to do this on one iteration per plate (of both dada & LULU as the pattern will hold across all according to library (plate) rather than anything else)
#BiocManager::install("decontam")
library(decontam)
library(ggplot2)
library(dplyr)
## Looking at plate 1
df <- p1nopool.lulu.index # Put sample_data into a ggplot-friendly data.frame
df <- df[order(df$totseq),]
df$Index <- seq(nrow(df))
df$alpha <- case_when(
  df$sampletype == "Sample" ~ 0.8,
  df$sampletype == "Positive" ~ 0.8,
  (df$sampletype != "Sample" & df$sampletype != "Positive"  ~ 1))

## plotting library sizes per replicate on a log scale - looking for nasty (i.e. large and sample like) values of control samples (except positive samples)
ggplot(data=df, aes(x=Index, y=log(totseq), color=sampletype)) + geom_point(aes(alpha = alpha)) + #ylim(0,200000) +
  facet_grid(cols = vars(rep)) +
  ggtitle("Log library size per sample (points) per PCR replicate (panels)")

ggplot(data=df, aes(x=Index, y=totseq, color=sampletype)) + geom_point(aes(alpha = alpha)) + #ylim(0,200000) +
  facet_grid(cols = vars(rep)) +
  ggtitle("Log library size per sample (points) per PCR replicate (panels)")

## For plate 1 we can see a problem with replicate 2 PCR - many of the control samples encroach in the sample library sizes and many sample library sizes are small

## Looking at plate 2
df <- p2nopool.lulu.index # Put sample_data into a ggplot-friendly data.frame
df <- df[order(df$totseq),]
df$Index <- seq(nrow(df))
df$alpha <- case_when(
  df$sampletype == "Sample" ~ 0.8,
  df$sampletype == "Positive" ~ 0.8,
  (df$sampletype != "Sample" & df$sampletype != "Positive"  ~ 1))

## plotting library sizes per replicate on a log scale - looking for nasty (i.e. large and sample like) values of control samples (except positive samples)
ggplot(data=df, aes(x=Index, y=log(totseq), color=sampletype)) + geom_point(aes(alpha = alpha)) + #ylim(0,200000) +
  facet_grid(cols = vars(rep))+
  ggtitle("Log library size per sample (points) per PCR replicate (panels)")

## For plate 2 we can see a problem with replicate 4 PCR - Consistently very much smaller library sizes than the other PCR replicates
##checking again without log scale
ggplot(data=df, aes(x=Index, y=totseq, color=sampletype)) + geom_point(aes(alpha = alpha)) + #ylim(0,200000) +
  facet_grid(cols = vars(rep))+
  ggtitle("Library size per sample (points) per PCR replicate (panels)")

## OK so:
## from plate 1 we should drop PCR replicate 2 due to inconsistencies with other replicates between samples and a seeming contamination issue.
## from plate 2 we shoudl drop PCR replicate 4. There was clearly some systematic error here that resulted in poor library representation.

# a quick function for these:
droprep2p1 <- function(x){
  x[!grepl("r2", rownames(x)),]
}
droprep4p2 <- function(x){
  x[!grepl("r4", rownames(x)),]
}
### Need to drop these replicates from plate1 and plate2 from both OTU matrices and the .index dataframes - putting them in lists for lapply
p1indexinfo <- list(p1seqtab.nochim.index, p1seqtabPP.nochim.index , p1seqtabpsPP.nochim.index, p1nopool.lulu.index, p1pool.lulu.index, p1pspool.lulu.index)
names(p1indexinfo) <- c("p1seqtab.nochim", "p1seqtabPP.nochim" , "p1seqtabpsPP.nochim", "p1nopool.lulu", "p1pool.lulu", "p1pspool.lulu")

p2indexinfo <- list(p2seqtab.nochim.index, p2seqtabPP.nochim.index , p2seqtabpsPP.nochim.index, p2nopool.lulu.index, p2pool.lulu.index, p2pspool.lulu.index)
names(p2indexinfo) <- c("p2seqtab.nochim", "p2seqtabPP.nochim" , "p2seqtabpsPP.nochim", "p2nopool.lulu", "p2pool.lulu", "p2pspool.lulu")

p1data <- list(P1seqtab.nochim, P1seqtabPP.nochim , P1seqtabpsPP.nochim, p1nopool.lulu, p1pool.lulu, p1pspool.lulu)
names(p1data) <- c("p1seqtab.nochim", "p1seqtabPP.nochim" , "p1seqtabpsPP.nochim", "p1nopool.lulu", "p1pool.lulu", "p1pspool.lulu") 

p2data <- list(P2seqtab.nochim, P2seqtabPP.nochim , P2seqtabpsPP.nochim, p2nopool.lulu, p2pool.lulu, p2pspool.lulu)
names(p2data) <- c("p2seqtab.nochim", "p2seqtabPP.nochim" , "p2seqtabpsPP.nochim", "p2nopool.lulu", "p2pool.lulu", "p2pspool.lulu")


## making single dataframe with all samples so indexing functions later can call it regardless of plate

sample.index.info <- rbind(p1indexinfo[[1]], p2indexinfo[[1]])
sample.index.info$fullminrep <- substr(sample.index.info$full,1,nchar(sample.index.info$full)-2)
v <- subset(sample.index.info, sampletype == "ExtractBlank")
unique(v$sample)


### droping the identified bad replicates from each plate
p1indexinfo1 <- lapply(p1indexinfo, droprep2p1)
p1data1 <-  lapply(p1data, droprep2p1)

p2indexinfo1 <- lapply(p2indexinfo, droprep4p2)
p2data1 <-  lapply(p2data, droprep4p2)

##exporting this info for results summary
write.csv(p1indexinfo1$p1seqtab.nochim, "plate1postrepremoval.csv")
write.csv(p2indexinfo1$p2seqtab.nochim, "plate2postrepremoval.csv")

## Now we have curated lists of data and index info where spurious PCR replicates are removed. We can now discount OTUs that appear in control samples
## and discard all control samples from data matrices, then move on to analysing the sample data.

## First lets assess the rate of tagjumping across our entire libraries for each library made (plate 1 & Plate 2)
# PLate 1 rates of tag jumping:
lapply(p1indexinfo1, function(x) sum(x$totseq))
lapply(p1indexinfo1, function(x) sum(subset(x, sampletype == "Tagcatch")$totseq)/sum(x$totseq))
# PLate 2 rates of tag jumping:
lapply(p2indexinfo1, function(x) sum(x$totseq))
lapply(p2indexinfo1, function(x) sum(subset(x, sampletype == "Tagcatch")$totseq)/sum(x$totseq))

# So tag jump rates between 5 and 7 reads per 100,000 - nothing to be concerned about
# so 1 per 20,000 reads.. 

## We can control for NTCs first. Here we simply take the largest value in each OTU that occurs in any NTC and subtract that from all samples in the plate.
## Function to do so:

ntc.change <- function(x){
  mind <- apply(x[grep("N_", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

p1data1.ntc <- lapply(p1data1, ntc.change)
p2data1.ntc <- lapply(p2data1, ntc.change)

## Now to control for extraction blanks - function to 1) batch samples into their cycling group, then extraction batch 2) select the highest number of sequences in each OTU in the extraction bkank & 
## stubtract that number from the samples OTUs (if it exists there) & and make it so the lowest number in the OTUs is zero (prevent negative values)

setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data")

sample_data <- read.csv("sample_data.csv")
str(sample_data)

# Here x is the list of matrices, y is the sample metadata (extract blank, experiment etc) - function to group by cycling group then extract blank
ntc.to.blankcontrol <-  function (x){
  z = merge(sample.index.info, sample_data, by.x='fullminrep', by.y="newname2", all.x=TRUE)
  z$full.name =  paste(substr(z$name,start=1,stop=2), z$full, sep = "")
  row.names(z) <- z$full.name
  z = merge(x, z, by =  'row.names', all.x=TRUE)
  rownames(z) <- z$Row.names 
  z = split(z, z$cyclinggroup)
  z = lapply(z, function(x) split(x, x$ExtractBlank))
  dropnames <- colnames(z[[1]][[1]][, c(which(nchar(colnames(z[[1]][[1]]))< 20))])
  z <- lapply(z, function(x) lapply(x, function(x) x[!(names(x) %in% dropnames)]))
  z = lapply(z, function(x) lapply(x, as.matrix))
  return(z)
}

p1data1.ntc.exgroups <- lapply(p1data1.ntc, ntc.to.blankcontrol)
p2data1.ntc.exgroups <- lapply(p2data1.ntc, ntc.to.blankcontrol)

p1data1.ntc.exgroups <-  unlist(p1data1.ntc.exgroups, recursive=FALSE)
p2data1.ntc.exgroups <-  unlist(p2data1.ntc.exgroups, recursive=FALSE)
## Sweep through cycling groups and extract blank groups and remove the largest value in each OTU that occurs in any Extract Control and subtract that from all samples in the cycling group:
## Function to do so:

exb.change <- function(x){
  mind <- apply(x[grep("E_", rownames(x)), ], 2, function(y) max(y, na.rm = TRUE))
  x1 <- sweep(x, 2, mind)
  x1 <- pmax(x1,0)
  return(x1)
}

## Subtracting Extract controls
## some warnings here but just for control samples with no exblanks - ignore
p1data1.ntc.excon <- lapply(p1data1.ntc.exgroups, function(k) lapply(k, exb.change)) ## some warnings here but just for control samples with no exblanks - ignore
p2data1.ntc.excon <- lapply(p2data1.ntc.exgroups, function(k) lapply(k, exb.change))

p1data1.ntc.excon1 <-  unlist(p1data1.ntc.excon, recursive=FALSE)
p2data1.ntc.excon1 <-  unlist(p2data1.ntc.excon, recursive=FALSE)

## Joining list elements after extract blank control and removing anything that isnt a sample
# Function to do so:
combine_matrices_by_prefix <- function(input_list, prefix_patterns) {
  # Initialize a list to store the combined matrices for each prefix
  combined_matrices <- list()
  
  # Iterate over each prefix pattern
  for (pattern in prefix_patterns) {
    # Find list names matching the pattern
    matching_names <- grep(pattern, names(input_list), value = TRUE)
    
    # Initialize a matrix to store the combined data for the current prefix
    combined_matrix <- NULL
    
    # Iterate over matching names for the current prefix and bind them together into a single matrix
    for (name in matching_names) {
      if (is.null(combined_matrix)) {
        combined_matrix <- input_list[[name]]
      } else {
        combined_matrix <- rbind(combined_matrix, input_list[[name]])
      }
    }
    
    # Store the combined matrix for the current prefix in the list
    combined_matrices[[pattern]] <- combined_matrix
  }
  
  # Return the list of combined matrices
  return(combined_matrices)
}

prefix_patterns_p1 <- c("p1seqtab.nochim", "p1seqtabPP.nochim", "p1seqtabpsPP.nochim", 
                     "p1nopool.lulu", "p1pool.lulu", "p1pspool.lulu")

prefix_patterns_p2 <- c("p2seqtab.nochim", "p2seqtabPP.nochim", "p2seqtabpsPP.nochim", 
                        "p2nopool.lulu", "p2pool.lulu", "p2pspool.lulu")

p1.ntc.excon <- combine_matrices_by_prefix(p1data1.ntc.excon1, prefix_patterns_p1)
p2.ntc.excon <- combine_matrices_by_prefix(p2data1.ntc.excon1, prefix_patterns_p2) # note that this removes all blank wells controls 

## removing all control samples so that just samples remain. 

p1.ntc.excon.samples <- lapply(p1.ntc.excon, function(x) x[grep("S_", rownames(x)), ])
p2.ntc.excon.samples <- lapply(p2.ntc.excon, function(x) x[grep("S_", rownames(x)), ])

View(p1.ntc.excon.samples$p1seqtab.nochim)
View(p2.ntc.excon.samples$p2seqtab.nochim)

## making lists where elements are individual samples and performing replicate control (i.e. only keeping OTUs in 1/3, 2/3 or 3/3 replicates)
## firstly a function to make sample-wise lists from current lists
controlledblanks.to.samplelist <-  function (x){
  z = merge(sample.index.info, sample_data, by.x='fullminrep', by.y="newname2", all.x=TRUE)
  z$full.name =  paste(substr(z$name,start=1,stop=2), z$full, sep = "")
  row.names(z) <- z$full.name
  z = merge(x, z, by =  'row.names', all.x=TRUE)
  z <- split(z, z$fullminrep)
  dropnames <- colnames(z[[1]][, c(which(nchar(colnames(z[[1]]))< 20))])
  z <- lapply(z, function(x) x[!(names(x) %in% dropnames)])
  z = lapply(z, as.matrix)
  return(z)
}

p1.ntc.excon.samplelists <- lapply(p1.ntc.excon.samples, controlledblanks.to.samplelist)
p2.ntc.excon.samplelists <- lapply(p2.ntc.excon.samples, controlledblanks.to.samplelist)



## some functions to throw out OTUs (turn to zero) if the OTU is not present in 2 or 3 replicates
rep.groups2 <- function(x){
  r2 <- apply(x, 2, function(c) replace(c, sum(c!=0)<2, 0))
  return(r2)
}
rep.groups3 <- function(x){
  r3 <- apply(x, 2, function(c) replace(c, sum(c!=0)<3, 0))
  return(r3)
}

## 2/3 replicate screening
p1.ntc.excon.samplelistsR2 <- lapply(p1.ntc.excon.samplelists, function(x) lapply(x, rep.groups2))
p2.ntc.excon.samplelistsR2 <- lapply(p2.ntc.excon.samplelists, function(x) lapply(x, rep.groups2))

## 3/3 replicate screening
p1.ntc.excon.samplelistsR3 <- lapply(p1.ntc.excon.samplelists, function(x) lapply(x, rep.groups3))
p2.ntc.excon.samplelistsR3 <- lapply(p2.ntc.excon.samplelists, function(x) lapply(x, rep.groups3))

## All replicate control has been performed. Now to combine replicates into single row per sample and join plates for final datasets
## Joining matrices together between plates now
library(dada2)
p1.ntc.excon.samplelists <- lapply(p1.ntc.excon.samplelists, function(x) lapply(x, as.matrix))
p2.ntc.excon.samplelists <- lapply(p2.ntc.excon.samplelists, function(x) lapply(x, as.matrix))
p1.ntc.excon.samplelistsR2 <- lapply(p1.ntc.excon.samplelistsR2, function(x) lapply(x, as.matrix))
p2.ntc.excon.samplelistsR2 <- lapply(p2.ntc.excon.samplelistsR2, function(x) lapply(x, as.matrix))
p1.ntc.excon.samplelistsR3 <- lapply(p1.ntc.excon.samplelistsR3, function(x) lapply(x, as.matrix))
p2.ntc.excon.samplelistsR3 <- lapply(p2.ntc.excon.samplelistsR3, function(x) lapply(x, as.matrix))

to.one.matrix <- function(x){
  lah <- do.call(rbind.data.frame, x)
  rownames(lah) <- names(x)
  colnames(lah) <- names(x[[1]])
  lah <- as.matrix(lah)
  return(lah)
}
## Making single matrices per plate with each sample on one row by adding within OTUs per sample
make.to.single <-  function (x){
  y <- lapply(x, function(w) lapply(w, function(w) colSums(w)))
  z <- lapply(y, function(w) to.one.matrix(w))
  return(z)
}


p1.ntc.excon.R1 <- make.to.single(p1.ntc.excon.samplelists)
p1.ntc.excon.R2 <- make.to.single(p1.ntc.excon.samplelistsR2)
p1.ntc.excon.R3 <- make.to.single(p1.ntc.excon.samplelistsR3)

p2.ntc.excon.R1 <- make.to.single(p2.ntc.excon.samplelists)
p2.ntc.excon.R2 <- make.to.single(p2.ntc.excon.samplelistsR2)
p2.ntc.excon.R3 <- make.to.single(p2.ntc.excon.samplelistsR3)

## Joining matrices from plates to a single matrix for each of the data and lulu iterations

join.plates <- function(x, y) {
  
  combined_matrices <- list()
  
  for (k in seq_along(x)) {
    combined_matrix <- mergeSequenceTables(x[[k]], y[[k]])
    combined_matrices[[k]] <- combined_matrix
  }
  
  return(combined_matrices)
}

R1joined <- join.plates(p1.ntc.excon.R1, p2.ntc.excon.R1)
R2joined <- join.plates(p1.ntc.excon.R2, p2.ntc.excon.R2)
R3joined <- join.plates(p1.ntc.excon.R3, p2.ntc.excon.R3)
View(R1joined[[6]])
write.csv(R2joined[[6]], "testing2rep.csv")
View(R2joined[[6]])
View(R3joined[[2]])

vecnames <- c("dada.nopool.nochim", "dada.pooled.nochim" , "dada.pspool.nochim", "dada.nopool.nc.lulu", "dada.pooled.nc.lulu", "dada.pspool.nc.lulu") 
names(R1joined) <- vecnames
names(R2joined) <- vecnames
names(R3joined) <- vecnames

## Now getting sequences in order to both assign species according to our reference library and also to blast the other sequences

uniquesToFasta(as.matrix(R1joined[[1]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.nopool.fasta", ids=paste0("OTU", seq(length(getSequences(R1joined[[1]])))))
uniquesToFasta(as.matrix(R1joined[[2]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.pool.fasta", ids=paste0("OTU", seq(length(getSequences(R1joined[[2]])))))
uniquesToFasta(as.matrix(R1joined[[3]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.pspool.fasta", ids=paste0("OTU", seq(length(getSequences(R1joined[[3]])))))
uniquesToFasta(as.matrix(R1joined[[4]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1lulu.nopool.fasta", ids=paste0("OTU", seq(length(getSequences(R1joined[[4]])))))
uniquesToFasta(as.matrix(R1joined[[5]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1lulu.pool.fasta", ids=paste0("OTU", seq(length(getSequences(R1joined[[5]])))))
uniquesToFasta(as.matrix(R1joined[[6]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1lulu.pspool.fasta", ids=paste0("OTU", seq(length(getSequences(R1joined[[6]])))))

uniquesToFasta(as.matrix(R2joined[[1]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.nopool.fasta", ids=paste0("OTU", seq(length(getSequences(R2joined[[1]])))))
uniquesToFasta(as.matrix(R2joined[[2]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.pool.fasta", ids=paste0("OTU", seq(length(getSequences(R2joined[[2]])))))
uniquesToFasta(as.matrix(R2joined[[3]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.pspool.fasta", ids=paste0("OTU", seq(length(getSequences(R2joined[[3]])))))
uniquesToFasta(as.matrix(R2joined[[4]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2lulu.nopool.fasta", ids=paste0("OTU", seq(length(getSequences(R2joined[[4]])))))
uniquesToFasta(as.matrix(R2joined[[5]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2lulu.pool.fasta", ids=paste0("OTU", seq(length(getSequences(R2joined[[5]])))))
uniquesToFasta(as.matrix(R2joined[[6]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2lulu.pspool.fasta", ids=paste0("OTU", seq(length(getSequences(R2joined[[6]])))))

uniquesToFasta(as.matrix(R3joined[[1]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R3dada.nopool.fasta", ids=paste0("OTU", seq(length(getSequences(R3joined[[1]])))))
uniquesToFasta(as.matrix(R3joined[[2]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R3dada.pool.fasta", ids=paste0("OTU", seq(length(getSequences(R3joined[[2]])))))
uniquesToFasta(as.matrix(R3joined[[3]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R3dada.pspool.fasta", ids=paste0("OTU", seq(length(getSequences(R3joined[[3]])))))
uniquesToFasta(as.matrix(R3joined[[4]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R3lulu.nopool.fasta", ids=paste0("OTU", seq(length(getSequences(R3joined[[4]])))))
uniquesToFasta(as.matrix(R3joined[[5]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R3lulu.pool.fasta", ids=paste0("OTU", seq(length(getSequences(R3joined[[5]])))))
uniquesToFasta(as.matrix(R3joined[[6]]), fout="/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R3lulu.pspool.fasta", ids=paste0("OTU", seq(length(getSequences(R3joined[[6]])))))


R1otupool.spc.d1 <- assignSpecies(getSequences(R1joined[[1]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R1otupool.spc.d2 <- assignSpecies(getSequences(R1joined[[2]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R1otupool.spc.d3 <- assignSpecies(getSequences(R1joined[[3]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R1otupool.spc.l4 <- assignSpecies(getSequences(R1joined[[4]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R1otupool.spc.l5 <- assignSpecies(getSequences(R1joined[[5]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R1otupool.spc.l6 <- assignSpecies(getSequences(R1joined[[6]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )

R2otupool.spc.d1 <- assignSpecies(getSequences(R2joined[[1]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R2otupool.spc.d2 <- assignSpecies(getSequences(R2joined[[2]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R2otupool.spc.d3 <- assignSpecies(getSequences(R2joined[[3]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R2otupool.spc.l4 <- assignSpecies(getSequences(R2joined[[4]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R2otupool.spc.l5 <- assignSpecies(getSequences(R2joined[[5]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R2otupool.spc.l6 <- assignSpecies(getSequences(R2joined[[6]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )

R3otupool.spc.d1 <- assignSpecies(getSequences(R3joined[[1]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R3otupool.spc.d2 <- assignSpecies(getSequences(R3joined[[2]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R3otupool.spc.d3 <- assignSpecies(getSequences(R3joined[[3]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R3otupool.spc.l4 <- assignSpecies(getSequences(R3joined[[4]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R3otupool.spc.l5 <- assignSpecies(getSequences(R3joined[[5]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )
R3otupool.spc.l6 <- assignSpecies(getSequences(R3joined[[6]]), "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/X204SC24022146-Z01-F001/01.RawData/FASTA_POTURD.fasta", allowMultiple=TRUE )


### getting a single matrix from list of samples (with single row)
R1SppAssList <- list(R1otupool.spc.d1, R1otupool.spc.d2,R1otupool.spc.d3, R1otupool.spc.l4, R1otupool.spc.l5, R1otupool.spc.l6)
names(R1SppAssList)<- vecnames

R2SppAssList <- list(R2otupool.spc.d1, R2otupool.spc.d2,R2otupool.spc.d3, R2otupool.spc.l4, R2otupool.spc.l5, R2otupool.spc.l6)
names(R2SppAssList)<- vecnames


getunassignedprep <- function(x) {
  x <- as.data.frame(x)
  x$OTU <- paste0("OTU", seq(nrow(x)))
  return(x)
}

R1SppAssListprep <- lapply(R1SppAssList, getunassignedprep)
R2SppAssListprep <- lapply(R2SppAssList, getunassignedprep)

## getting assigned and unassigned
speciesassigned <- function(x){
  x1 <- as.data.frame(x)
  x1$sequence <- rownames(x1)
  x1 <- x1[!is.na(x1$Genus),]
  x1$abundance <- x1$Genus
  x1[is.na(x1)] <- 0
  return(x1)
}

no.speciesassigned.to.fasta <- function(x){
  x1 <- as.data.frame(x)
  x1$sequence <- rownames(x1)
  x1 <- x1[is.na(x1$Genus),]
  x1$abundance <- x1$Genus
  x1[is.na(x1)] <- 0
  return(x1)
}

R1SppAssList <- lapply(R1SppAssListprep, speciesassigned)
R1SppUnAssList <- lapply(R1SppAssListprep, no.speciesassigned.to.fasta)

R2SppAssList <- lapply(R2SppAssListprep, speciesassigned)
R2SppUnAssList <- lapply(R2SppAssListprep, no.speciesassigned.to.fasta)

R1blast.d1 <- uniquesToFasta(R1SppUnAssList$dada.nopool.nochim, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.nopool.nochim.nospec1.fasta", ids=R1SppUnAssList$dada.nopool.nochim$OTU)
R1blast.d2 <- uniquesToFasta(R1SppUnAssList$dada.pooled.nochim, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.pool.nochim.nospec1.fasta", ids=R1SppUnAssList$dada.pooled.nochim$OTU)
R1blast.d3 <- uniquesToFasta(R1SppUnAssList$dada.pspool.nochim, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.pspool.nochim.nospec1.fasta", ids=R1SppUnAssList$dada.pspool.nochim$OTU)
R1blast.l4 <- uniquesToFasta(R1SppUnAssList$dada.nopool.nc.lulu, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.nopool.nc.lulu.nospec1.fasta", ids=R1SppUnAssList$dada.nopool.nc.lulu$OTU)
R1blast.l5 <- uniquesToFasta(R1SppUnAssList$dada.pooled.nc.lulu, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.pool.nc.lulu.nospec1.fasta", ids=R1SppUnAssList$dada.pooled.nc.lulu$OTU)
R1blast.l6 <- uniquesToFasta(R1SppUnAssList$dada.pspool.nc.lulu, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R1dada.pspool.nc.lulu.nospec1.fasta", ids=R1SppUnAssList$dada.pspool.nc.lulu$OTU)

R2blast.d1 <- uniquesToFasta(R2SppUnAssList$dada.nopool.nochim, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.nopool.nochim.nospec1.fasta", ids=R2SppUnAssList$dada.nopool.nochim$OTU)
R2blast.d2 <- uniquesToFasta(R2SppUnAssList$dada.pooled.nochim, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.pool.nochim.nospec1.fasta", ids=R2SppUnAssList$dada.pooled.nochim$OTU)
R2blast.d3 <- uniquesToFasta(R2SppUnAssList$dada.pspool.nochim, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.pspool.nochim.nospec1.fasta", ids=R2SppUnAssList$dada.pspool.nochim$OTU)
R2blast.l4 <- uniquesToFasta(R2SppUnAssList$dada.nopool.nc.lulu, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.nopool.nc.lulu.nospec1.fasta", ids=R2SppUnAssList$dada.nopool.nc.lulu$OTU)
R2blast.l5 <- uniquesToFasta(R2SppUnAssList$dada.pooled.nc.lulu, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.pool.nc.lulu.nospec1.fasta", ids=R2SppUnAssList$dada.pooled.nc.lulu$OTU)
R2blast.l6 <- uniquesToFasta(R2SppUnAssList$dada.pspool.nc.lulu, fout= "/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data/R2dada.pspool.nc.lulu.nospec1.fasta", ids=R2SppUnAssList$dada.pspool.nc.lulu$OTU)

## Blasted sequences not assigned with reference library (used blast nt_euk database), imported to MEGAN, used naive LCA algorithm 
## with a 0.95 min identity cut off to find the lowest common ancestor from BLAST matches

## exported taxID of LCA outcome and output below:
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data")

R1blast.d1.m <-  read.csv("R1dada.nopool.nochim.nospec1-ex.txt", header =FALSE)
R1blast.d2.m <-  read.csv("R1dada.pool.nochim.nospec1-ex.txt", header =FALSE)
R1blast.d3.m <-  read.csv("R1dada.pspool.nochim.nospec1-ex.txt", header =FALSE)
R1blast.l4.m <-  read.csv("R1dada.nopool.nc.lulu.nospec1-ex.txt", header =FALSE)
R1blast.l5.m <-  read.csv("R1dada.pool.nc.lulu.nospec1-ex.txt", header =FALSE)
R1blast.l6.m <-  read.csv("R1dada.pspool.nc.lulu.nospec1-ex.txt", header =FALSE)

R2blast.d1.m <-  read.csv("R2dada.nopool.nochim.nospec1-ex.txt", header =FALSE)
R2blast.d2.m <-  read.csv("R2dada.pool.nochim.nospec1-ex.txt", header =FALSE)
R2blast.d3.m <-  read.csv("R2dada.pspool.nochim.nospec1-ex.txt", header =FALSE)
R2blast.l4.m <-  read.csv("R2dada.nopool.nc.lulu.nospec1-ex.txt", header =FALSE)
R2blast.l5.m <-  read.csv("R2dada.pool.nc.lulu.nospec1-ex.txt", header =FALSE)
R2blast.l6.m <-  read.csv("R2dada.pspool.nc.lulu.nospec1-ex.txt", header =FALSE)

setwd("/Users/glennd/Downloads") ## where my locally constrcuted DB is
## Used locally constructured sql database for taxonomy
#install.packages("taxonomizr")
#install.packages("taxize")
#install.packages("taxizedb")
#remotes::install_github("ropensci/taxizedb")
library("taxonomizr")
library("taxize")
library("taxizedb")
####
#### ONLY RUN BELOW ONCE TO CONSTRUCT TAXONOMY DB
#prepareDatabase('accessionTaxa.sql') ## database prepared in Downloads directory
####
####
setwd("/Users/glennd/Downloads")
#function to prepare taxID tables and return taxonomy
colrowtax <- function(x){
  lahs <- x
  colnames(lahs) <- c("OTU","taxID")
  lahs1 <- lahs
  lahs<-getTaxonomy(lahs[,2],'accessionTaxa.sql')
  rownames(lahs) <- lahs1$OTU
  return(lahs)
}

##getting list of megan results
R1all.blast.megan <- list(R1blast.d1.m, R1blast.d2.m , R1blast.d3.m, R1blast.l4.m, R1blast.l5.m, R1blast.l6.m )
names(R1all.blast.megan) <- vecnames

R2all.blast.megan <- list(R2blast.d1.m, R2blast.d2.m , R2blast.d3.m, R2blast.l4.m, R2blast.l5.m, R2blast.l6.m )
names(R2all.blast.megan) <-  vecnames

R1all.blast.megan.tax <-  lapply(R1all.blast.megan, colrowtax)
R2all.blast.megan.tax <-  lapply(R2all.blast.megan, colrowtax) 
head(R1all.blast.megan.tax$dada.nopool.nochim)

## reading in our RefLib TaxIDs to get their taxonomy to parse with Genbank Taxonomy and Megan results


View(R1SppAssList$dada.nopool.nochim)


## Changed (ref-taxIDs.csv) this with nex TaxIDs as per Bobs curation update April 8, 2024
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data")
reflbids <- read.csv("ref-taxIDs1.csv")
reflbids <-unique(reflbids)
reflbids <-reflbids[,c(1,2)]
str(reflbids)

setwd("/Users/glennd/Downloads")
reflbids1 <-  colrowtax(reflbids)
head(reflbids1)
reflbids1 <-  as.data.frame(reflbids1)
## preparing to parse RifLib and Genbank dataframes - ## No OTU (row names) numbers are matched to POTUs in reflbids1 table
str(SppAssList$dada.nopool.nochim)
reflbids1$Genus <- rownames(reflbids1)
View(reflbids1)

R1SppAssList <- lapply(R1SppAssList, function(x) merge(x, reflbids1, by = "Genus", all.x = TRUE))
R2SppAssList <- lapply(R2SppAssList, function(x) merge(x, reflbids1, by = "Genus", all.x = TRUE))

library(rlist)
setwd("/Users/glennd/Documents/GitHub/legenDNAry")
list.save(R1SppAssList, 'Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R1.rds')
list.save(R2SppAssList, 'Raw_data/Reference_library_filtering/OTU-to-RefIDs-List_R2.rds')

View(R1SppAssList$dada.pspool.nc.lulu)
View(R2SppAssList$dada.pspool.nc.lulu)

library(textshape)
R1SppAssList <- lapply(R1SppAssList, function(x) textshape::column_to_rownames(x, loc = "sequence"))
R2SppAssList <- lapply(R2SppAssList, function(x) textshape::column_to_rownames(x, loc = "sequence"))

#now just keeping relevant columns (sequences as rownames, OTU and taxonomy)
str(R1SppAssList$dada.nopool.nochim)
R1SppAssList <- lapply(R1SppAssList, function(x) x <- x[,c(3,5:11)])
R2SppAssList <- lapply(R2SppAssList, function(x) x <- x[,c(3,5:11)])

str(R2SppAssList$dada.nopool.nochim)
##making same structure to parse with genbank assigned taxonomies
## Genbanktax
R1all.blast.megan.tax <- lapply(R1all.blast.megan.tax, function(x) as.data.frame(x))
R2all.blast.megan.tax <- lapply(R2all.blast.megan.tax, function(x) as.data.frame(x))
head(R1all.blast.megan.tax$R1blast.dnp.m, n = 20)
library(tibble)
R1all.blast.megan.tax <- lapply(R1all.blast.megan.tax, function(x) tibble::rownames_to_column(x, var = "OTU"))
R2all.blast.megan.tax <- lapply(R2all.blast.megan.tax, function(x) tibble::rownames_to_column(x, var = "OTU"))
str(R2all.blast.megan.tax$dada.nopool.nochim)
## need to add sequences back to these tables megan taxonmy tables
#View(SppUnAssList$dada.nopool.nochim)
R1SppUnAssList <- lapply(R1SppUnAssList, function(x) x[,c(3,4)])
R2SppUnAssList <- lapply(R2SppUnAssList, function(x) x[,c(3,4)])
# Merge function
merge_function <- function(x, y) {
  merge(x, y, by = "OTU", all.x = TRUE)
}
# Merge lists
R1SppUnAssList <- Map(merge_function, R1all.blast.megan.tax, R1SppUnAssList)
R2SppUnAssList <- Map(merge_function, R2all.blast.megan.tax, R2SppUnAssList)
library(textshape)
R1SppUnAssList <- lapply(R1SppUnAssList, function(x) textshape::column_to_rownames(x, loc = "sequence"))
R2SppUnAssList <- lapply(R2SppUnAssList, function(x) textshape::column_to_rownames(x, loc = "sequence"))
##both lists are the same structure etc
head(R1SppAssList$dada.nopool.nochim)
head(R1SppUnAssList$dada.nopool.nochim)
## adding ID source column
R1SppAssList <- mapply(cbind, R1SppAssList, "idsource"="RefLib", SIMPLIFY=F)
R2SppAssList <- mapply(cbind, R2SppAssList, "idsource"="RefLib", SIMPLIFY=F)
R1SppUnAssList <- mapply(cbind, R1SppUnAssList, "idsource"="GenBank", SIMPLIFY=F)
R2SppUnAssList <- mapply(cbind, R2SppUnAssList, "idsource"="GenBank", SIMPLIFY=F)
## merging all taxonomy data for each iteration to all.tax
library(purrr)
## all.taxs.source contains all the taxonomy information and whether it is from a local reference or from Genbank
R1all.taxs.source <- purrr::map2(R1SppAssList,R1SppUnAssList,rbind)
R2all.taxs.source <- purrr::map2(R2SppAssList,R2SppUnAssList,rbind) 
R1all.taxs.source.otu <- lapply(R1all.taxs.source, function(x) textshape::column_to_rownames(x, loc = "OTU"))
R2all.taxs.source.otu <- lapply(R2all.taxs.source, function(x) textshape::column_to_rownames(x, loc = "OTU"))

View(all.taxs.source.otu$dada.nopool.nochim)

## all.taxs is just the relevant columns for a phyloseq object
R1all.tax <- lapply(R1all.taxs.source.otu, function(x) as.matrix(x[,c(1:7)]))
R2all.tax <- lapply(R2all.taxs.source.otu, function(x) as.matrix(x[,c(1:7)]))

View(R1all.tax$dada.nopool.nochim)
## Now all dataframes, taxonomy and sequence information can be made into phyloseq objects
View(sample_data)
row.names(sample_data) <- sample_data$newname2
## Getting sequence data in right format:
setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data")
library(Biostrings)
R1dnp_seq <- readDNAStringSet("R1dada.nopool.fasta")
R1dp_seq <- readDNAStringSet("R1dada.pool.fasta")
R1dpsp_seq <- readDNAStringSet("R1dada.pspool.fasta")
R1dlnp_seq <- readDNAStringSet("R1lulu.nopool.fasta")
R1dlp_seq <- readDNAStringSet("R1lulu.pool.fasta")
R1dlpsp_seq <- readDNAStringSet("R1lulu.pspool.fasta")

R2dnp_seq <- readDNAStringSet("R2dada.nopool.fasta")
R2dp_seq <- readDNAStringSet("R2dada.pool.fasta")
R2dpsp_seq <- readDNAStringSet("R2dada.pspool.fasta")
R2dlnp_seq <- readDNAStringSet("R2lulu.nopool.fasta")
R2dlp_seq <- readDNAStringSet("R2lulu.pool.fasta")
R2dlpsp_seq <- readDNAStringSet("R2lulu.pspool.fasta")

#test <- read.fasta("dada.nopool.fasta")
#row.names(dnp_seq) <- dnp_seq$seq.name

R1seqlist <- list(R1dnp_seq, R1dp_seq, R1dpsp_seq, R1dlnp_seq, R1dlp_seq, R1dlpsp_seq)
R2seqlist <- list(R2dnp_seq, R2dp_seq, R2dpsp_seq, R2dlnp_seq, R2dlp_seq, R2dlpsp_seq)
View(R1seqlist)
## OTU data is in 6 dataframe in R1joined, R2joined & R3 joined. Tax data is in all.tax, sample_data is sample_data 
## and sequence data is in the xx_seq files correspponding to each of the six dada/lulu iterations
library(phyloseq)
make.phylo <- function (x, y, z){ # x = OTU data to this point, y = taxonomic data, z = sample data, k = sequence data
  k1 <- Biostrings::DNAStringSet(colnames(x))
  k2 <- as.matrix(x)
  metadata(k1)$name <- paste0("OTU", seq(length(getSequences(x))))
  colnames(k2) <- paste0("OTU", seq(length(getSequences(k2))))
  test1 <- phyloseq(otu_table(k2, taxa_are_rows = FALSE), tax_table(as.matrix(y)), sample_data(z))
  dna1 <- Biostrings::DNAStringSet(colnames(x))
  ## naming sequences as OTUs
  names(dna1) <- paste0("OTU", seq(length(getSequences(x))))
  wanted <- merge_phyloseq(test1, dna1)
  return(wanted)
}

sample_data <- read.csv("sample_data_v2.csv")
rownames(sample_data) <- sample_data$newname2
R1dnp <- make.phylo(R1joined[[1]], R1all.tax[[1]], sample_data)
R1dp <- make.phylo(R1joined[[2]], R1all.tax[[2]], sample_data)
R1dpsp <- make.phylo(R1joined[[3]], R1all.tax[[3]], sample_data)
R1lnp <- make.phylo(R1joined[[4]], R1all.tax[[4]], sample_data)
R1lp <- make.phylo(R1joined[[5]], R1all.tax[[5]], sample_data)
R1lpsp <- make.phylo(R1joined[[6]], R1all.tax[[6]], sample_data)

R1phylo <- list(R1dnp, R1dp, R1dpsp, R1lnp, R1lp, R1lpsp )
names(R1phylo) <- c("R1.dada.nopool", "R1.dada.pool" , "R1.dada.pspool", "R1.lulu.nopool", "R1.lulu.pool", "R1.lulu.pspool" )

R2dnp <- make.phylo(R2joined[[1]], R2all.tax[[1]], sample_data)
R2dp <- make.phylo(R2joined[[2]], R2all.tax[[2]], sample_data)
R2dpsp <- make.phylo(R2joined[[3]], R2all.tax[[3]], sample_data)
R2lnp <- make.phylo(R2joined[[4]], R2all.tax[[4]], sample_data)
R2lp <- make.phylo(R2joined[[5]], R2all.tax[[5]], sample_data)
R2lpsp <- make.phylo(R2joined[[6]], R2all.tax[[6]], sample_data)

R2phylo <- list(R2dnp, R2dp, R2dpsp, R2lnp, R2lp, R2lpsp )
names(R2phylo) <- c("R2.dada.nopool", "R2.dada.pool" , "R2.dada.pspool", "R2.lulu.nopool", "R2.lulu.pool", "R2.lulu.pspool" )

## Now all data is grouped into the six bioinfomatic filters (dada-nopool, dada.pool, dada.pspool, lulu.nopool, lulu.pool, lulu.pspool)
## and available by OTU prevelance across PCR replicates R1 = in one or more replicates, R2 = in two or more replicates 

## These are all in the lists R1phylo, R2phylo & R3phylo and so functions can be applied across entire lists or just to one element.

## First thing to do it to tidy up the phyloseq objects by removing OTUs with no BLAST matches and those that aren't plants

R1phylo.plants <- lapply(R1phylo, function(x) subset_taxa(x, phylum =="Streptophyta"))
R2phylo.plants <- lapply(R2phylo, function(x) subset_taxa(x, phylum =="Streptophyta"))

## Now Tidying up the phyloseqs object with phyloseq_validate (removing taxa with zero reads)
library(microViz)
R1phylo.plants <- lapply(R1phylo.plants, function (x) phyloseq_validate(x, remove_undetected = TRUE))
R2phylo.plants <- lapply(R2phylo.plants, function (x) phyloseq_validate(x, remove_undetected = TRUE))


## Now fixing up taxonomy tables so that entires go all the way out to species (no NAs and each different OTU of the same Taxa is numbered)

R1phylo.plants.tax <- lapply(R1phylo.plants, function (x) tax_fix(x,
  min_length = 4,
  unknowns = c(""),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified"))

R2phylo.plants.tax <- lapply(R2phylo.plants, function (x) tax_fix(x,
                                                                  min_length = 4,
                                                                  unknowns = c(""),
                                                                  sep = " ", anon_unique = TRUE,
                                                                  suffix_rank = "classified"))

## adding another level of OTU taxonomy to distinguish between different OTUs
## <- lapply(all.blast.megan.tax, function(x) tibble::rownames_to_column(x, var = "OTU"))
remotes::install_github("mikemc/speedyseq")
library(speedyseq)

R1phylo.plants.tax <- lapply(R1phylo.plants.tax, function (x) {
  mutate_tax_table(x, speciesOTU = paste0(species, sep = '_', seq(nrow(tax_table(x)))))
})
R2phylo.plants.tax <- lapply(R2phylo.plants.tax, function (x) {
  mutate_tax_table(x, speciesOTU = paste0(species, sep = '_', seq(nrow(tax_table(x)))))
})

## We can now checkout library sizes to see if there are any samples that need to go...
tab <- otu_table(R1phylo.plants.tax$R1.dada.nopool)
tab <- rowSums(tab)
df <- sort(tab)
df
head(as.data.frame(df))
boxplot(df)
p <- ggplot(as.data.frame(df), aes(x = 1, y = df)) + 
  geom_violin(trim=FALSE)
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)

## We can drop samples with >10000 sequences in total (Probably 20,000 is better)
View(sample_data(R1phylo.plants.tax$R1.dada.nopool))
R1phylo.plants.tax <- lapply(R1phylo.plants.tax, function (x) subset_samples(x, newname2 != "C24_P2"))
R1phylo.plants.tax <- lapply(R1phylo.plants.tax, function (x) subset_samples(x, newname2 != "C183_P1"))
R2phylo.plants.tax <- lapply(R2phylo.plants.tax, function (x) subset_samples(x, newname2 != "C24_P2"))
R2phylo.plants.tax <- lapply(R2phylo.plants.tax, function (x) subset_samples(x, newname2 != "C183_P1"))

#################################################################################################################
##################################################################################################
###########################################################################################
### Note - here R environment saved as "CescPRdata.RData"--
#save.image("CescPRdata_v2.RData")
# setwd("/Users/glennd/Documents/Cesc-PR-eDNA/Soil_eDNA_data")

# If loading from here all libraries used above need reloading...
library(dada2)
library(textshape)
library(tibble)
library(purrr)
library(Biostrings)
library(phyloseq)
library(microViz)
library(speedyseq)
library(phylosmith)
library(vegan)
library(metagMisc)
library(ggplot2)
library(ggExtra)
library(ggplotify)
library(ggpubr)
library(iNEXT)

# Load data from above
load("Raw_data/CescPRdata_v2.RData") # - remove first # to load environment.. 

## Looking at rarefaction curves to assess sequencing coverage after removing small libraries > 10000 sequences
tab <- R1phylo.plants.tax[[5]]
ttab <- otu_table(tab)
class(ttab) <- "matrix" 
par(mfrow =c(2,2))
raremin <- min(rowSums(ttab))
rarecurve(ttab[1:25,], step = 100, sample = raremin, col = "blue", label = FALSE)
rarecurve(ttab[26:51,], step = 100, sample = raremin, col = "blue", label = FALSE)
rarecurve(ttab[52:77,], step = 100, sample = raremin, col = "blue", label = FALSE)
rarecurve(ttab[76:101,], step = 100, sample = raremin, col = "blue", label = FALSE)


## Using R1 (OTU found in 1/3 or more replicate PCRs), 
## with lulu filtering and pooling.. so this is R1phylo.plants.tax[[5]]
## if you View(R1phylo.plants.tax) you can see each of the bioinformatic iterations and,
## R1 = OTU found in 1/3 or more replicate PCRs
## R2 = OTU found in 2/3 or more replicate PCRs
## R3 = OTU found in all (3/3) or more replicate PCRs
View(R1phylo.plants.tax)
#########################################
## pulling out the DILUTION EXPERIMENT
tab <- R1phylo.plants.tax[[5]]
tabinc <- subset_samples(tab, dilution == "yes")
tabinc <- phyloseq_validate(tabinc, remove_undetected = TRUE)
tabinc
## reordering samples to make plots more informative
remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)
samorder <-sample.data.frame(tabinc)
samorder <- arrange(samorder, realsample, paircombo)

n <- rownames(samorder)

##  Visualising paired-dilution relative compositions
a<- tabinc %>%
  comp_barplot(
    tax_level = "speciesOTU",
    sample_order = n ,
    #label = "realsample", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    #taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other OTUs", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) +
  coord_flip()
a
## Seems a Consistent difference in th incubation experiment samples only...
## lets see if the number of taxa detected differs between treatments

otucounts <- phyloseq_ntaxa_by_tax(tabinc, TaxRank = "speciesOTU")
otucounts1 <- phyloseq_ntaxa_by_tax(tabinc, TaxRank = "phylum")
otusummary <- otucounts %>%
  group_by(paircombo) %>%
  dplyr::summarise(
    OTUs = sum(N.OTU),
    .groups = "drop"
  )
dilu.richness <- merge(otucounts1, otusummary, by = "paircombo", all.x = TRUE)
View(dilu.richness)
## in order to plot paired samples, removing second repeat of leS2A3 (ddiltue) - it has an equal number of OTUs (22) anyway
dilu.richness <- dilu.richness[-c(31), ]
library(ggpubr)
library(dplyr)
library(tidyr)
df <- dilu.richness
df1 <- df %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = OTUs)
df <- mdilu
df2 <- df %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = fab)

ggpaired(df1, cond1 = "clean", cond2 = "ddilute",
         fill = "condition", palette = "npg", line.color = "gray")
colnames(df2) <- c("realsample", "clean.lib", "ddilute.lib")
## important to remember that this is all confounded by library size -- we can add that to this paired table.. and normalize OTUs according to it. 
fab <- otu_table(tabinc)
fab <- rowSums(fab)
fab <- sort(fab)
mdilu <- merge(dilu.richness, data.frame(fab, newname2 = names(fab)), by = "newname2", all.x = TRUE)
str(mdilu)
richlib <- merge(df1, df2, by = "realsample", all.x = TRUE)
## normalizing OTUs by librarysize by OTUs/librarysize - OTUs per 10,000 reads
richlib$cleanNO <- 10000*(richlib$clean/richlib$clean.lib)
richlib$diluteNO <- 10000*(richlib$ddilute/richlib$ddilute.lib)

## boxplot on the normalized OTUs
ggpaired(richlib, cond1 = "cleanNO", cond2 = "diluteNO",
         fill = "condition", palette = "npg", line.color = "gray", legend.title="", tickslab= FALSE) + 
  theme(axis.text.x=element_blank())+
  xlab("Extraction Protocol") +
  ylab("OTU's per 10000 reads") +
  scale_fill_hue(labels = c("Cleaned", "Phenol.Dilu")) 

## seeing if we can T-test this difference -- normality in groups first 
shapiro.test(richlib$cleanNO)
shapiro.test(richlib$diluteNO)
## looks like both are not normal ( P < 0.05 = reject the null hypothesis that the data are normally distributed)
## we can use non-parametric tests
wilcox.test(richlib$cleanNO, richlib$diluteNO, paired = TRUE, alternative = "two.sided")

## We can also check this with a diversity index which uses the distribution of OTUs across libraries and sample sizes within libraies to calculate it
## (The shannon index) - we use shannon here because chao is really sensitive to the lack of singletons (a problem with sequence data)
plot_richness(tabinc, x="DNAtreat",  measures=c( "Shannon"))
plot_richness(tabinc, x="DNAtreat", measures="Shannon", color = "DNAtreat")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

alpha.diversity <- estimate_richness(tabinc, measures=c("Shannon"))
chaodilu <- merge(dilu.richness, data.frame(alpha.diversity, newname2 = rownames(alpha.diversity)), by = "newname2", all.x = TRUE)
df3 <- chaodilu %>% 
  pivot_wider(id_cols = realsample, names_from = DNAtreat, values_from = Shannon)
colnames(df3) <- c("realsample", "clean.shan", "dilute.shan")
richlib <- merge(richlib, df3, by = "realsample", all.x = TRUE)
richlib

shapiro.test(richlib$clean.shan)
shapiro.test(richlib$dilute.shan)
## data are normal with shannon indices so we can use a paired sample t test

t.test(richlib$clean.shan, richlib$dilute.shan, paired = TRUE, alternative = "two.sided")
## plot of differences
ggpaired(richlib, cond1 = "clean.shan", cond2 = "dilute.shan",
         fill = "condition", palette = "npg", line.color = "gray", legend.title="", tickslab= FALSE) + 
  theme(axis.text.x=element_blank())+
  xlab("Extraction Protocol") +
  ylab("Shannon Diversity Index") +
  scale_fill_hue(labels = c("Cleaned", "Phenol.Dilu")) 

## ok so with the raw normalised data we have no evidence of a difference in alpha diversity (meaning cleaning the DNA made no difference)..but 
## tests with shannon index indicate significantly greater A diversity with the unclean DNA. FML...
## so buying those expensive kits and running through an entire column extract protocol was a waste of time... in fact counter productive
## excellent ....

###################################################################################
###################################################################
################################################
## taking a quick look at the single samples Vs the mixes from those experiments
## i.e. where we took DNA samples from the individual soil samples before we made the metasample for extraction
tab <- lenient
tabinc <- subset_samples(tab, homo != "")
tabinc <- phyloseq_validate(tabinc, remove_undetected = TRUE)
tabinc
## reordering samples to make plots more informative
remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)
samorder <-sample.data.frame(tabinc)
library(dplyr)
samorder <- arrange(samorder, realsample1, homo)

n <- rownames(samorder)
n
##  Visualising paired-dilution relative compositions (top 20 taxa by total seq count)
a<- tabinc %>%
  comp_barplot(
    tax_level = "speciesOTU",
    sample_order = n ,
    #label = "realsample", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    #taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other OTUs", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) +
  coord_flip()
a
### Even from this we can see that the three individual soil samples we took per station are completely different before combined into a metasample
## we can visualize this with an ordination appropriate for compositional data (PCA with CLR tranformed data)
tabinc %>%
  tax_transform("clr", rank = "speciesOTU") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "homo", shape = "realsample1", size = 3) +
  scale_colour_brewer(palette = "Dark2")
## pretty clear here as well there is substantial varition betwen subsamples and without creating metasamples there might be sampling atrefacts
## (look at how different the ID1 samples are and how 2 ID1 samples are more similar to C8 samples then any other ID1 sample....)
## lets look at OTU overlap between sub-samples and metasamples
## to be fair there really should be some sort of rarefying here to account for uneven sample sizes.....
library(BiocManager)
#BiocManager::install("microbiome")
library(microbiome)
#remotes::install_github("Russel88/MicEco")
library(MicEco)
## one for each sample
c12 <- subset_samples(tab, realsample1 == "C12")
c12 <- subset_samples(c12, DNAtreat == "clean")
c12<- phyloseq_validate(c12, remove_undetected = TRUE)
rarfun  <- function(x) {
  rfy <- min(rowSums(otu_table(x)))
  sar.rats.rarefy <- rarefy_even_depth(x, sample.size=rfy, replace=FALSE, rngseed = 1)
  return(sar.rats.rarefy)
}
c12r <- rarfun(c12)
ps_venn(
  c12, 
  group = "realsample"
)
c12r <- rarfun(c12)
ps_venn(
  c12r, 
  group = "realsample"
)
c18 <- subset_samples(tab, realsample1 == "C18")
c18 <- phyloseq_validate(c18, remove_undetected = TRUE)
c18r <- rarfun(c18)
ps_venn(
  c18r, 
  group = "realsample"
)
c8 <- subset_samples(tab, realsample1 == "C8")
c8 <- subset_samples(c8, DNAtreat == "clean")
c8 <- phyloseq_validate(c8, remove_undetected = TRUE)
c8r <- rarfun(c8)

ps_venn(
  c8r, 
  group = "realsample"
)
ID1 <- subset_samples(tab, realsample1 == "ID1")
ps_venn(
  ID1, 
  group = "realsample"
)

#### taking a quick look at the same with the small grid samples...
tab <- R1phylo.plants.tax[[5]]
tabinc <- subset_samples(tab, experiment == "bigplot/smallss")
tabinc <- phyloseq_validate(tabinc, remove_undetected = TRUE)
tabinc
rfy <- min(rowSums(otu_table(tabinc)))
Gr <- rarfun(tabinc)
heatmap(otu_table(Gr))
plot_heatmap(Gr)
Gr %>%
  tax_transform("clr", rank = "species") %>%
  comp_heatmap(colors = heat_palette(sym = TRUE), name = "CLR", sample_names_show = TRUE)
Gr %>%
  tax_transform("compositional", rank = "speciesOTU") %>%
  comp_heatmap(sample_names_show = TRUE)
Grfilt %>%
  tax_transform("compositional", rank = "speciesOTU") %>%
  comp_heatmap(sample_names_show = TRUE)
tabincfilt %>%
  tax_transform("compositional", rank = "species") %>%
  comp_heatmap(sample_names_show = TRUE)
Grfilt <- Gr %>%
  tax_filter(min_prevalence = 5) %>%
  tax_fix()
rowSums(otu_table(tabincfilt))
k1 <- phyloseq_ntaxa_by_tax(
tabinc,
  TaxRank = "speciesOTU",
  relative = F,
  add_meta_data = F
)
table(k$Sample, k$N.OTU)
table(k1$Sample, k1$N.OTU)
## reordering samples to make plots more informative
remotes::install_github("kstagaman/phyloseqCompanion")
library(phyloseqCompanion)
samorder <-sample.data.frame(tabinc)
samorder <- arrange(samorder, realsample)

n <- rownames(samorder)
n
##  Visualising paired-dilution relative compositions (top 20 taxa by total seq count)
a<- tabinc %>%
  comp_barplot(
    tax_level = "speciesOTU",
    sample_order = n ,
    #label = "realsample", # name an alternative variable to label axis
    n_taxa = 20, # give more taxa unique colours
    #taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "), # remove underscores
    other_name = "Other OTUs", # set custom name for the "other" category
    merge_other = FALSE, # split the "Other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) +
  coord_flip()
a
### holy hell these are all crazy different in composition as well (remember this was a 2x2m grid)
## lets look at a venn diagram with just the 10cm grid (the G5 samples)
g5 <- subset_samples(tabinc, realsample1 == "G5")
ps_venn(
  g5, 
  group = "realsample"
)
## remarkable... 2/3 sample has => unique taxa as it has shared taxa between any other sample..
## lets look at the small grid V medium grid samples
tabinc  %>%
  tax_transform("clr", rank = "speciesOTU") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "factorlevel", size = 4) +
  scale_colour_brewer(palette = "Dark2")
## the small grid samples are barely the most similar to each other in multivariate space compared to the larger grid samples..
## we cn also see how grid size stacks up across the entire sample set
tabb <- subset_samples(tab, gridsize != "")
tabb
tabb <- phyloseq_validate(tabb, remove_undetected = TRUE)
tabb  %>%
  tax_transform("clr", rank = "speciesOTU") %>%
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>%
  ord_plot(color = "gridname", shape = "gridsize", size = 3) +
  stat_ellipse(aes(colour = gridname)) +
  scale_colour_brewer(palette = "Dark2")

######### Some visualizations of the whole dataset

remotes::install_github("schuyler-smith/phylosmith")
remotes::install_github("schuyler-smith/phyloschuyler")
library(phylosmith)
library(phyloseq)
library(vegan)
## note this function - here you can get "core taxa".. for these parameters all taxa in >1% of samples (frequency), over a relative abundance of 0.001 per sample
k2 <-phylosmith::taxa_core(tab, treatment = NULL, subset = NULL,frequency = 0.01, abundance_threshold = 0.001)

## check out metacoder for full dataset visualizations
library(metacoder)
metpr <- parse_phyloseq(tab)
metpr$data$tax_abund <- calc_taxon_abund(metpr, "otu_table")
metpr$data$tax_data <- calc_obs_props(metpr, "otu_table")
metpr$data$tax_abund <- calc_taxon_abund(metpr, "tax_data")

metpr$data$tax_occ <- calc_n_samples(metpr, "tax_abund")

set.seed(1) # This makes the plot appear the same each time it is run 
heat_tree(metpr, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford")


