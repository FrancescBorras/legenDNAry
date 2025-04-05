############ 
### MAKE SPECIES TABLE FOR PUBLICATION
############
library(readxl)
library(googlesheets4)

### Read abundance data (LFDP 2023 summary)
abun <- readRDS("Raw_data/LFDP2023-extract-v2-20240427.RDA")
# Drop a 'species' with blank name
abun <- abun[abun$spcode != "",]

### Read the taxonomic data sheet 
sp <- read_excel("Raw_data/LFDP-SPcodes.xlsx")

# Filter taxonomic data sheet to species in the 2023 LFDP 
sp <- sp[sp$`SPECIES CODE` %in% abun$spcode,]

### Read the clean sequence sheet
url <- "https://docs.google.com/spreadsheets/d/1Ihe1dXkENK4zZ0lBsN1_1tDnaIvVD3u9oiNHjtMRHzo/edit?gid=1431127478#gid=1431127478"
gs4_deauth()
otu <- read_sheet(url, sheet="Leaf OTUs cleaned")

### Merge the abundances to the species codes
sp$LFDP2023_Abundance <- abun$total_abund[match(sp$`SPECIES CODE`, abun$spcode)]
sp$LFDP2023_BasalArea <- abun$total_ba[match(sp$`SPECIES CODE`, abun$spcode)]

### Fix some different species codes
otu$mnemonic[otu$mnemonic=="BEIPEN alpha1"] <- "BEIPEN"
otu$mnemonic[otu$mnemonic=="CALANT"] <- "CALCAL"
otu$mnemonic[otu$mnemonic=="MANIND1"] <- "MANIND"
otu$mnemonic[otu$mnemonic=="GENAME1"] <- "GENAME"
otu$mnemonic[otu$mnemonic=="MYRLEA"] <- "MYRLEP"
otu$mnemonic[otu$mnemonic=="FERN1"] <- "CYAARB"

otu$mnemonic[otu$mnemonic=="FERN1"] <- "CYAARB"


# See which still don't match up
otu$mnemonic[!otu$mnemonic %in% sp$`SPECIES CODE`]
sp$`SPECIES CODE`[!sp$`SPECIES CODE` %in% otu$mnemonic]

### Merge the sequences to the species codes
sp$seq <- otu$highestOTU[match(sp$`SPECIES CODE`, otu$mnemonic)]

### Drop/Edit some columns
sp$LFDP2023 <- NULL
sp$`Life form` <- NULL
sp$`BLAST Comments` <- NULL
sp$gOTU_2 <- NULL
names(sp)[which(names(sp)=="SPECIES CODE")] <- "SPCODE"
names(sp)[which(names(sp)=="POTU")] <- "ASV"
names(sp)[which(names(sp)=="gOTU")] <- "OTU"
names(sp)[which(names(sp)=="seq")] <- "Sequence"

### Change 0 basal area herbaciouos species to NA (DBH is not recorded)
sp$total_ba[sp$total_ba==0] <- NA

### Save
writexl::write_xlsx(sp, "/Users/au529793/Desktop/PR-Taxonomy_Table.xlsx")


sp <- read_xlsx("/Users/au529793/Desktop/PR-Taxonomy_Table.xlsx")



length(unique(sp$OTU[!is.na(sp$OTU)]))

mean(!is.na(sp$OTU))

sum(sp$total_abund[!is.na(sp$OTU)])/sum(sp$total_abund)
sum(sp$total_abund[!is.na(sp$OTU)])/sum(sp$total_abund)






