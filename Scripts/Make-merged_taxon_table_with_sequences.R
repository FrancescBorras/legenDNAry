### GENERATE MERGED TABLE WITH POTUS, SEQUENCES, GOTUS, ETC.

library(googlesheets4)
gs4_deauth()

# Read in Cesc's DNA table
leaf.otus <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1Ihe1dXkENK4zZ0lBsN1_1tDnaIvVD3u9oiNHjtMRHzo/edit?usp=sharing", sheet="Leaf OTUs")

# Read in the species list from the LFDP census data
codes <- readxl::read_xlsx("Raw_data/LFDP-SPcodes.xlsx")

# Which codes from Cesc's table do not match the LFDP codes?
leaf.otus$mnemonic[!leaf.otus$mnemonic %in% codes$`SPECIES CODE`]

# Change some of these that should be included...
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "FERN1"] <- "CYAARB"
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "GENAME1"] <- "GENAME"
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "MANIND1"] <- "MANIND"
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "CALANT"] <- "CALCAL"
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "MYRLEA"] <- "MYRLEP"
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "BANANA"] <- "HELCAR"
leaf.otus$mnemonic[leaf.otus$mnemonic %in% "BEIPEN alpha1"] <- "BEIPEN"

# Exclude species in LFDP table that were not included in the 2023 census
codes <- codes[codes$LFDP2023==1,]

# Which codes from LFDP are not in Cesc's table?
sort(codes$`SPECIES CODE`[!codes$`SPECIES CODE` %in% leaf.otus$mnemonic])

# Drop unknown species from census data
codes <- codes[!codes$`SPECIES CODE` %in% "UNKSPP",]

# Merge the two tables based on the LFDP census table
compiled <- cbind(codes, leaf.otus[match(codes$`SPECIES CODE`, leaf.otus$mnemonic),])

# Fill collected column
compiled$collected[is.na(compiled$collected)] <- "no"

# Drop some irrelevant columns
compiled <- compiled[,!colnames(compiled) %in% c("LFDP2023","mnemonic","scientific name","To sort")]

# Save data
writexl::write_xlsx(compiled, 
                    "Processed_data/merged_taxon_table_with_sequences-20241122.xlsx")







### READ BACK IN THIS COMPILED DATA FOR SOME SUMMARY NUMBERS
compiled <- readxl::read_xlsx("Processed_data/merged_taxon_table_with_sequences-20241122.xlsx")

# How many species were sequenced for ref library?
length(unique(compiled$`SPECIES CODE`[!is.na(compiled$highestOTU)])) # 107

# How many total tree species in the 2023 census?
length(sort(unique(compiled$`SPECIES CODE`))) # 141

# How many OTUs (gOTUs in the end) are there in the final eDNA dataset?
length(unique(compiled$gOTU[!is.na(compiled$gOTU)])) # 78

# How many OTUs (POTUs in the end) are there in the final eDNA dataset?
# length(unique(compiled$POTU[!is.na(compiled$POTU)])) # 82

# (Note: how many unique POTUs were there that remained in the final dataset?)
# length(unique(compiled$POTU[!is.na(compiled$gOTU)])) # 83

# How many gOTUs are there in the final dataset with unique haplotypes?
sum(table(compiled$gOTU[!is.na(compiled$gOTU) & !is.na(compiled$highestOTU)])==1) # 62

sum(table(compiled$gOTU[!is.na(compiled$gOTU) & !is.na(compiled$highestOTU)])>1) # 62

# How many gOTUs are there in the final dataset without unique haplotypes?
107 - 62 # 45

sum(table(compiled$gOTU)==1) # 56


compiled$gOTU

x <- table(compiled$gOTU)
x[x==1] <- 0
sum(x)
length(x[x>0])

141 - (69 + 56)
125 - (69 + 56)
125+16

125-107

compiled$`SPECIES CODE`[!compiled$`SPECIES CODE` %in% tree23$`SPECIES CODE`]
tree23$`SPECIES CODE`[!tree23$`SPECIES CODE` %in% compiled$`SPECIES CODE`]



sum(is.na(compiled$gOTU))
sum(!is.na(compiled$gOTU))



# How many gOTUs are there in the final dataset with unique haplotypes?
sum(table(compiled$gOTU[!is.na(compiled$gOTU)])==1) # 56


# What % of total final gOTUs are unique haplotypes?
56/107 # 52%

# What % of total final gOTUs are not unique haplotypes?
(107 - 56)/107 # 48%








