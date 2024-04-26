### SUMMARIZE LFDP DATA ROUND 2

# Get total abundance and basal area per species at the full plot and in LU clusters


library(EDIutils)
raw6 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "325c43057e0dd4e1cd6a13fa5125a76d")
census <- readr::read_csv(file = raw6)

# 0. LOAD THE LFDP CENSUS DATA AS 'census'
head(census)
census$Mnemonic <- as.factor(census$Mnemonic)

# 1. COMPUTE BASAL AREA OF EACH STEM
census$BA <- pi * (census$DBH/2000)^2

table(census$Mnemonic)
tapply(census$BA, census$Mnemonic, sum, )

head(census$Status)


### NEED TO ADD LAND USE PART...
