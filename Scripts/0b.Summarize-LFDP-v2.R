### SUMMARIZE LFDP DATA ROUND 2

### This script computes total abundance and basal area per species (of live stems) for the full plot and in the two main land-use history parts of the plot.

# # 0. LOAD THE LFDP CENSUS DATA AS 'census'
library(EDIutils)
raw6 <- read_data_entity(packageId = "knb-lter-luq.119.1545979", entityId = "325c43057e0dd4e1cd6a13fa5125a76d")
census <- readr::read_csv(file = raw6)

# 1. LOAD THE LFDP CENSUS DATA AS 'census'
head(census)
census$Mnemonic <- as.factor(census$Mnemonic)

# 2. DROP DEAD STEMS
census <- census[census$Status %in% "alive",]

# 3. COMPUTE BASAL AREA OF EACH STEM
census$BA <- pi * (census$DBH/2000)^2

# 4. DOWNLOAD THE PHYSICAL ENVIRONMENT DATA
env <- readr::read_csv(EDIutils::read_data_entity(packageId = "knb-lter-luq.47.381050",
                                                  entityId = "21c7cce729a96de675aaa9281526eb4a"))
# 5. ADD LAND USE CATEGORY TO CENSUS DATA
census$LU <- ifelse(env$CoverClass==4, 'h', 'l')[match(census$Quadrat, env$Quadrat)]

# 6. COMPUTE TOTAL ABUNDANC AND BASAL AREA IN WHOLE PLOT AND LAND-USE HISTORY
df <- data.frame(spcode=names(table(census$Mnemonic)),
                 total_abund=as.numeric(table(census$Mnemonic)),
                 low_LU_abund=as.numeric(table(census$Mnemonic[census$LU=='l'])),
                 high_LU_abund=as.numeric(table(census$Mnemonic[census$LU=='h'])),
                 total_ba=as.numeric(tapply(census$BA, census$Mnemonic, sum, na.rm=T)),
                 low_ba=as.numeric(tapply(census$BA[census$LU=='l'], 
                                          census$Mnemonic[census$LU=='l'], sum, na.rm=T)),
                 high_ba=as.numeric(tapply(census$BA[census$LU=='h'], 
                                           census$Mnemonic[census$LU=='h'], sum, na.rm=T)))


# 7. SAVE RESULTS AND SEND OUTPUT TO: <robert.muscarella@ebc.uu.se>
# saveRDS(df, "Raw_data/LFDP2016-extract-v2-20240427.RDA")
# saveRDS(df, "Raw_data/LFDP2023-extract-v2-20240427.RDA")


