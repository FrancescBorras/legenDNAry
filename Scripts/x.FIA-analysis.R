





### Compare eDNA with FIA data...


library(rFIA)

# Get FIA data
# Somehow condense to same gOTUs...
# Compare composition of FIA plots with DNA data from LFDP
# Map how similarity of the DNA sample declines across the island


# Download FIA data for Puerto Rico 
fia <- getFIA(states='PR', dir='Raw_data/prFIA/')

# Load data
prfia <- readFIA("Raw_data/prFIA") #fix the download bj
head(prfia)

lapply(prfia, head)
prfia$PLOT

prfia$TREE


#Only estimates for the most recent inventory year
MR_PR_Trees <- clipFIA(PR_Trees)#mostRecent = TRUE) #subset of the most recent dara (MR)
MR_PR_Trees_tpa <- tpa(MR_PR_Trees)
head(MR_PR_Trees_tpa)

#All Inventory Years Available
PR_Trees_tpa<-tpa(PR_Trees)
head(PR_Trees_tpa)

# Gorup data by species
PR_Trees_sp <- tpa(PR_Trees , bySpecies = TRUE)
PR_Trees_sp_total <- tpa(MR_PR_Trees , bySpecies = TRUE, totals = TRUE, treeDomain = TRUE, areaDomain = TRUE)
head(PR_Trees_sp)
View(PR_Trees_sp)

#Abundnce Data are saved in Results/FIA
FIA_abbundance <- write.csv(PR_Trees_sp, "Data/FIA/FIA_abundnace.csv")





