install.packages("mobsim")

library(mobsim)
comm_rand <- sim_poisson_community(s_pool = 30, n_sim = 300)
comm_agg <- sim_thomas_community(s_pool = 30, n_sim = 300, sigma = 0.05, mother_points = 1)

plot(comm_rand)

par(mfrow = c(1,2))
samples_rand <- sample_quadrats(comm_rand, avoid_overlap = TRUE)
samples_agg <- sample_quadrats(comm_agg, avoid_overlap = TRUE)



# 1. Simluate community
    # - Different spatial aggregations
    # - Different SAD (how many rare species?)

# 2. Sample community
    
# 3. Rarefy sample to reflect detectability in eDNA
    # - Directly proportional to abundance
    # - Random

# 4. Compute diversity / composition comparisons between sample and "True" data


1. Simluate community
- Different spatial aggregations
- Different SAD (how many rare species?)
2. Sample community
3. Rarefy sample to reflect detectability in eDNA
- Directly proportional to abundance
- Random
4. Compute diversity / composition comparisons between sample and "True" data

What is correspondence between e.g. species richness of stem and dna sample when you have stems reflected in DNA as proportional to their actual abundance but you have different levels of rarity.