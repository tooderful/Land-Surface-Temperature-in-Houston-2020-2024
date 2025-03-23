packs <- c('tidyverse', 'sf', 'spdep', 'spData', 'sfdep', 'viridis', 'gridExtra', 'spatialreg')
lapply(packs, require, character.only = TRUE)

# read in the census block shapefile
tx.sf <- st_read('harrisCountyProcessedData/blocks/Blocks.shp', quiet=TRUE)
st_crs(tx.sf)
tx.sf <- st_transform(tx.sf, crs = 6588) # EPSG code for NAD1983(2011) Texas South Central (ftUS)

# -------------------- HARRIS W/ REDLINING -------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
harrisR_data = read_csv('harrisDataFR.csv', show_col_types = FALSE)
harrisR_data$Landscape_Shape_Index[is.na(harrisR_data$Landscape_Shape_Index)] <- 0
harrisR_data <- harrisR_data %>%
  mutate(CTBKEY = as.character(CTBKEY))
# Join the data frames based on the CTBKEY column
harrisR_sf <- tx.sf %>%
  inner_join(harrisR_data, by = "CTBKEY")
plot(harrisR_sf)

# Defining Spatial Neighbors
harrisR.nb <- poly2nb(harrisR_sf, queen=TRUE) # Queen continuity
### Want row standardized first-order, adjacency matrix
# # but style='W' does not allow us to get a symmetric matrix
harrisR.rsw <- nb2mat(harrisR.nb, style='B', zero.policy = TRUE) # Row standardization

# creates a distance matrix between all census blocks
dist.mat <- st_distance(harrisR_sf)
# determines which areas have no neighbors
no.neighbors <- which(rowSums(harrisR.rsw) == 0) # 412

# determine which areas are closest to these areas
nearest.area <- rep(NA, length(no.neighbors))
for(i in 1:length(no.neighbors)){
  distances <- dist.mat[no.neighbors[i],]
  nearest.area[i] <- which(distances == distances[order(distances)][2])
  
  # Update W
  harrisR.rsw[no.neighbors[i], nearest.area[i]] <- 1
  harrisR.rsw[nearest.area[i], no.neighbors[i]] <- 1   # Keep W symmetric
}

colnames(harrisR.rsw) <- rownames(harrisR.rsw)
isSymmetric.matrix(harrisR.rsw)

harrisR.rsw <- mat2listw(harrisR.rsw, style='B', zero.policy=TRUE)
# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran <- moran.test(harrisR_data$LST, harrisR.rsw, 
                     alternative = 'greater') #  Test for clustering --> the default
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 0.704355166, while the expected value is -0.001046025. Because the Moran's I value 
# is larger than the E(I), there is evidence of some significant positive spatial 
# auto-correlation and clustering.

# Spatial Regression
sar <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + FctImp + NDVI +    
                  PerCover_Forest + grade,listw=harrisR.rsw,data=harrisR_sf) # NOT OK --> system is computationally singular
sar8 <- spautolm(LST ~ FctImp + grade,listw=harrisR.rsw,data=harrisR_sf) # one continuous and grade OK
# FctImp is 21.0375874
# gradeB is -408.7690590
# gradeC is -446.8715241
# gradeD is -44.4651061
sar9 <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + FctImp + NDVI +    
                  PerCover_Forest,listw=harrisR.rsw,data=harrisR_sf) # NOT OK - no grade --> system is computationally singular
sar10 <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + NDVI +    
                  PerCover_Forest + grade,listw=harrisR.rsw,data=harrisR_sf) # NOT OK - no FctImp --> system is computationally singular
sar11 <- spautolm(LST ~ NDVI + Landscape_Shape_Index + FctImp + grade,listw=harrisR.rsw,data=harrisR_sf) # OK - no Patch_Density and PerCover_Forest
# NDVI is -1.101465e-01
# LSI is -2.615683e-03
# FctImp is 1.784450e+01
# gradeB is -3.996837e+02
# gradeC is -4.710507e+02
# gradeD is -7.068902e+01
sar12 <- spautolm(LST ~ Landscape_Shape_Index + FctImp + NDVI +    
                  PerCover_Forest + grade,listw=harrisR.rsw,data=harrisR_sf) # OK - no Patch_Density
# LSI is -2.530612e-03
# FctImp is 1.758736e+01
# NDVI is -1.054684e-01
# PerCover_Forest is -2.307867e+03
# gradeB is -4.119354e+02
# gradeC is -4.816313e+02
# gradeD is -8.136159e+01
sar13 <- spautolm(LST ~ Patch_Density + FctImp + NDVI +    
                  PerCover_Forest + grade,listw=harrisR.rsw,data=harrisR_sf) # NOT OK - no LSI
sar14 <- spautolm(LST ~ Landscape_Shape_Index + FctImp + NDVI + PerCover_Forest,
                  listw=harrisR.rsw,data=harrisR_sf) # OK - no PD and no grade
# LSI is -2.474076e-03
# FctImp is 1.755700e+01
# NDVI is -1.024587e-01
# PerCover_Forest is -1.964537e+03

# -------------------- HARRIS W/O REDLINING ------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
harris_data = read_csv('harrisDataF.csv', show_col_types = FALSE)
harris_data$Landscape_Shape_Index[is.na(harris_data$Landscape_Shape_Index)] <- 0
harris_data <- harris_data %>%
  mutate(CTBKEY = as.character(CTBKEY))
# Join the data frames based on the CTBKEY column
harris_sf <- tx.sf %>%
  inner_join(harris_data, by = "CTBKEY") # 14371 census blocks

# Defining Spatial Neighbors
harris.nb <- poly2nb(harris_sf, queen=TRUE) # Queen continuity
### Want row standardized first-order, adjacency matrix
# # but style='W' does not allow us to get a symmetric matrix
harris.rsw <- nb2mat(harris.nb, style='B', zero.policy = TRUE) # Row standardization

# creates a distance matrix between all census blocks
dist.mat2 <- st_distance(harris_sf)
# determines which areas have no neighbors
no.neighbors2 <- which(rowSums(harris.rsw) == 0) # 1388

# determine which areas are closest to these areas
nearest.area2 <- rep(NA, length(no.neighbors2))
for(i in 1:length(no.neighbors2)){
  distances2 <- dist.mat2[no.neighbors2[i],]
  nearest.area2[i] <- which(distances2 == distances2[order(distances2)][2])
  
  # Update W
  harris.rsw[no.neighbors2[i], nearest.area2[i]] <- 1
  harris.rsw[nearest.area2[i], no.neighbors2[i]] <- 1   # Keep W symmetric
}

colnames(harris.rsw) <- rownames(harris.rsw)
isSymmetric.matrix(harris.rsw)

harris.rsw <- mat2listw(harris.rsw, style='B', zero.policy=TRUE)

# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran2 <- moran.test(harris_data$LST, harris.rsw, 
                     alternative = 'greater') #  Test for clustering --> the default
gmoran2
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 9.403676e-01, while the expected value is -6.958942e-05. Because the Moran's I value 
# is significantly larger than the E(I), there is evidence of some significant 
# positive spatial auto-correlation and clustering.

# Spatial Regression
sar <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + FctImp + NDVI +    
                  PerCover_Forest,listw=harris.rsw,data=harris_sf) # NOT OK --> system is computationally singular
sar2 <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + NDVI +    
                  PerCover_Forest,listw=harris.rsw,data=harris_sf) # NOT OK - no FctImp --> system is computationally singular
sar3 <- spautolm(LST ~ NDVI + Landscape_Shape_Index + FctImp,listw=harris.rsw,data=harris_sf) # OK - no Patch Density and PerCoverForest
# NDVI is -3.340581e-01
# LSI is 3.668341e-05
# FctImp is 2.075663e+01
sar4 <- spautolm(LST ~ Landscape_Shape_Index + FctImp + NDVI +    
                  PerCover_Forest,listw=harris.rsw,data=harris_sf) # OK - no patch density
# LSI is 2.682187e-05
# FctImp is 1.950768e+01
# NDVI is -3.310703e-01
# PerCover_Forest is -1.441637e+03
sar5 <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + FctImp + NDVI,listw=harris.rsw,data=harris_sf) # NOT OK - no PerCoverForest --> system is computationally singular
sar6 <- spautolm(LST ~ Patch_Density + FctImp + NDVI +    
                  PerCover_Forest,listw=harris.rsw,data=harris_sf) # NOT OK - no LSI --> system is computationally singular
sar7 <- spautolm(LST ~ Patch_Density + Landscape_Shape_Index + FctImp +    
                  PerCover_Forest,listw=harris.rsw,data=harris_sf) # NOT OK - no NDVI --> system is computationally singular
