packs <- c('tidyverse', 'sf', 'spdep', 'spData', 'sfdep', 'viridis', 'gridExtra', 'spatialreg')
lapply(packs, require, character.only = TRUE)

# read in the census block shapefile
la.sf <- st_read('laCountyProcessedData/2020_Census_Blocks/2020_Census_Blocks.shp', quiet=TRUE)
st_crs(la.sf)
la.sf <- st_transform(la.sf, crs = 6424) # EPSG code for NAD1983(2011) California V (ftUS)

# -------------------- LA W/ REDLINING -----------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
laR_data = read_csv('laDataFR.csv', show_col_types = FALSE)
laR_data$LandscapeShapeIndex[is.na(laR_data$LandscapeShapeIndex)] <- 0
laR_data <- laR_data %>%
  mutate(CTCB20 = as.character(CTCB20))
# Join the data frames based on the CTBKEY column
laR_sf <- la.sf %>%
  inner_join(laR_data, by = "CTCB20")
laR_sf2 <- laR_sf %>%
  select(OBJECTID, CTCB20, LST, NDVI, FctImp, LandscapeShapeIndex, PatchDensity, PerCover_Forest, grade, geometry)
laR_sf2 = na.omit(laR_sf2)
plot(laR_sf)

# Defining Spatial Neighbors
laR.nb <- poly2nb(laR_sf2, queen=TRUE) # Queen continuity
### Want row standardized first-order, adjacency matrix
# # but style='W' does not allow us to get a symmetric matrix
laR.rsw <- nb2mat(laR.nb, style='B', zero.policy = TRUE) # Row standardization

# creates a distance matrix between all census blocks
dist.mat <- st_distance(laR_sf2)
# determines which areas have no neighbors
no.neighbors <- which(rowSums(laR.rsw) == 0) # 2251

# determine which areas are closest to these areas
nearest.area <- rep(NA, length(no.neighbors))
for(i in 1:length(no.neighbors)){
  distances <- dist.mat[no.neighbors[i],]
  nearest.area[i] <- which(distances == distances[order(distances)][2])
  
  # Update W
  laR.rsw[no.neighbors[i], nearest.area[i]] <- 1
  laR.rsw[nearest.area[i], no.neighbors[i]] <- 1   # Keep W symmetric
}

colnames(laR.rsw) <- rownames(laR.rsw)
isSymmetric.matrix(laR.rsw)

laR.rsw <- mat2listw(laR.rsw, style='B', zero.policy=TRUE)

# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran <- moran.test(laR_sf2$LST, laR.rsw, 
                      alternative = 'greater') #  Test for clustering --> the default
gmoran
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 1.1116605279, while the expected value is -0.0001276976. Because the Moran's I value 
# is significantly larger than the E(I), there is evidence of some significant positive 
# spatial auto-correlation and clustering.

# Spatial Regression
sar <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp + NDVI +    
                  PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK --> system is computationally singular
sar2 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + NDVI +    
                  PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no fctImp --> system is computationally singular
sar3 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no NDVI --> system is computationally singular
sar4 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp + NDVI +  
                  grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no PerCover_Forest --> system is computationally singular
sar5 <- spautolm(LST ~ LandscapeShapeIndex + FctImp + NDVI +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # OK - no PatchDensity
# LSI is 3.093447e-05
# FctImp is 1.226140e+01
# NDVI is -1.523219e-01
# PerCover_Forest is -1.936843e+03
# gradeB is 2.295651e+02
# gradeC is 3.479789e+02
# gradeD is 3.824976e+02
sar6 <- spautolm(LST ~ PatchDensity + FctImp + NDVI +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no LSI --> system is computationally singular
sar7 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp + NDVI +    
                   PerCover_Forest,listw=laR.rsw,data=laR_sf2) # NOT OK - no grade --> system is computationally singular
sar8 <- spautolm(LST ~ LandscapeShapeIndex + FctImp + NDVI + PerCover_Forest,
                 listw=laR.rsw,data=laR_sf2) # OK - no grade and PatchDensity
# LSI is 6.350280e-05
# FctImp is 1.339360e+01
# NDVI is -1.936543e-01
# PerCover_Forest is -1.855657e+03

# -------------------- LA W/O REDLINING ----------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
la_data = read_csv('laDataF.csv', show_col_types = FALSE)
la_data$LandscapeShapeIndex[is.na(la_data$LandscapeShapeIndex)] <- 0
la_data <- la_data %>%
  mutate(CTCB20 = as.character(CTCB20))
# Join the data frames based on the CTBKEY column
la_sf <- la.sf %>%
  inner_join(la_data, by = "CTCB20")

# Defining Spatial Neighbors
la.nb <- poly2nb(la_sf, queen=TRUE) # Queen continuity
### Want row standardized first-order, adjacency matrix
# # but style='W' does not allow us to get a symmetric matrix
la.rsw <- nb2mat(la.nb, style='B', zero.policy = TRUE) # Row standardization

# creates a distance matrix between all census blocks
dist.mat2 <- st_distance(la_sf)
# determines which areas have no neighbors
no.neighbors2 <- which(rowSums(la.rsw) == 0) # Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 'which': vector memory limit of 16.0 Gb reached, see mem.maxVSize()

# UNABLE TO RUN DUE TO LIMIT OF VECTOR MEMORY ----------------------------------
# determine which areas are closest to these areas
nearest.area2 <- rep(NA, length(no.neighbors2))
for(i in 1:length(no.neighbors)){
  distances2 <- dist.mat2[no.neighbors2[i],]
  nearest.area2[i] <- which(distances2 == distances2[order(distances2)][2])
  
  # Update W
  la.rsw[no.neighbors2[i], nearest.area2[i]] <- 1
  la.rsw[nearest.area2[i], no.neighbors2[i]] <- 1   # Keep W symmetric
}

colnames(la.rsw) <- rownames(la.rsw)
isSymmetric.matrix(la.rsw)

la.rsw <- mat2listw(la.rsw, style='B', zero.policy=TRUE)

# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran2 <- moran.test(la_sf$LST, la.rsw, 
                      alternative = 'greater') #  Test for clustering --> the default
gmoran2
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 141.7929, while the expected value is -5.208876e-05. Because the Moran's I value 
# is significantly larger than the E(I), there is evidence of some significant positive 
# spatial auto-correlation and clustering.

# Spatial Regression
sar <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp + NDVI +    
                  PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK --> system is computationally singular
sar2 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + NDVI +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no fctImp --> system is computationally singular
sar3 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no NDVI --> system is computationally singular
sar4 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp + NDVI +  
                   grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no PerCover_Forest --> system is computationally singular
sar5 <- spautolm(LST ~ LandscapeShapeIndex + FctImp + NDVI +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # OK - no PatchDensity
# LSI is 3.093447e-05
# FctImp is 1.226140e+01
# NDVI is -1.523219e-01
# PerCover_Forest is -1.936843e+03
# gradeB is 2.295651e+02
# gradeC is 3.479789e+02
# gradeD is 3.824976e+02
sar6 <- spautolm(LST ~ PatchDensity + FctImp + NDVI +    
                   PerCover_Forest + grade,listw=laR.rsw,data=laR_sf2) # NOT OK - no LSI --> system is computationally singular
sar7 <- spautolm(LST ~ PatchDensity + LandscapeShapeIndex + FctImp + NDVI +    
                   PerCover_Forest,listw=laR.rsw,data=laR_sf2) # NOT OK - no grade --> system is computationally singular

