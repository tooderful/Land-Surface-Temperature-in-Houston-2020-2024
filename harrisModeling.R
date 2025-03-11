packs <- c('tidyverse', 'sf', 'spdep', 'spData', 'sfdep', 'viridis', 'gridExtra')
lapply(packs, require, character.only = TRUE)

# read in the census block shapefile
tx.sf <- st_read('harrisCountyProcessedData/blocks/Blocks.shp', quiet=TRUE)
st_crs(tx.sf)
tx.sf <- st_transform(tx.sf, crs = 6588) # EPSG code for NAD1983(2011) Texas South Central (ftUS)

# -------------------- HARRIS W/ REDLINING -------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
harrisR_data = read_csv('harrisDataFR.csv', show_col_types = FALSE)
harrisR_data <- harrisR_data %>%
  mutate(CTBKEY = as.character(CTBKEY))
# Join the data frames based on the CTBKEY column
harrisR_sf <- tx.sf %>%
  inner_join(harrisR_data, by = "CTBKEY")
plot(harrisR_sf)

# Defining Spatial Neighbors
harrisR.nb <- poly2nb(harrisR_sf, queen=TRUE) # Queen continuity
harrisR.rsw <- nb2listw(harrisR.nb, style='W', zero.policy = TRUE) # Row standardization

# CURRENTLY INCORRECT, GIVING ZERO ISLANDS??
isolated_polygons <- lengths(harrisR.nb) == 0
harrisR_sf$island <- isolated_polygons
harrisR_sf$island <- as.logical(harrisR_sf$island)
print(sum(harrisR_sf$island))

# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran <- moran.test(harrisR_data$LST, harrisR.rsw, 
                     alternative = 'greater') #  Test for clustering --> the default
# z-value is the Moran I statistic standard deviate or the Moran I statistics - expectation / squareroot of variance
z=(0.432614830 - (-0.001838235))/sqrt(0.002383118)
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 8.899592, while the expected value is -0.001838235. Because the Moran's I value 
# is larger than the E(I), there is evidence of some significant positive spatial 
# auto-correlation and clustering.

# -------------------- HARRIS W/O REDLINING ------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
harris_data = read_csv('harrisDataF.csv', show_col_types = FALSE)
harris_data <- harris_data %>%
  mutate(CTBKEY = as.character(CTBKEY))
# Join the data frames based on the CTBKEY column
harris_sf <- tx.sf %>%
  inner_join(harris_data, by = "CTBKEY")

# Defining Spatial Neighbors
harris.nb <- poly2nb(harris_sf, queen=TRUE) # Queen continuity
harris.rsw <- nb2listw(harris.nb, style='W', zero.policy = TRUE) # Row standardization

# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran2 <- moran.test(harris_data$LST, harris.rsw, 
                     alternative = 'greater') #  Test for clustering --> the default
# z-value is the Moran I statistic standard deviate or the Moran I statistics - expectation / squareroot of variance
z=(8.278723e-01 - (-7.702973e-05))/sqrt(5.993047e-05)
gmoran2
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 106.9498, while the expected value is -7.702973e-05. Because the Moran's I value 
# is significantly larger than the E(I), there is evidence of some significant 
# positive spatial auto-correlation and clustering.
