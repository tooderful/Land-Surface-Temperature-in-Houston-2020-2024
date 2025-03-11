packs <- c('tidyverse', 'sf', 'spdep', 'spData', 'sfdep', 'viridis', 'gridExtra')
lapply(packs, require, character.only = TRUE)

# read in the census block shapefile
la.sf <- st_read('laCountyProcessedData/2020_Census_Blocks/2020_Census_Blocks.shp', quiet=TRUE)
st_crs(la.sf)
la.sf <- st_transform(la.sf, crs = 6424) # EPSG code for NAD1983(2011) California V (ftUS)

# -------------------- LA W/ REDLINING -----------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
laR_data = read_csv('laDataFR.csv', show_col_types = FALSE)
laR_data <- laR_data %>%
  mutate(CTCB20 = as.character(CTCB20))
# Join the data frames based on the CTBKEY column
laR_sf <- la.sf %>%
  inner_join(laR_data, by = "CTCB20")
plot(laR_sf)

# Defining Spatial Neighbors
laR.nb <- poly2nb(laR_sf, queen=TRUE) # Queen continuity
laR.rsw <- nb2listw(laR.nb, style='W', zero.policy = TRUE) # Row standardization

# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran3 <- moran.test(laR_data$LST, laR.rsw, 
                      alternative = 'greater') #  Test for clustering --> the default
# z-value is the Moran I statistic standard deviate or the Moran I statistics - expectation / squareroot of variance
z=(0.7122937778 - (-0.0001699524))/sqrt(0.0001985806 )
gmoran3
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 50.55852, while the expected value is -0.0001699524. Because the Moran's I value 
# is significantly larger than the E(I), there is evidence of some significant positive 
# spatial auto-correlation and clustering.

# -------------------- LA W/O REDLINING ----------------------------------------
# filter the census blocks for those that are in the harrisFR.csv
la_data = read_csv('laDataF.csv', show_col_types = FALSE)
la_data <- la_data %>%
  mutate(CTCB20 = as.character(CTCB20))
# Join the data frames based on the CTBKEY column
la_sf <- la.sf %>%
  inner_join(la_data, by = "CTCB20")

# Defining Spatial Neighbors
la.nb <- poly2nb(la_sf, queen=TRUE) # Queen continuity
la.rsw <- nb2listw(la.nb, style='W', zero.policy = TRUE) # Row standardization


# Moran's I statistic - testing spatial autocorrelation
# The null hypothesis is that there is no spatial autocorrelation.
# The alternative hypothesis is that there is some spatial autocorrelation, specifically positive spatial autocorrelation or clustering in which the Moran's I values are significantly larger than E(I).
# The significance level we establish is 0.05.
# P-value approach
gmoran4 <- moran.test(la_data$LST, la.rsw, 
                      alternative = 'greater') #  Test for clustering --> the default
# z-value is the Moran I statistic standard deviate or the Moran I statistics - expectation / squareroot of variance
z=(9.420728e-01 - (-5.208876e-05))/sqrt(4.414767e-05 )
gmoran4
# analytical approach for the Global Moran's I test - we reject the null in favor 
# of the alternative (p-value = 2.2e-16 < \alpha = 0.05). The Moran's I equals 
# 141.7929, while the expected value is -5.208876e-05. Because the Moran's I value 
# is significantly larger than the E(I), there is evidence of some significant positive 
# spatial auto-correlation and clustering.
