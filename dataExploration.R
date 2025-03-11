library(readr)

# -------------------- HARRIS W/ REDLINING -------------------------------------
# Read in the data
harrisR_data = read_csv('harrisDataFR.csv', show_col_types = FALSE)
# Look at the first 6 observations
head(harrisR_data)
# Check the dimension
dim(harrisR_data)
summary(harrisR_data)
names(harrisR_data)

# check how many NA values in LSI column (meaning there is no greenspace)
sum(is.na(harrisR_data$Landscape_Shape_Index)) # 945
# change those NA values to 0
harrisR_data$Landscape_Shape_Index[is.na(harrisR_data$Landscape_Shape_Index)] <- 0

# Fit the multiple linear regression model
harrisR_model = lm(formula = LST ~ Patch_Density + 
                        Landscape_Shape_Index + FctImp + NDVI +    
                        PerCover_Forest + grade,data = harrisR_data)
summary(harrisR_model)
# check for linearity and equal variance assumption
plot(resid(harrisR_model)~fitted(harrisR_model), main="Residuals vs. Predicted Plot")
# No trend → equal variance, Linearity → they are close to the zero line

# Get the model residuals
model_residuals = harrisR_model$residuals
# Plot the result
hist(model_residuals)
# The histogram is skewed to the left; hence we can not conclude the 
# normality with enough confidence. Instead of the histogram, let’s look at the 
# residuals along the normal Q-Q plot. If there is normality, then the values should follow a straight line.

# Plot the residuals
qqnorm(model_residuals)
# Plot the Q-Q line
qqline(model_residuals)
# From the plot, we can observe that the residuals roughly lie in a straight line. 
# Then we can assume that the residuals of the model follow a normal distribution. 

# -------------------- LA W/ REDLINING -----------------------------------------
# Read in the data
laR_data = read_csv('laDataFR.csv', show_col_types = FALSE)
# Look at the first 6 observations
head(laR_data)
# Check the dimension
dim(laR_data)
summary(laR_data)
names(laR_data)

# check how many NA values in LSI column (meaning there is no greenspace)
sum(is.na(laR_data$LandscapeShapeIndex)) # 9828
# change those NA values to 0
laR_data$LandscapeShapeIndex[is.na(laR_data$LandscapeShapeIndex)] <- 0

# Fit the multiple linear regression model
laR_model = lm(formula = LST ~ PatchDensity + 
                     LandscapeShapeIndex + FctImp + NDVI +    
                     PerCover_Forest + grade,data = laR_data)
summary(laR_model)

# check for linearity and equal variance assumption
plot(resid(laR_model)~fitted(laR_model), main="Residuals vs. Predicted Plot")
# No trend → equal variance, Linearity → they are close to the zero line

# Get the model residuals
model_residuals2 = laR_model$residuals
# Plot the result
hist(model_residuals2)
# The histogram is ever so slightly skewed to the right; hence we can not conclude the 
# normality with enough confidence. Instead of the histogram, let’s look at the 
# residuals along the normal Q-Q plot. If there is normality, then the values should follow a straight line.

# Plot the residuals
qqnorm(model_residuals2)
# Plot the Q-Q line
qqline(model_residuals2)
# From the plot, we can observe that the residuals do not roughly lie in a straight line,
# as there are long tails.
# Then we can assume that the residuals of the model do not follow a normal distribution.

# -------------------- HARRIS W/O REDLINING -------------------------------------
# Read in the data
harris_data = read_csv('harrisDataF.csv', show_col_types = FALSE)
# Look at the first 6 observations
head(harris_data)
# Check the dimension
dim(harris_data)
summary(harris_data)
names(harris_data)

# check how many NA values in LSI column (meaning there is no greenspace)
sum(is.na(harris_data$Landscape_Shape_Index)) # 12545
# change those NA values to 0
harris_data$Landscape_Shape_Index[is.na(harris_data$Landscape_Shape_Index)] <- 0

# Fit the multiple linear regression model
harris_model = lm(formula = LST ~ Patch_Density + 
                     Landscape_Shape_Index + FctImp + NDVI +    
                     PerCover_Forest,data = harris_data)
summary(harris_model)

# check for linearity and equal variance assumption
plot(resid(harris_model)~fitted(harris_model), main="Residuals vs. Predicted Plot")
# No trend → equal variance, Linearity → they are close to the zero line

# Get the model residuals
model_residuals3 = harris_model$residuals
# Plot the result
hist(model_residuals3)
# The histogram is slightly skewed to the left; hence we can not conclude the 
# normality with enough confidence. Instead of the histogram, let’s look at the 
# residuals along the normal Q-Q plot. If there is normality, then the values should follow a straight line.

# Plot the residuals
qqnorm(model_residuals3)
# Plot the Q-Q line
qqline(model_residuals3)
# From the plot, we can observe that the residuals do not roughly lie in a straight line. 
# Then we can assume that the residuals of the model do not follow a normal distribution.

# -------------------- LA W/O REDLINING ----------------------------------------
# Read in the data
la_data = read_csv('laDataF.csv', show_col_types = FALSE)
# Look at the first 6 observations
head(la_data)
# Check the dimension
dim(la_data)
summary(la_data)
names(la_data)

# check how many NA values in LSI column (meaning there is no greenspace)
sum(is.na(la_data$LandscapeShapeIndex)) # 21521
# change those NA values to 0
la_data$LandscapeShapeIndex[is.na(la_data$LandscapeShapeIndex)] <- 0

# Fit the multiple linear regression model
la_model = lm(formula = LST ~ PatchDensity + 
                 LandscapeShapeIndex + FctImp + NDVI +    
                 PerCover_Forest,data = la_data)
summary(la_model)

# check for linearity and equal variance assumption
plot(resid(la_model)~fitted(la_model), main="Residuals vs. Predicted Plot")
# has a weird trend → NOT equal variance, NO linearity → they are not close to the zero line

# Get the model residuals
model_residuals4 = laR_model$residuals
# Plot the result
hist(model_residuals4)
# The histogram is ever so slightly skewed to the right; hence we can not conclude the 
# normality with enough confidence. Instead of the histogram, let’s look at the 
# residuals along the normal Q-Q plot. If there is normality, then the values should follow a straight line.

# Plot the residuals
qqnorm(model_residuals4)
# Plot the Q-Q line
qqline(model_residuals4)
# From the plot, we can observe that the residuals do not roughly lie in a straight line,
# as there are long tails.
# Then we can assume that the residuals of the model do not follow a normal distribution.
