# Regression Diagnostics

# Suspected outliers and possibly cases with high leverage should be studied individually
# to decide whether or not they should be included. Systematic features in residual plots, 
# e.g. curvature or apparent non-constant variance—require to modify the structure of the 
# model to match the data more closely. 
# But distinct problems can also interact: if the errors have a skewed distribution, 
# apparent outliers may be produced in the direction of the skew. Transforming the 
# response to make the errors less skewed can solve this problem. Similarly, properly 
# modeling a nonlinear relationship may bring apparently outlying observations in line 
# with the rest of the data.

# Load required packages
utilis <- c('car', 'gvlma', 'nnet', 'mice', 'splines')
invisible(lapply(utilis, require, character.only = T))

# Load the list of imputed dataset
if (exists("imp_samp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_sample.rds' file
  imp_folder <- dirname(imp_path)
  imp_samp <- readRDS(imp_path)
}

# original dataset (with missing data)
data <- mice::complete(imp_samp, action = 30)
# data1 <- data[-c(1751, 336, 785, 2371, 1297, 1645), ] # influential point

# Define the main additive model including prenatal and postnatal stress, sex, age, 
# ethnicity, and maternal covariates (BMI, smoking and alcohol consumption)
intern  <- lm( intern_score_13_z ~ prenatal_stress_z + postnatal_stress_z + sex + age_child + 
                ethnicity +  m_bmi_before_pregnancy + m_smoking + m_drinking, data = data) 
fat.per <- update( intern, tot_fat_percent_13_z ~ .) 
fat.and <- update( intern, andr_fat_mass_13_z ~ .)

intern.pre  <- update( intern, . ~ . - postnatal_stress_z) 
fat.per.pre <- update( fat.per, . ~ . - postnatal_stress_z) 
fat.and.pre <- update( fat.and, . ~ . - postnatal_stress_z)
intern.pos  <- update( intern, . ~ . - prenatal_stress_z) 
fat.per.pos <- update( fat.per, . ~ . - prenatal_stress_z) 
fat.and.pos <- update( fat.and, . ~ . - prenatal_stress_z)

# Try transformations of the outcome
intern.sqrt <- update( intern, sqrt(.) ~ .)

# Evaluate the interaction between stress periods for each model
intern.int  <- update( intern,  . ~ . + prenatal_stress_z:postnatal_stress_z)
fat.per.int <- update( fat.per, . ~ . + prenatal_stress_z:postnatal_stress_z)
fat.and.int <- update( fat.and, . ~ . + prenatal_stress_z:postnatal_stress_z)

# Evaluate possible nonlinear effects of the two stress periods for each model 
intern.nlin <- update( intern, . ~ . - prenatal_stress_z - postnatal_stress_z 
                       + splines::ns(prenatal_stress_z, 3) + splines::ns(postnatal_stress_z, 3))
fat.per.nlin <- update( intern.nlin, tot_fat_percent_13_z ~ .)
fat.and.nlin <- update( intern.nlin, andr_fat_mass_13_z ~ .)

# Evaluate non-linearity for the age covariate. 
intern.nlin.AGE <- update( intern, . ~ . - age_child + splines::ns(age_child, 3))
fat.per.nlin.AGE <- update( intern.nlin.AGE, tot_fat_percent_13_z ~ .)
fat.and.nlin.AGE <- update( intern.nlin.AGE, andr_fat_mass_13_z ~ .)
# Evaluate non-linearity for maternal BMI
intern.nlin.BMI <- update( intern, . ~ . - m_bmi_before_pregnancy + splines::ns(m_bmi_before_pregnancy, 3))
fat.per.nlin.BMI <- update( intern.nlin.BMI, tot_fat_percent_13_z ~ .)
fat.and.nlin.BMI <- update( intern.nlin.BMI, andr_fat_mass_13_z ~ .)

anova(intern.int, intern)
anova(fat.per.int, fat.per)
anova(fat.and.int, fat.and)

anova(intern.pre, intern)
anova(intern.pos, intern)
anova(fat.per.pre, fat.per)
anova(fat.per.pos, fat.per)
anova(fat.and.pre, fat.and)
anova(fat.and.pos, fat.and)

anova(intern.nlin, intern)
anova(fat.per.nlin, fat.per)
anova(fat.and.nlin, fat.and)

anova(intern.nlin.AGE, intern)
anova(fat.per.nlin.AGE, fat.per)
anova(fat.and.nlin.AGE, fat.and)

anova(intern.nlin.BMI, intern)
anova(fat.per.nlin.BMI, fat.per)
anova(fat.and.nlin.BMI, fat.and)

# grp <- nnet::multinom(risk_groups_rec ~ prenatal_stress_z + postnatal_stress_z + 
#                  sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, 
#                  model = T, data = data)

### ======================= GENERAL FIT DIAGNOSTICS ======================== ###

# RESIDUALS VS. PREDICTORS & VS. FITTED VALUES
# Useful for revealing problems but less so for determining the exact nature of the problem. 
car::residualPlots(intern, main = "Internalizing")
car::residualPlots(intern.pre, main = "Internalizing (prenatal)")
car::residualPlots(intern.pos, main = "Internalizing (postnatal)")
car::residualPlots(fat.per, main = "Fat mass")
car::residualPlots(fat.per.pre, main = "Fat mass (prenatal)")
car::residualPlots(fat.per.pre, main = "Fat mass (postnatal)")
# If a linear model is correctly specified, the Pearson residuals are independent 
# of the fitted values and predictors, so these graphs should be null plots, with 
# no systematic features. One non-null plot (curved general trend) is enough to 
# suggest that the specified model does not match the data, generally implying a 
# failure of one or more assumptions. In the residual plot for a factor the boxes 
# should all have about the same center and spread.
# A lack-of-fit test is also computed for each numeric predictor (i.e. a t-test for 
# (regressor)^2 added to the model. Significant p-values indicate lack-of-fit, 
# confirming the nonlinear pattern visible in the graph. 
# For the plot of residuals vs fitted values, the Tukey’s test for non-additivity 
# is obtained by adding the squares of the fitted values to the model and refitting. 
# The test confirms the impression of curvature (the fitted model is not adequate).

# MARGINAL MODEL PLOTS
# A variation on the basic residual plot displaying the marginal relationship between 
# the response and each predictor (whereas the plot vs fitted values displays the 
# conditional distribution of the response given the fit o f the model).
# A lowess smooth is fit to the data (in blue) and to the fitted values (red dashed)
# and the two curves should match on each of the plots. Any mismatches are evidence 
# that the model does not fit the data well.
car::marginalModelPlots(intern, main = "Internalizing") # sd = T to check variance assumptions
car::marginalModelPlots(fat.per, main = "Fat mass")

### ============================= UNUSUAL DATA ============================= ###

# Unusual data in regression include outliers, high-leverage points, and influential 
# observations. 

# --------------------------------- Leverage --------------------------------- #
# Observations that are relatively far from the center of the regressor space, 
# taking account of the correlational pattern among the regressors, have potentially 
# greater influence on the least-squares regression coefficients: they have high 
# leverage. The most common measures of leverage are the h_i, or hat-values.
# One way of examining the hat-values and other individual-observation diagnostic 
# statistics is to construct index plots, graphing the statistics against the 
# corresponding observation indices.
car::influenceIndexPlot(intern, main = "Internalizing")
car::influenceIndexPlot(fat.per, main = "Fat mass")
# produces index plots of Cook’s distances, Studentized residuals, corresponding 
# Bonferroni p values for outlier testing, and the hat-values.
# Ids that stand out among hat-values, have regressor values are unusual relative 
# to the other ids. 

# Influence Plot is an an alternative to index plots of diagnostic statistics:
car::influencePlot(intern, id.method="identify", main="Influence Plot - internalizing", 
                   sub="Circle size is proportial to Cook's Distance" )
car::influencePlot(fat.per, id.method="identify", main="Influence Plot - fat mass", 
                   sub="Circle size is proportial to Cook's Distance" )

# --------------------------------- Outliers --------------------------------- #
# Regression outliers are y values that are unusual conditional on the values of 
# the predictors. Under a mean-shift outlier model: the t-statistic for testing 
# the null hypothesis of no mean-shift has n —k —2 df if errors are normally 
# distributed == Studentized_residual.

# The generic qqPlot function plots Studentized residuals against the corresponding 
# quantiles of t(n-k—2). By default, it generates a 95% pointwise confidence envelope 
# for the Studentized residuals, using a parametric version of the bootstrap.
# Also returns the observations with the largest absolute Studentized residuals.
car::qqPlot(intern, main="QQ Plot - internalizing") 
car::qqPlot(fat.per, main="QQ Plot - fat mass")
# Obs that stray outside of the confidence envelope are problematic. When the distribution 
# of Studentized residuals looks heavy-tailed compared to the reference t distribution: 
# perhaps a method of robust regression would be more appropriate for the data.

# Bonferroni p-value for most extreme obs (i.e. highest Studentized residuals)
car::outlierTest(intern) 
car::outlierTest(fat.per)

# ------------------------- Influential Observations ------------------------- #
# ADDED-VARIABLE PLOTS
# or partial-regression plots, are useful for studying the impact of observations 
# on regression coefficients. They contrast residuals from two auxiliary regressions:
# 1. Regress y on all the regressors excluding x. Res = the part of y that is not 
#    explained by all the regressors except x.
# 2. Regress x on the other regressors. Res = the part o f x that is not explained 
#    by the other regressors, that remains when we condition on them.
# Obs  that are farthest form the mean / with the largest absolute residuals are 
# indicated. Points far to the left or right represent observations for which the 
# value of x is unusual given the values of the other regressors. There points may 
# influence the slopes when are isolated. 
car::avPlots(intern, main = "Internalizing")
car::avPlots(fat.per, main = "Fat mass") 
# Note, these plots are misleading when diagnosing other sorts o f problems, such 
# as non-linearity. But are a useful diagnostic for finding potentially jointly 
# influential points: sets of points that are out of line with the rest of the data, 
# at the extreme left or right of the horizontal axis. 

# Leverage plots generalize added-variable plots to terms with more than 1 df 
# (factor or polynomial regressors). For terms with 1 df, leverage plots are very 
# similar to added-variable plots, except that slope is always equal to 1, not to 
# the corresponding regression coefficient. # car::leveragePlots(fit1) # car::leveragePlots(fit2)

# ---------------------------- Influence measures  --------------------------- #
# An observation that is both outlying and has high leverage exerts influence on 
# the regression coefficients: if the observation is removed, the coefficients 
# change considerably. 
# The most common summary measure of influence is Cook's distance. If any noteworthy 
# D are apparent, it is prudent to remove the corresponding cases one at the time,
# refit the regression, and see how the results change.

# Cook's D plot: identify D values > 4/(n-k-1)
plot(intern, which=4, cook.levels = 4/((nrow(data) - length(intern$coefficients) - 2)))
# Compare combination of exclusions
fit1.1 <- update(intern, subset = ! rownames(data) %in% c(2452))
compareCoefs(intern, fit1.1)

# Cook's D plot: identify D values > 4/(n-k-1)
plot(fat.per, which=4, cook.levels = 4/((nrow(data) - length(fat.per$coefficients) - 2)))
# Compare combination of exclusions
fit2.1 <- update(fat.per, subset = ! rownames(data) %in% c(501))
compareCoefs(fat.per, fit2.1)

# Rather than summarizing influence by looking at all coefficients simultaneously, 
# we could look at individual differences:
# dfb <- dfbetas(intern) 

# plot(dfb[, "prenatal_stress_z"], type = "h")
# labs <- ifelse(abs(dfb[, "prenatal_stress_z"]) > 0.2, row.names(data), "")
# text(dfb[, "prenatal_stress_z"],labs, pos=2, col='blue',cex=0.6)

# plot(dfb[ , c("prenatal_stress_z", "postnatal_stress_z")])
# identify(dfb[,"prenatal_stress_z"], dfb[,"postnatal_stress_z"], rownames(data))
# The negative relationship between the values for the two regressors reflects the 
# positive correlation of the regressors themselves. Observations on the top left
# make x coefficient smaller and y coefficient larger, and vice-versa. 

### ====================== NORMALITY OF RESIDUALS ========================== ###

# Departures from the assumption of normally distributed errors is probably the 
# most difficult problem to diagnose. Residuals can have substantially different 
# variances, can be strongly correlated, and tend to behave more like a normal 
# sample than do the original errors (supemormality property).
# A quantile-comparison plot of Studentized residuals against the t distribution 
# is useful in drawing attention to the tail behavior of the residuals, possibly 
# revealing heavy-tailed or skewed distributions. # qqPlot(fit, main="QQ Plot")
# Distribution of studentized residuals
plot(density(rstudent(intern)), main="Distribution of Studentized Residuals - Internalizing")
plot(density(rstudent(fat.per)), main="Distribution of Studentized Residuals - Fat mass")

# BOX-COX TRANSFORMATIONS of the outcome
# summary(pi <- car::powerTransform(fit2, family ="yjPower")) # ....
# INVERSE RESPONSE PLOTS
# car::inverseResponsePlot(fit2, id.n=4)

# PREDICTOR TRANSFORMATIONS ...

### =========================== HOMOSCEDASTICITY =========================== ###

# Non-constant error variance test
car::ncvTest(intern) 
car::ncvTest(fat.per)
# plot studentized residuals vs. fitted values
car::spreadLevelPlot(intern)
car::spreadLevelPlot(fat.per)

### ============================ COLLINEARITY ============================== ###
vif(intern) # variance inflation factors
sqrt(vif(intern)) > 2 # problem?
vif(fat.per) # variance inflation factors
sqrt(vif(fat.per)) > 2 # problem?

### ============================ NONLINEARITY ============================== ###

# Component + residual plot
car::crPlots(intern) 
car::crPlots(fat.per)
# Ceres plots (generalization of component + residual plot)
car::ceresPlots(intern)
car::ceresPlots(fat.per)

# Test for Autocorrelated Errors
car::durbinWatsonTest(intern)
car::durbinWatsonTest(fat.per)


# Global test of model assumptions #############################################

gvmodel1 <- gvlma::gvlma(intern)
gvmodel2 <- gvlma::gvlma(fat.per)
summary(gvmodel1)
summary(gvmodel2)
