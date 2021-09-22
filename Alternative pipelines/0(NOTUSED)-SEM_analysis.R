################################################################################
####################### NOT USED IN THE FINAL VERSION ##########################
################################################################################

# Hola,
# The following code runs the final SEM analysis for the project "Early-life stress as a risk factor 
# for poor mental and physical health in children: A population-based study" looking at the 
# association between ELS and psycho-cardio-metabolic multi-morbidity at age 10.  

# All it requires is just one file: the imputed dataset we built using the "Imputation.R"
# script, that you can find in this repository: https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

#### ---------------------------- Dependencies ---------------------------- ####

library(lavaan) 
library(semPlot)
library(car)

# check if the path to the dataset is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

#------------------------------------------------------------------------------#
# Load the (complete) dataset
datarisk <- readRDS(paste(pathtodata, 'ELS_PCM_imputed.rds', sep = ""))

#------------------------------------------------------------------------------#
# Check collinearity 

regpost <- lm(postrisk ~ post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization, 
          data = datarisk)
regpre <- lm(prerisk ~ pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress, 
             data = datarisk)
reg <- lm(els ~ prerisk + postrisk, data = datarisk)

# print out model results
anova(reg)
summary(regore)

# test of multicollinearity
vif(reg)

#------------------------------------------------------------------------------#
# Check summed scores

datarisk$presum <- rowSums(datarisk[c("pre_life_events","pre_contextual_risk", 
                                      "pre_personal_stress", "pre_interpersonal_stress")], na.rm = F)
datarisk$postsum <- rowSums(datarisk[c("post_life_events","post_contextual_risk", 
                                      "post_parental_risk", "post_interpersonal_risk", 
                                      "post_direct_victimization")], na.rm = F)
datarisk$elssum <- rowSums(datarisk[c("presum","postsum")],  na.rm = F)
datarisk$multisum <- rowSums(datarisk[c("intern_score_z", "fat_mass_z")], na.rm = F)

# simple model: 
sum_model <- "
# Latent overall risk (pre and postnatal contributions constrained to be equal)
#    els =~ same*presum + same*postsum
# Latent outcome: psycho-cardio-metabolic-risk (mental and physical contributions constrained to be equal)
    pcmr =~ same*intern_score_z + same*fat_mass_z
# fix variance of all latent factors to unity (to standardize)
#    els ~~ 1*els
    pcmr ~~ 1*pcmr
# Regressions
    pcmr ~ els
    intern_score_z ~ els
    fat_mass_z ~ els"

# Fit the model
fitsum = sem(sum_model, datarisk, std.lv = TRUE) # fix the variances of all the latent variables to unity 
# (factor loadings of the first indicator will no longer be fixed to 1).

summary(fitsum, fit.measures = TRUE, standardized = TRUE) 

##----------------------------------------------------------------------------##
out_model = "
pcmr =~ same*intern_score_z + same*fat_mass_z
pcmr ~~ 1*pcmr"
fitout = cfa(out_model, datarisk)
summary(fitout, fit.measures = TRUE, standardized = TRUE)

##----------------------------------------------------------------------------##
## ------------------------------- ANALYSIS --------------------------------- ##
##----------------------------------------------------------------------------##

# In the model definition, I first specify the correlations and latent variables, 
# then I estimate the association latent overall risk exposure and and latent pcm 
# outcome via linear regression. In then adjust for confounding from age and gender
# (adj_model1) and maternal drinking, smoking and BMI during pregnancy (adj_model2).

# Basic (un-adjusted) model description
basic_model <- "
# Latent prenatal risk
    prerisk =~ NA*pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress
# Latent postnatal risk
    postrisk =~ NA*post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization
# Latent overall risk (pre and postnatal contributions constrained to be equal)
    els =~ samev*prerisk + samev*postrisk
# Latent outcome: psycho-cardio-metabolic-risk (mental and physical contributions constrained to be equal)
    pcmr =~ samev2*intern_score_z + samev2*fat_mass_z
# fix variance of all latent factors to unity (to standardize)
    postrisk ~~ 1*postrisk
    prerisk ~~ 1*prerisk
    els ~~ 1*els
    pcmr ~~ 1*pcmr
# Residual correlations
    pre_life_events ~~ post_life_events
    pre_contextual_risk ~~ post_contextual_risk
    pre_personal_stress ~~ post_parental_risk
    pre_interpersonal_stress ~~ post_interpersonal_risk
# Regressions
    pcmr ~ els
    intern_score_z ~ els
    fat_mass_z ~ els"

# Fit the model
fitit = sem(basic_model, datarisk, 
            std.lv = TRUE) # fix the variances of all the latent variables to unity 
                           # (factor loadings of the first indicator will no longer be fixed to 1).

# Inspect the parameter’ table, which provides a summary of what parameters are free 
# in the model (i.e., have to be estimated), and what parameters were requested by the user in the model syntax
parTable(fitit)
# ‘user’ refers to a parameter request explicitly in the syntax
# non-zero values for the ‘free’ column denote parameters that are freely estimated.
varTable(fitit)

# Inspect the results
summary(fitit, fit.measures = TRUE, standardized = TRUE) # OR parameterestimates(fitit), to get the unstandardized parameters
# OR:
stres = standardizedsolution(fitit)
stres[ stres$op == "=~",] # just the loadings (with CIs)
mean(abs(stres[ stres$op == "=~", "est.std"])) # mean absolute loading 
stres[ stres$op == "~~" & stres$lhs != stres$rhs, ] # Just correlations
stres[ stres$op == "~~" & stres$lhs == stres$rhs, ] # Just residual variances

# Free parameters
coef(fitit)

# Fit measures
# We select a subset of fit indices, following in Schuurmans et al, approach
fitmeasures(fitit, c('chisq', 'df', 'pvalue', 'cfi', 'tli', 
                     'rmsea', 'rmsea.ci.lower', 'rmsea.ci.upper', 'srmr'))

# Other model features
inspect(fitit, 'r2') # Variance explained by the latent factor in each indicator (i.e. square of loadings)
inspect(fitit, "sampstat") # Indicators' covariance matrix
std_fit <- inspect(fitit, "std") # Standardized model matrices, "est" for unstandardized
std_fit$lambda # Standardized loadings
std_fit$psi # Latent variable correlation matrix

## -------------------------------- Plotting -------------------------------- ##
semPaths(fitit, what = "std", 
         layout = "tree2", 
         nCharNodes = 8, # number of characters that are not omitted
         sizeMan = 7, 
         sizeMan2 = 5)

## ----------------------------- Adjusted models ---------------------------- ##

# adding confoundinfg terms to the model description
adj_model1 <- paste(basic_model, "+ sex + age_child")
adj_model2 <- paste(adj_model1, "+ m_bmi_berore_pregnancy + m_smoking + m_drinking")

fititbetter = sem(adj_model1, datarisk)
fititbest = sem(adj_model2, datarisk)
## Warning message:
##   In lav_object_post_check(object) :
##   lavaan WARNING: some estimated lv variances are negative

summary(fititbetter, standardized = TRUE)
summary(fititbest, standardized = TRUE)

## ------------------- Compare basic and adjusted models -------------------- ##

# Store all model fits into a list
fits <- list()
fits$fit1 <- fitit; fits$fit2 <- fititbetter; fits$fit3 <- fititbest

# Combine fitmeausures into a table (tip: sapply is a function designed to loop through lists)
round(sapply(fits, 
      function(X) fitmeasures(X,  c('chisq', 'df', 'pvalue', 'cfi', 'tli', 'rmsea', 'rmsea.ci.lower', 'rmsea.ci.upper', 'srmr'))), 
      3) # Round them to 3 decimals for better visualization

# Test differences between models
anova(fits$fit1, fits$fit2, fits$fit3)

## ---------------------------- Model improvement --------------------------- ##

# Because the fit isn't great, we first inspect the residuals
# i.e. to what extent is the correlation between two items in the data is not captured 

# Standardized residuals correlation matrix
resid(fitit, type = "cor") # OR resid(fitit) for unstandardized residuals covariance matrix

# MODICIFICATION INDICES
mod_ind <- modificationindices(fitit) # Returns a list of all the possible modifications to the model 
# i.e. all the possible cross-loadings and correlated residuals.
sorted_mi = mod_ind[order(mod_ind$mi, decreasing=TRUE), ] # Sort them by size
head(sorted_mi, 10) # Top 10 modification indices
sorted_mi[sorted_mi$mi > 10, ] # All modification indices that are above 10

## ------------------------- Save latent variables -------------------------- ##

# Create a dataframe with the extracted latent variables
pfit <- data.frame(predict(fitit)) 

# Explore their relation 
plot(pfit$prerisk, pfit$postrisk)
plot(pfit$els, pfit$pcmr)

# add predicted values to the dataset
datarisk <- cbind(datarisk, pfit)
datarisk2 <- cbind(datarisk, pfit$prerisk, pfit$postrisk)
names(datarisk2)[124:125] = c('prerisk', 'postrisk')

################################################################################
# Model using extracted variables 

extr_model <- "
# Latent overall risk
    els =~ prerisk + postrisk
# Latent outcome (psycho-cardio-metabolic-risk)
    pcmr =~ intern_score_z + fat_mass_z
# Regressions
    pcmr ~ els"  # pcmr ~ intern_score_z  # pcmr ~ fat_mass_z are other alternatives

fititagain = sem(extr_model, datarisk2)
summary(fititagain, standardized = TRUE)

semPaths(fititagain, "std", 
         layout = "tree2", 
         nCharNodes = 8, # number of characters that are not omitted
         sizeMan = 7, 
         sizeMan2 = 5)

################################################################################
# Model only prenatal 
pre_model <- "
# Latent prenatal risk
    prerisk =~ pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress
# Latent outcome (psycho-cardio-metabolic-risk)
    pcmr =~ intern_score_z + fat_mass_z
# Regressions
    pcmr ~ prerisk"  # pcmr ~ intern_score_z  # pcmr ~ fat_mass_z are other alternatives

fititpre = sem(pre_model, datarisk)
summary(fititpre, standardized = TRUE)

semPaths(fititpre, "std", 
         layout = "tree2", 
         nCharNodes = 8, # number of characters that are not omitted
         sizeMan = 7, 
         sizeMan2 = 5)       

################################################################################
# Model only post
post_model <- "
# Latent postnatal risk
    postrisk =~ post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization
# Latent outcome (psycho-cardio-metabolic-risk)
    pcmr =~ intern_score_z + fat_mass_z
# Regressions
    pcmr ~ postrisk"  # pcmr ~ intern_score_z  # pcmr ~ fat_mass_z are other alternatives

fititpost = sem(post_model, datarisk)
summary(fititpost, standardized = TRUE)

semPaths(fititpost, "std", 
         layout = "tree2", 
         nCharNodes = 8, # number of characters that are not omitted
         sizeMan = 7, 
         sizeMan2 = 5)       

################################################################################
# Model mediation 
med_model <- "
# Latent prenatal risk
    prerisk =~ NA*pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress
# Latent postnatal risk
    postrisk =~ NA*post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization
# Latent outcome (psycho-cardio-metabolic-risk)
    pcmr =~ samev*intern_score_z + samev*fat_mass_z
# fix variance of all latent factors to unity (to standardize)
    postrisk ~~ 1*postrisk
    prerisk ~~ 1*prerisk
    pcmr ~~ 1*pcmr 
# Residual correlations
    pre_life_events ~~ post_life_events
    pre_contextual_risk ~~ post_contextual_risk
    pre_personal_stress ~~ post_parental_risk
    pre_interpersonal_stress ~~ post_interpersonal_risk
# direct effect
    pcmr ~ c*prerisk
# mediator
    postrisk ~ a*prerisk
    pcmr ~ b*postrisk
# indirect effect (a*b)
    ab := a*b 
# total effect
    total := c + (a*b)"  # pcmr ~ intern_score_z  # pcmr ~ fat_mass_z are other alternatives

med_model <- "
# Latent outcome (psycho-cardio-metabolic-risk)
    pcmr =~ samev*intern_score_z + samev*fat_mass_z
# fix variance of all latent factors to unity (to standardize)
    pcmr ~~ 1*pcmr 
# direct effect
    pcmr ~ c*presum
# mediator
    postsum ~ a*presum
    pcmr ~ b*postsum
# indirect effect (a*b)
    ab := a*b 
# total effect
    total := c + (a*b)" 

fitmed = sem(med_model, datarisk)
summary(fitmed, standardized = TRUE)

fitmeasures(fitmed, c('chisq', 'df', 'pvalue', 'cfi', 'tli', 
                     'rmsea', 'rmsea.ci.lower', 'rmsea.ci.upper', 'srmr'))

semPaths(fitmed, "std", 
         layout = "tree2", 
         nCharNodes = 8, # number of characters that are not omitted
         sizeMan = 7, 
         sizeMan2 = 5)       

