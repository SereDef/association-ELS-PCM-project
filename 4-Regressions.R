
# Final step in our journey,
# The following code runs a set of analyses for the project "Early-life stress as a risk factor 
# for poor mental and physical health in children: A population-based study" looking at the 
# association between ELS and psycho-cardio-metabolic multi-morbidity at age 10.  

# All it requires is just one file: the imputed dataset we built using the "Imputation.R"
# script, that you can find in this repository: https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

# As usual, here are the packages we need 
library(mice);
library(nnet);
library(openxlsx)

# check if the path to the dataset is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

#------------------------------------------------------------------------------#
# In the previous step (3-Imputation.R), we created several complete versions of the data 
# by replacing the missing values with plausible data values. Now we need to estimate the 
# parameters of interest from each imputed dataset and pool the estimates into one.

# Load the list of imputed dataset
imp <- readRDS(paste(pathtodata,'imputation_list.rds', sep = ""))

# Define a customized function that pools regression results and builds a dataframe 
# that is easy to read and save
modeltable <- function(fit, logm = FALSE) {
  p_fit <- mice::pool(fit) # pool results 
  mod <- summary(p_fit) # extract relevant information
  mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
  if (logm == FALSE) {
    mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)) # adjusted R2
  } else {
    levels(mod$y.level) <- c("intern", "fatmas", "multim") # make group comparisons easier to read
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
  }
  return(mod)
}

################################################################################
################################################################################

# First, we assess the effect of stress on internalizing and fat mass separately 
# in two multiple linear regressions. Each model has a minimally adjusted (sex and age)
# and fully adjusted version (sex, age, maternal BMI, smoking, drinking)

# There we go:
#------------------------------ INTERNALIZING ---------------------------------#
fit_is <- with(imp, lm(intern_score_z ~ prenatal_stress_z*postnatal_stress_z + sex + age_child))
mod1 <- modeltable(fit_is)

# Fully adjusted model
fit_is_full <- with(imp, lm(intern_score_z ~ prenatal_stress_z*postnatal_stress_z + sex + age_child + 
                              m_bmi_berore_pregnancy + m_smoking + m_drinking))
mod2 <- modeltable(fit_is_full)


#------------------------------- FAT MASS -------------------------------------#
fit_fm <- with(imp, lm(fat_mass_z ~ prenatal_stress_z*postnatal_stress_z + sex + age_child))
mod3 <- modeltable(fit_fm)

fit_fm_full <- with(imp, lm(fat_mass_z ~ prenatal_stress_z*postnatal_stress_z + sex + age_child + 
                              m_bmi_berore_pregnancy + m_smoking + m_drinking))
mod4 <- modeltable(fit_fm_full)

################################################################################
################################################################################
################################################################################

# Now let's examine the multimorbidity aspect using four groups defined in script 2.
# This includes a minimally adjusted and fully adjusted Multinomial Logistic Regression
# estimated using a Neural Network. 

# NOTE: the ratio of number of cases is important! 
# It is always easy to predict the larger category because probabilistically more likely
# rule of thumb is no higher that 1 to 5 ratio. 

# Anywho, let's run the model! 

fit_grp <- with(imp, nnet::multinom(risk_groups ~ prenatal_stress_z*postnatal_stress_z + 
                                      sex + age_child, model = T));
mod5 <- modeltable(fit_grp, logm = T)

# Fully adjusted model 
fit_grp_full <- with(imp, nnet::multinom(risk_groups ~ prenatal_stress_z*postnatal_stress_z + 
                                sex + age_child + m_bmi_berore_pregnancy + m_smoking + m_drinking, model = T));
mod6 <- modeltable(fit_grp_full, logm = T)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Finally,
# Let's have a look at the contribution of individual domains of stress to the probability 
# of belonging to each group. Again fist a minumally adjusted and then a fully adjusted 
# model is run. 

fit_grp_dm <- with(imp, nnet::multinom(risk_groups ~ pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress +
                                 post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                 sex + age_child, model = T))
mod7 <- modeltable(fit_grp_dm, logm = T)

# Fully adjusted model
fit_grp_dm_full <- with(imp, nnet::multinom(risk_groups ~ pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress +
                                   post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                   sex + age_child + m_bmi_berore_pregnancy + m_smoking + m_drinking, model = T))
mod8 <- modeltable(fit_grp_dm_full, logm = T)

################################################################################

# Export the outputs of summary statistics into an xlsx file with one model per sheet

modls <- list("1.intern_min" = mod1, "2.intern_ful" = mod2, 
              "3.fatmas_min" = mod3, "4.fatmas_ful" = mod4, 
              "5.riskgrp_min" = mod5, "6.riskgrp_ful" = mod6, 
              "7.domains_min" = mod7, "8.domains_ful" = mod8)

write.xlsx(modls, file = paste0(pathtodata, "Results.xlsx"))


################################################################################
############################### THE END ########################################
################################################################################