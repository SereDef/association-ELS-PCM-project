# Hi again, 
# the following code runs the imputation of missing data, necessary for the analysis
# of the association between ELS and psycho-cardio-metabolic multi-morbidity in children. 
# Ok, let's get started!

#### ---------------------------- Dependencies ---------------------------- ####

PATH_RESULTS <- '' # ATTENTION! DEFINE HERE THE PATH WHERE YOU WANT ALL RESULTS TO BE STORED.

# Load libraries
library(mice); library(miceadds)

#-------------------------------------------------------------------------------
flowchart <- function(df, return_selected_sample = F) {
  fc <- list(initial_sample = nrow(df))
  # enough prenatal variables
  step1 <- df[df$pre_percent_missing < 50.0,]
  loss <- nrow(step1) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_pren = loss, after_pren_selection = nrow(step1))
  # enough postnatal variables
  step2 <- step1[step1$post_percent_missing < 50.0,]
  loss <- nrow(step2) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_post = loss, after_post_selection = nrow(step2))
  # no selection on the outcome
  step3 <- step2[step2$twin == 0,] 
  loss <- nrow(step3) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_twins = loss, after_twin_selection = nrow(step3))
  #  exclude siblings
  finalsample <- step3[step3$sibling == 0,]
  loss <- nrow(finalsample) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_siblings = loss, final_sample = nrow(finalsample))
  
  print(fc)
  
  if (return_selected_sample == T) { return(finalsample) }
}

#-------------------------------------------------------------------------------
# define a quick function for passive imputation of domain scores
passive_imp_formula <- function(domain, add = "", type_nr = 1) {
  conc <- paste(domain, collapse = " + ")
  if (type_nr > 1) { div <- type_nr } else { div <- length(domain) }
  str <- paste0("~I( (", conc, ") / ", div, ")")
  if (add != "") {
    str <- paste0("~I( (", conc, " + ", add, ") / ", div+1, ")")
  }
  return(str)
}

#-------------------------------------------------------------------------------
# define a function that randomly shuffles internalizing and fat mass values and 
# returns the new size of the "randomly multimorbid" group.
permute <- function(df) { 
  # create empty dataset for permutation
  perm <- data.frame(1:nrow(df))
  
  perm$new_int <- sample(df$int)
  perm$new_fat <- sample(df$fat)
  perm$new_groups <- ifelse(perm$new_int == 0 & perm$new_fat == 0, 0, 
                            ifelse(perm$new_int == 1 & perm$new_fat == 0, 1, 
                                   ifelse(perm$new_int == 0 & perm$new_fat == 1, 2, 3)))
  new_n <- unname(summary(as.factor(perm$new_groups))[4])
  return(new_n)
}
#-------------------------------------------------------------------------------
# define a function that construncts the groups used as outcome for the third set 
# of models 
construct_grp <- function(int_var, fm_var, df, cutoff = 0.8, permute = T) {
  df$int = ifelse(df[, int_var] > quantile(df[, int_var], probs = cutoff, na.rm = T), 1, 0) 
  df$fat = ifelse(df[, fm_var]  > quantile(df[, fm_var],  probs = cutoff, na.rm = T), 1, 0) 
  df$risk_groups = rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    if ( is.na(df$int[i]) | is.na(df$fat[i]) )  { df$risk_groups[i] = NA
    } else if (df$int[i] == 0 & df$fat[i] == 0) { df$risk_groups[i] = 0   # Healthy
    } else if (df$int[i] == 1 & df$fat[i] == 0) { df$risk_groups[i] = 1   # High internalizing  only
    } else if (df$int[i] == 0 & df$fat[i] == 1) { df$risk_groups[i] = 2   # High fat mass only
    } else {                                      df$risk_groups[i] = 3 } # Multimorbid
  }
  # # Let's first factor that bad boy 
  df$risk_groups = factor(df$risk_groups, 
                          levels = c(0:3), 
                          labels = c("healthy", "internalizing_only", "cardiometabolic_only", "multimorbid"))
  message(paste("\nCombining:", int_var, "and", fm_var, "\n"))
  print(summary(df$risk_groups))
  corz = cor(df[, c(int_var, fm_var)], use = 'complete.obs')
  plot(df[, int_var], df[, fm_var], main = paste("Corr =", round(corz[1,2], 2)), xlab = int_var, ylab = fm_var, 
       col = c("darkgreen", "blue", "darkgoldenrod2", "red")[df$risk_groups])
  if (permute == T) {
    count <- 0 
    iterations <- 1000
    set.seed(310896)
    for (i in 1:iterations) {
      origN <- unname(summary(df$risk_groups)[4])
      randN <- permute(df)
      if (randN > origN) {
        count <- count + 1
      }
    }
    pval = format(round(count / iterations, 10), nsmall = 3)
    if (pval == 0.000) { pval = "< .001"}
    cat("Permutation p-value:", pval)
  }
  return(df$risk_groups)
}

#-------------------------------------------------------------------------------
# Organize variable names into domains to specify them later more easily
pre_LE <- c('family_member_died_pre',
            'friend_relative_died_pre',
            'family_member_ill_pre', 
            'friend_relative_ill_pre',
            'sick_or_accident_pre',
            'moved_pre',
            'blood_loss',
            'pregnancy_worried',
            'baby_worried',
            'burglary_or_car_theft_pre',
            'work_problems_pre',
            'abortion_pre',
            'married_pre',
            'unemployed_pre')
pre_CR <- c('income_reduced_pre',
            'homeless_pregnancy',
            'major_financial_problems_pre',
            'housing_adequacy_pre',
            'housing_basic_living_pre',
            'housing_defects_pre',
            'm_education_pre',
            'p_education_pre')
pre_PR <- c('m_criminal_record_pre',   #  without 'early_pregnancy': same variable already in postnatal PR score
            'p_criminal_record_pre',
            'm_attempted_suicide_pre',
            'm_depression_pre',
            'm_anxiety_pre',
            'm_interpersonal_sensitivity_pre',
            'p_depression_pre',
            'p_anxiety_pre',
            'p_interpersonal_sensitivity_pre') 
pre_IR <- c('divorce_pre',
            'p_rejected_child_pre',
            'p_went_away_pre',
            'conflict_in_family_pre',
            'argued_fam_friends_pre',
            'conflict_family_violence_pre',
            'marital_status_pregnancy',
            'family_affection',
            'family_size_pregnancy',
            'family_problems',
            'family_support',
            'social_network_emotional',
            'social_network_practical')
post_LE <- c('sick_or_accident_18m','sick_or_accident_30m','sick_or_accident_3y','sick_or_accident_4y','sick_or_accident_5y','sick_or_accident_6y','sick_or_accident_9y',
             'family_member_ill_8m', 'family_member_ill_21m', 'family_member_ill_3y', 'family_member_ill_4y', 'family_member_ill_5y','family_member_ill_6y','family_member_ill_9y',
             'smbd_important_died_8m','smbd_important_died_21m','smbd_important_died_3y','smbd_important_died_4y','smbd_important_died_5y','smbd_important_died_6y','smbd_important_died_9y',
             'separated_from_parent_18m','separated_from_parent_30m','separated_from_parent_3y','separated_from_parent_4y','separated_from_parent_5y','separated_from_parent_6y','separated_from_parent_8y',
             'moved_18m','moved_30m','moved_3y','moved_4y','moved_5y','moved_6y','moved_9y',
             'pet_died_18m','pet_died_30m','pet_died_3y','pet_died_4y','pet_died_5y','pet_died_6y','pet_died_9y',
             'started_nursery_18m', 'started_nursery_30m', 'started_nursery_3y', 'started_nursery_4y', 'started_nursery_5y', 
             'acquired_new_parent_18m', 'acquired_new_parent_30m', 'acquired_new_parent_3y', 'acquired_new_parent_4y', 'acquired_new_parent_5y', 'acquired_new_parent_6y', 'acquired_new_parent_8y',
             'change_carer_18m', 'change_carer_30m', 'change_carer_3y', 'change_carer_4y', 'change_carer_5y', 'change_carer_6y', 'change_carer_8y',
             'friend_relative_ill_8m', 'friend_relative_ill_21m', 'friend_relative_ill_3y', 'friend_relative_ill_4y', 'friend_relative_ill_5y', 'friend_relative_ill_6y', 'friend_relative_ill_9y',
             'partner_died_8m', 'partner_died_21m', 'partner_died_3y', 'partner_died_4y', 'partner_died_5y', 'partner_died_6y', 'partner_died_9y',
             'burglary_or_car_theft_21m', 'burglary_or_car_theft_3y', 'burglary_or_car_theft_4y', 'burglary_or_car_theft_5y', 'burglary_or_car_theft_6y', 'burglary_or_car_theft_9y',
             'separated_from_smbd_18m', 'separated_from_smbd_30m', 'separated_from_smbd_3y', 'separated_from_smbd_4y', 'separated_from_smbd_5y', 'separated_from_smbd_6y', 'separated_from_smbd_8y', 
             'lost_best_friend_8y',
             'new_sibling_18m', 'new_sibling_30m', 'new_sibling_3y', 'new_sibling_4y', 'new_sibling_5y', 'new_sibling_6y', 'new_sibling_9y',
             'ch_had_fright_18m', 'ch_had_fright_30m', 'ch_had_fright_3y', 'ch_had_fright_4y', 'ch_had_fright_5y', 'ch_had_fright_6y', 'ch_had_fright_8y')
post_CR <- c('homeless_childhood_8m', 'homeless_childhood_21m', 'homeless_childhood_3y', 'homeless_childhood_4y', 'homeless_childhood_5y', 'homeless_childhood_6y', 'homeless_childhood_9y',
             'major_financial_problems_8m', 'major_financial_problems_21m', 'major_financial_problems_3y', 'major_financial_problems_4y', 'major_financial_problems_5y', 'major_financial_problems_6y', 'major_financial_problems_9y',
             'income_reduced_8m', 'income_reduced_21m', 'income_reduced_3y', 'income_reduced_4y', 'income_reduced_5y', 'income_reduced_6y', 'income_reduced_9y',
             'unemployed_8m','unemployed_21m', 'unemployed_3y', 'unemployed_4y', 'unemployed_5y', 'unemployed_6y', 'unemployed_9y',
             'housing_adequacy_2y', 'housing_adequacy_4y',
             'housing_basic_living_2y', 'housing_basic_living_4y',
             'housing_defects_2y', 'housing_defects_4y',
             'm_education',
             'p_education',
             'neighbourhood_problems_21m', 'neighbourhood_problems_3y')
post_PR <- c('work_problems_8m','work_problems_21m','work_problems_3y','work_problems_4y','work_problems_5y','work_problems_6y','work_problems_9y',
             'criminal_record_parent_8m', 'criminal_record_parent_21m', 'criminal_record_parent_3y', 'criminal_record_parent_4y', 'criminal_record_parent_5y', 'criminal_record_parent_6y', 'criminal_record_parent_9y',
             'miscarriage_or_abortion_8m', 'miscarriage_or_abortion_21m', 'miscarriage_or_abortion_3y', 'miscarriage_or_abortion_4y', 'miscarriage_or_abortion_5y', 'miscarriage_or_abortion_6y', 'miscarriage_or_abortion_9y',
             'm_attempted_suicide_8m', 'm_attempted_suicide_21m', 'm_attempted_suicide_3y', 'm_attempted_suicide_4y', 'm_attempted_suicide_5y', 'm_attempted_suicide_6y', 'm_attempted_suicide_9y',
             'm_age',
             'p_age',
             'm_depression_8m', 'm_depression_21m', 'm_depression_3y', 'm_depression_4y', 'm_depression_5y', 'm_depression_6y', 'm_depression_9y',
             'p_depression_8m', 'p_depression_21m', 'p_depression_3y', 'p_depression_4y', 'p_depression_5y', 'p_depression_6y', 'p_depression_9y',
             'm_anxiety_8m', 'm_anxiety_21m', 'm_anxiety_3y', 'm_anxiety_5y', 'm_anxiety_6y',
             'p_anxiety_8m', 'p_anxiety_21m', 'p_anxiety_3y', 'p_anxiety_4y', 'p_anxiety_5y', 'p_anxiety_6y', 'p_anxiety_9y')
post_IR <- c( 'divorce_8m', 'divorce_21m', 'divorce_3y', 'divorce_4y', 'divorce_5y', 'divorce_6y', 'divorce_9y',
              'p_rejected_child_8m', 'p_rejected_child_21m', 'p_rejected_child_3y', 'p_rejected_child_4y', 'p_rejected_child_5y', 'p_rejected_child_6y', 'p_rejected_child_9y',
              'p_went_away_8m', 'p_went_away_21m', 'p_went_away_3y', 'p_went_away_4y', 'p_went_away_5y', 'p_went_away_6y', 'p_went_away_9y',
              'conflict_in_family_21m', 'conflict_in_family_3y', 'conflict_in_family_4y', 'conflict_in_family_5y', 'conflict_in_family_6y', 'conflict_in_family_9y',
              'conflict_family_violence_8m', 'conflict_family_violence_21m', 'conflict_family_violence_3y', 'conflict_family_violence_4y', 'conflict_family_violence_5y', 'conflict_family_violence_6y', 'conflict_family_violence_9y',
              'm_new_partner_8m', 'm_new_partner_21m', 'm_new_partner_3y', 'm_new_partner_4y', 'm_new_partner_5y', 'm_new_partner_6y', 'm_new_partner_9y',
              'argued_fam_friends_8m', 'argued_fam_friends_21m', 'argued_fam_friends_3y', 'argued_fam_friends_4y', 'argued_fam_friends_5y', 'argued_fam_friends_6y', 'argued_fam_friends_9y')
post_DV <- c('bullying_8y',
             'physical_violence_18m', 'physical_violence_30m', 'physical_violence_3y', 'physical_violence_4y', 'physical_violence_5y', 'physical_violence_6y', 'physical_violence_8y',
             'sexual_abuse_18m', 'sexual_abuse_30m', 'sexual_abuse_3y', 'sexual_abuse_4y', 'sexual_abuse_5y', 'sexual_abuse_6y', 'sexual_abuse_8y',
             'p_cruelty_physical_8m', 'p_cruelty_physical_21m', 'p_cruelty_physical_3y', 'p_cruelty_physical_4y', 'p_cruelty_physical_5y', 'p_cruelty_physical_6y', 'p_cruelty_physical_9y',
             'm_cruelty_physical_8m', 'm_cruelty_physical_21m', 'm_cruelty_physical_3y', 'm_cruelty_physical_4y', 'm_cruelty_physical_5y', 'm_cruelty_physical_6y', 'm_cruelty_physical_9y',
             'p_cruelty_emotional_8m', 'p_cruelty_emotional_21m', 'p_cruelty_emotional_3y', 'p_cruelty_emotional_4y', 'p_cruelty_emotional_5y', 'p_cruelty_emotional_6y', 'p_cruelty_emotional_9y',
             'm_cruelty_emotional_21m', 'm_cruelty_emotional_3y', 'm_cruelty_emotional_4y', 'm_cruelty_emotional_5y', 'm_cruelty_emotional_6y', 'm_cruelty_emotional_9y')
domains <- c('pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk', 
             'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')
outcomes_13y <- c('intern_score_13', 'total_fat_13', 'tot_fat_percent_13', 'andr_fat_mass_13')
outcomes_10y <- c('intern_score_10', 'total_fat_10', 'tot_fat_percent_10')
risk_grps <-c("risk_groups_tot", "risk_groups_andr", "risk_groups_perc", "risk_groups_tot_REC", "risk_groups_andr_REC", "risk_groups_perc_REC")
covars   <- c('sex', 'age_child', 'm_bmi_before_pregnancy', 'm_smoking', 'm_drinking')
auxil    <- c('m_dep_cont_pregnancy', 'p_dep_cont_pregnancy', # for postnatal only
              'm_bmi_7yrs', 'm_dep_cont_childhood', 'p_dep_cont_childhood', # for prenatal only
              'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_age_cont')
exclusion_criteria <- c('pre_percent_missing', 'post_percent_missing', 'twin', 'sibling')

################################################################################
################################################################################

# Load the dataframes with prenatal, postnatal stress and the outcome dataset
pre  <- readRDS(file.path(PATH_RESULTS, "prenatal_stress.rds"))
post <- readRDS(file.path(PATH_RESULTS, "postnatal_stress.rds"))
out  <- readRDS(file.path(PATH_RESULTS, "PCMout_cov_aux.rds"))

ELS <- cbind(pre, post, out)

# For the sake of time efficiency (and my mental health) let's select only those 
# variables that are needed for imputation nd subsequent sample selection. Once I 
# am at it, I also order them by domain. This is important because mice is sensitive
# to the order of the variables in the set (even though this may be a version-specific issue)
vars <- c('IDC', 
          # all variables for prenatal risk
          pre_LE, pre_CR, pre_PR, pre_IR,
          # all variables for postnatal risk
          post_LE, post_CR, post_PR, post_IR, post_DV,
          # all domain scores
          domains,
          # cumulative prenatal and postnatal stress exposure
          'prenatal_stress', "postnatal_stress",
          # outcome variables and covariates + additional auxiliary variables for imputation
          outcomes_13y, outcomes_10y, covars, auxil, exclusion_criteria)

ELS <- ELS[, vars]

# check out the sample
flowchart(ELS)

################################################################################
##---------------------------- Imputation model ----------------------------- ##
################################################################################

# We started with a dry run to specify the default arguments.
imp0 <- mice(ELS, maxit = 0, defaultMethod = rep('pmm', ncol(ELS)))
# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- imp0$method
# We use passive imputation for the domain scores. This means that the indicator items  
# are imputed first, and then, using these complete items, mean domain scores are 
# derived by the formula specified below.
meth['pre_life_events']         <- passive_imp_formula(pre_LE)
meth['pre_contextual_risk']     <- passive_imp_formula(pre_CR)
meth['pre_parental_risk']       <- passive_imp_formula(pre_PR)
meth['pre_interpersonal_risk']  <- passive_imp_formula(pre_IR)
meth['post_life_events']          <- passive_imp_formula(post_LE, type_nr = 17)
meth['post_contextual_risk']      <- passive_imp_formula(post_CR, type_nr = 10)
meth['post_parental_risk']        <- passive_imp_formula(post_PR, type_nr = 10)
meth['post_interpersonal_risk']   <- passive_imp_formula(post_IR, type_nr = 7)
meth['post_direct_victimization'] <- passive_imp_formula(post_DV, type_nr = 7)

# We also use passive imputation for the period specific cumulative ELS scores.
meth['prenatal_stress'] <- "~I( pre_life_events + pre_contextual_risk + pre_parental_risk + pre_interpersonal_risk )"
meth['postnatal_stress'] <- "~I( post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization )"

# We are going to need a different set of predictors for the different variables we impute 
# so let's define them using the predictormatrix, that gives these instructions to mice
predictormatrix <- imp0$predictorMatrix

# Do not impute nor use IDC, any of the exclusion criteria
predictormatrix[, c("IDC", exclusion_criteria) ] <- 
  predictormatrix[c("IDC", exclusion_criteria),] <- 0
# Leave the domain and total ELS scores to the passive imputation
predictormatrix[c(domains, 'prenatal_stress', 'postnatal_stress'), ] <- 0
# Do not use cumulative pre and postnatal scores as predictors in the items imputation
# nor the outocme timepoint
predictormatrix[c(pre_LE, pre_CR, pre_PR, pre_IR,# all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV), # all variables for postnatal risk
                c('prenatal_stress', 'postnatal_stress', outcomes_13y, outcomes_10y)]  <- 0

               ### Impute auxiliary variables and covariates ###
# To prevent multicollinearity, auxiliary variables are imputed given the domain scores 
# and the outcomes, but not the single items.
predictormatrix[c(outcomes_10y, covars, auxil),
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  'prenatal_stress', 'postnatal_stress',
                  outcomes_10y, covars, auxil)] <- 0

                            ### Impute outcomes ###
# Each primary outcome is imputed based on total ELS scores, auxiliary variables, 
# including the previous assessment of the same outcome, and the covariates. 
predictormatrix['intern_score_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains, outcomes_13y, outcomes_10y[!outcomes_10y == 'intern_score_10'])] <- 0
predictormatrix['total_fat_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains, outcomes_13y, outcomes_10y[!outcomes_10y == 'total_fat_10'])] <- 0
predictormatrix['andr_fat_mass_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains, outcomes_13y, outcomes_10y[!outcomes_10y == 'total_fat_10'])] <- 0
predictormatrix['tot_fat_percent_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains, outcomes_13y, outcomes_10y[!outcomes_10y == 'tot_fat_percent_10'])] <- 0

# All single items are then imputed such that we impute the items within a domain score given 
# the other items in that domain, the remaining domain scores and the auxiliary variables.
# Auxiliary variables were selected because they are either related to missingness or to 
# the domain scores themselves. When information is available prenatally and postnatally, 
# we use the different period to minimize bias (E.G., for the imputation of prenatal items
# we used BMI of the mother when the child was 5, and for imputation of postnatal items
# we used BMI of the mother during pregnancy).
# This method is more efficient (lower computational load) and the multicollinearity issue is 
# also resolved. The technique has been found to reduce standard error substantially compared to 
# complete-case analysis (Plumpton et al., 2010), and it outperforms other existing techniques (Eekhout et al., 2018). 

# So, let's adjust the predictormatrix such that ‘redundant’ items were not used as a predictor.

                                 ### PRENATAL ###
# LE domain 
predictormatrix[c(pre_LE), # LE
                c(pre_CR, pre_PR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # m_education is auxiliary for prenatal variables
                  post_PR, post_IR[-grep('divorce', post_IR)], post_DV,                # divorce is auxiliary for prenatal variables 
                  'pre_life_events', 'm_bmi_before_pregnancy', auxil[1:2])] <- 0
# CR domain
predictormatrix[c(pre_CR),
                c(pre_LE, pre_PR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # m_education is auxiliary for prenatal variables
                  post_PR, post_IR[-grep('divorce', post_IR)], post_DV,                # divorce is auxiliary for prenatal variables
                  'pre_contextual_risk', 'm_bmi_before_pregnancy', auxil[1:2])] <- 0
# PR domain 
predictormatrix[c(pre_PR),
                c(pre_LE, pre_CR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # m_education is auxiliary for prenatal variables
                  post_PR, post_IR[-grep('divorce', post_IR)], post_DV,                # divorce is auxiliary for prenatal variables
                  'pre_parental_risk', 'm_bmi_before_pregnancy', auxil[1:2])] <- 0
# IR domain
predictormatrix[c(pre_IR),
                c(pre_LE, pre_CR, pre_PR, post_LE, post_CR[!post_CR == 'm_education'], # m_education is auxiliary for prenatal variables
                  post_PR, post_IR[-grep('divorce', post_IR)], post_DV,                # divorce is auxiliary for prenatal variables
                  'pre_interpersonal_risk', 'm_bmi_before_pregnancy', auxil[1:2])] <- 0

                                ### POSTNATAL ###
# LE domain 
predictormatrix[c(post_LE),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],          # m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], # marital_status_pregnancy is auxiliary for postnatal variables
                  post_CR, post_PR, post_IR, post_DV,
                  'post_life_events', 'm_bmi_before_pregnancy', auxil[3:5])] <- 0       
# CR domain
predictormatrix[c(post_CR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],          # m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], # marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_PR, post_IR, post_DV,
                  'post_contextual_risk', 'm_bmi_before_pregnancy', auxil[3:5])] <- 0
# PR domain 
predictormatrix[c(post_PR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],          # m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], # marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_IR, post_DV,
                  'post_parental_risk', 'm_bmi_before_pregnancy', auxil[3:5])] <- 0
# IR domain
predictormatrix[c(post_IR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],          # m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], # marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_PR, post_DV,
                  'post_interpersonal_risk', 'm_bmi_before_pregnancy', auxil[3:5])] <- 0
# DV domain 
predictormatrix[c(post_DV),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],          # m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], # marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_PR, post_IR,
                  'post_direct_victimization', 'm_bmi_before_pregnancy', auxil[3:5])] <- 0

# OPTIONAL : Quickly check the matrix to make sure it looks legit
# pheatmap::pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain low computational load,
# we do 60 iterations using 30 imputed datasets (following Isabel's approach)
imp <- mice(ELS, m = 30, # nr of imputed datasets
            maxit = 60, #nr of iteration taken to impute missing values
            seed = 310896, # set a seed for the random number generation in case i need to generate the same dataset again
            method = meth,
            visitSequence = VisSeq, 
            predictorMatrix = predictormatrix)

################################################################################
################################################################################

# Apply exclusion criteria to select the sample for analysis: exclude participants that do not
# have enough prenatal and postnatal information (i.e. < 50% of the items in each period).
pren50cutoff <- miceadds::subset_datlist( imp, subset = imp$data$pre_percent_missing < 50.0,  toclass = 'mids')
post50cutoff <- miceadds::subset_datlist( pren50cutoff, subset = pren50cutoff$data$post_percent_missing < 50.0 )
# Not specifying toclass argument in the last call transforms mids object into a datalist object.

# Standardize prenatal, postnatal stress and the primary outcomes to obtain standard betas from regression
sdatlist <- miceadds::scale_datlist(post50cutoff, orig_var = c('prenatal_stress', 'postnatal_stress', 'intern_score_13',
                                                               'total_fat_13', 'andr_fat_mass_13', 'tot_fat_percent_13'), 
                                    trafo_var = paste0( c('prenatal_stress', 'postnatal_stress', 'intern_score_13', 
                                                          'total_fat_13', 'andr_fat_mass_13', 'tot_fat_percent_13'), "_z") )
# Reconvert back to mids object
post50cutoff <- miceadds::datlist2mids(sdatlist)

no_twins <- miceadds::subset_datlist( post50cutoff, subset = post50cutoff$data$twin == 0,  toclass = 'mids') 
no_sibls <- miceadds::subset_datlist( no_twins, subset = !is.na() & no_twins$data$sibling == 0, toclass = 'mids')

################################################################################
# Finally, we need to construct and add the risk group variables in  each imputed dataset. 
# to do so we will use the function construct_grp defined in 0-Functions.R
# We will be doing this only on the selected sample to save some time. The same code can 
# be adapted to use the pre-selection "imputation" variable. 

# convert imputations into a long dataframe
long <- mice::complete(no_sibls, action = "long", include = TRUE) 
# initialize empty long dataframe for return
return_dat <- data.frame() 
#  apply compute_grp to each imputed dataset
for (n_imp in unique(long$.imp)) {
  dset = long[long$.imp == n_imp, ] # one iteration only
  message(paste("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Imputation nr: ", n_imp, 
                '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'))
  # Contruct the risk group varible for each fat measure
  dset$risk_groups_tot  <- construct_grp('intern_score_13', 'total_fat_13', dset)
  dset$risk_groups_perc <- construct_grp('intern_score_13', 'tot_fat_percent_13', dset)
  dset$risk_groups_andr <- construct_grp('intern_score_13', 'andr_fat_mass_13', dset)
  # For later analyses also save a version where multimorbid is the reference. 
  dset$risk_groups_tot_REC  <- relevel(dset$risk_groups_tot, 'multimorbid')
  dset$risk_groups_perc_REC <- relevel(dset$risk_groups_perc, 'multimorbid')
  dset$risk_groups_andr_REC <- relevel(dset$risk_groups_andr, 'multimorbid')
  
  # Append results to long-type dataframe
  return_dat <- rbind.data.frame(return_dat, dset) 
}
finalset <- mice::as.mids(return_dat) # reconvert long dataframe to mids object

################################################################################
#### ------------------------- complete and save -------------------------- ####
################################################################################

# I save the mids object (i.e. list of imputed datasets)
saveRDS(imp, file.path(PATH_RESULTS,'imputation_list_full.rds'))
saveRDS(finalset, file.path(PATH_RESULTS,'imputation_list_sample.rds'))

# I also save the last imputed dataset for sanity checks
full_imputed <- complete(imp, 30) 
ELSPCM_imputed <- complete(finalset, 30) 
ELSPCM_original <- complete(finalset, 0) 
saveRDS(full_imputed,   file.path(PATH_RESULTS, 'imputed30_full.rds'))
saveRDS(ELSPCM_imputed, file.path(PATH_RESULTS, 'imputed30_sample.rds'))
write.csv(ELSPCM_original, file.path(PATH_RESULTS, 'original_sample.csv'))


################################################################################
