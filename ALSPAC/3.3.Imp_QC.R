# As usual, here are the packages we need 
library(mice);

# Load the list of imputed dataset
if (exists("imp_samp") == F) { 
  imp_path <- file.choose() # ATTENTION! Choose the 'imputation_list_sample.rds' file
  dirname <- dirname(imp_path)
  ifelse(!dir.exists(file.path(dirname, 'imp-QC')), dir.create(file.path(dirname, 'imp-QC')), F)
  imp_samp <- readRDS(imp_path)
  imp_full <- readRDS(file.path(dirname, "imputation_list_full.rds"))
}

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
# Assess the validity of a mean solution vs the single imputations. 
################################################################################
# for dichotomous variables, replace NA with the median across imputations
# for continuous (numeric) variables, replace NA with the mean across imputations
replance_mean <- function(orig) {
  for (var in colnames(orig)) {
    misslist <- orig$IDC[which(is.na(orig[, var]))]
    cat(paste(var, length(misslist), sep = " - "))
    cat('\n')
    if (var %in% c(pre_LE, pre_CR, pre_PR, pre_IR, post_LE, post_CR, post_PR, post_IR, post_DV, 
                   'sex', 'ethnicity', risk_grps)) {
      pool_mean <- with(impdat, by(impdat, IDC, function(x) c(median(x[,var]), sd(x[,var]))))
    } else { 
      pool_mean <- with(impdat, by(impdat, IDC, function(x) c(mean(x[,var]), sd(x[,var]))))
    }
    
    for (cld in orig$IDC) {
      child_m = pool_mean[dimnames(pool_mean)$IDC == cld]
      if (child_m[[1]][2] > 0.000) {
        orig[orig$IDC == cld, var] <- child_m[[1]][1] }
    }
  } 
  return(orig)
}
# ------------------------------------------------------------------------------
# Retrieve original (pre-imputation) dataset
original_set <- imp_samp$data
# Stack imputed datasets in long format, exclude the original data
impdat <- complete(imp_samp, action="long", include = FALSE)

# Resolve problem with categorical variables 
for (cat_var in c('sex', 'ethnicity', risk_grps)) {
  original_set[, cat_var] <- as.numeric(original_set[, cat_var])
  impdat[, cat_var] <- as.numeric(impdat[, cat_var])
}

# Compute pre imp means and sds
orig_means <- sapply(original_set[,-1], mean, na.rm = T, USE.NAMES = T)
orig_sds <- sapply(original_set[,-1], sd, na.rm = T, USE.NAMES = T)

# Create mean of imputation dataset
# BEWARE: with ALSPAC dataset, this next step takes ~ 7 hours
dataset <- replance_mean(original_set)

best = data.frame(1:31)
for (var in colnames(original_set)) {
  m = orig_means[var]
  m_avg = mean(dataset[, var], na.rm = T)
  pool_imp <- with(impdat, by(impdat, .imp, function(x) abs(mean(as.numeric(x[, var])) - m)))
  pool_imp[['31']] <- abs(m_avg - m)
  best[, var] = pool_imp
}
best = best[, -1]
cols <- rainbow(ncol(best))
i = 1

pdf(file.path(dirname, "imp-QC", "mean-diff-iterations.pdf"))
for (var in colnames(original_set)) {
  plot(best[,var], type = 'l', col = 'red', lwd=3, 
       xlim = c(0, 32), ylim = c(0, 0.5), main = paste(var, '- NA:', sum(is.na(original_set[, var]))),  
       xlab="Iteration", ylab="Absolute mean difference", xaxt="n")
  axis(side = 1, at = seq(1, 31, by = 1), las=2)
  segments(1, best[1, var], 30, best[30, var], col = 'black')
  i = i+1
}
dev.off()

# EVALUATE THE CONVERGENCE OF THE IMPUTATION
# Extract the means of the imputed values for each iteration # imp$chainMean
# The number of chains is equal to the number of imputed datasets. A chain refers 
# to the chain of regression models that is used to generate the imputed values. 
# The length of each chain is equal to the number of iterations.

# The convergence can be visualised by plotting the means in a convergence plot. 
pdf(file.path(dirname, "imp-QC", "covergence-plots-full.pdf")) # this works only for full imp
plot(imp_full)
dev.off()
# The variance between the imputation chains should be equal to the variance within 
# the chains, which indicates healthy convergence.

# COMPARE IMPUTED AND OBSERVED VALUES
# Next, we inspect and compare the density of the incomplete and imputed data
# For loops do not work with this wacky mice function so just run the loop below 
# and paste it: lil help
# for (var in c(pre_LE, pre_CR, pre_PR, pre_IR, post_LE, post_CR, post_PR, post_IR, post_DV, covars, auxil)) {
#    cat(paste('densityplot(imp, ~', var, ')', '\n')) }

pdf(file.path(dirname,"imp-QC","imp-vs-obs.pdf"))
densityplot(imp, ~ family_member_died_pre ) 
densityplot(imp, ~ friend_relative_died_pre ) 
densityplot(imp, ~ family_member_ill_pre ) 
densityplot(imp, ~ friend_relative_ill_pre ) 
densityplot(imp, ~ sick_or_accident_pre ) 
densityplot(imp, ~ moved_pre ) 
densityplot(imp, ~ blood_loss ) 
densityplot(imp, ~ pregnancy_worried ) 
densityplot(imp, ~ baby_worried ) 
densityplot(imp, ~ burglary_or_car_theft_pre ) 
densityplot(imp, ~ work_problems_pre ) 
densityplot(imp, ~ abortion_pre ) 
densityplot(imp, ~ married_pre ) 
densityplot(imp, ~ unemployed_pre ) 
densityplot(imp, ~ income_reduced_pre ) 
densityplot(imp, ~ homeless_pregnancy ) 
densityplot(imp, ~ major_financial_problems_pre ) 
densityplot(imp, ~ housing_adequacy_pre ) 
densityplot(imp, ~ housing_basic_living_pre ) 
densityplot(imp, ~ housing_defects_pre ) 
densityplot(imp, ~ m_education_pre ) 
densityplot(imp, ~ p_education_pre ) 
densityplot(imp, ~ m_criminal_record_pre ) 
densityplot(imp, ~ p_criminal_record_pre ) 
densityplot(imp, ~ m_attempted_suicide_pre ) 
densityplot(imp, ~ m_depression_pre ) 
densityplot(imp, ~ m_anxiety_pre ) 
densityplot(imp, ~ m_interpersonal_sensitivity_pre ) 
densityplot(imp, ~ p_depression_pre ) 
densityplot(imp, ~ p_anxiety_pre ) 
densityplot(imp, ~ p_interpersonal_sensitivity_pre ) 
densityplot(imp, ~ divorce_pre ) 
densityplot(imp, ~ p_rejected_child_pre ) 
densityplot(imp, ~ p_went_away_pre ) 
densityplot(imp, ~ conflict_in_family_pre ) 
densityplot(imp, ~ argued_fam_friends_pre ) 
densityplot(imp, ~ conflict_family_violence_pre ) 
densityplot(imp, ~ marital_status_pregnancy ) 
densityplot(imp, ~ family_affection ) 
densityplot(imp, ~ family_size_pregnancy ) 
densityplot(imp, ~ family_problems ) 
densityplot(imp, ~ family_support ) 
densityplot(imp, ~ social_network_emotional ) 
densityplot(imp, ~ social_network_practical ) 
densityplot(imp, ~ sick_or_accident_18m ) 
densityplot(imp, ~ sick_or_accident_30m ) 
densityplot(imp, ~ sick_or_accident_3y ) 
densityplot(imp, ~ sick_or_accident_4y ) 
densityplot(imp, ~ sick_or_accident_5y ) 
densityplot(imp, ~ sick_or_accident_6y ) 
densityplot(imp, ~ sick_or_accident_9y ) 
densityplot(imp, ~ family_member_ill_8m ) 
densityplot(imp, ~ family_member_ill_21m ) 
densityplot(imp, ~ family_member_ill_3y ) 
densityplot(imp, ~ family_member_ill_4y ) 
densityplot(imp, ~ family_member_ill_5y ) 
densityplot(imp, ~ family_member_ill_6y ) 
densityplot(imp, ~ family_member_ill_9y ) 
densityplot(imp, ~ smbd_important_died_8m ) 
densityplot(imp, ~ smbd_important_died_21m ) 
densityplot(imp, ~ smbd_important_died_3y ) 
densityplot(imp, ~ smbd_important_died_4y ) 
densityplot(imp, ~ smbd_important_died_5y ) 
densityplot(imp, ~ smbd_important_died_6y ) 
densityplot(imp, ~ smbd_important_died_9y ) 
densityplot(imp, ~ separated_from_parent_18m ) 
densityplot(imp, ~ separated_from_parent_30m ) 
densityplot(imp, ~ separated_from_parent_3y ) 
densityplot(imp, ~ separated_from_parent_4y ) 
densityplot(imp, ~ separated_from_parent_5y ) 
densityplot(imp, ~ separated_from_parent_6y ) 
densityplot(imp, ~ separated_from_parent_8y ) 
densityplot(imp, ~ moved_18m ) 
densityplot(imp, ~ moved_30m ) 
densityplot(imp, ~ moved_3y ) 
densityplot(imp, ~ moved_4y ) 
densityplot(imp, ~ moved_5y ) 
densityplot(imp, ~ moved_6y ) 
densityplot(imp, ~ moved_9y ) 
densityplot(imp, ~ pet_died_18m ) 
densityplot(imp, ~ pet_died_30m ) 
densityplot(imp, ~ pet_died_3y ) 
densityplot(imp, ~ pet_died_4y ) 
densityplot(imp, ~ pet_died_5y ) 
densityplot(imp, ~ pet_died_6y ) 
densityplot(imp, ~ pet_died_9y ) 
densityplot(imp, ~ started_nursery_18m ) 
densityplot(imp, ~ started_nursery_30m ) 
densityplot(imp, ~ started_nursery_3y ) 
densityplot(imp, ~ started_nursery_4y ) 
densityplot(imp, ~ started_nursery_5y ) 
densityplot(imp, ~ acquired_new_parent_18m ) 
densityplot(imp, ~ acquired_new_parent_30m ) 
densityplot(imp, ~ acquired_new_parent_3y ) 
densityplot(imp, ~ acquired_new_parent_4y ) 
densityplot(imp, ~ acquired_new_parent_5y ) 
densityplot(imp, ~ acquired_new_parent_6y ) 
densityplot(imp, ~ acquired_new_parent_8y ) 
densityplot(imp, ~ change_carer_18m ) 
densityplot(imp, ~ change_carer_30m ) 
densityplot(imp, ~ change_carer_3y ) 
densityplot(imp, ~ change_carer_4y ) 
densityplot(imp, ~ change_carer_5y ) 
densityplot(imp, ~ change_carer_6y ) 
densityplot(imp, ~ change_carer_8y ) 
densityplot(imp, ~ friend_relative_ill_8m ) 
densityplot(imp, ~ friend_relative_ill_21m ) 
densityplot(imp, ~ friend_relative_ill_3y ) 
densityplot(imp, ~ friend_relative_ill_4y ) 
densityplot(imp, ~ friend_relative_ill_5y ) 
densityplot(imp, ~ friend_relative_ill_6y ) 
densityplot(imp, ~ friend_relative_ill_9y ) 
densityplot(imp, ~ partner_died_8m ) 
densityplot(imp, ~ partner_died_21m ) 
densityplot(imp, ~ partner_died_3y ) 
densityplot(imp, ~ partner_died_4y ) 
densityplot(imp, ~ partner_died_5y ) 
densityplot(imp, ~ partner_died_6y ) 
densityplot(imp, ~ partner_died_9y ) 
densityplot(imp, ~ burglary_or_car_theft_21m ) 
densityplot(imp, ~ burglary_or_car_theft_3y ) 
densityplot(imp, ~ burglary_or_car_theft_4y ) 
densityplot(imp, ~ burglary_or_car_theft_5y ) 
densityplot(imp, ~ burglary_or_car_theft_6y ) 
densityplot(imp, ~ burglary_or_car_theft_9y ) 
densityplot(imp, ~ separated_from_smbd_18m ) 
densityplot(imp, ~ separated_from_smbd_30m ) 
densityplot(imp, ~ separated_from_smbd_3y ) 
densityplot(imp, ~ separated_from_smbd_4y ) 
densityplot(imp, ~ separated_from_smbd_5y ) 
densityplot(imp, ~ separated_from_smbd_6y ) 
densityplot(imp, ~ separated_from_smbd_8y ) 
densityplot(imp, ~ lost_best_friend_8y ) 
densityplot(imp, ~ new_sibling_18m ) 
densityplot(imp, ~ new_sibling_30m ) 
densityplot(imp, ~ new_sibling_3y ) 
densityplot(imp, ~ new_sibling_4y ) 
densityplot(imp, ~ new_sibling_5y ) 
densityplot(imp, ~ new_sibling_6y ) 
densityplot(imp, ~ new_sibling_9y ) 
densityplot(imp, ~ ch_had_fright_18m ) 
densityplot(imp, ~ ch_had_fright_30m ) 
densityplot(imp, ~ ch_had_fright_3y ) 
densityplot(imp, ~ ch_had_fright_4y ) 
densityplot(imp, ~ ch_had_fright_5y ) 
densityplot(imp, ~ ch_had_fright_6y ) 
densityplot(imp, ~ ch_had_fright_8y ) 
densityplot(imp, ~ homeless_childhood_8m ) 
densityplot(imp, ~ homeless_childhood_21m ) 
densityplot(imp, ~ homeless_childhood_3y ) 
densityplot(imp, ~ homeless_childhood_4y ) 
densityplot(imp, ~ homeless_childhood_5y ) 
densityplot(imp, ~ homeless_childhood_6y ) 
densityplot(imp, ~ homeless_childhood_9y ) 
densityplot(imp, ~ major_financial_problems_8m ) 
densityplot(imp, ~ major_financial_problems_21m ) 
densityplot(imp, ~ major_financial_problems_3y ) 
densityplot(imp, ~ major_financial_problems_4y ) 
densityplot(imp, ~ major_financial_problems_5y ) 
densityplot(imp, ~ major_financial_problems_6y ) 
densityplot(imp, ~ major_financial_problems_9y ) 
densityplot(imp, ~ income_reduced_8m ) 
densityplot(imp, ~ income_reduced_21m ) 
densityplot(imp, ~ income_reduced_3y ) 
densityplot(imp, ~ income_reduced_4y ) 
densityplot(imp, ~ income_reduced_5y ) 
densityplot(imp, ~ income_reduced_6y ) 
densityplot(imp, ~ income_reduced_9y ) 
densityplot(imp, ~ unemployed_8m ) 
densityplot(imp, ~ unemployed_21m ) 
densityplot(imp, ~ unemployed_3y ) 
densityplot(imp, ~ unemployed_4y ) 
densityplot(imp, ~ unemployed_5y ) 
densityplot(imp, ~ unemployed_6y ) 
densityplot(imp, ~ unemployed_9y ) 
densityplot(imp, ~ housing_adequacy_2y ) 
densityplot(imp, ~ housing_adequacy_4y ) 
densityplot(imp, ~ housing_basic_living_2y ) 
densityplot(imp, ~ housing_basic_living_4y ) 
densityplot(imp, ~ housing_defects_2y ) 
densityplot(imp, ~ housing_defects_4y ) 
densityplot(imp, ~ m_education ) 
densityplot(imp, ~ p_education ) 
densityplot(imp, ~ neighbourhood_problems_21m ) 
densityplot(imp, ~ neighbourhood_problems_3y ) 
densityplot(imp, ~ work_problems_8m ) 
densityplot(imp, ~ work_problems_21m ) 
densityplot(imp, ~ work_problems_3y ) 
densityplot(imp, ~ work_problems_4y ) 
densityplot(imp, ~ work_problems_5y ) 
densityplot(imp, ~ work_problems_6y ) 
densityplot(imp, ~ work_problems_9y ) 
densityplot(imp, ~ criminal_record_parent_8m ) 
densityplot(imp, ~ criminal_record_parent_21m ) 
densityplot(imp, ~ criminal_record_parent_3y ) 
densityplot(imp, ~ criminal_record_parent_4y ) 
densityplot(imp, ~ criminal_record_parent_5y ) 
densityplot(imp, ~ criminal_record_parent_6y ) 
densityplot(imp, ~ criminal_record_parent_9y ) 
densityplot(imp, ~ miscarriage_or_abortion_8m ) 
densityplot(imp, ~ miscarriage_or_abortion_21m ) 
densityplot(imp, ~ miscarriage_or_abortion_3y ) 
densityplot(imp, ~ miscarriage_or_abortion_4y ) 
densityplot(imp, ~ miscarriage_or_abortion_5y ) 
densityplot(imp, ~ miscarriage_or_abortion_6y ) 
densityplot(imp, ~ miscarriage_or_abortion_9y ) 
densityplot(imp, ~ m_attempted_suicide_8m ) 
densityplot(imp, ~ m_attempted_suicide_21m ) 
densityplot(imp, ~ m_attempted_suicide_3y ) 
densityplot(imp, ~ m_attempted_suicide_4y ) 
densityplot(imp, ~ m_attempted_suicide_5y ) 
densityplot(imp, ~ m_attempted_suicide_6y ) 
densityplot(imp, ~ m_attempted_suicide_9y ) 
densityplot(imp, ~ m_age ) 
densityplot(imp, ~ p_age ) 
densityplot(imp, ~ m_depression_8m ) 
densityplot(imp, ~ m_depression_21m ) 
densityplot(imp, ~ m_depression_3y ) 
densityplot(imp, ~ m_depression_4y ) 
densityplot(imp, ~ m_depression_5y ) 
densityplot(imp, ~ m_depression_6y ) 
densityplot(imp, ~ m_depression_9y ) 
densityplot(imp, ~ p_depression_8m ) 
densityplot(imp, ~ p_depression_21m ) 
densityplot(imp, ~ p_depression_3y ) 
densityplot(imp, ~ p_depression_4y ) 
densityplot(imp, ~ p_depression_5y ) 
densityplot(imp, ~ p_depression_6y ) 
densityplot(imp, ~ p_depression_9y ) 
densityplot(imp, ~ m_anxiety_8m ) 
densityplot(imp, ~ m_anxiety_21m ) 
densityplot(imp, ~ m_anxiety_3y ) 
densityplot(imp, ~ m_anxiety_5y ) 
densityplot(imp, ~ m_anxiety_6y ) 
densityplot(imp, ~ p_anxiety_8m ) 
densityplot(imp, ~ p_anxiety_21m ) 
densityplot(imp, ~ p_anxiety_3y ) 
densityplot(imp, ~ p_anxiety_4y ) 
densityplot(imp, ~ p_anxiety_5y ) 
densityplot(imp, ~ p_anxiety_6y ) 
densityplot(imp, ~ p_anxiety_9y ) 
densityplot(imp, ~ divorce_8m ) 
densityplot(imp, ~ divorce_21m ) 
densityplot(imp, ~ divorce_3y ) 
densityplot(imp, ~ divorce_4y ) 
densityplot(imp, ~ divorce_5y ) 
densityplot(imp, ~ divorce_6y ) 
densityplot(imp, ~ divorce_9y ) 
densityplot(imp, ~ p_rejected_child_8m ) 
densityplot(imp, ~ p_rejected_child_21m ) 
densityplot(imp, ~ p_rejected_child_3y ) 
densityplot(imp, ~ p_rejected_child_4y ) 
densityplot(imp, ~ p_rejected_child_5y ) 
densityplot(imp, ~ p_rejected_child_6y ) 
densityplot(imp, ~ p_rejected_child_9y ) 
densityplot(imp, ~ p_went_away_8m ) 
densityplot(imp, ~ p_went_away_21m ) 
densityplot(imp, ~ p_went_away_3y ) 
densityplot(imp, ~ p_went_away_4y ) 
densityplot(imp, ~ p_went_away_5y ) 
densityplot(imp, ~ p_went_away_6y ) 
densityplot(imp, ~ p_went_away_9y ) 
densityplot(imp, ~ conflict_in_family_21m ) 
densityplot(imp, ~ conflict_in_family_3y ) 
densityplot(imp, ~ conflict_in_family_4y ) 
densityplot(imp, ~ conflict_in_family_5y ) 
densityplot(imp, ~ conflict_in_family_6y ) 
densityplot(imp, ~ conflict_in_family_9y ) 
densityplot(imp, ~ conflict_family_violence_8m ) 
densityplot(imp, ~ conflict_family_violence_21m ) 
densityplot(imp, ~ conflict_family_violence_3y ) 
densityplot(imp, ~ conflict_family_violence_4y ) 
densityplot(imp, ~ conflict_family_violence_5y ) 
densityplot(imp, ~ conflict_family_violence_6y ) 
densityplot(imp, ~ conflict_family_violence_9y ) 
densityplot(imp, ~ m_new_partner_8m ) 
densityplot(imp, ~ m_new_partner_21m ) 
densityplot(imp, ~ m_new_partner_3y ) 
densityplot(imp, ~ m_new_partner_4y ) 
densityplot(imp, ~ m_new_partner_5y ) 
densityplot(imp, ~ m_new_partner_6y ) 
densityplot(imp, ~ m_new_partner_9y ) 
densityplot(imp, ~ argued_fam_friends_8m ) 
densityplot(imp, ~ argued_fam_friends_21m ) 
densityplot(imp, ~ argued_fam_friends_3y ) 
densityplot(imp, ~ argued_fam_friends_4y ) 
densityplot(imp, ~ argued_fam_friends_5y ) 
densityplot(imp, ~ argued_fam_friends_6y ) 
densityplot(imp, ~ argued_fam_friends_9y ) 
densityplot(imp, ~ bullying_8y ) 
densityplot(imp, ~ physical_violence_18m ) 
densityplot(imp, ~ physical_violence_30m ) 
densityplot(imp, ~ physical_violence_3y ) 
densityplot(imp, ~ physical_violence_4y ) 
densityplot(imp, ~ physical_violence_5y ) 
densityplot(imp, ~ physical_violence_6y ) 
densityplot(imp, ~ physical_violence_8y ) 
densityplot(imp, ~ sexual_abuse_18m ) 
densityplot(imp, ~ sexual_abuse_30m ) 
densityplot(imp, ~ sexual_abuse_3y ) 
densityplot(imp, ~ sexual_abuse_4y ) 
densityplot(imp, ~ sexual_abuse_5y ) 
densityplot(imp, ~ sexual_abuse_6y ) 
densityplot(imp, ~ sexual_abuse_8y ) 
densityplot(imp, ~ p_cruelty_physical_8m ) 
densityplot(imp, ~ p_cruelty_physical_21m ) 
densityplot(imp, ~ p_cruelty_physical_3y ) 
densityplot(imp, ~ p_cruelty_physical_4y ) 
densityplot(imp, ~ p_cruelty_physical_5y ) 
densityplot(imp, ~ p_cruelty_physical_6y ) 
densityplot(imp, ~ p_cruelty_physical_9y ) 
densityplot(imp, ~ m_cruelty_physical_8m ) 
densityplot(imp, ~ m_cruelty_physical_21m ) 
densityplot(imp, ~ m_cruelty_physical_3y ) 
densityplot(imp, ~ m_cruelty_physical_4y ) 
densityplot(imp, ~ m_cruelty_physical_5y ) 
densityplot(imp, ~ m_cruelty_physical_6y ) 
densityplot(imp, ~ m_cruelty_physical_9y ) 
densityplot(imp, ~ p_cruelty_emotional_8m ) 
densityplot(imp, ~ p_cruelty_emotional_21m ) 
densityplot(imp, ~ p_cruelty_emotional_3y ) 
densityplot(imp, ~ p_cruelty_emotional_4y ) 
densityplot(imp, ~ p_cruelty_emotional_5y ) 
densityplot(imp, ~ p_cruelty_emotional_6y ) 
densityplot(imp, ~ p_cruelty_emotional_9y ) 
densityplot(imp, ~ m_cruelty_emotional_21m ) 
densityplot(imp, ~ m_cruelty_emotional_3y ) 
densityplot(imp, ~ m_cruelty_emotional_4y ) 
densityplot(imp, ~ m_cruelty_emotional_5y ) 
densityplot(imp, ~ m_cruelty_emotional_6y ) 
densityplot(imp, ~ m_cruelty_emotional_9y ) 

densityplot(imp_samp, ~ sex ) 
densityplot(imp_samp, ~ age_child ) 
densityplot(imp_samp, ~ m_bmi_before_pregnancy ) 
densityplot(imp_samp, ~ m_smoking ) 
densityplot(imp_samp, ~ m_drinking ) 
densityplot(imp_samp, ~ ethnicity ) 
densityplot(imp_samp, ~ prenatal_stress_z ) 
densityplot(imp_samp, ~ postnatal_stress_z ) 
densityplot(imp_samp, ~ intern_score_13_z )
densityplot(imp_samp, ~ tot_fat_percent_13_z )
densityplot(imp_samp, ~ andr_fat_mass_13_z ) 
densityplot(imp_samp, ~ total_fat_13_z )
densityplot(imp_samp, ~ risk_groups_perc ) 
densityplot(imp_samp, ~ risk_groups_andr ) 
densityplot(imp_samp, ~ risk_groups_tot )
densityplot(imp_samp, ~ risk_groups_perc_REC ) 
densityplot(imp_samp, ~ risk_groups_andr_REC )
densityplot(imp_samp, ~ risk_groups_tot_REC ) 
dev.off()

# Make a missing data indicator (name it miss) for domain scores and check the 
# relation in the imputed data. To do so, plot the imputed values against their 
# respective calculated values.

# same story: For loops do not work with this wacky mice function so just run the 
# loop below and paste it: lil help
# passive_imp_formula <- function(domain) {
#   conc <- paste(domain[[1]], collapse = " + ")
#   str <- paste0("~I (", conc, " / ", length(domain[[1]]))
#   return(str)
# }
# dnames <- c('pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk',
#             'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')
# ds <- list(pre_LE, pre_CR, pre_PR, pre_IR,
#         post_LE, post_CR, post_PR, post_IR, post_DV)
# 
# for (i in 1:9) {
#   cat(paste0('miss <- is.na(original_set[, "', dnames[i], '" ])', '\n',
#             'xyplot(imp, ', dnames[i], ' ', passive_imp_formula(ds[i]), '),', '\n',
#             'na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "',dnames[i], '") \n')) }

pdf(file.path(dirname,"imp-QC","main-&-domains.pdf"))
stripplot(imp_samp, prenatal_stress_z ) 
stripplot(imp_samp, postnatal_stress_z ) 
stripplot(imp_samp, intern_score_13_z ) 
stripplot(imp_samp, tot_fat_percent_13_z ) 
stripplot(imp_samp, andr_fat_mass_13_z ) 
stripplot(imp_samp, total_fat_13_z ) 
miss <- is.na(original_set[, "pre_life_events" ])
xyplot(imp, pre_life_events ~I (family_member_died_pre + friend_relative_died_pre + family_member_ill_pre + friend_relative_ill_pre + sick_or_accident_pre + moved_pre + blood_loss + pregnancy_worried + baby_worried + burglary_or_car_theft_pre + work_problems_pre + abortion_pre + married_pre + unemployed_pre / 14),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_life_events") 
miss <- is.na(original_set[, "pre_contextual_risk" ])
xyplot(imp, pre_contextual_risk ~I (income_reduced_pre + homeless_pregnancy + major_financial_problems_pre + housing_adequacy_pre + housing_basic_living_pre + housing_defects_pre + m_education_pre + p_education_pre / 8),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_contextual_risk") 
miss <- is.na(original_set[, "pre_parental_risk" ])
xyplot(imp, pre_parental_risk ~I (m_criminal_record_pre + p_criminal_record_pre + m_attempted_suicide_pre + m_depression_pre + m_anxiety_pre + m_interpersonal_sensitivity_pre + p_depression_pre + p_anxiety_pre + p_interpersonal_sensitivity_pre / 9),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_parental_risk") 
miss <- is.na(original_set[, "pre_interpersonal_risk" ])
xyplot(imp, pre_interpersonal_risk ~I (divorce_pre + p_rejected_child_pre + p_went_away_pre + conflict_in_family_pre + argued_fam_friends_pre + conflict_family_violence_pre + marital_status_pregnancy + family_affection + family_size_pregnancy + family_problems + family_support + social_network_emotional + social_network_practical / 13),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_interpersonal_risk") 
miss <- is.na(original_set[, "post_life_events" ])
xyplot(imp, post_life_events ~I (sick_or_accident_18m + sick_or_accident_30m + sick_or_accident_3y + sick_or_accident_4y + sick_or_accident_5y + sick_or_accident_6y + sick_or_accident_9y + family_member_ill_8m + family_member_ill_21m + family_member_ill_3y + family_member_ill_4y + family_member_ill_5y + family_member_ill_6y + family_member_ill_9y + smbd_important_died_8m + smbd_important_died_21m + smbd_important_died_3y + smbd_important_died_4y + smbd_important_died_5y + smbd_important_died_6y + smbd_important_died_9y + separated_from_parent_18m + separated_from_parent_30m + separated_from_parent_3y + separated_from_parent_4y + separated_from_parent_5y + separated_from_parent_6y + separated_from_parent_8y + moved_18m + moved_30m + moved_3y + moved_4y + moved_5y + moved_6y + moved_9y + pet_died_18m + pet_died_30m + pet_died_3y + pet_died_4y + pet_died_5y + pet_died_6y + pet_died_9y + started_nursery_18m + started_nursery_30m + started_nursery_3y + started_nursery_4y + started_nursery_5y + acquired_new_parent_18m + acquired_new_parent_30m + acquired_new_parent_3y + acquired_new_parent_4y + acquired_new_parent_5y + acquired_new_parent_6y + acquired_new_parent_8y + 
                                   change_carer_18m + change_carer_30m + change_carer_3y + change_carer_4y + change_carer_5y + change_carer_6y + change_carer_8y + friend_relative_ill_8m + friend_relative_ill_21m + friend_relative_ill_3y + friend_relative_ill_4y + friend_relative_ill_5y + friend_relative_ill_6y + friend_relative_ill_9y + partner_died_8m + partner_died_21m + partner_died_3y + partner_died_4y + partner_died_5y + partner_died_6y + partner_died_9y + burglary_or_car_theft_21m + burglary_or_car_theft_3y + burglary_or_car_theft_4y + burglary_or_car_theft_5y + burglary_or_car_theft_6y + burglary_or_car_theft_9y + separated_from_smbd_18m + separated_from_smbd_30m + separated_from_smbd_3y + separated_from_smbd_4y + separated_from_smbd_5y + separated_from_smbd_6y + separated_from_smbd_8y + lost_best_friend_8y + new_sibling_18m + new_sibling_30m + new_sibling_3y + new_sibling_4y + new_sibling_5y + new_sibling_6y + new_sibling_9y + ch_had_fright_18m + ch_had_fright_30m + ch_had_fright_3y + ch_had_fright_4y + ch_had_fright_5y + ch_had_fright_6y + ch_had_fright_8y /17),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_life_events") 
miss <- is.na(original_set[, "post_contextual_risk" ])
xyplot(imp, post_contextual_risk ~I (homeless_childhood_8m + homeless_childhood_21m + homeless_childhood_3y + homeless_childhood_4y + homeless_childhood_5y + homeless_childhood_6y + homeless_childhood_9y + major_financial_problems_8m + major_financial_problems_21m + major_financial_problems_3y + major_financial_problems_4y + major_financial_problems_5y + major_financial_problems_6y + major_financial_problems_9y + income_reduced_8m + income_reduced_21m + income_reduced_3y + income_reduced_4y + income_reduced_5y + income_reduced_6y + income_reduced_9y + unemployed_8m + unemployed_21m + unemployed_3y + unemployed_4y + unemployed_5y + unemployed_6y + unemployed_9y + housing_adequacy_2y + housing_adequacy_4y + housing_basic_living_2y + housing_basic_living_4y + housing_defects_2y + housing_defects_4y + m_education + p_education + neighbourhood_problems_21m + neighbourhood_problems_3y / 10),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_contextual_risk") 
miss <- is.na(original_set[, "post_parental_risk" ])
xyplot(imp, post_parental_risk ~I (work_problems_8m + work_problems_21m + work_problems_3y + work_problems_4y + work_problems_5y + work_problems_6y + work_problems_9y + criminal_record_parent_8m + criminal_record_parent_21m + criminal_record_parent_3y + criminal_record_parent_4y + criminal_record_parent_5y + criminal_record_parent_6y + criminal_record_parent_9y + miscarriage_or_abortion_8m + miscarriage_or_abortion_21m + miscarriage_or_abortion_3y + miscarriage_or_abortion_4y + miscarriage_or_abortion_5y + miscarriage_or_abortion_6y + miscarriage_or_abortion_9y + m_attempted_suicide_8m + m_attempted_suicide_21m + m_attempted_suicide_3y + m_attempted_suicide_4y + m_attempted_suicide_5y + m_attempted_suicide_6y + m_attempted_suicide_9y + m_age + p_age + m_depression_8m + m_depression_21m + m_depression_3y + m_depression_4y + m_depression_5y + m_depression_6y + m_depression_9y + p_depression_8m + p_depression_21m + p_depression_3y + p_depression_4y + p_depression_5y + p_depression_6y + p_depression_9y + m_anxiety_8m + m_anxiety_21m + m_anxiety_3y + m_anxiety_5y + m_anxiety_6y + p_anxiety_8m + p_anxiety_21m + p_anxiety_3y + p_anxiety_4y + p_anxiety_5y + p_anxiety_6y + p_anxiety_9y / 10),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_parental_risk") 
miss <- is.na(original_set[, "post_interpersonal_risk" ])
xyplot(imp, post_interpersonal_risk ~I (divorce_8m + divorce_21m + divorce_3y + divorce_4y + divorce_5y + divorce_6y + divorce_9y + p_rejected_child_8m + p_rejected_child_21m + p_rejected_child_3y + p_rejected_child_4y + p_rejected_child_5y + p_rejected_child_6y + p_rejected_child_9y + p_went_away_8m + p_went_away_21m + p_went_away_3y + p_went_away_4y + p_went_away_5y + p_went_away_6y + p_went_away_9y + conflict_in_family_21m + conflict_in_family_3y + conflict_in_family_4y + conflict_in_family_5y + conflict_in_family_6y + conflict_in_family_9y + conflict_family_violence_8m + conflict_family_violence_21m + conflict_family_violence_3y + conflict_family_violence_4y + conflict_family_violence_5y + conflict_family_violence_6y + conflict_family_violence_9y + m_new_partner_8m + m_new_partner_21m + m_new_partner_3y + m_new_partner_4y + m_new_partner_5y + m_new_partner_6y + m_new_partner_9y + argued_fam_friends_8m + argued_fam_friends_21m + argued_fam_friends_3y + argued_fam_friends_4y + argued_fam_friends_5y + argued_fam_friends_6y + argued_fam_friends_9y / 7),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_interpersonal_risk") 
miss <- is.na(original_set[, "post_direct_victimization" ])
xyplot(imp, post_direct_victimization ~I (bullying_8y + physical_violence_18m + physical_violence_30m + physical_violence_3y + physical_violence_4y + physical_violence_5y + physical_violence_6y + physical_violence_8y + sexual_abuse_18m + sexual_abuse_30m + sexual_abuse_3y + sexual_abuse_4y + sexual_abuse_5y + sexual_abuse_6y + sexual_abuse_8y + p_cruelty_physical_8m + p_cruelty_physical_21m + p_cruelty_physical_3y + p_cruelty_physical_4y + p_cruelty_physical_5y + p_cruelty_physical_6y + p_cruelty_physical_9y + m_cruelty_physical_8m + m_cruelty_physical_21m + m_cruelty_physical_3y + m_cruelty_physical_4y + m_cruelty_physical_5y + m_cruelty_physical_6y + m_cruelty_physical_9y + p_cruelty_emotional_8m + p_cruelty_emotional_21m + p_cruelty_emotional_3y + p_cruelty_emotional_4y + p_cruelty_emotional_5y + p_cruelty_emotional_6y + p_cruelty_emotional_9y + m_cruelty_emotional_21m + m_cruelty_emotional_3y + m_cruelty_emotional_4y + m_cruelty_emotional_5y + m_cruelty_emotional_6y + m_cruelty_emotional_9y / 7),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_direct_victimization")
dev.off()


