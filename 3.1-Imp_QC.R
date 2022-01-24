# As usual, here are the packages we need 
library(mice);

# Load the list of imputed dataset
if (exists("imp_samp") == F) { 
  imp_path <- file.choose() # ATTENTION! Choose the 'imputation_list_sample.rds' file
  dirname <- dirname(imp_path)
  ifelse(!dir.exists(file.path(dirname, 'imp-QC')), dir.create(file.path(dirname, 'imp-QC')), F)
  imp_samp <- readRDS(imp_path)
  #imp_full <- readRDS(file.path(dirname, "imputation_list_full.rds"))
}

# Organize variable names into domains to specify them later more easily
pre_LE <- c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 
            'health', 'unemployed', 'work_study_problems', 'moved_house', 'blood_loss', 'examination', 
            'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery')
pre_CR <- c('financial_problems', 'trouble_pay_pregnancy', 'income_reduced', 'housing_defects', 'housing_adequacy', 
            'housing_basic_living', 'm_education_pregnancy', 'p_education_pregnancy')
pre_PR <- c('m_depression_pregnancy', 'm_anxiety_pregnancy', 'm_interp_sensitivity_pregnancy', 
            'p_depression_pregnancy', 'p_anxiety_pregnancy', 'p_interp_sensitivity_pregnancy', 'm_violence_people', 
            'm_violence_property', 'm_criminal_record', 'p_criminal_record') # without age, as the same variable is already in postnatal PR score
pre_IR <- c('difficulties_contacts','difficulties_partner','difficulties_family_friend','marital_status_pregnancy',
            'divorce_pregnancy','family_support','family_acceptance','family_affection','family_acception','family_trust',
            'family_painful_feelings','family_decisions','family_conflict','family_decisions_problems',
            'family_plans','family_talk_sadness', 'family_talk_worries', 'family_size_pregnancy')
post_LE <- c('sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died',
             'school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary')
post_CR <- c('material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once',
             'income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education')
post_PR <- c('tension_at_work','m_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression',
             'm_anxiety','p_anxiety','m_age','p_age')
post_IR <- c('marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member',
             'conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend')
post_DV <- c('m_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment',
             'sexual_behavior','rumors_or_gossip')
domains <- c('pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk', 
             'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')
outcomes_13y <- c('intern_score_13', 'total_fat_13', 'andr_fat_mass_13', 'tot_fat_percent_13')
outcomes_09y <- c('intern_score_09', 'total_fat_09', 'andr_fat_mass_09', 'tot_fat_percent_09')
risk_grps <-c("risk_groups_tot", "risk_groups_andr", "risk_groups_perc", "risk_groups_tot_REC", "risk_groups_andr_REC", "risk_groups_perc_REC")
covars   <- c('sex', 'age_child', 'm_bmi_before_pregnancy', 'm_smoking', 'm_drinking')
auxil    <- c('m_bmi_pregnancy','m_dep_cont_pregnancy', 'p_dep_cont_pregnancy', # for postnatal only 
              'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs',               # for prenatal only 
              'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_age_cont')
exclusion_criteria <- c('pre_percent_missing', 'post_percent_missing', 'twin', 'mother')

# ------------------------------------------------------------------------------
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
orig_means <- sapply(original_set, mean, na.rm = T, USE.NAMES = T)
orig_sds <- sapply(original_set, sd, na.rm = T, USE.NAMES = T)

# Create mean of imputation dataset
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
# for (var in colnames(original_set)) { cat(paste('densityplot(imp, ~', var, ')', '\n')) }

pdf(file.path(dirname,"imp-QC","imp-vs-obs.pdf"))
densityplot(imp_samp, ~ family_member_died ) 
densityplot(imp_samp, ~ friend_relative_died ) 
densityplot(imp_samp, ~ family_member_ill_pregnancy ) 
densityplot(imp_samp, ~ admitted_to_hospital ) 
densityplot(imp_samp, ~ health ) 
densityplot(imp_samp, ~ unemployed ) 
densityplot(imp_samp, ~ work_study_problems ) 
densityplot(imp_samp, ~ moved_house ) 
densityplot(imp_samp, ~ blood_loss ) 
densityplot(imp_samp, ~ examination ) 
densityplot(imp_samp, ~ baby_worried ) 
densityplot(imp_samp, ~ pregnancy_worried ) 
densityplot(imp_samp, ~ obstetric_care ) 
densityplot(imp_samp, ~ pregnancy_planned ) 
densityplot(imp_samp, ~ victim_robbery ) 
densityplot(imp_samp, ~ financial_problems ) 
densityplot(imp_samp, ~ trouble_pay_pregnancy ) 
densityplot(imp_samp, ~ income_reduced ) 
densityplot(imp_samp, ~ housing_defects ) 
densityplot(imp_samp, ~ housing_adequacy ) 
densityplot(imp_samp, ~ housing_basic_living ) 
densityplot(imp_samp, ~ m_education_pregnancy ) 
densityplot(imp_samp, ~ p_education_pregnancy ) 
densityplot(imp_samp, ~ m_depression_pregnancy ) 
densityplot(imp_samp, ~ m_anxiety_pregnancy ) 
densityplot(imp_samp, ~ m_interp_sensitivity_pregnancy ) 
densityplot(imp_samp, ~ p_depression_pregnancy ) 
densityplot(imp_samp, ~ p_anxiety_pregnancy ) 
densityplot(imp_samp, ~ p_interp_sensitivity_pregnancy ) 
densityplot(imp_samp, ~ m_violence_people ) 
densityplot(imp_samp, ~ m_violence_property ) 
densityplot(imp_samp, ~ m_criminal_record ) 
densityplot(imp_samp, ~ p_criminal_record ) 
densityplot(imp_samp, ~ difficulties_contacts ) 
densityplot(imp_samp, ~ difficulties_partner ) 
densityplot(imp_samp, ~ difficulties_family_friend ) 
densityplot(imp_samp, ~ marital_status_pregnancy ) 
densityplot(imp_samp, ~ divorce_pregnancy ) 
densityplot(imp_samp, ~ family_support ) 
densityplot(imp_samp, ~ family_acceptance ) 
densityplot(imp_samp, ~ family_affection ) 
densityplot(imp_samp, ~ family_acception ) 
densityplot(imp_samp, ~ family_trust ) 
densityplot(imp_samp, ~ family_painful_feelings ) 
densityplot(imp_samp, ~ family_decisions ) 
densityplot(imp_samp, ~ family_conflict ) 
densityplot(imp_samp, ~ family_decisions_problems ) 
densityplot(imp_samp, ~ family_plans ) 
densityplot(imp_samp, ~ family_talk_sadness ) 
densityplot(imp_samp, ~ family_talk_worries ) 
densityplot(imp_samp, ~ family_size_pregnancy ) 
densityplot(imp_samp, ~ sick_or_accident ) 
densityplot(imp_samp, ~ family_member_ill ) 
densityplot(imp_samp, ~ smbd_important_ill ) 
densityplot(imp_samp, ~ parent_died ) 
densityplot(imp_samp, ~ smbd_important_died ) 
densityplot(imp_samp, ~ pet_died ) 
densityplot(imp_samp, ~ school_workload ) 
densityplot(imp_samp, ~ repeated_grade ) 
densityplot(imp_samp, ~ lost_smth_important ) 
densityplot(imp_samp, ~ moved ) 
densityplot(imp_samp, ~ changed_school ) 
densityplot(imp_samp, ~ friend_moved ) 
densityplot(imp_samp, ~ fire_or_burglary ) 
densityplot(imp_samp, ~ material_deprivation ) 
densityplot(imp_samp, ~ financial_difficulties ) 
densityplot(imp_samp, ~ neiborhood_problems ) 
densityplot(imp_samp, ~ trouble_pay_childhood ) 
densityplot(imp_samp, ~ income_once ) 
densityplot(imp_samp, ~ income_chronic ) 
densityplot(imp_samp, ~ unemployed_once ) 
densityplot(imp_samp, ~ unemployed_chronic ) 
densityplot(imp_samp, ~ m_education ) 
densityplot(imp_samp, ~ p_education ) 
densityplot(imp_samp, ~ tension_at_work ) 
densityplot(imp_samp, ~ m_interpersonal_sensitivity ) 
densityplot(imp_samp, ~ p_interpersonal_sensitivity ) 
densityplot(imp_samp, ~ m_depression ) 
densityplot(imp_samp, ~ p_depression ) 
densityplot(imp_samp, ~ m_anxiety ) 
densityplot(imp_samp, ~ p_anxiety ) 
densityplot(imp_samp, ~ m_age ) 
densityplot(imp_samp, ~ p_age ) 
densityplot(imp_samp, ~ marital_problems ) 
densityplot(imp_samp, ~ marital_status ) 
densityplot(imp_samp, ~ family_size ) 
densityplot(imp_samp, ~ m_fad_5yrs ) 
densityplot(imp_samp, ~ m_fad_9yrs ) 
densityplot(imp_samp, ~ p_fad_9yrs ) 
densityplot(imp_samp, ~ conflict_family_member ) 
densityplot(imp_samp, ~ conflict_smbd_else ) 
densityplot(imp_samp, ~ conflict_in_family ) 
densityplot(imp_samp, ~ divorce_childhood ) 
densityplot(imp_samp, ~ argument_friend ) 
densityplot(imp_samp, ~ m_harsh_parent ) 
densityplot(imp_samp, ~ p_harsh_parent ) 
densityplot(imp_samp, ~ bullying ) 
densityplot(imp_samp, ~ physical_violence ) 
densityplot(imp_samp, ~ physical_threats ) 
densityplot(imp_samp, ~ sexual_harrasment ) 
densityplot(imp_samp, ~ sexual_behavior ) 
densityplot(imp_samp, ~ rumors_or_gossip ) 
densityplot(imp_samp, ~ pre_life_events ) 
densityplot(imp_samp, ~ pre_contextual_risk ) 
densityplot(imp_samp, ~ pre_parental_risk ) 
densityplot(imp_samp, ~ pre_interpersonal_risk ) 
densityplot(imp_samp, ~ post_life_events ) 
densityplot(imp_samp, ~ post_contextual_risk ) 
densityplot(imp_samp, ~ post_parental_risk ) 
densityplot(imp_samp, ~ post_interpersonal_risk ) 
densityplot(imp_samp, ~ post_direct_victimization ) 
densityplot(imp_samp, ~ prenatal_stress ) 
densityplot(imp_samp, ~ postnatal_stress ) 
densityplot(imp_samp, ~ intern_score_13 ) 
densityplot(imp_samp, ~ tot_fat_percent_13 )
densityplot(imp_samp, ~ andr_fat_mass_13 ) 
densityplot(imp_samp, ~ total_fat_13 ) 
densityplot(imp_samp, ~ intern_score_09 )
densityplot(imp_samp, ~ tot_fat_percent_09 ) 
densityplot(imp_samp, ~ andr_fat_mass_09 ) 
densityplot(imp_samp, ~ total_fat_09 ) 
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
#             'xyplot(imp_samp, ', dnames[i], ' ', passive_imp_formula(ds[i]), '),', '\n', 
#             'na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "',dnames[i], '") \n')) }

pdf(file.path(dirname,"imp-QC","main-&-domains.pdf"))
stripplot(imp_samp, prenatal_stress_z ) 
stripplot(imp_samp, postnatal_stress_z ) 
stripplot(imp_samp, intern_score_13_z ) 
stripplot(imp_samp, tot_fat_percent_13_z ) 
stripplot(imp_samp, andr_fat_mass_13_z ) 
stripplot(imp_samp, total_fat_13_z ) 
miss <- is.na(original_set[, "pre_life_events" ])
xyplot(imp_samp, pre_life_events ~I (family_member_died + friend_relative_died + family_member_ill_pregnancy + admitted_to_hospital + health + unemployed + work_study_problems + moved_house + blood_loss + examination + baby_worried + pregnancy_worried + obstetric_care + pregnancy_planned + victim_robbery / 15),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_life_events") 
miss <- is.na(original_set[, "pre_contextual_risk" ])
xyplot(imp_samp, pre_contextual_risk ~I (financial_problems + trouble_pay_pregnancy + income_reduced + housing_defects + housing_adequacy + housing_basic_living + m_education_pregnancy + p_education_pregnancy / 8),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_contextual_risk") 
miss <- is.na(original_set[, "pre_parental_risk" ])
xyplot(imp_samp, pre_parental_risk ~I (m_depression_pregnancy + m_anxiety_pregnancy + m_interp_sensitivity_pregnancy + p_depression_pregnancy + p_anxiety_pregnancy + p_interp_sensitivity_pregnancy + m_violence_people + m_violence_property + m_criminal_record + p_criminal_record / 10),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_parental_risk") 
miss <- is.na(original_set[, "pre_interpersonal_risk" ])
xyplot(imp_samp, pre_interpersonal_risk ~I (difficulties_contacts + difficulties_partner + difficulties_family_friend + marital_status_pregnancy + divorce_pregnancy + family_support + family_acceptance + family_affection + family_acception + family_trust + family_painful_feelings + family_decisions + family_conflict + family_decisions_problems + family_plans + family_talk_sadness + family_talk_worries + family_size_pregnancy / 18),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_interpersonal_risk") 
miss <- is.na(original_set[, "post_life_events" ])
xyplot(imp_samp, post_life_events ~I (sick_or_accident + family_member_ill + smbd_important_ill + parent_died + smbd_important_died + pet_died + school_workload + repeated_grade + lost_smth_important + moved + changed_school + friend_moved + fire_or_burglary / 13),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_life_events") 
miss <- is.na(original_set[, "post_contextual_risk" ])
xyplot(imp_samp, post_contextual_risk ~I (material_deprivation + financial_difficulties + neiborhood_problems + trouble_pay_childhood + income_once + income_chronic + unemployed_once + unemployed_chronic + m_education + p_education / 10),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_contextual_risk") 
miss <- is.na(original_set[, "post_parental_risk" ])
xyplot(imp_samp, post_parental_risk ~I (tension_at_work + m_interpersonal_sensitivity + p_interpersonal_sensitivity + m_depression + p_depression + m_anxiety + p_anxiety + m_age + p_age / 9),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_parental_risk") 
miss <- is.na(original_set[, "post_interpersonal_risk" ])
xyplot(imp_samp, post_interpersonal_risk ~I (marital_problems + marital_status + family_size + m_fad_5yrs + m_fad_9yrs + p_fad_9yrs + conflict_family_member + conflict_smbd_else + conflict_in_family + divorce_childhood + argument_friend / 11),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_interpersonal_risk") 
miss <- is.na(original_set[, "post_direct_victimization" ])
xyplot(imp_samp, post_direct_victimization ~I (m_harsh_parent + p_harsh_parent + bullying + physical_violence + physical_threats + sexual_harrasment + sexual_behavior + rumors_or_gossip / 8),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_direct_victimization") 
dev.off()


