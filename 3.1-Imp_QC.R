# As usual, here are the packages we need 
library(mice);

# Load the list of imputed dataset
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_full.rds' file
  imp <- readRDS(imp_path)
}

# for dichotomous variables, replace NA with the median across imputations
replance_med <- function(orig, variables) {
  for (var in variables) {
    misslist <- orig$IDC[which(is.na(orig[, var]))]
    cat(paste(var, length(misslist), sep = " - "))
    cat('\n')
    
    pool_mean <- with(impdat, by(impdat, IDC, function(x) c(median(x[,var]),
                                                            sd(x[,var]))))
    for (cld in orig$IDC) {
      child_m = pool_mean[dimnames(pool_mean)$IDC == cld]
      if (child_m[[1]][2] > 0.000) {
        orig[orig$IDC == cld, var] <- as.integer(child_m[[1]][1])
      }
    }
  }
  return(orig)
}

# for continuous (numeric) variables, replace NA with the mean across imputations
replance_mean <- function(orig, variables) {
  for (var in variables) {
    misslist <- orig$IDC[which(is.na(orig[, var]))]
    cat(paste(var, length(misslist), sep = " - "))
    cat('\n')
    
    pool_mean <- with(impdat, by(impdat, IDC, function(x) c(mean(x[,var]),
                                                            sd(x[,var]))))
    for (cld in orig$IDC) {
      child_m = pool_mean[dimnames(pool_mean)$IDC == cld]
      if (child_m[[1]][2] > 0.000) {
        orig[orig$IDC == cld, var] <- child_m[[1]][1]
      }
    }
  }
  return(orig)
}

# structure variables into domains 
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
outcomes <- c('intern_score_z', 'fmi_z', 'risk_groups') # fat_mass_z
covars   <- c('sex', 'age_child', 'm_bmi_berore_pregnancy', 'm_smoking', 'm_drinking')
auxil    <- c('m_bmi_pregnancy','m_dep_cont_pregnancy', 'p_dep_cont_pregnancy', # for postnatal only 
              'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs',               # for prenatal only 
              'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_age_cont')
exclusion_criteria <- c('pre_percent_missing', 'post_percent_missing', 'twin', 'mother')

# retrieve original (pre-imputation) dataset
original_set <- imp$data
original_set$sex <- as.numeric(as.character(original_set$sex))
original_set$ethnicity <- as.numeric(as.character(original_set$ethnicity))
original_set$risk_groups <- as.numeric(as.character(original_set$risk_groups))

# Compute pre imp means and sds
orig_means <- sapply(original_set, mean, na.rm = T, USE.NAMES = T)
orig_sds <- sapply(original_set, sd, na.rm = T, USE.NAMES = T)

# Stack imputed datasets in long format, exclude the original data
impdat <- complete(imp, action="long", include = FALSE)
# Create mean of imputation dataset
dataset <- replance_med(original_set, c(pre_LE, pre_CR, pre_PR, pre_IR,
                                       post_LE, post_CR, post_PR, post_IR, post_DV))
# dataset_final <- replance_mean(dataset, c(covars[-1], auxil[-7]))

best = data.frame(1:31)

for (var in c(pre_LE, pre_CR, pre_PR, pre_IR,
              post_LE, post_CR, post_PR, post_IR, post_DV, covars[-1], auxil[-7])) {
  m = orig_means[var]
  m_avg = mean(dataset[, var], na.rm = T)
  pool_imp <- with(impdat, by(impdat, .imp, function(x) abs(mean(as.numeric(x[, var])) - m)))
  pool_imp[['31']] <- abs(m_avg - m)
  
  best[, var] = pool_imp

}

#m = c(seq(1.8, 0.8, length.out = 31))
#adj = as.matrix(best[, -1]) * m

#best[, -1] = adj
best = best[, -1]
cols <- rainbow(ncol(best))
i = 1

pdf("mean-diff-iterations.pdf")
for (var in c(pre_LE, pre_CR, pre_PR, pre_IR,
              post_LE, post_CR, post_PR, post_IR, post_DV, covars[-1], auxil[-c(7,10)])) {
  plot(best[,var], type = 'l', col = 'red', lwd=3, 
       xlim = c(0, 32), ylim = c(0, 0.5), main = var,  xlab="Iteration",
       ylab="Absolute mean difference", xaxt="n")
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
pdf("covergence-plots.pdf")
plot(imp, c(pre_LE, pre_CR, pre_PR, pre_IR,
            post_LE, post_CR, post_PR, post_IR, post_DV, covars, auxil))
dev.off()
# The variance between the imputation chains should be equal to the variance within 
# the chains, which indicates healthy convergence.

# COMPARE IMPUTED AND OBSERVED VALUES
# Next, we inspect and compare the density of the incomplete and imputed data
# For loops do not work with this wacky mice function so just run the loop below 
# and paste it: lil help
# for (var in c(pre_LE, pre_CR, pre_PR, pre_IR, post_LE, post_CR, post_PR, post_IR, post_DV, covars, auxil)) {
#   cat(paste('densityplot(imp, ~', var, ')', '\n')) }

pdf("imp-vs-obs.pdf")
par(mfrow=c(2,2))    
densityplot(imp, ~ family_member_died ) 
densityplot(imp, ~ friend_relative_died ) 
densityplot(imp, ~ family_member_ill_pregnancy ) 
densityplot(imp, ~ admitted_to_hospital ) 
densityplot(imp, ~ health ) 
densityplot(imp, ~ unemployed ) 
densityplot(imp, ~ work_study_problems ) 
densityplot(imp, ~ moved_house ) 
densityplot(imp, ~ blood_loss ) 
densityplot(imp, ~ examination ) 
densityplot(imp, ~ baby_worried ) 
densityplot(imp, ~ pregnancy_worried ) 
densityplot(imp, ~ obstetric_care ) 
densityplot(imp, ~ pregnancy_planned ) 
densityplot(imp, ~ victim_robbery ) 
densityplot(imp, ~ financial_problems ) 
densityplot(imp, ~ trouble_pay_pregnancy ) 
densityplot(imp, ~ income_reduced ) 
densityplot(imp, ~ housing_defects ) 
densityplot(imp, ~ housing_adequacy ) 
densityplot(imp, ~ housing_basic_living ) 
densityplot(imp, ~ m_education_pregnancy ) 
densityplot(imp, ~ p_education_pregnancy ) 
densityplot(imp, ~ m_depression_pregnancy ) 
densityplot(imp, ~ m_anxiety_pregnancy ) 
densityplot(imp, ~ m_interp_sensitivity_pregnancy ) 
densityplot(imp, ~ p_depression_pregnancy ) 
densityplot(imp, ~ p_anxiety_pregnancy ) 
densityplot(imp, ~ p_interp_sensitivity_pregnancy ) 
densityplot(imp, ~ m_violence_people ) 
densityplot(imp, ~ m_violence_property ) 
densityplot(imp, ~ m_criminal_record ) 
densityplot(imp, ~ p_criminal_record ) 
densityplot(imp, ~ difficulties_contacts ) 
densityplot(imp, ~ difficulties_partner ) 
densityplot(imp, ~ difficulties_family_friend ) 
densityplot(imp, ~ marital_status_pregnancy ) 
densityplot(imp, ~ divorce_pregnancy ) 
densityplot(imp, ~ family_support ) 
densityplot(imp, ~ family_acceptance ) 
densityplot(imp, ~ family_affection ) 
densityplot(imp, ~ family_acception ) 
densityplot(imp, ~ family_trust ) 
densityplot(imp, ~ family_painful_feelings ) 
densityplot(imp, ~ family_decisions ) 
densityplot(imp, ~ family_conflict ) 
densityplot(imp, ~ family_decisions_problems ) 
densityplot(imp, ~ family_plans ) 
densityplot(imp, ~ family_talk_sadness ) 
densityplot(imp, ~ family_talk_worries ) 
densityplot(imp, ~ family_size_pregnancy ) 
densityplot(imp, ~ sick_or_accident ) 
densityplot(imp, ~ family_member_ill ) 
densityplot(imp, ~ smbd_important_ill ) 
densityplot(imp, ~ parent_died ) 
densityplot(imp, ~ smbd_important_died ) 
densityplot(imp, ~ pet_died ) 
densityplot(imp, ~ school_workload ) 
densityplot(imp, ~ repeated_grade ) 
densityplot(imp, ~ lost_smth_important ) 
densityplot(imp, ~ moved ) 
densityplot(imp, ~ changed_school ) 
densityplot(imp, ~ friend_moved ) 
densityplot(imp, ~ fire_or_burglary ) 
densityplot(imp, ~ material_deprivation ) 
densityplot(imp, ~ financial_difficulties ) 
densityplot(imp, ~ neiborhood_problems ) 
densityplot(imp, ~ trouble_pay_childhood ) 
densityplot(imp, ~ income_once ) 
densityplot(imp, ~ income_chronic ) 
densityplot(imp, ~ unemployed_once ) 
densityplot(imp, ~ unemployed_chronic ) 
densityplot(imp, ~ m_education ) 
densityplot(imp, ~ p_education ) 
densityplot(imp, ~ tension_at_work ) 
densityplot(imp, ~ m_interpersonal_sensitivity ) 
densityplot(imp, ~ p_interpersonal_sensitivity ) 
densityplot(imp, ~ m_depression ) 
densityplot(imp, ~ p_depression ) 
densityplot(imp, ~ m_anxiety ) 
densityplot(imp, ~ p_anxiety ) 
densityplot(imp, ~ m_age ) 
densityplot(imp, ~ p_age ) 
densityplot(imp, ~ marital_problems ) 
densityplot(imp, ~ marital_status ) 
densityplot(imp, ~ family_size ) 
densityplot(imp, ~ m_fad_5yrs ) 
densityplot(imp, ~ m_fad_9yrs ) 
densityplot(imp, ~ p_fad_9yrs ) 
densityplot(imp, ~ conflict_family_member ) 
densityplot(imp, ~ conflict_smbd_else ) 
densityplot(imp, ~ conflict_in_family ) 
densityplot(imp, ~ divorce_childhood ) 
densityplot(imp, ~ argument_friend ) 
densityplot(imp, ~ m_harsh_parent ) 
densityplot(imp, ~ p_harsh_parent ) 
densityplot(imp, ~ bullying ) 
densityplot(imp, ~ physical_violence ) 
densityplot(imp, ~ physical_threats ) 
densityplot(imp, ~ sexual_harrasment ) 
densityplot(imp, ~ sexual_behavior ) 
densityplot(imp, ~ rumors_or_gossip ) 
densityplot(imp, ~ sex ) 
densityplot(imp, ~ age_child ) 
densityplot(imp, ~ m_bmi_berore_pregnancy ) 
densityplot(imp, ~ m_smoking ) 
densityplot(imp, ~ m_drinking ) 
densityplot(imp, ~ m_bmi_pregnancy ) 
densityplot(imp, ~ m_dep_cont_pregnancy ) 
densityplot(imp, ~ p_dep_cont_pregnancy ) 
densityplot(imp, ~ m_bmi_5yrs ) 
densityplot(imp, ~ m_dep_cont_3yrs ) 
densityplot(imp, ~ p_dep_cont_3yrs ) 
densityplot(imp, ~ ethnicity ) 
densityplot(imp, ~ parity ) 
densityplot(imp, ~ gest_age_birth ) 
densityplot(imp, ~ gest_weight ) 
densityplot(imp, ~ m_age_cont )
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

pdf("domain-passive-imp.pdf")
miss <- is.na(original_set[, "pre_life_events" ])
xyplot(imp, pre_life_events ~I (family_member_died + friend_relative_died + family_member_ill_pregnancy + admitted_to_hospital + health + unemployed + work_study_problems + moved_house + blood_loss + examination + baby_worried + pregnancy_worried + obstetric_care + pregnancy_planned + victim_robbery / 15),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_life_events") 
miss <- is.na(original_set[, "pre_contextual_risk" ])
xyplot(imp, pre_contextual_risk ~I (financial_problems + trouble_pay_pregnancy + income_reduced + housing_defects + housing_adequacy + housing_basic_living + m_education_pregnancy + p_education_pregnancy / 8),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_contextual_risk") 
miss <- is.na(original_set[, "pre_parental_risk" ])
xyplot(imp, pre_parental_risk ~I (m_depression_pregnancy + m_anxiety_pregnancy + m_interp_sensitivity_pregnancy + p_depression_pregnancy + p_anxiety_pregnancy + p_interp_sensitivity_pregnancy + m_violence_people + m_violence_property + m_criminal_record + p_criminal_record / 10),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_parental_risk") 
miss <- is.na(original_set[, "pre_interpersonal_risk" ])
xyplot(imp, pre_interpersonal_risk ~I (difficulties_contacts + difficulties_partner + difficulties_family_friend + marital_status_pregnancy + divorce_pregnancy + family_support + family_acceptance + family_affection + family_acception + family_trust + family_painful_feelings + family_decisions + family_conflict + family_decisions_problems + family_plans + family_talk_sadness + family_talk_worries + family_size_pregnancy / 18),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "pre_interpersonal_risk") 
miss <- is.na(original_set[, "post_life_events" ])
xyplot(imp, post_life_events ~I (sick_or_accident + family_member_ill + smbd_important_ill + parent_died + smbd_important_died + pet_died + school_workload + repeated_grade + lost_smth_important + moved + changed_school + friend_moved + fire_or_burglary / 13),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_life_events") 
miss <- is.na(original_set[, "post_contextual_risk" ])
xyplot(imp, post_contextual_risk ~I (material_deprivation + financial_difficulties + neiborhood_problems + trouble_pay_childhood + income_once + income_chronic + unemployed_once + unemployed_chronic + m_education + p_education / 10),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_contextual_risk") 
miss <- is.na(original_set[, "post_parental_risk" ])
xyplot(imp, post_parental_risk ~I (tension_at_work + m_interpersonal_sensitivity + p_interpersonal_sensitivity + m_depression + p_depression + m_anxiety + p_anxiety + m_age + p_age / 9),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_parental_risk") 
miss <- is.na(original_set[, "post_interpersonal_risk" ])
xyplot(imp, post_interpersonal_risk ~I (marital_problems + marital_status + family_size + m_fad_5yrs + m_fad_9yrs + p_fad_9yrs + conflict_family_member + conflict_smbd_else + conflict_in_family + divorce_childhood + argument_friend / 11),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_interpersonal_risk") 
miss <- is.na(original_set[, "post_direct_victimization" ])
xyplot(imp, post_direct_victimization ~I (m_harsh_parent + p_harsh_parent + bullying + physical_violence + physical_threats + sexual_harrasment + sexual_behavior + rumors_or_gossip / 8),
       na.groups = miss, cex = c(1, 1), pch = c(1, 20), ylab = "Domain Imputed", xlab = "Domain Calculated", main = "post_direct_victimization") 
dev.off()


