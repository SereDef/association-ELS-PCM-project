
# Hi again, 
# the following script is building the final dataset that will be used for the analysis 
# of the association between ELS and psycho-cardio-metabolic multi-morbidity in children. 
# I will first exclude participants based on predifined criteria & data availability, 
# then construct the main outcome variable (risk_group) and cumulative exposures. 
# Missing values will be then imputed (see Imputation.R)

# This time you will only need three files: the "PCM_allvars.rds" we build using the 
# "PCM_outcomes_covs_aux.R" script; "prenatal_stress.rds" and "postnatal_stress.rds", 
# for which, if you want to know more, check out https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

#### ---------------------------- Dependencies ---------------------------- ####

# Point to useful libraries
library(pheatmap) # optional 

# This will come in handy for exclusion
'%notin%' <- Negate('%in%')

# check if the path to the datasets is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

################################################################################
# Load datasets
pre_risk <- readRDS(paste(pathtodata, 'prenatal_stress.rds', sep = ""))
post_risk <- readRDS(paste(pathtodata, 'postnatal_stress.rds', sep = ""))
outcome <- readRDS(paste(pathtodata,'PCM_allvars.rds', sep = ""))

# merge them
ELS_PCM <- Reduce(function(x,y) merge(x = x, y = y, by = c('IDC', 'IDM'), all.x = TRUE),
       list(pre_risk, post_risk, outcome)) 

################################################################################
## -------------------- Exclude participants (flowchart) -------------------- ##
################################################################################

initial_sample <- nrow(ELS_PCM)

## First exclusion step: remove participants whose missing value frequency is too high. 
## (i.e > 50% missing) for each developmental period. 
 
  # Exclude children with high missingness in the prenatal period
prenatal_miss50 <- ELS_PCM[ELS_PCM$pre_percent_missing < 50,]
after_prenatal <- nrow(prenatal_miss50)
sub1 <- after_prenatal - initial_sample
 
  # Exclude children with high missingness in the postnatal period
postnatal_miss50 <- prenatal_miss50[prenatal_miss50$post_percent_missing < 50,]
after_postnatal <- nrow(postnatal_miss50)
sub2 <- after_postnatal - after_prenatal

## Second step: excluded children with missing internalizing or { CMR scores }. 
 
  # Exclude children with missing internalizing score
intern_miss <- postnatal_miss50[!is.na(postnatal_miss50$intern_score_z),] 
after_intern <- nrow(intern_miss)
sub3 <- after_intern - after_postnatal
 
  # Exclude children with missing { CMR score } fat mass*
cmr_miss <- intern_miss[!is.na(intern_miss$fat_mass_z),] 
after_cmr <- nrow(cmr_miss)
sub4 <- after_cmr - after_intern

## Third step: exclude all twins and select the sibling with better data 
 
  # Exclude twins
no_twins <- cmr_miss[cmr_miss$twin == 0, ]
after_twins <- nrow(no_twins)
sub5 <- after_twins - after_cmr
 
  # Select only one sibling (based on data availability or randomly).
# First, I determine a list of mothers that have more than one child in the set.
# NOTE: duplicated() is the best option I could find to determine which mother IDs
# recur more than once (very non-elegant, tbh, but using table() gets even uglier)
siblings_id = data.frame(no_twins$mother[duplicated(no_twins$mother)])
# duplicated() funtion does not allow to ignore NAs so I remove them manually.
# I also transform the numeric vector into a dataframe because of indexing problems.
siblings_id = siblings_id[!is.na(siblings_id),]; siblings_id = data.frame(siblings_id);
# Second, I create an empty vector to fill with the IDC of the sibling(s) with more 
# missing items or with a randomly picked sibling in case they have the same nr of missing. 
worse_sibling = rep(NA, dim(siblings_id)[1] + 1) # I will need the "+1" for triplets! 
# Loop through the mother IDs I previously identified and link them to IDCs
for (i in 1:dim(siblings_id)[1]) {
  siblings = no_twins[no_twins$mother == siblings_id[i,1], ] # identify the couples of siblings
  # For some reason when I run the line above 2 rows of NAs are created too, go figure. 
  # Let's get rid of them:
  siblings = siblings[rowSums(is.na(siblings)) != ncol(siblings), ]
  # There is one mother with 3 siblings, let's select the "best" one and get rid of the other two
  if (dim(siblings)[1] > 2) {
    nmiss = c( sum(is.na(siblings[1,])), sum(is.na(siblings[2,])), sum(is.na(siblings[3,])) )
    if (which.min(nmiss) == 1) { worse_sibling[i] = siblings[2,'IDC']
    worse_sibling[dim(siblings_id)[1] + 1] = siblings[3,'IDC']
    } else if (which.min(nmiss) == 2) { worse_sibling[i] = siblings[1,'IDC']
    worse_sibling[dim(siblings_id)[1] + 1] = siblings[3,'IDC']
    } else { worse_sibling[i] = siblings[1,'IDC'] 
    worse_sibling[dim(siblings_id)[1] + 1] = siblings[2,'IDC'] }
  }
  # otherwise, select the "worse" sibling (with more missing) and add to it the the black list
  if ( sum(is.na(siblings[1,])) > sum(is.na(siblings[2,])) ) {
    worse_sibling[i] = siblings[1,'IDC']
  } else if ( sum(is.na(siblings[1,])) == sum(is.na(siblings[2,])) ) {
      worse_sibling[i] = siblings[sample(1:2, 1),'IDC']
  } else { worse_sibling[i] = siblings[2,'IDC'] }
}

# Now we are finally ready to exclude siblings
final <- no_twins[no_twins$IDC %notin% worse_sibling, ]
after_siblings <- nrow(final)
sub6 <- after_siblings - after_twins

# Flowchart
flowchart <- list(initial_sample, sub1, after_prenatal, sub2, after_postnatal, sub3, 
                 after_intern, sub4, after_cmr, sub5, after_twins, sub6, after_siblings)

# Rename final dataset:
ELS_PCM <- final
cat(paste("Well, congrats! Your final dataset includes", after_siblings ,"participants."))

################################################################################
#### ------------------- Construct RISK GROUPS variable ------------------- ####
################################################################################

# Compute groups 
ELS_PCM$int = ifelse(ELS_PCM$intern_score_z > quantile(ELS_PCM$intern_score_z, probs = 0.8), 1, 0) # 590 risk, 2780 no risk
ELS_PCM$fat = ifelse(ELS_PCM$fat_mass_z > quantile(ELS_PCM$fat_mass_z, probs = 0.8), 1, 0) # 674 risk, 2696 no risk

ELS_PCM$risk_groups = rep(NA, after_siblings)
for (i in 1:after_siblings) {
  if (ELS_PCM$int[i] == 0 & ELS_PCM$fat[i] == 0) { ELS_PCM$risk_groups[i] = 0 # healthy
  } else if (ELS_PCM$int[i] == 1 & ELS_PCM$fat[i] == 0) { ELS_PCM$risk_groups[i] = 1 # High intenalizing  only
  } else if (ELS_PCM$int[i] == 0 & ELS_PCM$fat[i] == 1) { ELS_PCM$risk_groups[i] = 2 # High fat mass only
  } else { ELS_PCM$risk_groups[i] = 3 } # multimorbid
}

ELS_PCM$risk_groups = as.factor(ELS_PCM$risk_groups)
# summary(ELS_PCM$risk_groups)  ##    0    1    2    3 
                                ## 2272  425  516  158
# # Let's first factor that bad boy 
# imp$risk_groups = factor(groups$risk_groups, 
#                          levels = c(0:3), 
#                          labels = c("healthy", "internalizing_only", "cardiometabolic_only", "multimorbid"))

# #display groups 
# attach(groups); plot(intern_score_z, fat_mass_z, 
#                      col=c("cornflowerblue","darkgrey", "darkgoldenrod2","chartreuse4")[risk_groups]); detach(groups)
# 
# #display groups in relation to stress 
# attach(groups); plot(presum, postsum, 
#                      col=c("cornflowerblue","darkgrey", "darkgoldenrod2","chartreuse4")[risk_groups]); detach(groups)

################################################################################
#### -------------- Construct CUMULATIVE STRESS variables ----------------- ####
################################################################################

# compute sum scores for prenatal and postnatal stress exposure
ELS_PCM$prenatal_stress <- rowSums(ELS_PCM[c("pre_life_events","pre_contextual_risk", 
                                      "pre_personal_stress", "pre_interpersonal_stress")], na.rm = F)
ELS_PCM$postnatal_stress <- rowSums(ELS_PCM[c("post_life_events","post_contextual_risk", 
                                       "post_parental_risk", "post_interpersonal_risk", 
                                       "post_direct_victimization")], na.rm = F)

#------------------------------------------------------------------------------#
          ## OPTIONAL : check out missing pattern per stress domain ##
#------------------------------------------------------------------------------#
# Let's have a look at the pattern of missing in our domains after exclusion process.

# domainonly <- ELS_PCM[, c('pre_life_events', 'pre_contextual_risk', 'pre_personal_stress', 'pre_interpersonal_stress', 
#               'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')]
# summary(domainonly)

# Heatmap of missing values together with table
# missingpattern <- md.pairs(domainonly) # outputs four tables: ($rr) how many datapoints are observed
# ($rm) observed and missing, ($mr) missing versus observed and ($mm) is missing vs missing
# pheatmap(as.matrix(missingpattern$mm), display_numbers = T, number_format = "%.0f")

# missingtable = md.pattern(domainonly)

#------------------------------------------------------------------------------#
# For the sake of time efficiency (and my mental health) let us reatain only those 
# variables that are useful for imputation in the final dataset. Once I am at it,
# I also order them by domain. This is important because, it turns out that mice 
# is sensitive to the order of the variables in the set (even though this may be 
# a version-specific issue)

ELS_PCM_essentials = ELS_PCM[, c('IDC', 
                      # all variables for prenatal risk
                'family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems', 'moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", "m_education_pregnancy", # CR
                "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record",  # PS, without age, as the same variable is already in postnatal PR score
                "difficulties_contacts","difficulties_partner","difficulties_family_friend","marital_status_pregnancy","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                      # all variables for postnatal risk
                'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                'material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education', # CR
                'tension_at_work','m_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                      # all domain scores
                'pre_life_events', 'pre_contextual_risk', 'pre_personal_stress', 'pre_interpersonal_stress', 
                'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization',
                      # cumulative prenatal and postnatal stress exposure
                'prenatal_stress', "postnatal_stress",
                      # outcome variables and covariates
                'intern_score_z', 'fat_mass_z', 'risk_groups',
                'sex', 'age_child', 'm_bmi_berore_pregnancy', 'm_smoking', 'm_drinking',
                      # additional auxiliary variables for imputation
                'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_bmi_pregnancy', 'm_bmi_5yrs',
                'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs', 'm_age_cont')]

#------------------------------------------------------------------------------#
                      ## OPTIONAL : check out the dataset ##
#------------------------------------------------------------------------------#
# # Examine missing frequency per variable 
# percent_missing <- function(var) { sum(is.na(var)) / length(var) * 100 }
# missing_frequency = apply(ELS_PCM_essentials, 2, percent_missing)
# print(sort(missing_frequency, decreasing = T))
# 
# # Examine correlations 
# test_essentials = subset(ELS_PCM_essentials[, 2:113], select=-c(parent_died, p_age))
# correlations = round(cor(test_essentials, use = "complete.obs"),2)
# pheatmap(correlations, cluster_rows = F, cluster_cols = F)

################################################################################
#### --------------------------- save and run ----------------------------- ####
################################################################################

# Save the dataset in an .rds file, in the directory where the raw data are stored
saveRDS(ELS_PCM_essentials, paste(pathtodata,'ELSPCM_dataset.rds', sep =""))
