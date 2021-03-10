
# Hi again, 
# the following code is building the final dataset used for the analysis of the association 
# between ELS and psycho-cardio-metabolic multi-morbidity in children. I will first 
# exclude participants based on data availability, then impute the missing values and
# save a final (complete) dataset. 

# This time you will only need three files: the "PCM_allvars.rds" we build using the 
# "PCM_outcomes_covs_aux.R" script; "prenatal_stress.rds" and "postnatal_stress.rds", 
# for which, if you want to know more, check out https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

#### ---------------------------- Dependencies ---------------------------- ####

# Point to the necessary libraries
library(pheatmap); library(mice);

# This will come in handy for exclusion
'%notin%' <- Negate('%in%')

# check if the path to the datasets is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

#-------------------------------------------------------------------------------
# Load datasets
pre_risk <- readRDS(paste(pathtodata, 'prenatal_stress.rds', sep = ""))
post_risk <- readRDS(paste(pathtodata, 'postnatal_stress.rds', sep = ""))
outcome <- readRDS(paste(pathtodata,'PCM_allvars.rds', sep = ""))

# merge them
ELS_PCM <- Reduce(function(x,y) merge(x = x, y = y, by = c('IDC', 'IDM'), all.x = TRUE),
       list(pre_risk, post_risk, outcome)) 

##----------------------------------------------------------------------------##
## -------------------- Exclude participants (flowchart) -------------------- ##
##----------------------------------------------------------------------------##

initial_sample <- 9901

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

#-------------------------------------------------------------------------------
#                                    OPTIONAL                                  #
# Let's have a look at the pattern of missing in our domains after exclusion process.

# domainonly <- ELS_PCM[, c('pre_life_events', 'pre_contextual_risk', 'pre_personal_stress', 'pre_interpersonal_stress', 
#                                'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')]
# summary(domainonly)

# Heatmap of missing values together with table
# missingpattern <- md.pairs(domainonly) # outputs four tables: ($rr) how many datapoints are observed
# ($rm) observed and missing, ($mr) missing versus observed and ($mm) is missing vs missing
# pheatmap(as.matrix(missingpattern$mm), display_numbers = T, number_format = "%.0f")

# missingtable = md.pattern(domainonly)

# library(VIM); aggr_plot <- aggr(domainonly_final, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
#               labels=names(domainonly_final), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

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
                'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education', # CR
                'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                      # all domain scores
                'pre_life_events', 'pre_contextual_risk', 'pre_personal_stress', 'pre_interpersonal_stress', 
                'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization',
                      # outcome variables and covariates
                'intern_score_z', 'fat_mass_z',
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

#------------------------------------------------------------------------------#
##---------------------------- Imputation model ----------------------------- ##
#------------------------------------------------------------------------------#

# We started with a dry run to specify the default arguments.
imp0 <- mice(ELS_PCM_essentials, maxit = 0, 
             defaultMethod = rep('pmm',4)) # set the imputation method to predictive mean matching (PMM)* 
             #remove.collinear = F) # Because maternal age is measured twice in prenatal and postnatal 

# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- make.method(ELS_PCM_essentials)
# We use passive imputation for the domain scores. This means that the indicator items  
# are imputed first, and then, using these complete items, mean domain scores are 
# derived by the formula specified below.
meth['pre_life_events'] <- "~I( (family_member_died + friend_relative_died + family_member_ill_pregnancy + admitted_to_hospital + health + unemployed + work_study_problems + moved_house + blood_loss + examination + baby_worried + pregnancy_worried + obstetric_care + pregnancy_planned + victim_robbery) / 15)" 
meth['pre_contextual_risk'] <- "~I( (financial_problems + trouble_pay_pregnancy + income_reduced + housing_defects + housing_adequacy + housing_basic_living + m_education_pregnancy) / 7)"
meth['pre_personal_stress'] <- "~I( (m_age + m_depression_pregnancy + m_anxiety_pregnancy + m_interp_sensitivity_pregnancy + m_violence_people + m_violence_property + m_criminal_record) / 7)"
meth['pre_interpersonal_stress'] <- "~I( (difficulties_contacts + difficulties_partner + difficulties_family_friend + marital_status_pregnancy + divorce_pregnancy + family_support + family_acceptance + family_affection + family_acception + family_trust + family_painful_feelings + family_decisions + family_conflict + family_decisions_problems + family_plans + family_talk_sadness + family_talk_worries + family_size_pregnancy) / 18)"
meth['post_life_events'] <- "~I( (sick_or_accident + family_member_ill + smbd_important_ill + parent_died + smbd_important_died + pet_died + school_workload + repeated_grade + lost_smth_important + moved + changed_school + friend_moved + fire_or_burglary) / 13)"
meth['post_contextual_risk'] <- "~I( (tension_at_work + material_deprivation + financial_difficulties + neiborhood_problems + trouble_pay_childhood + income_once + income_chronic + unemployed_once + unemployed_chronic + m_education + p_education) / 11)"
meth['post_parental_risk'] <- "~I( (m_age + p_age + m_interpersonal_sensitivity + m_anxiety + m_depression + p_interpersonal_sensitivity + p_depression + p_anxiety) / 8)"
meth['post_interpersonal_risk'] <- "~I( (conflict_family_member + conflict_smbd_else + conflict_in_family + divorce_childhood + argument_friend + marital_problems + marital_status + family_size + m_fad_5yrs + m_fad_9yrs + p_fad_9yrs) / 11)"
meth['post_direct_victimization'] <- "~I( (physical_violence + physical_threats + sexual_harrasment + sexual_behavior + rumors_or_gossip + m_harsh_parent + p_harsh_parent + bullying) / 8)"

# We are going to need a different set of predictors for the different variables we impute 
# so let's define them using the predictormatrix, that gives these instructions to 
# mice
predictormatrix <- imp0$predictorMatrix

# Do not use IDC as predictor:
predictormatrix[, "IDC"] <- predictormatrix["IDC",] <- 0
# Do not use age_child as a predictor (no reason to believe it is associated with missingness)
predictormatrix[, "age_child"] <- predictormatrix["age_child",] <- 0

               ### Impute auxiliary variables and covariates ###

# To prevent multicollinearity, we adjust the predictor matrix such that the 
# auxiliary variables would be imputed given the domain scores and the outcomes, 
# but not by the single items.
predictormatrix[c('sex', 'age_child', 'm_bmi_berore_pregnancy', 'm_smoking', 'm_drinking',
                  'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_bmi_pregnancy', 'm_bmi_5yrs',
                  'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs', 'm_age_cont'),
                c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems','moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", "m_education_pregnancy", # CR
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","marital_status_pregnancy","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic','m_education','p_education', # CR
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                  'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'sex', 'age_child', 'm_bmi_berore_pregnancy', 'm_smoking', 'm_drinking',
                  'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_bmi_pregnancy', 'm_bmi_5yrs',
                  'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs', 'm_age_cont' )] <- 0

# All single items are then imputed on the item level. We adjust the predictormatrix 
# such that we impute the items within a domain score given the other items in that
# domain score, the remaining domain scores, the outcomes and the auxiliary variables.
  # Auxiliary variables were selected because they are either related to missingness 
  # or to the domain scores themselves. They are used in the prediction model for 
  # every other variable. When information is available prenatally as well postnatally, 
  # we use the different period to minimize bias (E.G., for the imputation of prenatal
  # items we used BMI of the mother when the child was aged 5, and for imputation 
  # of postnatal items we used BMI of the mother during pregnancy).
# This method is more efficient as it does not use all available items, hence the 
# computational load is lower. The multicollinearity issue is also resolved. 
# The technique has been found to reduce standard error substantially compared to 
# complete-case analysis (Plumpton et al., 2010), and it outperforms other existing 
# techniques (Eekhout et al., 2018). 

# So, let's adjust the predictormatrices such that ‘redundant’ items were not used as a predictor.

                                  ### PRENATAL ###
# LE domain 
predictormatrix[c('family_member_died','friend_relative_died','family_member_ill_pregnancy','admitted_to_hospital','health','unemployed','work_study_problems','moved_house','blood_loss','examination','baby_worried','pregnancy_worried','obstetric_care','pregnancy_planned','victim_robbery'), # LE
                c("financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", "m_education_pregnancy", # CR
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","marital_status_pregnancy","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'p_education', # CR minus m_education that is auxiliary for prenatal variables
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age', 'p_age', # PR 
                  'marital_problems','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR minus marital_status that is auxiliary for prenatal variables 
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'pre_life_events',
                  'm_bmi_berore_pregnancy', 'm_bmi_pregnancy', 'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy')] <- 0
# CR domain
predictormatrix[c("financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", "m_education_pregnancy"), # CR
                c('family_member_died','friend_relative_died','family_member_ill_pregnancy','admitted_to_hospital','health','unemployed','work_study_problems','moved_house','blood_loss','examination','baby_worried','pregnancy_worried','obstetric_care','pregnancy_planned','victim_robbery', # LE
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","marital_status_pregnancy","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic','p_education', # CR minus m_education that is auxiliary for prenatal variables
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age', 'p_age', # PR 
                  'marital_problems','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR minus marital_status that is auxiliary for prenatal variables 
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'pre_contextual_risk',
                  'm_bmi_berore_pregnancy', 'm_bmi_pregnancy', 'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy')] <- 0

# PS domain 
# I did not allow m_education_pregnancy and dep here because of multicollinearity
predictormatrix[c("m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record"), # PS
                c('family_member_died','friend_relative_died','family_member_ill_pregnancy','admitted_to_hospital','health','unemployed','work_study_problems','moved_house','blood_loss','examination','baby_worried','pregnancy_worried','obstetric_care','pregnancy_planned','victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living","m_education_pregnancy", # CR
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","marital_status_pregnancy","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'p_education', # CR minus m_education that is auxiliary for prenatal variables
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR 
                  'marital_problems','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR  minus marital_status that is auxiliary for prenatal variables
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'pre_personal_stress',
                  'm_bmi_berore_pregnancy', 'm_bmi_pregnancy', 'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy')] <- 0

# IP domain
predictormatrix[c("difficulties_contacts","difficulties_partner","difficulties_family_friend","marital_status_pregnancy","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy"), # IS
                c('family_member_died','friend_relative_died','family_member_ill_pregnancy','admitted_to_hospital','health','unemployed','work_study_problems','moved_house','blood_loss','examination','baby_worried','pregnancy_worried','obstetric_care','pregnancy_planned','victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", "m_education_pregnancy", # CR
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record",  # PS
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'p_education', # CR minus m_education that is auxiliary for prenatal variables
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR 
                  'marital_problems','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR  minus marital_status that is auxiliary for prenatal variables
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'pre_interpersonal_stress',
                  'm_bmi_berore_pregnancy', 'm_bmi_pregnancy', 'm_dep_cont_pregnancy', 'p_dep_cont_pregnancy')] <- 0
  
                                  ### POSTNATAL ###
# LE domain 
predictormatrix[c('sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary'), # LE
                c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems','moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", # CR minus "m_education_pregnancy" that is auxiliary for postnatal variables
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS 
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", 
                   # |-> IS minus "marital_status_pregnancy" that is auxiliary for postnatal variables
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic','m_education','p_education', # CR
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                  'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'post_life_events',
                  'm_bmi_berore_pregnancy', 'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs')] <- 0          
# CR domain
predictormatrix[c('tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education'), # CR
                c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems','moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", # CR minus "m_education_pregnancy" that is auxiliary for postnatal variables
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS 
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                   # |-> IS minus "marital_status_pregnancy" that is auxiliary for postnatal variables
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                  'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'post_contextual_risk',
                  'm_bmi_berore_pregnancy', 'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs')] <- 0
  
# PR domain 
# I did not allow m_education_pregnancy and dep here because of multicollinearity
predictormatrix[c('m_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age'), # PR
                c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems','moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", # CR minus "m_education_pregnancy" that is auxiliary for postnatal variables
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS 
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                   # |-> IS minus "marital_status_pregnancy" that is auxiliary for postnatal variables
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education', # CR
                  'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'post_parental_risk',
                  'm_bmi_berore_pregnancy', 'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs')] <- 0

# IP domain
predictormatrix[c('marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend'), # IR
                c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems','moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", # CR minus "m_education_pregnancy" that is auxiliary for postnatal variables
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS 
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                  # |-> IS minus "marital_status_pregnancy" that is auxiliary for postnatal variables
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education', # CR
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                  'm_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip', # DV
                  'post_interpersonal_risk',
                  'm_bmi_berore_pregnancy', 'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs')] <- 0

# DV domain
predictormatrix[c('m_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment','sexual_behavior','rumors_or_gossip'), # DV
                c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 'health', 'unemployed', 'work_study_problems','moved_house', 'blood_loss', 'examination', 'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery', # LE
                  "financial_problems", "trouble_pay_pregnancy", "income_reduced", "housing_defects", "housing_adequacy", "housing_basic_living", # CR minus "m_education_pregnancy" that is auxiliary for postnatal variables
                  "m_depression_pregnancy", "m_anxiety_pregnancy", "m_interp_sensitivity_pregnancy", "m_violence_people", "m_violence_property", "m_criminal_record", # PS 
                  "difficulties_contacts","difficulties_partner","difficulties_family_friend","divorce_pregnancy","family_support","family_acceptance","family_affection","family_acception","family_trust","family_painful_feelings","family_decisions","family_conflict","family_decisions_problems","family_plans","family_talk_sadness", "family_talk_worries", "family_size_pregnancy", # IS
                   # |-> IS minus "marital_status_pregnancy" that is auxiliary for postnatal variables
                  'sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died','school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary', # LE
                  'tension_at_work','material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once','income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education', # CR
                  'm_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression','m_anxiety','p_anxiety','m_age','p_age', # PR
                  'marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member','conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend', # IR
                  'post_direct_victimization',
                  'm_bmi_berore_pregnancy', 'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs')] <- 0

# OPTIONAL :Quickly check the matrix to make sure it looks legit
# pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain
# low computational load, we do 60 iterations using 30 imputed datasets (following 
# Isabel's approach)
imputation <- mice(ELS_PCM_essentials, m = 30, # nr of imputed datasets
                   maxit = 60, #nr of iteration taken to impute missing values
                   seed = 310896, # set a seed for the random number generation in case i need to generate the same dataset again
                   method = meth,
                   visitSequence = VisSeq, 
                   predictorMatrix = predictormatrix)

### There were no logged events in the imputation. 
### Visual inspection of the convergence graphs showed convergence after 20 to 30 iterations.
# plot(imputation)

################### OPTIONAL CHECKS (beware: it takes time) ####################
# # Inspecting the distribution of observed and imputed values
# stripplot(imputation, pch = 20, cex = 1.2) # red dots are imputed values
# # A scatter plot is also useful to spot unusual patterns in two vars
# xyplot(imputation, pre_life_events ~ post_life_events | .imp, pch = 20, cex = 1.4)

################################################################################
#### ------------------------- complete and save -------------------------- ####
################################################################################

# I save the mids object (i.e. list of imputed datasets)
saveRDS(imputation, paste(pathtodata,'imputation_list.rds', sep = ""))

# I also save the last imputed dataset for sanity checks
ELS_PCM_imputed <- complete(imputation, 30) 
saveRDS(ELS_PCM_imputed, paste(pathtodata,'ELS_PCM_imputed.rds', sep = ""))

################################################################################
