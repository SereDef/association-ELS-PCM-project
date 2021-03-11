
# Hi again, 
# the following code runs the imputation of missing data, necessary for 
# the analysis of the association between ELS and psycho-cardio-metabolic multi-morbidity 
# in children. 

# All you need it the file we created in the 'Dataset_contruction.R' script: the 
# "ELSPCM_dataset.rds" 
# Ok, let's get started!

#### ---------------------------- Dependencies ---------------------------- ####

# Point to the necessary libraries
library(mice);

# This will come in handy for exclusion
'%notin%' <- Negate('%in%')

# check if the path to the datasets is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

################################################################################
# Load datasets
ELSPCM <- readRDS(paste(pathtodata, 'ELSPCM_dataset.rds', sep = ""))

#------------------------------------------------------------------------------#
##---------------------------- Imputation model ----------------------------- ##
#------------------------------------------------------------------------------#

# We started with a dry run to specify the default arguments.
imp0 <- mice(ELSPCM, maxit = 0, 
             defaultMethod = rep('pmm',4)) # set the imputation method to predictive mean matching (PMM)* 
             #remove.collinear = F) # Because maternal age is measured twice in prenatal and postnatal 

# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- make.method(ELSPCM)
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

# We also use passive imputation for the period specific cumulative ELS scores.
meth['prenatal_stress'] <- "~I( pre_life_events + pre_contextual_risk + pre_personal_stress + pre_interpersonal_stress )"
meth['postnatal_stress'] <- "~I( post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization )"

# We are going to need a different set of predictors for the different variables we impute 
# so let's define them using the predictormatrix, that gives these instructions to 
# mice
predictormatrix <- imp0$predictorMatrix

# Do not use IDC as predictor:
predictormatrix[, "IDC"] <- predictormatrix["IDC",] <- 0
# Do not use age_child as a predictor (no reason to believe it is associated with missingness)
predictormatrix[, "age_child"] <- predictormatrix["age_child",] <- 0
# Do not use outcome groups as a predictor
# predictormatrix[, "risk_groups"] <- predictormatrix["risk_groups",] <- 0

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
                  'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'pre_life_events', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'pre_contextual_risk', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'pre_personal_stress', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'pre_interpersonal_stress', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'post_life_events', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'post_contextual_risk', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'post_parental_risk', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'post_interpersonal_risk', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
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
                  'post_direct_victimization', 'prenatal_stress', 'postnatal_stress', 'risk_groups',
                  'm_bmi_berore_pregnancy', 'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs')] <- 0

# OPTIONAL :Quickly check the matrix to make sure it looks legit
# pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain
# low computational load, we do 60 iterations using 30 imputed datasets (following 
# Isabel's approach)
imputation <- mice(ELSPCM, m = 30, # nr of imputed datasets
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
saveRDS(ELS_PCM_imputed, paste(pathtodata,'ELSPCM_imputed.rds', sep = ""))

################################################################################
