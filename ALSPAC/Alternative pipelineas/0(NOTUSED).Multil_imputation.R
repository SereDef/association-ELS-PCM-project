
# Hi again, 
# the following code runs the imputation of missing data, necessary for 
# the analysis of the association between ELS and psycho-cardio-metabolic multi-morbidity 
# in children. 
# Ok, let's get started!

passive_imp_formula <- function(domain, add = "") {
  conc <- paste(domain, collapse = " + ")
  str <- paste0("~I( (", conc, ") / ", length(domain), ")")
  if (add != "") {
    str <- paste0("~I( (", conc, " + ", add, ") / ", length(domain)+1, ")")
  }
  return(str)
}

# Load libraries
library(mice);
library(miceadds)

# Load the dataframes with prenatal, postnatal stress and the outcome dataset
pre  <- readRDS(file.path(alspac_folder, "prenatal_stress.rsd"))
post <- readRDS(file.path(alspac_folder, "postnatal_stress.rsd"))
out  <- readRDS(file.path(alspac_folder, "PCMout_cov_aux.rsd"))

ELS <- cbind(pre, post, out)

################################################################################

# add NAs for missing time points so that there would be an equal number of time points for all variable types

ELS <-  cbind(ELS, STARTED_NURSERY_6Y = NA, STARTED_NURSERY_9Y =  NA, BURGLARY_OR_CAR_THEFT_8M = NA,
              HOUSING_ADEQUACY_8M = NA, HOUSING_ADEQUACY_21M = NA, HOUSING_ADEQUACY_5Y = NA, HOUSING_ADEQUACY_6Y = NA, HOUSING_ADEQUACY_9Y =  NA,
              HOUSING_BASIC_LIVING_8M = NA, HOUSING_BASIC_LIVING_21M = NA, HOUSING_BASIC_LIVING_5Y = NA, HOUSING_BASIC_LIVING_6Y = NA, HOUSING_BASIC_LIVING_9Y = NA,
              HOUSING_DEFECTS_8M = NA, HOUSING_DEFECTS_21M  = NA, housing_defects_2y = NA, housing_defects_4y  = NA, HOUSING_DEFECTS_5Y = NA, HOUSING_DEFECTS_6Y = NA, HOUSING_DEFECTS_9Y = NA,
              NEIGHBOURHOOD_PROBLEMS_8M = NA, NEIGHBOURHOOD_PROBLEMS_4Y = NA, NEIGHBOURHOOD_PROBLEMS_5Y = NA, NEIGHBOURHOOD_PROBLEMS_6Y = NA, NEIGHBOURHOOD_PROBLEMS_9Y = NA,
              M_ANXIETY_4Y = NA, M_ANXIETY_9Y = NA,
              M_CRUELTY_EMOTIONAL_8M = NA)


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
# all postnatal variables, when a timepoint in missing I make it up with a CAPITAL LETTERS NAME
post_LE <- c('sick_or_accident_18m','sick_or_accident_30m','sick_or_accident_3y','sick_or_accident_4y','sick_or_accident_5y','sick_or_accident_6y','sick_or_accident_9y',
             'family_member_ill_8m', 'family_member_ill_21m', 'family_member_ill_3y', 'family_member_ill_4y', 'family_member_ill_5y','family_member_ill_6y','family_member_ill_9y',
             'smbd_important_died_8m','smbd_important_died_21m','smbd_important_died_3y','smbd_important_died_4y','smbd_important_died_5y','smbd_important_died_6y','smbd_important_died_9y',
             'separated_from_parent_18m','separated_from_parent_30m','separated_from_parent_3y','separated_from_parent_4y','separated_from_parent_5y','separated_from_parent_6y','separated_from_parent_8y',
             'moved_18m','moved_30m','moved_3y','moved_4y','moved_5y','moved_6y','moved_9y',
             'pet_died_18m','pet_died_30m','pet_died_3y','pet_died_4y','pet_died_5y','pet_died_6y','pet_died_9y',
             'started_nursery_18m', 'started_nursery_30m', 'started_nursery_3y', 'started_nursery_4y', 'started_nursery_5y', 'STARTED_NURSERY_6Y','STARTED_NURSERY_9Y',
             'acquired_new_parent_18m', 'acquired_new_parent_30m', 'acquired_new_parent_3y', 'acquired_new_parent_4y', 'acquired_new_parent_5y', 'acquired_new_parent_6y', 'acquired_new_parent_8y',
             'change_carer_18m', 'change_carer_30m', 'change_carer_3y', 'change_carer_4y', 'change_carer_5y', 'change_carer_6y', 'change_carer_8y',
             'friend_relative_ill_8m', 'friend_relative_ill_21m', 'friend_relative_ill_3y', 'friend_relative_ill_4y', 'friend_relative_ill_5y', 'friend_relative_ill_6y', 'friend_relative_ill_9y',
             'partner_died_8m', 'partner_died_21m', 'partner_died_3y', 'partner_died_4y', 'partner_died_5y', 'partner_died_6y', 'partner_died_9y',
             'BURGLARY_OR_CAR_THEFT_8M', 'burglary_or_car_theft_21m', 'burglary_or_car_theft_3y', 'burglary_or_car_theft_4y', 'burglary_or_car_theft_5y', 'burglary_or_car_theft_6y', 'burglary_or_car_theft_9y',
             'separated_from_smbd_18m', 'separated_from_smbd_30m', 'separated_from_smbd_3y', 'separated_from_smbd_4y', 'separated_from_smbd_5y', 'separated_from_smbd_6y', 'separated_from_smbd_8y', 
             'lost_best_friend_8y',
             'new_sibling_18m', 'new_sibling_30m', 'new_sibling_3y', 'new_sibling_4y', 'new_sibling_5y', 'new_sibling_6y', 'new_sibling_9y',
             'ch_had_fright_18m', 'ch_had_fright_30m', 'ch_had_fright_3y', 'ch_had_fright_4y', 'ch_had_fright_5y', 'ch_had_fright_6y', 'ch_had_fright_8y')

post_CR <- c('homeless_childhood_8m', 'homeless_childhood_21m', 'homeless_childhood_3y', 'homeless_childhood_4y', 'homeless_childhood_5y', 'homeless_childhood_6y', 'homeless_childhood_9y',
             'major_financial_problems_8m', 'major_financial_problems_21m', 'major_financial_problems_3y', 'major_financial_problems_4y', 'major_financial_problems_5y', 'major_financial_problems_6y', 'major_financial_problems_9y',
             'income_reduced_8m', 'income_reduced_21m', 'income_reduced_3y', 'income_reduced_4y', 'income_reduced_5y', 'income_reduced_6y', 'income_reduced_9y',
             'unemployed_8m','unemployed_21m', 'unemployed_3y', 'unemployed_4y', 'unemployed_5y', 'unemployed_6y', 'unemployed_9y',
             'HOUSING_ADEQUACY_8M','HOUSING_ADEQUACY_21M','housing_adequacy_2y', 'housing_adequacy_4y','HOUSING_ADEQUACY_5Y','HOUSING_ADEQUACY_6Y','HOUSING_ADEQUACY_9Y',
             'HOUSING_BASIC_LIVING_8M','HOUSING_BASIC_LIVING_21M','housing_basic_living_2y', 'housing_basic_living_4y','HOUSING_BASIC_LIVING_5Y', 'HOUSING_BASIC_LIVING_6Y','HOUSING_BASIC_LIVING_9Y',
             'HOUSING_DEFECTS_8M','HOUSING_DEFECTS_21M','housing_defects_2y', 'housing_defects_4y','HOUSING_DEFECTS_5Y','HOUSING_DEFECTS_6Y','HOUSING_DEFECTS_9Y',
             'm_education',
             'p_education',
             'NEIGHBOURHOOD_PROBLEMS_8M','neighbourhood_problems_21m', 'neighbourhood_problems_3y','NEIGHBOURHOOD_PROBLEMS_4Y','NEIGHBOURHOOD_PROBLEMS_5Y','NEIGHBOURHOOD_PROBLEMS_6Y','NEIGHBOURHOOD_PROBLEMS_9Y')

post_PR <- c('work_problems_8m','work_problems_21m','work_problems_3y','work_problems_4y','work_problems_5y','work_problems_6y','work_problems_9y',
             'criminal_record_parent_8m', 'criminal_record_parent_21m', 'criminal_record_parent_3y', 'criminal_record_parent_4y', 'criminal_record_parent_5y', 'criminal_record_parent_6y', 'criminal_record_parent_9y',
             'miscarriage_or_abortion_8m', 'miscarriage_or_abortion_21m', 'miscarriage_or_abortion_3y', 'miscarriage_or_abortion_4y', 'miscarriage_or_abortion_5y', 'miscarriage_or_abortion_6y', 'miscarriage_or_abortion_9y',
             'm_attempted_suicide_8m', 'm_attempted_suicide_21m', 'm_attempted_suicide_3y', 'm_attempted_suicide_4y', 'm_attempted_suicide_5y', 'm_attempted_suicide_6y', 'm_attempted_suicide_9y',
             'm_age',
             'p_age',
             'm_depression_8m', 'm_depression_21m', 'm_depression_3y', 'm_depression_4y', 'm_depression_5y', 'm_depression_6y', 'm_depression_9y',
             'p_depression_8m', 'p_depression_21m', 'p_depression_3y', 'p_depression_4y', 'p_depression_5y', 'p_depression_6y', 'p_depression_9y',
             'm_anxiety_8m', 'm_anxiety_21m', 'm_anxiety_3y', 'M_ANXIETY_4Y', 'm_anxiety_5y', 'm_anxiety_6y','M_ANXIETY_9Y',
             'p_anxiety_8m', 'p_anxiety_21m', 'p_anxiety_3y', 'p_anxiety_4y', 'p_anxiety_5y', 'p_anxiety_6y', 'p_anxiety_9y')

post_IR <- c('divorce_8m', 'divorce_21m', 'divorce_3y', 'divorce_4y', 'divorce_5y', 'divorce_6y', 'divorce_9y',
             'p_rejected_child_8m', 'p_rejected_child_21m', 'p_rejected_child_3y', 'p_rejected_child_4y', 'p_rejected_child_5y', 'p_rejected_child_6y', 'p_rejected_child_9y',
             'p_went_away_8m', 'p_went_away_21m', 'p_went_away_3y', 'p_went_away_4y', 'p_went_away_5y', 'p_went_away_6y', 'p_went_away_9y',
             'conflict_in_family_8m', 'conflict_in_family_21m', 'conflict_in_family_3y', 'conflict_in_family_4y', 'conflict_in_family_5y', 'conflict_in_family_6y', 'conflict_in_family_9y',
             'conflict_family_violence_8m', 'conflict_family_violence_21m', 'conflict_family_violence_3y', 'conflict_family_violence_4y', 'conflict_family_violence_5y', 'conflict_family_violence_6y', 'conflict_family_violence_9y',
             'm_new_partner_8m', 'm_new_partner_21m', 'm_new_partner_3y', 'm_new_partner_4y', 'm_new_partner_5y', 'm_new_partner_6y', 'm_new_partner_9y',
             'argued_fam_friends_8m', 'argued_fam_friends_21m', 'argued_fam_friends_3y', 'argued_fam_friends_4y', 'argued_fam_friends_5y', 'argued_fam_friends_6y', 'argued_fam_friends_9y')

post_DV <- c('bullying_8y',
             'physical_violence_18m', 'physical_violence_30m', 'physical_violence_3y', 'physical_violence_4y', 'physical_violence_5y', 'physical_violence_6y', 'physical_violence_9y',
             'sexual_abuse_18m', 'sexual_abuse_30m', 'sexual_abuse_3y', 'sexual_abuse_4y', 'sexual_abuse_5y', 'sexual_abuse_6y', 'sexual_abuse_8y',
             'p_cruelty_physical_8m', 'p_cruelty_physical_21m', 'p_cruelty_physical_3y', 'p_cruelty_physical_4y', 'p_cruelty_physical_5y', 'p_cruelty_physical_6y', 'p_cruelty_physical_9y',
             'm_cruelty_physical_8m', 'm_cruelty_physical_21m', 'm_cruelty_physical_3y', 'm_cruelty_physical_4y', 'm_cruelty_physical_5y', 'm_cruelty_physical_6y', 'm_cruelty_physical_9y',
             'p_cruelty_emotional_8m', 'p_cruelty_emotional_21m', 'p_cruelty_emotional_3y', 'p_cruelty_emotional_4y', 'p_cruelty_emotional_5y', 'p_cruelty_emotional_6y', 'p_cruelty_emotional_9y',
             'M_CRUELTY_EMOTIONAL_8M', 'm_cruelty_emotional_21m', 'm_cruelty_emotional_3y', 'm_cruelty_emotional_4y', 'm_cruelty_emotional_5y', 'm_cruelty_emotional_6y', 'm_cruelty_emotional_9y')

outcomes <- c('intern_score_z', 'fat_mass_z', 'risk_groups')
covars   <- c('sex', 'age_child', 'm_bmi_before_pregnancy', 'm_smoking', 'm_drinking')
auxil    <- c('m_dep_cont_pregnancy', 'p_dep_cont_pregnancy', # for postnatal only
              'm_bmi_7yrs', 'm_dep_cont_childhood', 'p_dep_cont_childhood', # for prenatal only
              'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_age_cont')
exclusion_criteria <- c('pre_percent_missing', 'post_percent_missing', 'twin')

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
         'pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk', 
         'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization',
         # cumulative prenatal and postnatal stress exposure
         'prenatal_stress', "postnatal_stress",
         # outcome variables and covariates + additional auxiliary variables for imputation
         outcomes, covars, auxil, exclusion_criteria)

ELS <- ELS[, vars]
#------------------------------------------------------------------------------#

# make some fake data
# ELS <- data.frame(replicate(length(vars), sample(c(0:1, NA), 1000, rep=TRUE)))
# colnames(ELS) <- vars
#ELS[, which(grepl("[A-Z]", names(ELS)))] <- NA

post_LE_t <- c('sick_or_accident',
         'family_member_ill',
         'smbd_important_died',
         'separated_from_parent',
         'moved',
         'pet_died',
         'started_nursery',
         'acquired_new_parent',
         'change_carer',
         'friend_relative_ill',
         'partner_died',
         'burglary_or_car_theft',
         'separated_from_smbd',
         'lost_best_friend_8y',
         'new_sibling',
         'ch_had_fright')
post_CR_t <-c( 'homeless_childhood',
         'major_financial_problems' ,
         'income_reduced',
         'unemployed',
         'housing_adequacy',
         'housing_basic_living',
         'housing_defects',
         'm_education' ,
         'p_education' ,
         'neighbourhood_problems')
post_PR_t <-c('work_problems',
         'criminal_record_parent',
         'miscarriage_or_abortion',
         'm_attempted_suicide',
         'm_age',
         'p_age',
         'm_depression',
         'p_depression',
         'm_anxiety',
         'p_anxiety')
post_IR_t <- c('divorce',
         'p_rejected_child',
         'p_went_away',
         'conflict_in_family',
         'conflict_family_violence',
         'm_new_partner',
         'argued_fam_friends')
post_DV_t <- c('bullying_8y',
         'physical_violence',
         'sexual_abuse',
         'p_cruelty_physical',
         'm_cruelty_physical',
         'p_cruelty_emotional',
         'm_cruelty_emotional')

nms <- c(post_LE_t, post_CR_t, post_PR_t, post_IR_t, post_DV_t)

pos = list()

for (v in 1:length(nms)) {
  
  ps = which(grepl(paste0(nms[v], "_[0-9]"), names(ELS), ignore.case = T))
  if (length(ps) > 1) {
    
    pos[[v]] <- c(ps)
  }
}

pos <- pos[-which(sapply(pos, is.null))]

dl <- reshape(data = ELS, idvar = "IDC", timevar = "Time", direction="long", 
              varying = pos,
              v.names = nms,
              sep="_")

# The long dataset is not ordered yet by ID and Time. This can be done by using the order function.

ELSlong <- dl[order(dl$IDC, dl$Time), ]

#------------------------------------------------------------------------------#
##---------------------------- Imputation model ----------------------------- ##
#------------------------------------------------------------------------------#

# We started with a dry run to specify the default arguments.
imp0 <- mice(ELSlong, maxit = 0, remove.collinear = F)
# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- imp0$method
meth[c(which(names(meth) == 'sick_or_accident'):ncol(ELSlong))] <- '2l.bin'
# set the imputation method to Generalized Linear Mixed model (20-50 iterations 
# are recommended and may be hard to run in small datasets)
# We use passive imputation for the domain scores. This means that the indicator items  
# are imputed first, and then, using these complete items, mean domain scores are 
# derived by the formula specified below.

meth['pre_life_events']         <- passive_imp_formula(pre_LE)
meth['pre_contextual_risk']     <- passive_imp_formula(pre_CR)
meth['pre_parental_risk']       <- passive_imp_formula(pre_PR)
meth['pre_interpersonal_risk']  <- passive_imp_formula(pre_IR)

# below hasn't been changed to ALSPAC yet ####### NO to fix
meth['post_life_events']          <- passive_imp_formula(post_LE_t)
meth['post_contextual_risk']      <- passive_imp_formula(post_CR_t)
meth['post_parental_risk']        <- passive_imp_formula(post_PR_t)
meth['post_interpersonal_risk']   <- passive_imp_formula(post_IR_t)
meth['post_direct_victimization'] <- passive_imp_formula(post_DV_t)

# We also use passive imputation for the period specific cumulative ELS scores.
meth['prenatal_stress'] <- "~I( pre_life_events + pre_contextual_risk + pre_parental_risk + pre_interpersonal_risk )"
meth['postnatal_stress'] <- "~I( post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization )"

# We are going to need a different set of predictors for the different variables we impute 
# so let's define them using the predictormatrix, that gives these instructions to mice
predictormatrix <- imp0$predictorMatrix

# Do not use cumulative pre and postnatal scores as predictors
predictormatrix[, c('prenatal_stress', 'postnatal_stress') ]  <- 0
# Do not impute nor use IDC, any of the outcomes, exclusion criteria or age_child 
# as predictors (no reason to believe age at outcome is associated with missingness)
predictormatrix[, c("IDC", outcomes, "age_child", exclusion_criteria)] <- 0
predictormatrix[c("IDC", outcomes, "age_child", exclusion_criteria, 'Time'), ] <- 0


           ### Impute auxiliary variables and covariates ###

# To prevent multicollinearity, we adjust the predictor matrix such that the 
# auxiliary variables would be imputed given the domain scores and the outcomes, 
# but not by the single items.
predictormatrix[c(covars, auxil),
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE_t, post_CR_t, post_PR_t, post_IR_t, post_DV_t, # all variables for postnatal risk
                  covars, auxil, 'Time')] <- 0

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
predictormatrix[c(pre_LE), # LE
                c(pre_CR, pre_PR, pre_IR, post_LE_t, post_CR_t[!post_CR_t == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR_t, post_IR_t[!post_IR_t == 'divorce'], post_DV_t,             # divorce is auxiliary for prenatal variables 
                  'pre_life_events', 'm_bmi_before_pregnancy', auxil[1:2], 'Time')] <- 0
# CR domain
predictormatrix[c(pre_CR),
                c(pre_LE, pre_PR, pre_IR, post_LE_t, post_CR_t[!post_CR_t == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR_t, post_IR_t[!post_IR_t == 'divorce'], post_DV_t,             # divorce is auxiliary for prenatal variables
                  'pre_contextual_risk', 'm_bmi_before_pregnancy',  auxil[1:2], 'Time')] <- 0

# PR domain 
predictormatrix[c(pre_PR),
                c(pre_LE, pre_CR, pre_IR, post_LE_t, post_CR_t[!post_CR_t == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR_t, post_IR_t[!post_IR_t == 'divorce'], post_DV_t,             # divorce is auxiliary for prenatal variables
                  'pre_parental_risk', 'm_bmi_before_pregnancy',  auxil[1:2], 'Time')] <- 0

# IR domain
predictormatrix[c(pre_IR),
                c(pre_LE, pre_CR, pre_PR, post_LE_t, post_CR_t[!post_CR_t == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR_t, post_IR_t[!post_IR_t == 'divorce'], post_DV_t,             # divorce is auxiliary for prenatal variables
                  'pre_interpersonal_risk', 'm_bmi_before_pregnancy',  auxil[1:2], 'Time')] <- 0

                                ### POSTNATAL ###
# LE domain 
predictormatrix[c(post_LE_t),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'divorce_pre'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_CR_t, post_PR_t, post_IR_t, post_DV_t,
                  'post_life_events', 'm_bmi_before_pregnancy',  auxil[3:5])] <- 0       
# CR domain
predictormatrix[c(post_CR_t),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'divorce_pre'], #  divorce_pre is auxiliary for postnatal variables
                  post_LE_t, post_PR_t, post_IR_t, post_DV_t,
                  'post_contextual_risk', 'm_bmi_before_pregnancy',  auxil[3:5])] <- 0

# PR domain 
predictormatrix[c(post_PR_t),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'divorce_pre'], #  divorce_pre is auxiliary for postnatal variables
                  post_LE_t, post_CR_t, post_IR_t, post_DV_t,
                  'post_parental_risk', 'm_bmi_before_pregnancy',  auxil[3:5])] <- 0

# IR domain
predictormatrix[c(post_IR_t),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'divorce_pre'], #  divorce_pre is auxiliary for postnatal variables
                  post_LE_t, post_CR_t, post_PR_t, post_DV_t,
                  'post_interpersonal_risk', 'm_bmi_before_pregnancy',  auxil[3:5])] <- 0

# DV domain 
predictormatrix[c(post_DV_t),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pre'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'divorce_pre'], #  divorce_pre is auxiliary for postnatal variables
                  post_LE_t, post_CR_t, post_PR_t, post_IR_t,
                  'post_direct_victimization', 'm_bmi_before_pregnancy',  auxil[3:5])] <- 0

# Define the cluster variable and set Time as a first level independent variable
predictormatrix[c(post_LE_t[!post_LE_t == 'lost_best_friend_8y'], 
                  post_CR_t[!post_CR_t %in% c('m_education', 'p_education')], 
                  post_PR_t[!post_PR_t %in% c('m_age', 'p_age')], post_IR_t, 
                  post_DV_t[!post_DV_t == 'bullying_8y']), c('IDC', 'Time')] <- matrix(rep(c(-2, 2),each=44),nrow=44)

# This can help
mat_lvl1 <- function(dim) {
  m <- matrix(2, nrow = dim, ncol = dim)
  diag(m) = 0
  return(m)
}

# Define level 1 independent variables (imputation model with a fixed effect, random intercept and random slope.)
predictormatrix[c(post_LE_t[!post_LE_t == 'lost_best_friend_8y']), 
                c(post_LE_t[!post_LE_t == 'lost_best_friend_8y'])] <- mat_lvl1(length(post_LE_t)-1)
predictormatrix[c(post_CR_t[!post_CR_t %in% c('m_education', 'p_education')]), 
                c(post_CR_t[!post_CR_t %in% c('m_education', 'p_education')])] <- mat_lvl1(length(post_CR_t)-2)
predictormatrix[c(post_PR_t[!post_PR_t %in% c('m_age', 'p_age')]), 
                c(post_PR_t[!post_PR_t %in% c('m_age', 'p_age')])] <- mat_lvl1(length(post_PR_t)-2)
predictormatrix[c(post_IR_t), c(post_IR_t)] <- mat_lvl1(length(post_IR_t))
predictormatrix[c(post_DV_t[!post_DV_t == 'bullying_8y']), 
                c(post_DV_t[!post_DV_t == 'bullying_8y'])] <- mat_lvl1(length(post_DV_t)-1)

# pheatmap::pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain
# low computational load, we do 60 iterations using 30 imputed datasets (following 
# Isabel's approach)
imp <- mice(ELSlong, m = 2, # nr of imputed datasets
                   maxit = 2, #nr of iteration taken to impute missing values
                   seed = 310896, # set a seed for the random number generation in case i need to generate the same dataset again
                   method = meth,
                   visitSequence = VisSeq, 
                   predictorMatrix = predictormatrix)

warnings()
imp$loggedEvents


# # complete
# imp.data <- as.list(1:5)
# for(i in 1:5){
#   imp.data[[i]] <- complete(m.out, action=i)
# }
# 
# # reshape
# imp.data <- lapply(imp.data, melt, id=c("ID","GROUP"))
# 
# # analyse
# imp.fit <- lapply(imp.data, FUN=function(x){
#   lmer(value ~ as.numeric(variable)+(1|ID), data=x) 
# })
# imp.res <- sapply(imp.fit, fixef)

