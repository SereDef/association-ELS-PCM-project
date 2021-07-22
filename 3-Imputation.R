
# Hi again, 
# the following script is building the final dataset that will be used for the analysis 
# of the association between ELS and psycho-cardio-metabolic multi-morbidity in children. 
# I will first run multiple imputation of missing data, then exclude participants based 
# on predefined criteria & data availability, 

# You will only need the "ELSPCM_dataset.rds" file that we build using the 
# "2-Dataset_preprocessing.R" script.
# If you want to know more about the prenatal and postnatal score construction, check 
# out https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

source('0-Functions.R') # where flowchart(), select_sibling() and the domains are defined

# Load datasets
ELSPCM <- readRDS(paste0(pathtodata, 'ELSPCM_dataset.rds'))

################################################################################
# For the sake of time efficiency (and my mental health) let's select only those 
# variables that are needed for imputation and subsequent sample selection. 
# I will use the variable names defined in 0-Functions.R. 
# Once I am at it, I also order them by domain. This is important because mice is sensitive
# to the order of the variables in the set (even though this may be a version-specific issue)

ELSPCM_essentials <- ELSPCM[, c('IDC', 
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
                     outcomes, covars, auxil, exclusion_criteria)]

siblings_to_exclude <- flowchart(ELSPCM_essentials)

#------------------------------------------------------------------------------#
##---------------------------- Imputation model ----------------------------- ##
#------------------------------------------------------------------------------#

# We started with a dry run to specify the default arguments.
imp0 <- mice(ELSPCM_essentials, maxit = 0, 
             defaultMethod = rep('pmm',4)) # set the imputation method to predictive mean matching (PMM)* 
             # remove.collinear = F)

# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- make.method(ELSPCM_essentials)
# We use passive imputation for the domain scores. This means that the indicator items  
# are imputed first, and then, using these complete items, mean domain scores are 
# derived by the formula specified below.
meth['pre_life_events']           <- "~I( (family_member_died + friend_relative_died + family_member_ill_pregnancy + admitted_to_hospital + health + unemployed + work_study_problems + moved_house + blood_loss + examination + baby_worried + pregnancy_worried + obstetric_care + pregnancy_planned + victim_robbery) / 15)" 
meth['pre_contextual_risk']       <- "~I( (financial_problems + trouble_pay_pregnancy + income_reduced + housing_defects + housing_adequacy + housing_basic_living + m_education_pregnancy + p_education_pregnancy) / 8)"
meth['pre_parental_risk']         <- "~I( (m_age + m_depression_pregnancy + m_anxiety_pregnancy + m_interp_sensitivity_pregnancy + p_depression_pregnancy + p_anxiety_pregnancy + p_interp_sensitivity_pregnancy + m_violence_people + m_violence_property + m_criminal_record + p_criminal_record) / 11)"
meth['pre_interpersonal_risk']    <- "~I( (difficulties_contacts + difficulties_partner + difficulties_family_friend + marital_status_pregnancy + divorce_pregnancy + family_support + family_acceptance + family_affection + family_acception + family_trust + family_painful_feelings + family_decisions + family_conflict + family_decisions_problems + family_plans + family_talk_sadness + family_talk_worries + family_size_pregnancy) / 18)"
meth['post_life_events']          <- "~I( (sick_or_accident + family_member_ill + smbd_important_ill + parent_died + smbd_important_died + pet_died + school_workload + repeated_grade + lost_smth_important + moved + changed_school + friend_moved + fire_or_burglary) / 13)"
meth['post_contextual_risk']      <- "~I( (material_deprivation + financial_difficulties + neiborhood_problems + trouble_pay_childhood + income_once + income_chronic + unemployed_once + unemployed_chronic + m_education + p_education) / 10)"
meth['post_parental_risk']        <- "~I( (tension_at_work + m_age + p_age + m_interpersonal_sensitivity + m_anxiety + m_depression + p_interpersonal_sensitivity + p_depression + p_anxiety) / 9)"
meth['post_interpersonal_risk']   <- "~I( (conflict_family_member + conflict_smbd_else + conflict_in_family + divorce_childhood + argument_friend + marital_problems + marital_status + family_size + m_fad_5yrs + m_fad_9yrs + p_fad_9yrs) / 11)"
meth['post_direct_victimization'] <- "~I( (physical_violence + physical_threats + sexual_harrasment + sexual_behavior + rumors_or_gossip + m_harsh_parent + p_harsh_parent + bullying) / 8)"

# We also use passive imputation for the period specific cumulative ELS scores.
meth['prenatal_stress']  <- "~I( pre_life_events + pre_contextual_risk + pre_parental_risk + pre_interpersonal_risk )"
meth['postnatal_stress'] <- "~I( post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization )"

# We are going to need a different set of predictors for the different variables we impute 
# so let's define them using the predictormatrix, that gives these instructions to 
# mice
predictormatrix <- imp0$predictorMatrix

# Do not use cumulative pre and postnatal scores as predictors
predictormatrix[, c('prenatal_stress', 'postnatal_stress') ]  <- 0
# Do not impute nor use IDC, any of the outcomes, exclusion criteria or age_child 
# as predictors (no reason to believe age at outcome is associated with missingness)
predictormatrix[, c("IDC", outcomes, "age_child", exclusion_criteria)] <- 
  predictormatrix[c("IDC", outcomes, "age_child", exclusion_criteria),] <- 0

               ### Impute auxiliary variables and covariates ###

# To prevent multicollinearity, we adjust the predictor matrix such that the 
# auxiliary variables would be imputed given the domain scores and the outcomes, 
# but not by the single items.
predictormatrix[c(covars, auxil),
                c(pre_LE, pre_CR, pre_PR, pre_IR,# all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  covars, auxil )] <- 0

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
                c(pre_CR, pre_PR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables 
                  'pre_life_events', 'm_bmi_berore_pregnancy', auxil[1:3])] <- 0
# CR domain
predictormatrix[c(pre_CR),
                c(pre_LE, pre_PR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables
                  'pre_contextual_risk', 'm_bmi_berore_pregnancy', auxil[1:3])] <- 0

# PR domain 
predictormatrix[c(pre_PR),
                c(pre_LE, pre_CR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables
                  'pre_parental_risk', 'm_bmi_berore_pregnancy', auxil[1:3])] <- 0

# IR domain
predictormatrix[c(pre_IR),
                c(pre_LE, pre_CR, pre_PR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables
                  'pre_interpersonal_risk', 'm_bmi_berore_pregnancy', auxil[1:3])] <- 0
  
                                  ### POSTNATAL ###
# LE domain 
predictormatrix[c(post_LE),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_CR, post_PR, post_IR, post_DV,
                  'post_life_events', 'm_bmi_berore_pregnancy', auxil[4:6])] <- 0       
# CR domain
predictormatrix[c(post_CR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_PR, post_IR, post_DV,
                  'post_contextual_risk', 'm_bmi_berore_pregnancy', auxil[4:6])] <- 0
  
# PR domain 
predictormatrix[c(post_PR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_IR, post_DV,
                  'post_parental_risk', 'm_bmi_berore_pregnancy', auxil[4:6])] <- 0

# IR domain
predictormatrix[c(post_IR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_PR, post_DV,
                  'post_interpersonal_risk', 'm_bmi_berore_pregnancy', auxil[4:6])] <- 0

# DV domain
predictormatrix[c(post_DV),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_PR, post_IR,
                  'post_direct_victimization', 'm_bmi_berore_pregnancy', auxil[4:6])] <- 0
  
# OPTIONAL : Quickly check the matrix to make sure it looks legit
pheatmap::pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain
# low computational load, we do 60 iterations using 30 imputed datasets (following 
# Isabel's approach)
imputation <- mice(ELSPCM_essentials, m = 30, # nr of imputed datasets
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
#xyplot(imputation, pre_life_events ~ post_life_events | .imp, pch = 20, cex = 1.4)

pren50cutoff <- miceadds::subset_datlist( imputation,   subset = imputation$data$pre_percent_missing < 50.0,  toclass = 'mids')
post50cutoff <- miceadds::subset_datlist( pren50cutoff, subset = pren50cutoff$data$post_percent_missing < 50.0 )
# Not specifying toclass argument in the last call transforms mids object into a datalist object.

# Standardize prenatal and postnatal stress to obtain standard betas from regression
sdatlist <- miceadds::scale_datlist(post50cutoff, orig_var = c('prenatal_stress', 'postnatal_stress'), 
                                    trafo_var = paste0( c('prenatal_stress', 'postnatal_stress'), "_z") )

# Reconvert back to mids object
post50cutoff <- miceadds::datlist2mids(sdatlist)

out_int      <- miceadds::subset_datlist( post50cutoff, subset = !is.na(post50cutoff$data$intern_score_z),  toclass = 'mids')
out_fat      <- miceadds::subset_datlist( out_int,      subset = !is.na(out_int$data$fat_mass_z),  toclass = 'mids')
no_twins     <- miceadds::subset_datlist( out_fat,      subset = out_fat$data$twin == 0,  toclass = 'mids') 
finalset     <- miceadds::subset_datlist( no_twins,     subset = no_twins$data$IDC %notin% siblings_to_exclude, toclass = 'mids')

################################################################################
#### ------------------------- complete and save -------------------------- ####
################################################################################

# I save the mids object (i.e. list of imputed datasets)
saveRDS(post50cutoff, paste0(pathtodata,'imputation_list_ELS.rds'))
saveRDS(finalset, paste0(pathtodata,'imputation_list_ELSPCM.rds'))

# I also save the last imputed dataset for sanity checks
ELS_PCM_imputed <- complete(finalset, 30) 
saveRDS(ELS_PCM_imputed, paste0(pathtodata,'ELSPCM_imputed.rds'))

################################################################################
