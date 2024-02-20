
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

################################################################################
# Load datasets
pre_risk  <- readRDS(file.path(pathtoresults, 'prenatal_stress.rds'))
post_risk <- readRDS(file.path(pathtoresults, 'postnatal_stress.rds'))
outcome   <- readRDS(file.path(pathtoresults, 'PCM_allvars.rds'))

# merge them
ELSPCM <- Reduce(function(x,y) merge(x = x, y = y, by = c('IDC', 'IDM'), all = T),
                  list(pre_risk, post_risk, outcome) ) 

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
                     domains,
                     # cumulative prenatal and postnatal stress exposure
                    'prenatal_stress', 'postnatal_stress',
                     # outcome variables and covariates + additional auxiliary variables for imputation
                    outcomes_13y, outcomes_09y, covars, 'ethn_cont', auxil, exclusion_criteria)]

siblings_to_exclude <- flowchart(ELSPCM_essentials)

#------------------------------------------------------------------------------#
##---------------------------- Imputation model ----------------------------- ##
#------------------------------------------------------------------------------#

# We start with a dry run to specify the default arguments: i.e. setting the imputation
# method to predictive mean matching (PMM)
imp0 <- mice(ELSPCM_essentials, maxit = 0, defaultMethod = rep('pmm',4) ) 
# * PMM imputes a value randomly from a set of observed values whose predicted values 
#   are closest to the predicted value of the specified regression model. PMM has been 
#   said to perform quite well under circumstance where the categorical data is sparse 
#   (Van Buuren, 2018).

meth <- imp0$method
# We use passive imputation for the domain scores. This means that the indicator items  
# are imputed first, and then, using these complete items, mean domain scores are 
# derived by the formula specified using passive_imp_formula (see 0-Functions.R)
meth['pre_life_events']           <- passive_imp_formula(pre_LE)
meth['pre_contextual_risk']       <- passive_imp_formula(pre_CR)
meth['pre_parental_risk']         <- passive_imp_formula(pre_PR, add = "m_age")
meth['pre_interpersonal_risk']    <- passive_imp_formula(pre_IR)
meth['post_life_events']          <- passive_imp_formula(post_LE)
meth['post_contextual_risk']      <- passive_imp_formula(post_CR)
meth['post_parental_risk']        <- passive_imp_formula(post_PR)
meth['post_interpersonal_risk']   <- passive_imp_formula(post_IR)
meth['post_direct_victimization'] <- passive_imp_formula(post_DV)

# We also use passive imputation for the period specific cumulative ELS scores.
meth['prenatal_stress']  <- "~I( pre_life_events + pre_contextual_risk + pre_parental_risk + pre_interpersonal_risk )"
meth['postnatal_stress'] <- "~I( post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization )"

# We are going to need a different set of predictors for the different variables to impute
# so we define them by manually modifying the predictormatrix
predictormatrix <- imp0$predictorMatrix

# Do not impute nor use IDC, any of the exclusion criteria
predictormatrix[, c("IDC", exclusion_criteria) ] <- 
  predictormatrix[c("IDC", exclusion_criteria),] <- 0
# Leave the domain and total ELS scores to the passive imputation
predictormatrix[c(domains, 'prenatal_stress', 'postnatal_stress'), ] <- 0
# Do not use cumulative pre and postnatal scores as predictors in the items imputation
# nor the outcome timepoints
predictormatrix[c(pre_LE, pre_CR, pre_PR, pre_IR,# all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV), # all variables for postnatal risk
                c('prenatal_stress', 'postnatal_stress', outcomes_13y, outcomes_09y) ]  <- 0

               ### Impute auxiliary variables and covariates ###
# To prevent multicollinearity, auxiliary variables are imputed given the domain scores 
# and the outcomes, but not the single items.
predictormatrix[c(outcomes_09y, covars, 'ethn_cont', auxil),
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  'prenatal_stress', 'postnatal_stress',
                  outcomes_09y, covars, auxil)] <- 0

                           ### Impute outcomes ###
# Each primary outcome is imputed based on total ELS scores, auxiliary variables, 
# including the previous assessment of the same outcome, and the covariates. 
predictormatrix['intern_score_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains,
                  outcomes_13y, outcomes_09y[!outcomes_09y == 'intern_score_09'])] <- 0
predictormatrix['total_fat_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains,
                  outcomes_13y, outcomes_09y[!outcomes_09y == 'total_fat_09'])] <- 0
predictormatrix['andr_fat_mass_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains,
                  outcomes_13y, outcomes_09y[!outcomes_09y == 'andr_fat_mass_09'])] <- 0
predictormatrix['tot_fat_percent_13',
                c(pre_LE, pre_CR, pre_PR, pre_IR, # all variables for prenatal risk
                  post_LE, post_CR, post_PR, post_IR, post_DV, # all variables for postnatal risk
                  domains,
                  outcomes_13y, outcomes_09y[!outcomes_09y == 'tot_fat_percent_09'])] <- 0

# All single items are then imputed such that we impute the items within a domain score given 
# the other items in that domain, the remaining domain scores and the auxiliary variables.
  # Auxiliary variables were selected because they are either related to missingness or to 
  # the domain scores themselves. When information is available prenatally and postnatally, 
  # we use the different period to minimize bias (E.G., for the imputation of prenatal items
  # we used BMI of the mother when the child was 5, and for imputation of postnatal items
  # we used BMI of the mother during pregnancy).
# This method is more efficient (lower computational load) and the multicollinearity issue is 
# also resolved. The technique has been found to reduce standard error substantially compared to 
# complete-case analysis (Plumpton et al., 2016), and it outperforms other existing techniques (Eekhout et al., 2018). 

# So, let's adjust the predictormatrix such that ‘redundant’ items were not used as a predictor.

                                  ### PRENATAL ###
# LE domain 
predictormatrix[c(pre_LE), # LE
                c(pre_CR, pre_PR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables 
                  'pre_life_events', 'm_bmi_before_pregnancy', auxil[1:3])] <- 0
# CR domain
predictormatrix[c(pre_CR),
                c(pre_LE, pre_PR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables
                  'pre_contextual_risk', 'm_bmi_before_pregnancy', auxil[1:3])] <- 0
# PR domain 
predictormatrix[c(pre_PR),
                c(pre_LE, pre_CR, pre_IR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables
                  'pre_parental_risk', 'm_bmi_before_pregnancy', auxil[1:3])] <- 0
# IR domain
predictormatrix[c(pre_IR),
                c(pre_LE, pre_CR, pre_PR, post_LE, post_CR[!post_CR == 'm_education'], # because m_education is auxiliary for prenatal variables
                  post_PR, post_IR[!post_IR == 'marital_status'], post_DV,             # marital_status is auxiliary for prenatal variables
                  'pre_interpersonal_risk', 'm_bmi_before_pregnancy', auxil[1:3])] <- 0
  
                                  ### POSTNATAL ###
# LE domain 
predictormatrix[c(post_LE),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_CR, post_PR, post_IR, post_DV,
                  'post_life_events', 'm_bmi_before_pregnancy', auxil[4:6])] <- 0 
# CR domain
predictormatrix[c(post_CR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_PR, post_IR, post_DV,
                  'post_contextual_risk', 'm_bmi_before_pregnancy', auxil[4:6])] <- 0
# PR domain 
predictormatrix[c(post_PR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_IR, post_DV,
                  'post_parental_risk', 'm_bmi_before_pregnancy', auxil[4:6])] <- 0
# IR domain
predictormatrix[c(post_IR),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_PR, post_DV,
                  'post_interpersonal_risk', 'm_bmi_before_pregnancy', auxil[4:6])] <- 0
# DV domain
predictormatrix[c(post_DV),
                c(pre_LE, pre_CR[!pre_CR == 'm_education_pregnancy'],    #  m_education_pregnancy is auxiliary for postnatal variables
                  pre_PR, pre_IR[!pre_IR == 'marital_status_pregnancy'], #  marital_status_pregnancy is auxiliary for postnatal variables
                  post_LE, post_CR, post_PR, post_IR,
                  'post_direct_victimization', 'm_bmi_before_pregnancy', auxil[4:6])] <- 0
  
# OPTIONAL : Quickly check the matrix to make sure it looks legit
# pheatmap::pheatmap(predictormatrix, cluster_rows = F, cluster_cols = F)

# visit the sequence
VisSeq <- imp0$visitSequence

# Run the actual imputation. To ensure convergence among the variables but retain low computational load,
# we do 60 iterations using 30 imputed datasets (following Isabel's approach)
imputation <- mice(ELSPCM_essentials, m = 30, # nr of imputed datasets
                   maxit = 60, #nr of iteration taken to impute missing values
                   seed = 310896, # set a seed for the random number generation in case i need to generate the same dataset again
                   method = meth,
                   visitSequence = VisSeq, 
                   predictorMatrix = predictormatrix)

### There were no logged events in the imputation. 
### Visual inspection of the convergence graphs showed convergence after 20 to 30 iterations.

# ------------------------------------------------------------------------------
# Apply exclusion criteria to select the sample for analysis: exclude participants that do not
# have enough prenatal and postnatal information (i.e. < 50% of the items in each period).
pren50cutoff <- miceadds::subset_datlist( imputation,   subset = imputation$data$pre_percent_missing < 50.0,  toclass = 'mids')
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
no_sibls <- miceadds::subset_datlist( no_twins, subset = no_twins$data$IDC %notin% siblings_to_exclude, toclass = 'mids')

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
saveRDS(imputation, file.path(pathtoresults,'imputation_list_full.rds'))
saveRDS(finalset,   file.path(pathtoresults,'imputation_list_sample.rds'))

# I also save the last imputed dataset for sanity checks
full_imputed <- complete(imputation, 30) 
ELSPCM_imputed <- complete(finalset, 30) 
ELSPCM_original <- complete(finalset, 0) 
saveRDS(full_imputed,   file.path(pathtoresults, 'imputed30_full.rds'))
saveRDS(ELSPCM_imputed, file.path(pathtoresults, 'imputed30_sample.rds'))
write.csv(ELSPCM_original, file.path(pathtoresults, 'original_sample.csv'))

################################################################################
