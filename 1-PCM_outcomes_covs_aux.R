
# Hi there, 
# the following code is collecting and merging all necessary variables for the analysis
# of the association between ELS and psycho-cardio-metabolic multi-morbidity in children 
# This includes the outcomes of interest (internalizing problems and cardio-metabolic
# risk), the covariates that are going to be used as well as the auxiliary variables 
# used in the imputation of the final dataset. 
# It does not include data used to build the ELS score exposure, for which you can 
# check out https://github.com/SereDef/cumulative-ELS-score

## For this script you are going to need the following datasets from datamanagement
# CHILDCBCL9_10082016.sav, GR1093-E1_CBCL_18062020.sav, 
# CHILDFATMASS9_13092016.sav, CHILDDXATOTALBODY13_18102021.sav,
# CHILD-ALLGENERALDATA_07072020.sav, GEDRAGSGROEP_MaternalDrinking_22112016.sav, 
# MATERNALSMOKING_22112016.sav, GEDRAGSGROEP_MaternalDrinking_22112016.sav, 
# MOTHERANTHROPOMETRY_18022013.sav, GR1003-BSI D1_22112016.sav, and 
# BSI 3 years of age_GR1065 G1-GR1066 C1_22112016.sav, GR1004-BSI G1_22112016.sav

# -----------------------------------------------------------------------------#
##                        Ok, we are good to go!                              ##
# -----------------------------------------------------------------------------#

source('0-Functions.R') # where readquick() is defined
# When prompted with a window, please navigate to the directory where your input 
# files are stored. Choose any file in the directory to continue.
# Note: the code assumes that all (raw) data is stored in ONE folder.

################################################################################
#### ------------------ INTERNALIZING PROBLEMS ( @13 ) -------------------- ####
################################################################################

# The internalizing sub-scale of the Child Behavior Checklist (CBCL 6-18) (Achenbach, 1999)
# is an empirically based score derived from the widely used parent-report questionnaire. 
# It contains 32 items rated on a 3-point scale (0 = "Not true", 1 = "Somewhat or 
# sometimes true", 2 = "Very or often true") and referred to the past 6 months. 
# Hence the score ranges from 0 to 64. 

# Read in the dataset
cbcl09 <- readquick("CHILDCBCL9_10082016.sav") # 9901 obs. of 627 vars
cbcl13 <- readquick("GR1093-E1_CBCL_18062020.sav") # 9901 obs. of 331 vars

cbcl <- merge(cbcl09, cbcl13, by = "idc")

internalizing <- data.frame('IDC' = cbcl09$idc,
                            'age_child_cbcl'   = cbcl$agechild_gr1093,
                            'intern_score_09'  = cbcl$sum_int_9m, # weighted sum score internalizing scale (allowing 25% missing)
                            'intern_score_13'  = cbcl$sum_int_14)
# Calculation (based on SPSS script available on the Generation R V: drive): 
# if 24 out of 32 items are available (i.e. 75%), sum of the item
# scores * (32 / 32 - number of missing values). 

# One alternative could be to use the anxious/depressed empirical sub-scale (mother report)
# However I did give it a shot and it does not seem to perfom better than the internalizing one.

################################################################################
#### -------------------------- FAT MASS ( @13 ) -------------------------- ####
################################################################################

# Although the original plan for this project involved a more comprehensive measure
# of child metabolic syndrome (see CMR_score.R file for specifics), fat mass was
# the only metabolic risk variable that showed appreciable variance in such young 
# children and the highest correlation with internalizing (small but significant r =.12)
# It was also selected on the base of data availability both cross-sectionally and 
# for future longitudinal assessment. 

# Read in the datasets
dxa09  <- readquick("CHILDFATMASS9_13092016.sav") # 5862 obs of 28 vars
dxa13  <- readquick("CHILDDXATOTALBODY13_18102021.sav") # 5862 obs of 28 vars

dxa <- merge(dxa09, dxa13, by = 'idc')

# Select only the necessary measures
fat_mass <- data.frame('IDC' = dxa$idc,
             'age_child_dxa' = dxa$agechild13, # age at visit
          'andr_fat_mass_09' = dxa$fat_mass_androidchild9,
          'andr_fat_mass_13' = dxa$fat_mass_androidchild13,
              'total_fat_09' = dxa$fat_mass_totalchild9,
              'total_fat_13' = dxa$fat_mass_totalchild13,
        'tot_fat_percent_09' = dxa$avg_percent_fat*100, # Re-scale to match the age13 variable 
        'tot_fat_percent_13' = dxa$totalfatpercentagechild13 )

# fat_mass$fmi <- fat_mass$tot_fat / (fat_mass$height^2)
# hist(fat_mass$fmi)

################################################################################
# merge the two main (mental and physical) outcomes with child sex into one dataset
PCM_outcome <- merge(internalizing, fat_mass, by = 'IDC',  all = T)

################################################################################
#### ---------------------------- COVARIATES ------------------------------ ####
################################################################################

# Variables that will be used in the covariate models of this project are those 
# marked with ###. they include: 'sex', 'age_child', 'm_smoking', 'm_drinking' 
# and 'm_bmi_berore_pregnancy'.

# For the other demographic auxiliary variables (used for imputation): when they 
# were assessed both prenatally and postnatally, both measures are included.
# Auxiliary variables for this project include: 'ethnicity', 'm_age', parity', 
# 'm_smoking', 'gest_age_birth', 'gest_weight', 'm_bmi_pregnancy', 'm_bmi_5yrs',
# 'sex', 'm_dep_pregnancy', 'p_dep_pregnacy', "m_dep_3yrs", "p_dep_3yrs"

# ------------------------------------------------------------------------------
### AGE of the child
# Combine age of the child measured during first visit and at CBCL administration
# This value will serve as a covariate in the first adjusted model.
PCM_outcome$age_child <- (PCM_outcome$age_child_dxa + PCM_outcome$age_child_cbcl) / 2
    # OPTIONAL: check age difference between measurements
    # plot(PCM_outcome$age_child_dxa, PCM_outcome$age_child_cbcl)
    # summary(PCM_outcome$age_child_dxa - PCM_outcome$age_child_cbcl)

#-------------------------------------------------------------------------------
### MATERNAL SMOKING during pregnancy 
smokingv1 <- readquick("MATERNALSMOKING_22112016.sav") #  9778 obs of 11 variables

smoking <- data.frame('IDM' = smokingv1$idm, 
                'm_smoking' = smokingv1$smoke_all ) # (1) never a smoker; 
                                                    # (2) smoked until pregnancy was known (i.e., first trimester only); 
                                                    # (3) continued smoking during pregnancy.
#-------------------------------------------------------------------------------
### MATERNAL ALCOHOL CONSUMPTION during pregnancy
drinkingv1 <- readquick("GEDRAGSGROEP_MaternalDrinking_22112016.sav") # drinkingv2 <- readquick("MATERNALALCOHOL_22112016.sav") # old variable

drinking <- data.frame('IDM' = drinkingv1$idm, 
                'm_drinking' = drinkingv1$mdrink_updated) # (0) never; 
                                                          # (1) until pregnancy was known (i.e., first trimester only); 
                                                          # (2) continued during pregnancy occasionally;
                                                          # (3) continued during pregnancy frequently.
#-------------------------------------------------------------------------------
## Other variables
child_general <- readquick("CHILD-ALLGENERALDATA_07072020.sav") # 9901 obs of 122 

# Ethnicity recode – dichotomized into: dutch and non-dutch;
child_general$ethnicity <- ifelse(is.na(child_general$ethnfv2), NA,
                                        ifelse(child_general$ethnfv2 == 1, 0, 1))  # Dutch = 0, non-Dutch = 1
                           # American, western (300) Asian, western (500) European (700), Oceanie (800)
                           # Indonesian (2), Cape Verdian (3), Maroccan (4) Dutch Antilles (5) Surinamese 
                           # (6) Turkish (7) African (200), American, non western (400), Asian, non western (600)

general_cov_aux <- data.frame('IDC' = child_general$idc, 
                              'IDM' = child_general$idm,
                              'sex' = as.factor(child_general$gender),    ### 1 = boy; 2 = girl.
                        'ethnicity' = as.factor(child_general$ethnicity), ### 0 = Dutch, 1 = non-Dutch
           'm_bmi_before_pregnancy' = child_general$bmi_0,     ### self-reported Maternal BMI
                             'twin' = child_general$twin,      # (exclusion criteria)
                           'mother' = child_general$mother,    # mother id used to identify siblings (for exclusion)
                           'parity' = child_general$parity,    # parity (used for imputation)
                   'gest_age_birth' = child_general$gestbir,   # gestational age at birth (used for imputation)
                      'gest_weight' = child_general$weight,    # gestational weight (used for imputation)
                  'm_bmi_pregnancy' = child_general$bmi_1,     # maternal BMI during pregnancy (used for imputation)
                       'm_age_cont' = child_general$age_m_v2)  # maternal age at intake (used for imputation) 

#-------------------------------------------------------------------------------
# Maternal and paternal depression during pregnancy and at age 3 
bsi_pregnancy_m <- readquick('GR1003-BSI D1_22112016.sav') # 9778 obs of 261 vars
bsi_pregnancy_p <- readquick('GR1004-BSI G1_22112016.sav') # 9778 obs of 261 vars
bsi_3yrs        <- readquick('BSI 3 years of age_GR1065 G1-GR1066 C1_22112016.sav') # 9897 obs of 49 vars

# Depression during pregnancy
dep_pregnancy_m <- data.frame('IDM' = bsi_pregnancy_m$idm, 'm_dep_cont_pregnancy' = bsi_pregnancy_m$dep  )
dep_pregnancy_p <- data.frame('IDM' = bsi_pregnancy_p$idm, 'p_dep_cont_pregnancy' = bsi_pregnancy_p$dep_p)

# Merge it with the previous dataset
general_cov_aux <- Reduce(function(x, y) merge(x = x, y = y, by = 'IDM',  all.x = TRUE),
                                list(general_cov_aux, smoking, drinking, dep_pregnancy_m, dep_pregnancy_p))

# Depression @ 3y: Items 9, 16, 17, 18, 35, and 50
d <- data.frame(bsi_3yrs$g0100365, bsi_3yrs$g0100665, bsi_3yrs$g0100765, bsi_3yrs$g0100865, bsi_3yrs$g0101365, bsi_3yrs$g0102165, # mother report
                bsi_3yrs$c0100366, bsi_3yrs$c0100666, bsi_3yrs$c0100766, bsi_3yrs$c0100866, bsi_3yrs$c0101366, bsi_3yrs$c0102166) # father report
n_items_m <- rowSums(!is.na(d[,1:6])); n_items_p <- rowSums(!is.na(d[,7:12]))
bsi_3yrs$m_dep_cont_3yrs <- ifelse(n_items_m >= 5, (rowSums(d[,1:6]) / n_items_m)-1, NA)
bsi_3yrs$p_dep_cont_3yrs <- ifelse(n_items_p >= 5, (rowSums(d[,7:12])/ n_items_p)-1, NA)

dep_3yrs <- bsi_3yrs[, c('idc','m_dep_cont_3yrs', 'p_dep_cont_3yrs')]; colnames(dep_3yrs)[1] <- 'IDC'
# Merge it with the previous dataset
general_cov_aux <- merge(general_cov_aux, dep_3yrs, by = 'IDC',  all.x = T)

#-------------------------------------------------------------------------------
# Maternal BMI at age 5 (used for imputation)
m_anthropometry_5yrs <- readquick('MOTHERANTHROPOMETRY_18022013.sav')

m_bmi_5yrs <- data.frame('mother' = m_anthropometry_5yrs$mother,
                     'm_bmi_5yrs' = m_anthropometry_5yrs$bmimotherf5)

# Merge with the other general variables
general_cov_aux <- merge(general_cov_aux, m_bmi_5yrs, by = 'mother', all.x = T, incomparables = NA)

################################################################################
################################################################################

# merging outcome variables, covariates and auxiliary variables in one dataset
PCM_project = merge(PCM_outcome, general_cov_aux, by = 'IDC', all = T)

################################################################################
#### --------------------------- save and run ----------------------------- ####
################################################################################

# Save the dataset in an .rds file, in the directory where the raw data are stored
saveRDS(PCM_project, file.path(pathtoresults, 'PCM_allvars.rds'))

