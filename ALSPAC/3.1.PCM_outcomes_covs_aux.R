# Hi there, 
# the following code is collecting and merging all necessary variables for the analysis
# of the association between ELS and psycho-cardio-metabolic multi-morbidity in children 
# This includes the outcomes of interest (internalizing problems and cardio-metabolic
# risk), the covariates that are going to be used as well as the auxiliary variables 
# used in the imputation of the final dataset. 

# Outcome measures in ALSPAC including internalizing and fat mass variables, were recorded at 
# 10y, 13y, 15y, 17y and 22y. 

#### ---------------------------- Dependencies ---------------------------- ####

PATH_RESULTS <- '' # ATTENTION! DEFINE HERE THE PATH WHERE YOU WANT ALL RESULTS TO BE STORED.

# First, let's point to the necessary libraries
library(foreign)

# Check if the path to the data is already in memory, otherwise ask for it. 
if (exists("alspac_file") == F) { 
  alspac_file <- file.choose() 
  alspac_folder <- dirname(alspac_file) 
  # Read in the data
  alspac.table <- foreign::read.spss(alspac_file, use.value.label=TRUE, to.data.frame=TRUE) 
}

# Read in the data on parental depression 
dep <- readRDS(file.path(alspac_folder, "raw_parent_depr_anxiety.rds"))

# Change all names in alspac.table to lowercase
names(alspac.table)=tolower(names(alspac.table))

# Initiate a cov_out dataframe with id of the child as first column
cov_out <- data.frame("IDC" = paste(alspac.table$cidb2957, alspac.table$qlet, sep = "_"))

################################################################################
#### ---- INTERNALIZING PROBLEMS and FATMASS ( @ 10, 13, 15, 15 & 22 ) -_-- ####
################################################################################

select_age_outcome <- function(age, df) {
  
  if (age == 10) { vars = c('ku991a',    # DV: Age of study child at completion (months)
                            'ku707a',    # DV: SDQ emotional symptoms score (complete cases) 
                            'f9003c',    # Age (months) at F9 visit	
                            'f9dx135',   # Total Body - fat mass (g): F9
                            'f9dx010') } # DV: DXA weight (Kg): F9
  else if (age == 11) { vars = c('kw9991a',  # DV: Age of study child at completion (months)
                                 'kw6602a',  # DV: SDQ - Emotional symptoms score (complete cases)
                                 'fe003c',   # Age (months) at F11+ visit	
                                 'fedx135',  # Total Body - fat mass (g): DXA: F11
                                 'fedx016') } # DV: DXA weight (Kg): F11
  else if (age == 13) { vars = c('ta9991a',  # DV: Age of study child at completion (months)
                                 'ta7025a',  # DV: SDQ emotional symptoms score (prorated)
                                 'fg0011a',  # DV: Age of study child at attendance (months): TF2
                                 'fg3254',   # Total Body - fat mass (g): DXA: TF2
                                 'fg3207',   # DV: DXA weight (Kg): DXA: TF2
                                 'fg3257') } # Android - fat mass (g): DXA: TF2
  else if (age == 15) { vars = c('fh0011a',  # DV: Age of study child at attendance (months): TF3
                                 'fh6876',   # DV: Depression (self-report 6-band computer prediction, ICD-10 and DSM-IV): TF3
                                 'fh0011a',  # DV: Age of study child at attendance (months): TF3
                                 'fh2254',   # Total - fat mass (g): DXA: TF3
                                 'fh2208',   # Keyed Weight (kg): DXA: TF3	
                                 'fh2257') } # Android - fat mass (g): DXA: TF3
  else if (age == 17) { vars = c('tc9991a',  # DV: Age of study teenager at completion (months)
                                 'tc4025a',  # DV: SDQ emotional symptoms score (prorated)
                                 'fj003b',   # Age in years at clinic visit [F17]
                                 'fjdx135',  # Total: fat mass (g) [F17]	
                                 'fjmr022',  # M15: Weight (kgs) [F17]	
                                 'fjdx138') } # Android: fat mass (g) [F17]
  else if (age == 22) { vars = c('ypb9992',   # DV: Respondent age at completion (months)
                                 'ypb5180',   # DV: MFQ (moods and feelings questionnaire) score
                                 'fkar0010',  # Age at clinic visit (in months): F@24
                                 'fkdx1001',  # Total body fat mass (g): F@24
                                 'fkms1030',  # Weight (kg): F@24	
                                 'fkdx1041')} # Android fat mass (g): F@24
  
  div = 12; if (vars[3] == 'fj003b') { div = 1 }
  
  out.df <- data.frame( 
    paste0('age_child_sdq_', age) = as.numeric(levels(df[, vars[1]]))[df[, vars[1]]]/ 12,
     paste0('intern_score_', age) = as.numeric(levels(df[, vars[2]]))[df[, vars[2]]],
    paste0('age_child_dxa_', age) = as.numeric(levels(df[, vars[3]]))[df[, vars[3]]]/ div,
        paste0('total_fat_', age) = as.numeric(levels(df[, vars[4]]))[df[, vars[4]]],
       paste0('weight_dxa_', age) = as.numeric(levels(df[, vars[5]]))[df[, vars[5]]] ) 
  
  if (length(vars) > 5) {
    out.df[paste0('andr_fat_mass_', age)] <- as.numeric(levels(df[, vars[6]]))[df[, vars[6]]] 
  }
  return(out.df)       
}

out_10 <- select_age_outcome(10, alspac.table)
out_13 <- select_age_outcome(13, alspac.table)

# Body fat percentage calculation
out_10$tot_fat_percent_10 <- (out_10$total_fat_10 / (out_10$weight_dxa_10*1000)) * 100
out_13$tot_fat_percent_13 <- (out_13$total_fat_13 / (out_13$weight_dxa_13*1000)) * 100

cov_out = cbind(cov_out, out_10, out_13)

# ------------------------------------------------------------------------------
# Inspect the correlations between outcome variables
cor_outcome <- round(cor(cov_out[, -1], use = 'pairwise.complete.obs'), 2) 
write.csv(cor_outcome, file = file.path(PATH_RESULTS, "corr_mat_outcomes.csv"), row.names = T, quote = F)

################################################################################
#### ---------------------------- COVARIATES ------------------------------ ####
################################################################################

# Variables that will be used in the covariate models of this project are those 
# marked with ###. they include: 'sex', 'age_child', 'm_smoking', 'm_drinking' 
# and 'm_bmi_before_pregnancy'.

# For the other demographic auxiliary variables (used for imputation): when they 
# were assessed both prenatally and postnatally, both measures are included.
# Auxiliary variables for this project include: 'ethnicity', 'm_age_cont', parity', 
# 'gest_age_birth', 'gest_weight', 'm_bmi_pregnancy', 'm_bmi_7y',
# 'm_dep_pregnancy', 'p_dep_pregnacy', "m_dep_3y", "p_dep_3y"

# ------------------------------------------------------------------------------
### SEX of the child
cov_out$sex <- alspac.table$kz021 # 1 = Male; 2 = Female.

# ------------------------------------------------------------------------------
### AGE of the child
# Combine age of the child measured during internalising and fatmass measurement
# This value will serve as a covariate in the first adjusted model.

cov_out$age_child <- (cov_out$age_child_sdq_13 + cov_out$age_child_dxa_13) / 2

# OPTIONAL: check age difference between measurements
plot(cov_out$age_child_sdq_13, cov_out$age_child_dxa_13)
summary(cov_out$age_child_sdq_13 - cov_out$age_child_dxa_13)

#-------------------------------------------------------------------------------
### MATERNAL SMOKING during pregnancy 
# 3 categories (never, former and current smoker)

# Ever smoked: no = 0, yes = 1
cov_out$b650r <- ifelse(alspac.table$b650 == 'N', 0, ifelse(alspac.table$b650 == 'Y', 1, NA)) 
# CIGS smoked per day during pregnancy: none = 0, occasionally or >1 = 1
cov_out$c482r <- ifelse(alspac.table$c482 == 'None', 0, ifelse(alspac.table$c482 != 'DK', 1, NA)) 

cov_out$m_smoking <- ifelse(cov_out$b650r == 0 & cov_out$c482r == 0, 0,        # Never a smoker
                            ifelse(cov_out$b650r == 1 & cov_out$c482r == 0, 1, # Former smoker
                                   ifelse(cov_out$c482r == 1, 2, NA)))         # Current smoker

#-------------------------------------------------------------------------------
### MATERNAL ALCOHOL CONSUMPTION during pregnancy
# Combined maternal drinking during the first and the last trimester of pregnancy.

# Alcohol consumption in 1-3MTHS this PREG
cov_out$b721r <- ifelse(alspac.table$b721 == 'never', 0,   
                        ifelse(alspac.table$b721 == '<1 glass PWK', 1,
                               ifelse(alspac.table$b721 == '1+ glasses PWK', 2, 
                                      ifelse(alspac.table$b721 == '1-2 glasses PDAY', 3, 
                                             ifelse(alspac.table$b721 == '3-9 glasses PDAY', 4, 
                                                    ifelse(alspac.table$b721 == '10+ glasses PDAY', 5, NA))))))
# FREQ of alcohol use in last 2MTHS of PREG
cov_out$e220r <- ifelse(alspac.table$e220 == 'Not at all', 0, 
                        ifelse(alspac.table$e220 == '<1PWK', 1,
                               ifelse(alspac.table$e220 == 'At least 1PWK', 2, 
                                      ifelse(alspac.table$e220 == '1-2 glasses daily', 3, 
                                            ifelse(alspac.table$e220 == '3-9 glasses daily', 4, 
                                                   ifelse(alspac.table$e220 == '>9 glasses daily', 5, NA))))))

cov_out$m_drinking <- rowMeans(cov_out[, c('b721r', 'e220r')], na.rm = F)


#-------------------------------------------------------------------------------
## Other variables

# Ethnicity
cov_out$ethnicity <- ifelse(is.na(alspac.table$c800) & is.na(alspac.table$c801), NA, 
                            ifelse((alspac.table$c800 == 'White' & alspac.table$c801 == 'White') | 
                                     ((is.na(alspac.table$c800) | is.na(alspac.table$c801)) & 
                                     (alspac.table$c800 == 'White' | alspac.table$c801 == 'White')), 1, 0))

# Maternal weight (Kg)
cov_out$weight_pre <- as.numeric(levels(alspac.table$dw002))[alspac.table$dw002] # Pre-pregnancy weight (Kg)
cov_out$weight_7y  <- as.numeric(levels(alspac.table$m4220))[alspac.table$m4220] # kg 7y

# Maternal height (cm) 7y 1m (we treat height as a constant)
cov_out$height_7y <- as.numeric(levels(alspac.table$m4221))[alspac.table$m4221] / 100 # transform into meters

# BMI
cov_out$m_bmi_before_pregnancy <- cov_out$weight_pre / ((cov_out$height_7y)^2) # calculating pregnancy BMI 
cov_out$m_bmi_7yrs             <- cov_out$weight_7y  / ((cov_out$height_7y)^2) # calculating BMI at age 7

# Siblings identifier
cov_out$sibling <- ifelse(alspac.table$mult == 1, 1, ifelse(alspac.table$mult == 0, 0, NA))  # multiple pregnancies in ALSPAC (exclusion criteria)

# Twin identifier
cov_out$twin <- rep(0, nrow(cov_out))
cov_out$twin[which(duplicated(alspac.table$cidb2957) | duplicated(alspac.table$cidb2957, fromLast = T))] <- 1
                                  
# Auxiliary variables (for imputation)
cov_out$parity <- as.numeric(levels(alspac.table$b032))[alspac.table$b032]  # parity (used for imputation)

cov_out$gest_age_birth <- as.numeric(levels(alspac.table$bestgest))[alspac.table$bestgest] # gestational age at birth (used for imputation)
cov_out$gest_weight    <- as.numeric(levels(alspac.table$kz030))[alspac.table$kz030]       # gestational weight (used for imputation)
cov_out$m_age_cont     <- as.numeric(levels(alspac.table$mz028b))[alspac.table$mz028b]     # maternal age at intake (used for imputation) 

#-------------------------------------------------------------------------------
# Maternal and paternal depression during pregnancy and childhood are calculated 
# in the CCEI EPDS script. 

# Prenatal maternal depression
cov_out$m_dep_cont_pregnancy <- (dep$m_EPDS_total_18wg + dep$m_EPDS_total_32wg) / 2
# Postnatal maternal depression
post_dep_m <- sapply(dep[, c('f200', 'g290', 'h200a')], as.integer) #  NAs introduced by coercion from not depressed and very depressed
cov_out$m_dep_cont_childhood <- rowMeans(post_dep_m, na.rm = T)

# Prenatal paternal depression
cov_out$p_dep_cont_pregnancy <- dep$p_EPDS_total_18wg
# Postnatal paternal depression
post_dep_p <- sapply(dep[, c('pe290', 'pd200')], as.integer) 
cov_out$p_dep_cont_childhood <- rowMeans(post_dep_p, na.rm = T)

################################################################################
#### --------------------------- save and run ----------------------------- ####
################################################################################

# Save the dataset in the directory where you have the raw data
saveRDS(cov_out, file.path(PATH_RESULTS, "PCMout_cov_aux.rds"))

