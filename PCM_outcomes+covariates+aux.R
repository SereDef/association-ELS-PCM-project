
# Hi there, 
# the following code is building a dataset will all necessary variables for the analysis
# of the association between ELS and psycho-cardio-metabolic multi-morbidity in children 
# This includes the outcomes of interest (internalizing problems and cardio-metabolic
# risk), the covariates that are going to be used as well as the auxiliary variables 
# used in the imputation of the final dataset. 
# It does not include data used to build the ELS score exposure, for which you can 
# check out https://github.com/SereDef/cumulative-ELS-score

## For this script you are going to need the following datasets from datamanagement
# CHILDCBCL9_10082016.sav, CHILDFATMASS9_13092016.sav, CHILDBLOODPRESSURE9_21042016.sav,
# CHILDSERUM9_01082017.sav, CHILDGROWTH9_04062020.sav, CHILDOBESITY9_04062020.sav,
# CHILD-ALLGENERALDATA_07072020.sav, GEDRAGSGROEP_MaternalDrinking_22112016.sav, 
# MATERNALSMOKING_22112016.sav, GEDRAGSGROEP_MaternalDrinking_22112016.sav, 
# MOTHERANTHROPOMETRY_18022013.sav, GR1003-BSI D1_22112016.sav, and 
# BSI 3 years of age_GR1065 G1-GR1066 C1_22112016.sav, GR1004-BSI G1_22112016.sav

#### ---------------------------- Dependencies ---------------------------- ####

# First, let's point to the necessary libraries
library(foreign)
library(stats)

# Define the path to the data
pathtodata <- readline(prompt="Enter path to data: ")
# ATTENTION!!! You will be prompted with an "Enter path to data:" message 
# -> Enter the location of your datafiles. The code assumes that all (raw) data is 
# stored in ONE folder. Do not forget the final slash in your path, and, speaking of slashes, 
# beware of OS sensitive changes when you want to modify the structure of your dirs!

#### ------------------------------ FUNCTIONS ----------------------------- ####

# Read in the data and fix all SPSS weird missing codes into NAs
readquick <- function(filename, rootdir = pathtodata, exclude_col = "") { # only works for SPSS files
  dat <- read.spss(paste(rootdir, filename, sep=""), 
                         use.value.labels = F, to.data.frame = T)
  # Get rid of all capital letters in column names (so you don't have to worry)
  names(dat) <- tolower(names(dat))
  # Replace values of 777, 888 or 999 with NAs unless they are IDCs or IDMs 
  # If you do not want this to happen for any other column use the exclude_col argument. 
  for (i in 1:length(dat)) {
    if (colnames(dat)[i] == "idm" | colnames(dat)[i] == "idc" | colnames(dat)[i] == exclude_col) {
      dat[,i] <- dat[,i]
    } else {
      dat[,i] <- ifelse(dat[,i] == 777 | dat[,i] == 888 | dat[,i] == 999, NA, dat[,i]) }
  } 
  return(dat)
}
# ------------------------------------------------------------------------------
# Some cardio-metabolic variables need their SDSs to be adjusted based on sex 
# { and age ... but not done here as the relations were not linear } 
sex_adjust_scale <- function(x) {
  zscore = ifelse(PCM_outcome$sex == 1, 
                  yes = scale(x, center = mean(x[PCM_outcome$sex == 1], na.rm =T), 
                              scale = sd(x[PCM_outcome$sex == 1], na.rm =T)),
                  no = scale(x, center = mean(x[PCM_outcome$sex == 2], na.rm =T), 
                             scale = sd(x[PCM_outcome$sex == 2], na.rm =T)))
  return(zscore)
}

# -----------------------------------------------------------------------------#
##                        Ok, we are good to go now!                          ##
# -----------------------------------------------------------------------------#

################################################################################
#### ------------------ INTERNALIZING PROBLEMS ( @ 9 ) -------------------- ####
################################################################################

# Read in the dataset
cbcl <- readquick("CHILDCBCL9_10082016.sav") # 9901 obs. of 627 vars

# The internalizing sub-scale of the Child Behavior Checklist (CBCL 6-18) (Achenbach, 1999)
# is an empirically based score derived from the widely used parent-report questionnaire. 
# It contains 32 items rated on a 3-point scale (0 = "Not true", 1 = "Somewhat or 
# sometimes true", 2 = "Very or often true") and referred to the past 6 months. 
# Hence the score ranges from 0 to 64. 

# Internalizing scale @ 9 yrs # informant: MOTHER.
internalizing = cbcl[,c('idc',
                'agechild_cbcl9m',
                'cbcl5_9m',  # there is very little he/she enjoys
                'cbcl14_9m', # cries a lot
                'cbcl29_9m', # fears certain animals, situations or places other than school
                'cbcl30_9m', # fears going to school
                'cbcl31_9m', # fears he/she might think or do something bad
                'cbcl32_9m', # feels he/she has to be perfect
                'cbcl33_9m', # feels or complains that no one loves him/her
                'cbcl35_9m', # feels worthless or inferior
                'cbcl42_9m', # would rather be lone than with others
                'cbcl45_9m', # nervous, highstrung or tense
                'cbcl47_9m', # nightmares
                'cbcl49_9m', # Constipated, doesn't move bowels
                'cbcl50_9m', # too fearful or anxious
                'cbcl51_9m', # feels dizzy or lightheaded
                'cbcl52_9m', # feels too guilty
                'cbcl54_9m', # overtired without good reason
                'cbcl56a_9m', # physical problems without medical cause ...aches or pains(not stomach or headaches)
                'cbcl56b_9m', # ...headaches
                'cbcl56c_9m', # ...nausea, feels sick
                'cbcl56d_9m', # ...problems with eyes (not if corrected with glasses)
                'cbcl56e_9m', # ...rashes or other skin problems
                'cbcl56f_9m', # ...stomachaches
                'cbcl56g_9m', # ...vomiting, throwing up
                'cbcl65_9m', # refuses to talk
                'cbcl69_9m', # secretive, keeps things to self
                'cbcl71_9m', # self-conscious or easily embarrassed
                'cbcl75_9m', # too shy or timid
                'cbcl91_9m', # talks about killing self
                'cbcl102_9m', # underactive, slow moving, or lacks energy
                'cbcl103_9m', # unhappy, sad, or depressed
                'cbcl111_9m', # withdrawn, doesn't get involved with others
                'cbcl112_9m', # worries
                'nmisint_9m', # number of missing values in internalizing scale items
                'sum_int_9m')] # weighted sum score internalizing scale (allowing 25% missing)
                 # Calculation (based on SPSS script available on the V: drive): 
                 # if 24 out of 32 item sare available (i.e. 75%), sum of the item
                 # scores * (32 / 32 - nmisnt_9m). 
# Let's make it a bit more reader friendly 
colnames(internalizing)[35:36] <- c("n_missing_intern", "intern_score"); 

# One alternative could be to use the anxious/depressed empirical sub-scale (9y, mother report)
# However I did give it a shot and it does not seem to perfom better than the internalizing one.
# anx_dep = cbcl[,c('idc', 'agechild_cbcl9m', 'cbcl14_9m', # cries a lot
#  'cbcl29_9m', # fears certain animals, situations or places other than school
#  'cbcl30_9m', # fears going to school
#  'cbcl31_9m', # fears he/she might think or do something bad
#  'cbcl32_9m', # feels he/she has to be perfect
#  'cbcl33_9m', # feels or complains that no one loves him/her
#  'cbcl35_9m', # feels worthless or inferior
#  'cbcl45_9m', # nervous, highstrung or tense
#  'cbcl50_9m', # too fearful or anxious
#  'cbcl52_9m', # feels too guilty
#  'cbcl71_9m', # self-conscious or easily embarrassed
#  'cbcl91_9m', # talks about killing self
#  'cbcl112_9m', # worries
#  'nmisanx_9m',  # number of missing values in anxious/depressed scale items
#  'sum_anx_9m')] # weighted sum score anxious/depressed scale (allowing 25% missing)

# ------------------------------------------------------------------------------

################################################################################
#### ---------------- CARDIO-METABOLIC RISK SCORE ( @ 9 ) ----------------- ####
################################################################################

# Read in the datasets
fatmass <- readquick("CHILDFATMASS9_13092016.sav") # 5862 obs of 28 vars
bp <- readquick("CHILDBLOODPRESSURE9_21042016.sav") # 5862 obs of 28 vars
serum <- readquick("CHILDSERUM9_01082017.sav") # 5862 obs of 13 vars

# Select only the necessary measures
fat_mass = fatmass[,c('idc', 'agechild9_visit1', # this value (age) is the same for all other datasets
                      'fat_mass_androidchild9')]
blood_pressure = bp[,c('idc', 
                      'meansbp_ex1_child9',
                      'meandbp_ex1_child9')]
trig_hdlc_ins = serum[,c('idc', 
                       'triglyceriden9_clean', # or use the triglyceridenchild9_log?
                       'hdlcholesterolchild9_clean',
                       'insulinechild9_clean')] # or use the _log?

# Merge them into one dataset
cardio_metab <- Reduce(function(x,y) merge(x = x, y = y, by = 'idc',  all.x = TRUE),
                          list(fat_mass, blood_pressure, trig_hdlc_ins)) 
# Let's make it a bit more reader friendly 
colnames(cardio_metab)[3:8] <- c("fat_mass","sbp","dbp","triglycerides","hdl_cholesterol","insuline")


################################################################################
# merge the two main (mental and physical) outcomes into one dataset
PCM_outcome <- merge(internalizing, cardio_metab, by = 'idc', all.x = T)
################################################################################

# ------------------------------------------------------------------------------
# Add sex variable that I will need for correcting the cardio-metab z scores 
child_general <- readquick("CHILD-ALLGENERALDATA_07072020.sav") # 9901 obs of 122 
sex = child_general[,c('idc', 'gender')]; colnames(sex)[2] <- "sex"; # 1 = boy; 2 = girl.
PCM_outcome <- merge(sex, PCM_outcome, by = 'idc', all.x = T) 

# ------------------------------------------------------------------------------
# Take the standard deviation score of the variables of interest
PCM_outcome$intern_score_z = scale(PCM_outcome$intern_score)
  # these used the previously defined function that corrects SDSs for sex
PCM_outcome$fat_mass_z = sex_adjust_scale(PCM_outcome$fat_mass)
PCM_outcome$sbp_z = sex_adjust_scale(PCM_outcome$sbp)
PCM_outcome$dbp_z = sex_adjust_scale(PCM_outcome$dbp)
PCM_outcome$triglycerides_z = sex_adjust_scale(PCM_outcome$triglycerides)
PCM_outcome$hdl_cholesterol_z = sex_adjust_scale(PCM_outcome$hdl_cholesterol)
PCM_outcome$insuline_z = sex_adjust_scale(PCM_outcome$insuline)

# ------------------------------------------------------------------------------
# A continuous cumulative cardio-metabolic risk (CMR) index (Voortman et al., 2016) 
# is calculated as the sum of {age- and sex-specific} SD-scores (SDS) of: android fat 
# mass (BF); systolic or diastolic blood pressure (SBP/DBP); high-density lipoprotein-
# cholesterol (HDL-C); triglycerides (TG) and insulin levels. 
# Scores for HDL-C are multiplied by − 1 since higher HDL-C concentrations reflect 
# lower cardio-metabolic risk and scores for SBP and DBP are multiplied by 0.5 so that 
# each contributes half to the blood pressure component. 

# Thus, the cardio-metabolic risk factor score was calculated as: 
# BF% SDS + 0.5 × SBP SDS + 0.5 × DBP SDS + TG SDS + (−1 × HDL−C SDS) + insulin SDS. 

# PROVISIONAL (i.e. adjusted for sex only)
PCM_outcome$cmr_score = PCM_outcome$fat_mass_z +
  (0.5 * PCM_outcome$sbp_z) + (0.5 * PCM_outcome$dbp_z) +
  PCM_outcome$triglycerides_z +
  (-1 * PCM_outcome$hdl_cholesterol_z) +
  PCM_outcome$insuline_z

# A higher score reflects higher cardio-metabolic risk. When dichotomized, the score 
# can be seen as a proxy of child metabolic syndrome (van Gijssel et al., 2016). 
# For more specifics about the measurements see (Boyd et al., 2013; Voortman et al., 2016).
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Since the score did not seem to work very well, I also include BMI (@ 9y) and 
# it's categorization into underweight (1,2,3), normal weight, overweight and obese
# into the dataset
growth_9yrs <- readquick('CHILDGROWTH9_04062020.sav')
obesity_9yrs <- readquick('CHILDOBESITY9_04062020.sav') # Weight status (Cole)
# (-3) Thin, grade 3; (-2) Thin, grade 2; (-1) Thin, grade 1; (0) Normal weight; 
# (1) Overweight; (2) Obese;
bmi_9yrs = growth_9yrs[,c('idc', 'sdsbmiforage9childt')] # Body Mass Index (age 20w - 20y) The Netherlands, 2010 (Talma???)

# Not used # obesity_5yrs <- readquick('CHILDOBESITY5_08102014.sav') 
# obese_forever = merge(obesity_5yrs, obesity_9yrs, by = 'idc', all.x = T)

PCM_outcome <- Reduce(function(x,y) merge(x = x, y = y, by = 'idc',  all.x = TRUE),
                       list(PCM_outcome, bmi_9yrs, obesity_9yrs)) 
colnames(PCM_outcome)[53:54] <- c("bmi_z_4age","bmi_binned")

# ------------------------------------------------------------------------------
#                                     OPTIONAL                                 #
# Check if relation between these variables and child age is approximately linear
plot(PCM_outcome$agechild9_visit1, PCM_outcome$fat_mass,
     main = paste("FAT MASS - AGE (@9y): corr =", cor(PCM_outcome$fat_mass, PCM_outcome$agechild9_visit1, method = "pearson", use = "complete.obs")),
     xlab = "Age in years", ylab = "Android Fat Mass", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)
plot(PCM_outcome$agechild9_visit1, PCM_outcome$sbp,
     main = paste("SBP - AGE (@9y): corr =", cor(PCM_outcome$sbp, PCM_outcome$agechild9_visit1, method = "pearson", use = "complete.obs")),
     xlab = "Age in years", ylab = "Systolic Blood Pressure", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)
plot(PCM_outcome$agechild9_visit1, PCM_outcome$dbp,
     main = paste("DBP - AGE (@9y): corr =", cor(PCM_outcome$dbp, PCM_outcome$agechild9_visit1, method = "pearson", use = "complete.obs")),
     xlab = "Age in years", ylab = "Diastolic Blood Pressure", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)
plot(PCM_outcome$agechild9_visit1, PCM_outcome$triglycerides,
     main = paste("TRIGLYCERIDES - AGE (@9y): corr =", cor(PCM_outcome$triglycerides, PCM_outcome$agechild9_visit1, method = "pearson", use = "complete.obs")),
     xlab = "Age in years", ylab = "Triglycerides", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)
plot(PCM_outcome$agechild9_visit1, PCM_outcome$hdl_cholesterol,
     main = paste("HDL-CHOLESTEROL - AGE (@9y): corr =", cor(PCM_outcome$hdl_cholesterol, PCM_outcome$agechild9_visit1, method = "pearson", use = "complete.obs")),
     xlab = "Age in years", ylab = "HDL - Cholesterol", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)
plot(PCM_outcome$agechild9_visit1, PCM_outcome$insuline,
     main = paste("INSULINE - AGE (@9y): corr =", cor(PCM_outcome$insuline, PCM_outcome$agechild9_visit1, method = "pearson", use = "complete.obs")),
     xlab = "Age in years", ylab = "Insuline", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)

# ------------------------------------------------------------------------------
# Let's peak into the correlations between internalizing and cardio-metabolic vars
# Select only the variables of interest
relevant_vars = PCM_outcome[,45:54]

# Build a correlation matrix
test_cor = cor(relevant_vars, method = "pearson", use = "complete.obs")

# Represent the matrix I just built with an heatmap
# but first, default colors in R heatmap are horrifying so let's fix them
better_colors <- colorRampPalette(c("blue", "white", "red"))(n = 299)
heatmap(test_cor, col = better_colors, Rowv = NA, symm = T, scale = 'none')

# Scatterplots to explore the relations between Internalizing and other physical vars
                                   # CMR score
plot(relevant_vars$intern_score_z, relevant_vars$cmr_score, 
     main = paste("PCM Outcomes: corr =", round(test_cor[8,1],4)),
     xlab="Internalizing score (z score) ", ylab="Cardio-metabolic risk (z score)", 
     col=rgb(0,100,0,50,maxColorValue=255), pch=16)
                      # Individual components of the CMR score
plot(relevant_vars$intern_score_z, relevant_vars$fat_mass_z, 
     main = paste("INT - FAT MASS: corr =", round(test_cor[2,1],4)),
     xlab="Internalizing score (z score) ", ylab="Android Fat Mass (z score)", 
     col=rgb(100,50,50,100,maxColorValue=255), pch=16)
plot(relevant_vars$intern_score_z, relevant_vars$sbp_z, 
     main = paste("INT - SBP: corr =", round(test_cor[3,1],4)),
     xlab="Internalizing score (z score) ", ylab="Systolic Blood Pressure (z score)", 
     col=rgb(100,50,50,100,maxColorValue=255), pch=16)
plot(relevant_vars$intern_score_z, relevant_vars$dbp_z, 
     main = paste("INT - DBP: corr =", round(test_cor[4,1],4)),
     xlab="Internalizing score (z score) ", ylab="Diastolic Blood Pressure (z score)", 
     col=rgb(100,50,50,100,maxColorValue=255), pch=16)
plot(relevant_vars$intern_score_z, relevant_vars$triglycerides_z, 
     main = paste("INT - TRIGLYCERIDES: corr =", round(test_cor[5,1],4)),
     xlab="Internalizing score (z score) ", ylab="Triglycerides (z score)", 
     col=rgb(100,50,50,100,maxColorValue=255), pch=16)
plot(relevant_vars$intern_score_z, relevant_vars$hdl_cholesterol_z, 
     main = paste("INT - HDL-CHOLESTEROL: corr =", round(test_cor[6,1],4)),
     xlab="Internalizing score (z score) ", ylab="HDL - cholesterol (z score)", 
     col=rgb(100,50,50,100,maxColorValue=255), pch=16)
plot(relevant_vars$intern_score_z, relevant_vars$insuline_z, 
     main = paste("INT - INSULINE: corr =", round(test_cor[7,1],4)),
     xlab="Internalizing score (z score) ", ylab="Insuline (z score)", 
     col=rgb(100,50,50,100,maxColorValue=255), pch=16)
                                    # BMI
plot(relevant_vars$intern_score_z, relevant_vars$bmi_z_4age, 
     main = paste("INT - BMI (@9y): corr =", round(test_cor[9,1],4)),
     xlab="Internalizing score (z score) ", ylab="Body Mass Index - 9 yrs (z score)", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)
                                 # Obesity
plot(relevant_vars$intern_score_z, relevant_vars$bmi_binned, 
     main = paste("INT - Obesity (@9y): corr =", round(test_cor[10,1],4)),
     xlab="Internalizing score (z score) ", ylab="Obesity - 9 yrs", 
     col=rgb(1,100,150,150,maxColorValue=255), pch=16)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

################################################################################
#### ---------------------------- COVARIATES ------------------------------ ####
################################################################################

# Variables that will be used in the covariate models of this project are those 
# marked with ###. they include: 'sex', 'age_child', 'm_smoking', 'm_drinking' 
# and 'm_bmi_berore_pregnancy'.

# For the other demographic auxiliary variables (used for imputation): when they 
# were assessed both prenatally and postnatally, both measures are included.
# Auxiliary variables for this project include: 'ethnicity', 'm_age' (in the postnatal_ELS score), 
# 'parity', 'm_smoking', 'gest_age_birth', 'gest_weight', 'm_bmi_pregnancy', 'm_bmi_5yrs',
# 'sex', 'm_dep_pregnancy', 'p_dep_pregnacy', "m_dep_3yrs", "p_dep_3yrs"

# ------------------------------------------------------------------------------

### AGE of the child
# Combine age of the child measured during first visit and at CBCL administration
# This value will serve as a covariate in the first adjusted model.
PCM_outcome$age_child = (PCM_outcome$agechild9_visit1 + PCM_outcome$agechild_cbcl9m) / 2

#-------------------------------------------------------------------------------
### MATERNAL SMOKING during pregnancy 
smokingv1 <- readquick("MATERNALSMOKING_22112016.sav") #  9778 obs of 11 variables

smoking = smokingv1[,c('idm', 'smoke_all')] # (1) never a smoker; 
                                            # (2) smoked until pregnancy was known (i.e., first trimester only); 
                                            # (3) continued smoking during pregnancy.
colnames(smoking)[2] = "m_smoking"

#-------------------------------------------------------------------------------
### MATERNAL ALCOHOL CONSUMPTION during pregnancy
drinkingv1 <- readquick("GEDRAGSGROEP_MaternalDrinking_22112016.sav") #drinkingv2 <- readquick("MATERNALALCOHOL_22112016.sav") # old variable

drinking = drinkingv1[,c('idm', 'mdrink_updated')] # (0) never; 
                                                   # (1) until pregnancy was known (i.e., first trimester only); 
                                                   # (2) continued during pregnancy occasionally;
                                                   # (3) continued during pregnancy frequently.
colnames(drinking)[2] = "m_drinking"

#-------------------------------------------------------------------------------
## Others

# Ethnicity recode – dichotomized into: dutch, western and non-western;
for (i in 1:9901) {
  if (is.na(child_general$ethnfv2[i])) { child_general$ethnicity[i] <- NA
  } else if (child_general$ethnfv2[i] == 1) { # Dutch
    child_general$ethnicity[i] <- 0 
  } else if (child_general$ethnfv2[i] == 300 | child_general$ethnfv2[i] == 500 | child_general$ethnfv2[i] >= 700) { 
    # American, western (300) Asian, western (500) European (700), Oceanie (800)
    child_general$ethnicity[i] <- 1 
  } else { 
    child_general$ethnicity[i] <- 2 } 
  # Indonesian (2), Cape Verdian (3), Maroccan (4) Dutch Antilles (5) Surinamese 
  # (6) Turkish (7) African (200), American, non western (400), Asian, non western (600)
}

general_cov_auxiliary = child_general[,c('idc', 'idm', ### SEX was already included in the set. 
                               'ethnicity', ### ETHNICITY – dutch, western, non-western
                               'bmi_0',     ###	MATERNAL BMI – self-reported, before pregnancy
                               'twin',      # used for exclusion criteria 
                               'mother',    # mother id used to identify siblings (for exclusion)
                               'parity',    # parity (used for imputation)
                               'gestbir',   # gestational age at birth (used for imputation)
                               'weight',    # gestational weight (used for imputation)
                               'bmi_1')]    # maternal BMI during pregnancy (used for imputation)
# Again, let's try to keep it user friendly 
colnames(general_cov_auxiliary)[c(4,8:10)] = c("m_bmi_berore_pregnancy", "gest_age_birth", 
                                               "gest_weight", "m_bmi_pregnancy")
#-------------------------------------------------------------------------------
# Maternal BMI at age 5 (used for imputation)
m_anthropometry_5yrs = readquick('MOTHERANTHROPOMETRY_18022013.sav')
m_bmi_5yrs = m_anthropometry_5yrs[,c('mother','bmimotherf5')]; colnames(m_bmi_5yrs)[2] = "m_bmi_5yrs"

# Merge with the other general variables
general_cov_auxiliary = merge(general_cov_auxiliary, m_bmi_5yrs, by = 'mother',  all.x = TRUE)

#-------------------------------------------------------------------------------
# Maternal and paternal depression during pregnancy and at age 3 
bsi_pregnancy_m = readquick('GR1003-BSI D1_22112016.sav') # 9778 obs of 261 vars
bsi_pregnancy_p = readquick('GR1004-BSI G1_22112016.sav') # 9778 obs of 261 vars
bsi_3yrs = readquick('BSI 3 years of age_GR1065 G1-GR1066 C1_22112016.sav') # 9897 obs of 49 vars

# Depression during pregnancy
dep_pregnancy_m = bsi_pregnancy_m[, c('idm','dep')]; colnames(dep_pregnancy_m)[2] = c("m_dep_pregnancy")
dep_pregnancy_p = bsi_pregnancy_p[, c('idm','dep_p')]; colnames(dep_pregnancy_p)[2] = c("p_dep_pregnancy")

# Merge it with the previous dataset
general_cov_auxiliary <- Reduce(function(x,y) merge(x = x, y = y, by = 'idm',  all.x = TRUE),
                                list(general_cov_auxiliary, dep_pregnancy_m, dep_pregnancy_p))

# Depression @ 3y: Items 9, 16, 17, 18, 35, and 50
d <- data.frame(bsi_3yrs$g0100365, bsi_3yrs$g0100665, bsi_3yrs$g0100765, bsi_3yrs$g0100865, bsi_3yrs$g0101365, bsi_3yrs$g0102165, # mother report
                bsi_3yrs$c0100366, bsi_3yrs$c0100666, bsi_3yrs$c0100766, bsi_3yrs$c0100866, bsi_3yrs$c0101366, bsi_3yrs$c0102166) # father report
n_items_m <- rowSums(!is.na(d[,1:6])); n_items_p <- rowSums(!is.na(d[,7:12]))
bsi_3yrs$m_dep_3yrs <- ifelse(n_items_m >= 5, yes = (rowSums(d[,1:6])/n_items_m)-1, no = NA)
bsi_3yrs$p_dep_3yrs <- ifelse(n_items_p >= 5, yes = (rowSums(d[,7:12])/n_items_p)-1, no = NA)

dep_3yrs = bsi_3yrs[, c('idc','m_dep_3yrs', 'p_dep_3yrs')]
# Merge it with the previous dataset
general_cov_auxiliary = merge(general_cov_auxiliary, dep_3yrs, by = 'idc',  all.x = TRUE)

#-------------------------------------------------------------------------------
################################################################################
# Merge all covariates / auxiliary variables together
covariates_and_auxiliary <- Reduce(function(x,y) merge(x = x, y = y, by = 'idm',  all.x = TRUE),
                       list(general_cov_auxiliary, smoking, drinking)) 
################################################################################

################################################################################
# merging outcome variables, covariates and auxiliary variables in one dataset
PCM_project = merge(PCM_outcome, covariates_and_auxiliary, by = 'idc', all.x = T)

# A bit of a quick and dirty fix to make merging with ELS easier
colnames(PCM_project)[which(colnames(PCM_project) == 'idc')] <- toupper('idc')
colnames(PCM_project)[which(colnames(PCM_project) == 'idm')] <- toupper('idm')

################################################################################
#-------------------------------------------------------------------------------

################################################################################
#### --------------------------- save and run ----------------------------- ####
################################################################################

# Save the dataset in an .rds file, in the directory where the raw data are stored
saveRDS(PCM_project, paste(pathtodata,'PCM_allvars.rds'))
