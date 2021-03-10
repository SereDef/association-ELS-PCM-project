################################################################################
####################### NOT USED IN THE FINAL VERSION ##########################
################################################################################

# Hi, the following code is building a cardio-metabolic risk score at age 10. 
# The original score was created by Voortman et al., (2016) at age 6. 

## For this script you are going to need the following datasets from datamanagement
# CHILDFATMASS9_13092016.sav, CHILDBLOODPRESSURE9_21042016.sav,
# CHILDSERUM9_01082017.sav, CHILDGROWTH9_04062020.sav, CHILDOBESITY9_04062020.sav,
# CHILD-ALLGENERALDATA_07072020.sav 

#### ---------------------------- Dependencies ---------------------------- ####

# First, let's point to the necessary libraries
library(foreign)
library(stats)

# Defining the path to the data
# check if the path to the data is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

# ATTENTION!!! If prompted with an "Enter path to data:" message -> Enter the location
# of your datafiles. The code assumes that all (raw) data is stored in ONE folder.
# Do not forget the final slash in your path, and, speaking of slashes, beware of 
# OS sensitive changes when you want to modify the structure of your dirs!

#### ------------------------------ FUNCTIONS ----------------------------- ####

# Read in the data and fix all SPSS weird missing codes into NAs
readquick <- function(filename, rootdir = pathtodata, exclude_col = "") { 
  dat <- read.spss(paste(rootdir, filename, sep=""), # only works for SPSS files
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
# Some cardio-metabolic variables need their SDSs to be adjusted based on sex and age
# This function performs the correction using the residuals of a linear regression of 
# age on the component of interest, stratified by sex, as done by Voortman et al., (2016). 
# However, It also allows to perform correction for age only, sex only or both.

age_sex_adjust <- function(x, adjust_for_sex = TRUE, adjust_for_age = TRUE) {
  # ADJUST FOR AGE ONLY
  if (adjust_for_age == TRUE & adjust_for_sex == FALSE) {
    zscore = ifelse(is.na(cardio_metab[x]) == TRUE, yes = NA,
                    no = rstandard(lm (paste(x , "~ agechild9_visit1"), na.action=na.exclude, data = cardio_metab)))
  # ADJUST FOR SEX ONLY
  } else if (adjust_for_sex == TRUE & adjust_for_age == FALSE) {
    zscore = ifelse(cardio_metab$sex == 1, # = boy
                  yes = scale(cardio_metab[x], center = mean(cardio_metab[cardio_metab$sex == 1, x], na.rm =T),
                              scale = sd(cardio_metab[cardio_metab$sex == 1, x], na.rm =T)),
                  no = scale(cardio_metab[x], center = mean(cardio_metab[cardio_metab$sex == 2, x], na.rm =T),
                             scale = sd(cardio_metab[cardio_metab$sex == 2, x], na.rm =T)))
  # ADJUST FOR SEX & AGE
  } else { boys = cardio_metab[cardio_metab$sex == 1, ]; girls = cardio_metab[cardio_metab$sex == 2, ]; # Split dataset by sex
    zscore = ifelse(is.na(cardio_metab[x]) == TRUE, yes = NA,
                    no = ifelse(cardio_metab$sex == 1, # = boy
                                yes = rstandard(lm (paste(x , "~ agechild9_visit1"), na.action=na.exclude, data = boys)),
                                no = rstandard(lm (paste(x , "~ agechild9_visit1"), na.action=na.exclude, data = girls)) ))
  }
return(as.numeric(zscore))
}

# -----------------------------------------------------------------------------#
##                        Ok, we are good to go now!                          ##
# -----------------------------------------------------------------------------#

# First let's retrieve the variable indicating child sex 
child_general <- readquick("CHILD-ALLGENERALDATA_07072020.sav") # 9901 obs of 122 
sex = child_general[,c('idc', 'gender')]; colnames(sex)[2] <- "sex"; # 1 = boy; 2 = girl.


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
                       'triglyceridenchild9_log', # using the log, because of non-normality
                       'hdlcholesterolchild9_clean',
                       'insulinechild9_log')] # using the log, because of non-normality

# Merge them into one dataset
cardio_metab <- Reduce(function(x,y) merge(x = x, y = y, by = 'idc',  all.x = TRUE),
                          list(sex, fat_mass, blood_pressure, trig_hdlc_ins)) 
# Let's make it a bit more reader friendly 
colnames(cardio_metab)[4:9] <- c("fat_mass","sbp","dbp","triglycerides","hdl_cholesterol","insuline")

# ------------------------------------------------------------------------------
# Before we can use them in the analysis, the outcome variables need to be standardized. 
# so, here we take the standard deviation score using age_sex_adjust()

cardio_metab$fat_mass_z <- age_sex_adjust("fat_mass", adjust_for_sex = FALSE)
cardio_metab$sbp_z = age_sex_adjust("sbp", adjust_for_sex = FALSE)
cardio_metab$dbp_z = age_sex_adjust("dbp", adjust_for_sex = FALSE)
cardio_metab$triglycerides_z = age_sex_adjust("triglycerides", adjust_for_sex = FALSE)
cardio_metab$hdl_cholesterol_z = age_sex_adjust("hdl_cholesterol", adjust_for_sex = FALSE)
cardio_metab$insuline_z = age_sex_adjust("insuline", adjust_for_sex = FALSE)

# ------------------------------------------------------------------------------
# A continuous cumulative cardio-metabolic risk (CMR) index (Voortman et al., 2016)
# is calculated as the sum of age- and sex-specific SD-scores (SDS) of: android fat
# mass (BF); systolic or diastolic blood pressure (SBP/DBP); high-density lipoprotein-
# cholesterol (HDL-C); and the log-transformed triglycerides (TG) and insulin levels.
# The formula for the cardio-metabolic risk factor is:
# BF% SDS + 0.5 × SBP SDS + 0.5 × DBP SDS + TG SDS + (−1 × HDL−C SDS) + insulin SDS.

# Compute number of missing cardio-metabolic values for each subject
cardio_metab$cmr_miss <- rowSums(is.na(cardio_metab[tail(names(cardio_metab), 6)]))

# SBP and DBP are multiplied by 0.5 so that each contributes 1/2 to the blood pressure component.
cardio_metab$sbp_z_weighted = 0.5 * cardio_metab$sbp_z
cardio_metab$dbp_z_weighted = 0.5 * cardio_metab$dbp_z
# HDL-C is multiplied by − 1 since higher HDL-C concentrations reflect lower cardio-metabolic risk
cardio_metab$hdl_cholesterol_z_reversed = -1 * cardio_metab$hdl_cholesterol_z

# Ok, let's do this CRM score computation
cardio_metab$cmr_score <- ifelse(cardio_metab$cmr_miss <= 1,
                                yes = rowSums(cardio_metab[c("fat_mass_z",
                                              "sbp_z_weighted", "dbp_z_weighted",
                                              "triglycerides_z",
                                              "hdl_cholesterol_z_reversed",
                                              "insuline_z")], na.rm = TRUE),
                                no = NA)

# A higher score reflects higher cardio-metabolic risk. When dichotomized, the score
# can be seen as a proxy of child metabolic syndrome (van Gijssel et al., 2016).
# For more specifics about the measurements see (Boyd et al., 2013; Voortman et al., 2016).

# ------------------------------------------------------------------------------
                    ## OPTIONAL: check score consistency ##
# ------------------------------------------------------------------------------
# # Check the correlation between components and, once we are at it, let's peak into 
# # the correlations between internalizing and cardio-metabolic vars
# 
# # Build a correlation matrix
# test_cor = cor(cardio_metab, method = "pearson", use = "complete.obs")
# 
# # Represent the matrix I just built with an heatmap, but first, 
# # the default colors in R heatmap are horrifying so let's fix them
# better_colors <- colorRampPalette(c("blue", "white", "red"))(n = 299)
# heatmap(test_cor, col = better_colors, Rowv = NA, symm = T, scale = 'none')
# 
############# NEED INTERNALIZING TOO FOR THIS STEP! 
# # Scatterplot: Internalizing and CMR score
# plot(cardio_metab$intern_score_z, cardio_metab$cmr_score, 
#      main = paste("PCM Outcomes: corr =", round(test_cor[12,1],4)),
#      xlab="Internalizing score (z score) ", ylab="Cardio-metabolic risk (z score)", 
#      col=rgb(0,100,0,50,maxColorValue=255), pch=16)
# abline(lm(cmr_score ~ intern_score_z, data = cardio_metab, ), col = 'red')
# # For more exploratory plot-around, check the jupyter notebook in the multimorbidity folder

# ------------------------------------------------------------------------------
                          ## OPTIONAL: include BMI ##
# ------------------------------------------------------------------------------
# # Since the score did not seem to work very well, I also include BMI (@ 9y) and it's 
# # categorization into underweight (1,2,3), normal weight, overweight and obese
# growth_9yrs <- readquick('CHILDGROWTH9_04062020.sav')
# obesity_9yrs <- readquick('CHILDOBESITY9_04062020.sav') # Weight status (Cole)
#   # (-3) Thin, grade 3; (-2) Thin, grade 2; (-1) Thin, grade 1; (0) Normal weight; (1) Overweight; (2) Obese;
# bmi_9yrs = growth_9yrs[,c('idc', 'sdsbmiforage9childt')] # Body Mass Index (age 20w - 20y) The Netherlands, 2010 (Talma???)
# 
# # Not used # obesity_5yrs <- readquick('CHILDOBESITY5_08102014.sav') 
# # obese_forever = merge(obesity_5yrs, obesity_9yrs, by = 'idc', all.x = T)
# 
# cardio_metab <- Reduce(function(x,y) merge(x = x, y = y, by = 'idc',  all.x = TRUE),
#                        list(cardio_metab, bmi_9yrs, obesity_9yrs)) 
# colnames(cardio_metab)[21:22] <- c("bmi_z_4age","bmi_binned")

# ------------------------------------------------------------------------------

################################################################################
#### --------------------------- save and run ----------------------------- ####
################################################################################

# Save the dataset in an .rds file, in the directory where the raw data are stored
# saveRDS(cardio_metab, paste(pathtodata,'CMR_score_ageadj.rds', sep =""))

