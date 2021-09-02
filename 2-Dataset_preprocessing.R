
# Hi again, 
# the following script constructs the main outcome variable (risk_group) and cumulative 
# exposure that will be used for the analysis of the association between ELS and 
# psycho-cardio-metabolic multi-morbidity in children. 
# I will first merge all the previously created datafiles, then construct the main 
# outcome variable (risk_group) and cumulative exposures (prenatal_stress and postnatal_stress).
# Finally I will run a non-parametric permutation test to assess whether the Multimorbid
# group identified is greater in size than expected by chance alone. 
# Missing values will be then imputed in the next script (see 3-Imputation.R)

# This time you will only need three files: the "PCM_allvars.rds" we build using the 
# "PCM_outcomes_covs_aux.R" script; "prenatal_stress.rds" and "postnatal_stress.rds", 
# for which, if you want to know more, check out https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

source('0-Functions.R')

################################################################################
# Load datasets
pre_risk  <- readRDS(paste0(pathtodata, 'prenatal_stress.rds'))
post_risk <- readRDS(paste0(pathtodata, 'postnatal_stress.rds'))
outcome   <- readRDS(paste0(pathtodata, 'PCM_allvars.rds'))

# merge them
ELS_PCM <- Reduce(function(x,y) merge(x = x, y = y, by = c('IDC', 'IDM'), all = T),
                  list(pre_risk, post_risk, outcome) ) 

################################################################################
#### -------------- Construct CUMULATIVE STRESS variables ----------------- ####
################################################################################

# compute sum scores for prenatal and postnatal stress exposure
ELS_PCM$prenatal_stress <- rowSums(ELS_PCM[,c("pre_life_events", 
                                              "pre_contextual_risk", 
                                              "pre_parental_risk", 
                                              "pre_interpersonal_risk")], na.rm = F)

ELS_PCM$postnatal_stress <- rowSums(ELS_PCM[,c("post_life_events", 
                                               "post_contextual_risk", 
                                               "post_parental_risk", 
                                               "post_interpersonal_risk", 
                                               "post_direct_victimization")], na.rm = F)

################################################################################
#### ------------------- Construct RISK GROUPS variable ------------------- ####
################################################################################

# Compute >80th percentile internalizing and fat mass scores. 
ELS_PCM$int = ifelse(ELS_PCM$intern_score_z > quantile(ELS_PCM$intern_score_z, probs = 0.8, na.rm = T), 1, 0) #  867 risk, 4055 no risk
# ELS_PCM$fat = ifelse(ELS_PCM$fat_mass_z     > quantile(ELS_PCM$fat_mass_z,     probs = 0.8, na.rm = T), 1, 0) # 1128 risk, 4509 no risk
ELS_PCM$fat = ifelse(ELS_PCM$fmi_z > quantile(ELS_PCM$fmi_z, probs = 0.8, na.rm = T), 1, 0) # 1125 risk, 4498 no risk

# Compute groups
ELS_PCM$risk_groups = rep(NA, nrow(ELS_PCM))
for (i in 1:nrow(ELS_PCM)) {
  if (is.na(ELS_PCM$int[i]) | is.na(ELS_PCM$fat[i]))    { ELS_PCM$risk_groups[i] = NA
  } else if (ELS_PCM$int[i] == 0 & ELS_PCM$fat[i] == 0) { ELS_PCM$risk_groups[i] = 0   # Healthy
  } else if (ELS_PCM$int[i] == 1 & ELS_PCM$fat[i] == 0) { ELS_PCM$risk_groups[i] = 1   # High internalizing  only
  } else if (ELS_PCM$int[i] == 0 & ELS_PCM$fat[i] == 1) { ELS_PCM$risk_groups[i] = 2   # High fat mass only
  } else {                                                ELS_PCM$risk_groups[i] = 3 } # Multimorbid
}

# Let's factor that bad boy 
ELS_PCM$risk_groups = factor(ELS_PCM$risk_groups, 
                         levels = c(0:3), 
                         labels = c("healthy", "internalizing_only", "cardiometabolic_only", "multimorbid"))

summary(ELS_PCM$risk_groups)

#------------------------------------------------------------------------------#
# ------------------------- PERMUTATION TESTING -------------------------------#
#------------------------------------------------------------------------------#

count <- 0 
itarations <- 1000
set.seed(3100896)

for (i in 1:itarations) {
  origN <- unname(summary(ELS_PCM$risk_groups)[4])
  randN <- permute(ELS_PCM)
  if (randN > origN) {
    count <- count + 1
  }
}

pval = format(round(count / itarations, 10), nsmall = 3); cat("P-value: ", pval)

# Check the flowchart of sample selection but do not actually perform exclusion, 
# this will be done after imputation, but we are going to use these numbers as a 
# sanity check
finalsample <- flowchart(ELS_PCM, return_selected_sample = T)

summary(finalsample$risk_groups)

# Recompute permutation on the final sample
count <- 0 
itarations <- 1000
set.seed(3100896)

for (i in 1:itarations) {
  origN <- unname(summary(finalsample$risk_groups)[4])
  randN <- permute(finalsample)
  if (randN > origN) {
    count <- count + 1
  }
}

pval = format(round(count / itarations, 10), nsmall = 3); cat("P-value: ", pval)


################################################################################
################################################################################
################################################################################
# Create two more groups to account for the effect of tress on low fat mass

# # Compute groups 
# data$int = ifelse(data$intern_score_z > quantile(data$intern_score_z, probs = 0.8), 1, 0) # 590 risk, 2780 no risk
# data$fat = ifelse(data$fat_mass_z > quantile(data$fat_mass_z, probs = 0.8), 1, 0) # 674 risk, 2696 no risk
# data$lowfat = ifelse(data$fat_mass_z < quantile(data$fat_mass_z, probs = 0.2), 1, 0)
# 
# 
# data$groups = rep(NA, 3371)
# for (i in 1:3371) {
#   if (data$int[i] == 0 & data$fat[i] == 0 & data$lowfat[i] == 0) { data$groups[i] = 0 # healthy
#   } else if (data$int[i] == 0 & data$fat[i] == 0 & data$lowfat[i] == 1) { data$groups[i] = 1 # low fat only 
#   } else if (data$int[i] == 1 & data$fat[i] == 0 & data$lowfat[i] == 1) {  data$groups[i] = 2 # anorexia type (high internalizing low fat)
#   } else if (data$int[i] == 1 & data$fat[i] == 0 & data$lowfat[i] == 0) { data$groups[i] = 3 # High intenalizing  only
#   } else if (data$int[i] == 0 & data$fat[i] == 1 & data$lowfat[i] == 0) { data$groups[i] = 4 # High fat mass only
#   } else { data$groups[i] = 5 } # multimorbid
# }
# 
# data$groups = as.factor(data$groups)
# summary(data$groups)

# # Let's first factor that bad boy 
# data$groups = factor(data$groups,
#                      levels = c(0:5),
#                      labels = c("healthy", "low_fat", "anorexic_type", "internalizing_only", "cardiometabolic_only", "multimorbid"))
# 
# # display groups
# attach(data); plot(intern_score_z, fat_mass_z,
#                    col=c("cornflowerblue","black", "red", "darkgrey", "darkgoldenrod2","chartreuse4")[groups]); detach(data)
# 
# saveRDS(data, paste(pathtodata,"subgroup6.rds")



################################################################################
#### --------------------------- save and run ----------------------------- ####
################################################################################

# Save the dataset in an .rds file, in the directory where the raw data are stored
saveRDS(ELS_PCM, paste0(pathtodata,'ELSPCM_dataset.rds'))
