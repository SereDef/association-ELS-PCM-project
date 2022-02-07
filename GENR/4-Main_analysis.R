
############################# STATISTICAL ANALYSES #############################

# Load required packages
utilis <- c('nnet', 'CMAverse', 'mice', 'robustbase', 'openxlsx')
lapply(utilis, require, character.only = T);

# ============================================================================ #
# Load mids object
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_sample.rds' file 
  dirname <- dirname(imp_path)
  imp <- readRDS(imp_path)
}

# outlier exclusion
# imp <- miceadds::subset_datlist( imp, subset = ! imp$data$IDC %in% c(5105), toclass = 'mids')

# Transform into a datalist
dat <- miceadds::mids2datlist(imp)

# ============================================================================ #
# Define a customized function that pools regression results and builds a dataframe 
# that is easy to read and save
pool_fit <- function(fit, title, logm = FALSE, rev = FALSE, r_squared = TRUE) {
  p_fit <- mice::pool(fit) # pool results 
  mod <- summary(p_fit) # extract relevant information
  mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
  if (logm == F) {
    mod$lci <- round((mod$estimate - 1.96*mod$std.error), 4)
    mod$uci <- round((mod$estimate + 1.96*mod$std.error), 4)
    if (r_squared == T) {
      mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
      mod$rsq_adj <- c(pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)) # adjusted R2
    }
  } else {
    if (rev == T) { levels(mod$y.level) <- c("M:healthy", "M:intern", "M:fatmas") 
    } else { levels(mod$y.level) <- c("H:intern", "H:fatmas", "H:multim") } # make group comparisons easier to read
    mod$OR <- round(exp(-(mod$estimate)), 4)
    mod$lci <- round(exp((-(mod$estimate)) - 1.96*mod$std.error), 4)
    mod$uci <- round(exp((-(mod$estimate)) + 1.96*mod$std.error), 4)
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
    
  }
  print(mod)
  mod <- cbind(data.frame("model" = rep(title, nrow(mod))), mod)
  mod[nrow(mod)+1,] <- NA
  return(mod)
}

# ============================================================================ #
####################### MAIN ANALYSIS: STEPWISE APPROACH #######################
# ============================================================================ #

# Predict internalizing at 13 using covariates only 
fit_int_cov <- with(imp, lm(intern_score_13_z ~ 
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using prenatal stress only 
fit_int_pre <- with(imp, lm(intern_score_13_z ~ prenatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using postnatal stress only 
fit_int_pos <- with(imp, lm(intern_score_13_z ~ postnatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using an additive prenatal and postnatal model
fit_int_add <- with(imp, lm(intern_score_13_z ~ prenatal_stress_z + postnatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using prenatal by postnatal model
fit_int_int <- with(imp, lm(intern_score_13_z ~ prenatal_stress_z * postnatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
int_cov <- pool_fit(fit_int_cov, "Internalizing - base model")
int_pre <- pool_fit(fit_int_pre, "Internalizing - pren model")
int_pos <- pool_fit(fit_int_pos, "Internalizing - post model")
int_add <- pool_fit(fit_int_add, "Internalizing - pren+post model")
int_int <- pool_fit(fit_int_int, "Internalizing - pren*post model")
# combine them into a single frame
intern <- rbind(int_cov,int_pre,int_pos,int_add,int_int)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Predict fat mass percentage at 13 using covariates only 
fit_fat_cov <- with(imp, lm(tot_fat_percent_13_z ~ 
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal stress only 
fit_fat_pre <- with(imp, lm(tot_fat_percent_13_z ~ prenatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using postnatal stress only 
fit_fat_pos <- with(imp, lm(tot_fat_percent_13_z ~ postnatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using an additive prenatal and postnatal model
fit_fat_add  <- with(imp, lm(tot_fat_percent_13_z ~ prenatal_stress_z + postnatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal by postnatal model
fit_fat_int <- with(imp, lm(tot_fat_percent_13_z ~ prenatal_stress_z * postnatal_stress_z +
                    sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
fat_cov <- pool_fit(fit_fat_cov, "Fat mass % - base model")
fat_pre <- pool_fit(fit_fat_pre, "Fat mass % - pren model")
fat_pos <- pool_fit(fit_fat_pos, "Fat mass % - post model")
fat_add <- pool_fit(fit_fat_add, "Fat mass % - pren+post model")
fat_int <- pool_fit(fit_fat_int, "Fat mass % - pren*post model")
# combine them into a single frame
fatmas <- rbind(fat_cov,fat_pre,fat_pos,fat_add,fat_int)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 covariats only 
fit_grps_cov <- with(imp,  nnet::multinom(risk_groups_perc_REC ~
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre <- with(imp,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos <- with(imp,  nnet::multinom(risk_groups_perc_REC ~ postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add  <- with(imp,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z + postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int <- with(imp,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z * postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Pool the results
grps_cov <- pool_fit(fit_grps_cov, "Comorbidity - base model", logm =T, rev =T)
grps_pre <- pool_fit(fit_grps_pre, "Comorbidity - pren model", logm =T, rev =T)
grps_pos <- pool_fit(fit_grps_pos, "Comorbidity - post model", logm =T, rev =T)
grps_add <- pool_fit(fit_grps_add, "Comorbidity - pren+post model", logm =T, rev =T)
grps_int <- pool_fit(fit_grps_int, "Comorbidity - pren*post model", logm =T, rev =T)
# combine them into a single frame
comorb <- rbind(grps_cov,grps_pre,grps_pos,grps_add,grps_int)

# ============================================================================ #
####################### ADDITIONAL (1): STRATIFY BY SEX ########################
# ============================================================================ #

impM <- miceadds::subset_datlist(imp, subset = imp$data$sex == 1,  toclass = 'mids') 
impF <- miceadds::subset_datlist(imp, subset = imp$data$sex == 2,  toclass = 'mids') 

# ================================== MALES =================================== #
#  Predict internalizing at 13 using covariates only 
fit_int_cov.m <- with(impM, lm(intern_score_13_z ~ 
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_int_pre.m <- with(impM, lm(intern_score_13_z ~ prenatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
#  Predict internalizing at 13 using postnatal stress only 
fit_int_pos.m <- with(impM, lm(intern_score_13_z ~ postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using an additive prenatal and postnatal model
fit_int_add.m <- with(impM, lm(intern_score_13_z ~ prenatal_stress_z + postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using prenatal by postnatal model
fit_int_int.m <- with(impM, lm(intern_score_13_z ~ prenatal_stress_z * postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
int_cov.m <- pool_fit(fit_int_cov.m, "Male: internalizing - base model")
int_pre.m <- pool_fit(fit_int_pre.m, "Male: internalizing - pren model")
int_pos.m <- pool_fit(fit_int_pos.m, "Male: internalizing - post model")
int_add.m <- pool_fit(fit_int_add.m, "Male: internalizing - pren+post model")
int_int.m <- pool_fit(fit_int_int.m, "Male: internalizing - pren*post model")
# combine them into a single frame
intern.m <- rbind(int_cov.m,int_pre.m,int_pos.m,int_add.m,int_int.m)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict fat mass percentage at 13 using covariates only 
fit_fat_cov.m <- with(impM, lm(tot_fat_percent_13_z ~
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal stress only 
fit_fat_pre.m <- with(impM, lm(tot_fat_percent_13_z ~ prenatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using postnatal stress only 
fit_fat_pos.m <- with(impM, lm(tot_fat_percent_13_z ~ postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using an additive prenatal and postnatal model
fit_fat_add.m  <- with(impM, lm(tot_fat_percent_13_z ~ prenatal_stress_z + postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal by postnatal model
fit_fat_int.m <- with(impM, lm(tot_fat_percent_13_z ~ prenatal_stress_z * postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
fat_cov.m <- pool_fit(fit_fat_cov.m, "Male: fat mass % - base model")
fat_pre.m <- pool_fit(fit_fat_pre.m, "Male: fat mass % - pren model")
fat_pos.m <- pool_fit(fit_fat_pos.m, "Male: fat mass % - post model")
fat_add.m <- pool_fit(fit_fat_add.m, "Male: fat mass % - pren+post model")
fat_int.m <- pool_fit(fit_fat_int.m, "Male: fat mass % - pren*post model")
# combine them into a single frame
fatmas.m <- rbind(fat_cov.m,fat_pre.m,fat_pos.m,fat_add.m,fat_int.m)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 using covariates only 
fit_grps_cov.m <- with(impM,  nnet::multinom(risk_groups_perc_REC ~ 
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.m <- with(impM,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.m <- with(impM,  nnet::multinom(risk_groups_perc_REC ~ postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.m <- with(impM,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z + postnatal_stress_z +
                        age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int.m <- with(impM,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z * postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Pool the results
grps_cov.m <- pool_fit(fit_grps_cov.m, "Male: comorbidity - base model", logm =T, rev =T)
grps_pre.m <- pool_fit(fit_grps_pre.m, "Male: comorbidity - pren model", logm =T, rev =T)
grps_pos.m <- pool_fit(fit_grps_pos.m, "Male: comorbidity - post model", logm =T, rev =T)
grps_add.m <- pool_fit(fit_grps_add.m, "Male: comorbidity - pren+post model",logm =T, rev =T)
grps_int.m <- pool_fit(fit_grps_int.m, "Male: comorbidity - pren*post model",logm =T, rev =T)
# combine them into a single frame
comorb.m <- rbind(grps_cov.m,grps_pre.m,grps_pos.m,grps_add.m,grps_int.m)

# ================================= FEMALES ================================== #
#  Predict internalizing at 13 using covariates only 
fit_int_cov.f <- with(impF, lm(intern_score_13_z ~ 
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
#  Predict internalizing at 13 using prenatal stress only 
fit_int_pre.f <- with(impF, lm(intern_score_13_z ~ prenatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
#  Predict internalizing at 13 using postnatal stress only 
fit_int_pos.f <- with(impF, lm(intern_score_13_z ~ postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using an additive prenatal and postnatal model
fit_int_add.f <- with(impF, lm(intern_score_13_z ~ prenatal_stress_z + postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using prenatal by postnatal model
fit_int_int.f <- with(impF, lm(intern_score_13_z ~ prenatal_stress_z * postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
int_cov.f <- pool_fit(fit_int_cov.f, "Female: internalizing - base model")
int_pre.f <- pool_fit(fit_int_pre.f, "Female: internalizing - pren model")
int_pos.f <- pool_fit(fit_int_pos.f, "Female: internalizing - post model")
int_add.f <- pool_fit(fit_int_add.f, "Female: internalizing - pren+post model")
int_int.f <- pool_fit(fit_int_int.f, "Female: internalizing - pren*post model")
# combine them into a single frame
intern.f <- rbind(int_cov.f,int_pre.f,int_pos.f,int_add.f,int_int.f)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict fat mass percentage at 13 using covariates only 
fit_fat_cov.f <- with(impF, lm(tot_fat_percent_13_z ~
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal stress only 
fit_fat_pre.f <- with(impF, lm(tot_fat_percent_13_z ~ prenatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using postnatal stress only 
fit_fat_pos.f <- with(impF, lm(tot_fat_percent_13_z ~ postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using an additive prenatal and postnatal model
fit_fat_add.f  <- with(impF, lm(tot_fat_percent_13_z ~ prenatal_stress_z + postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal by postnatal model
fit_fat_int.f <- with(impF, lm(tot_fat_percent_13_z ~ prenatal_stress_z * postnatal_stress_z +
                      age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
fat_cov.f <- pool_fit(fit_fat_cov.f, "Female: fat mass % - base model")
fat_pre.f <- pool_fit(fit_fat_pre.f, "Female: fat mass % - pren model")
fat_pos.f <- pool_fit(fit_fat_pos.f, "Female: fat mass % - post model")
fat_add.f <- pool_fit(fit_fat_add.f, "Female: fat mass % - pren+post model")
fat_int.f <- pool_fit(fit_fat_int.f, "Female: fat mass % - pren*post model")
# combine them into a single frame
fatmas.f <- rbind(fat_cov.f,fat_pre.f,fat_pos.f,fat_add.f,fat_int.f)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 using covariates only 
fit_grps_cov.f <- with(impF,  nnet::multinom(risk_groups_perc_REC ~ 
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.f <- with(impF,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.f <- with(impF,  nnet::multinom(risk_groups_perc_REC ~ postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.f <- with(impF,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z + postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int.f <- with(impF,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z * postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Pool the results
grps_cov.f <- pool_fit(fit_grps_cov.f, "Female: comorbidity - base model", logm =T, rev =T)
grps_pre.f <- pool_fit(fit_grps_pre.f, "Female: comorbidity - pren model", logm =T, rev =T)
grps_pos.f <- pool_fit(fit_grps_pos.f, "Female: comorbidity - post model", logm =T, rev =T)
grps_add.f <- pool_fit(fit_grps_add.f, "Female: comorbidity - pren+post model",logm =T, rev =T)
grps_int.f <- pool_fit(fit_grps_int.f, "Female: comorbidity - pren*post model",logm =T, rev =T)
# combine them into a single frame
comorb.f <- rbind(grps_cov.f,grps_pre.f,grps_pos.f,grps_add.f,grps_int.f)

# ============================================================================ #
################# ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION ###################
# ============================================================================ #
# Let's have a look at the contribution of individual domains of stress to each outcome.

# Internalizing 
fit_int_dm <- with(imp, lm(intern_score_13_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_int_dm.m <- with(impM, lm(intern_score_13_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_int_dm.f <- with(impF, lm(intern_score_13_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Fat mass percentage
fit_fat_dm <- with(imp, lm(tot_fat_percent_13_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                   post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                   sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fat_dm.m <- with(impM, lm(tot_fat_percent_13_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                   post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                   age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fat_dm.f <- with(impF, lm(tot_fat_percent_13_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                   post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                   age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Risk groups
fit_grp_dm <- with(imp, nnet::multinom(risk_groups_perc_REC ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                   post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                   sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
fit_grp_dm.m <- with(impM, nnet::multinom(risk_groups_perc_REC ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
fit_grp_dm.f <- with(impF, nnet::multinom(risk_groups_perc_REC ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Pool the results
dm_int <- pool_fit(fit_int_dm, "Internalizing - domain model")
dm_fat <- pool_fit(fit_fat_dm, "Fat mass % - domain model")
dm_grp <- pool_fit(fit_grp_dm, "Comorbidity - domain model", logm = T, rev=T)
dm_int.m <- pool_fit(fit_int_dm.m, "Male: internalizing - domain model")
dm_fat.m <- pool_fit(fit_fat_dm.m, "Male: fat mass % - domain model")
dm_grp.m <- pool_fit(fit_grp_dm.m, "Male: comorbidity - domain model", logm = T, rev=T)
dm_int.f <- pool_fit(fit_int_dm.f, "Female: internalizing - domain model")
dm_fat.f <- pool_fit(fit_fat_dm.f, "Female: fat mass % - domain model")
dm_grp.f <- pool_fit(fit_grp_dm.f, "Female: comorbidity - domain model", logm = T, rev=T)
# combine them into a single frame
intern.dm <- rbind(dm_int, dm_int.m, dm_int.f)
fatmas.dm <- rbind(dm_fat, dm_fat.m, dm_fat.f)
comorb.dm <- rbind(dm_grp, dm_grp.m, dm_grp.f)

# ============================================================================ #
################## ADDITIONAL (3): CAUSAL MEDIATION ANALYSIS ###################
# ============================================================================ #

cma_loop <- function(outcome, datalist, ymod = "linear", allow_interaction = F, n_imputations = 30) {
  # initiate lists where model fits per imputation will be stored
  yregr_fit <- list(); mregr_fit <- list(); decom_fit <- list()
  drink_fit <- list(); smoke_fit <- list()
  # loop through imputed datasets
  for (i in 1:n_imputations) {
    message("\nImputation number: ", i) # update
    est <- CMAverse::cmest(data = dat[[i]], full = F,
           exposure = "prenatal_stress_z",  mediator = c("postnatal_stress_z"),
           outcome = outcome,  yreg = ymod, # outcome regression model (linear or multimod)
           yval = "multimorbid", # value of outcome at which causal effects on the risk/odds ratio scale are estimated
           EMint = allow_interaction, # exposure-mediator interaction in yreg (i.e. prenatal*postnatal). 
           model = "gformula", 
           # exposure-outcome confounder(s), exposure-mediator confounder(s) and mediator-outcome confounder(s) not affected by the exposure
           basec = c("sex", "age_child", "ethnicity", "m_bmi_before_pregnancy"), 
           # mediator-outcome confounder(s) affected by the exposure following temporal order
           postc = c("m_drinking", "m_smoking"), # regression model for each variable in postc (gformula): 
           postcreg = list(glm(m_drinking ~ prenatal_stress_z + ethnicity + m_bmi_before_pregnancy, family = gaussian(), data = dat[[i]]), 
                           glm(m_smoking ~ prenatal_stress_z + ethnicity + m_bmi_before_pregnancy, family = gaussian(), data = dat[[i]])), 
           mreg = list("linear"), # regression model for each mediator (i.e. postnatal ~ prenatal + covariates). 
           astar = 0, # control value for the exposure (i.e. prenatal stress = mean)
           a = 1, # the active value for the exposure (i.e. prenatal stress = + 1sd)
           mval = list(0), # control value for mediator (i.e. postnatal stress = mean).
           estimation = "imputation", # method for estimating causal effects. paramfunc is alternative. 
           inference = "bootstrap", # method for estimating standard errors of causal effects. delta is alternative. 
           nboot = 800) # Default is 200. # boot.ci.type # (percentile or bias-corrected and accelerated (BCa))
    s <- summary(est)
    
    yregr_fit[[i]] <- s[["reg.output"]][["yreg"]]
    mregr_fit[[i]] <- s[["reg.output"]][["mreg"]][[1]]
    drink_fit[[i]] <- s[["reg.output"]][["postcreg"]][[1]]
    smoke_fit[[i]] <- s[["reg.output"]][["postcreg"]][[2]]
    decom_fit[[i]] <- s[["summarydf"]] 
  }
  
  setClass("Totfit", slots = list(yreg ='list', mreg ='list', c_dr ='list', c_sm ='list', decomp ='list'))
  totfit <- new("Totfit", yreg = yregr_fit, mreg = mregr_fit, c_dr = drink_fit, c_sm = smoke_fit, decomp = decom_fit)
  
  return(totfit)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pool_cma_regs <- function(fit, title) {
  p = mice::pool(mice::as.mira(fit)) 
  s = summary(p)
  s$sign <- ifelse(s$p.value < 0.05, '*', '') # add a column to highlight significant terms
  print(s)
  s <- cbind(data.frame("model" = rep(title, nrow(s))), s)
  s[nrow(s)+1,] <- NA
  return(s)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rubin’s Rules # https://bookdown.org/mwheymans/bookmi/rubins-rules.html
pool_cma_decomp <- function(fit, n_parameters = 7) {
  m = length(fit@decomp); # number of imputations
  n = nrow(dat[[1]]) # sample size
  p = n_parameters
  
  # Extract estimates and SEs for all imputations
  est <- sapply(fit@decomp, "[", , "Estimate") 
  se <-  lapply(fit@decomp, "[", , "Std.error") 
  
  # Estimate mean
  d <- Reduce("+", fit@decomp) / m
  # add column for effect names
  d <- cbind("Effect" = rownames(fit@decomp[[1]]), d)
  
  # Within imputation variance
  vw <- rowMeans(sapply(se, function(x) x^2))
  # Between imputation variance
  vb = c()
  for (i in 1:p) {
    pooled_est = d$Estimate[i]
    ith_est = est[i, ]
    est_var <- sapply(ith_est, function(x) (x - pooled_est)^2)
    mean_var <- Reduce("+", est_var) / (m - 1)
    vb[i] <- mean_var
  }
  # Total variance
  vt <- vw + vb + (vb/m)
  # Pooled SE
  d$SE_p <- sqrt(vt)
  
  # Fraction of Missing information: proportion of variation in parameter due to the missing data.
  lambda <- (vb + (vb/m)) / vt
  # Pooled degrees of freedom
  df_old <- (m - 1) / lambda^2
  df_obs <- (((n - p) + 1)/((n - p) + 3)) * (n - p)*(1 - lambda)
  df_adj <- (df_old*df_obs) / (df_old+df_obs)
  
  # Wald test pooled
  d$W_p <- (d$Estimate^2) / d$SE_p
  t <- pt(d$W_p, df_adj, lower.tail=F)
  
  # Pooled CIs
  d$CIL_p <- d$Estimate - (t * d$SE_p)
  d$CIU_p <- d$Estimate + (t * d$SE_p)
  # Pooled p-value
  d$Pv_p <- round(t, 4)
  d$sign <- ifelse(d$P.val < 0.05, '*', '')
  
  print(d)
  return(d)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit the CMA model to all imputed datasets
fit_int_cma <- cma_loop('intern_score_13_z', dat)
fit_fat_cma <- cma_loop('tot_fat_percent_13_z', dat)
# fit_fata_cma <- cma_loop('andr_fat_mass_13_z', dat)
# fit_grps_cma <- cma_loop('risk_groups_perc_REC', dat, ymod = "multinomial")

yreg_int_cma <- pool_cma_regs(fit_int_cma@yreg, "Inernalizing - yreg")
yreg_fat_cma <- pool_cma_regs(fit_fat_cma@yreg, "Fat mass % - yreg")
medreg_cma   <- pool_cma_regs(fit_int_cma@mreg, "Mediator regression")
drink_cma    <- pool_cma_regs(fit_int_cma@c_dr, "Drink regression")
smoke_cma    <- pool_cma_regs(fit_int_cma@c_sm, "Smoke regression")
mediation <- rbind(yreg_int_cma, yreg_fat_cma, medreg_cma, drink_cma, smoke_cma)

intern.cma <- pool_cma_decomp(fit_int_cma)
fatmas.cma <- pool_cma_decomp(fit_fat_cma)

# ============================================================================ #

# ============================================================================ #
###################### SENSITIVITY (1): ANDROID FAT MASS #######################
# ============================================================================ #
# Predict android fat mass at 13 using covariates only 
fit_fata_cov <- with(imp, lm(andr_fat_mass_13_z ~ 
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_cov.m <- with(impM, lm(andr_fat_mass_13_z ~ 
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_cov.f <- with(impF, lm(andr_fat_mass_13_z ~ 
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict android fat mass at 13 using prenatal stress only 
fit_fata_pre <- with(imp, lm(andr_fat_mass_13_z ~ prenatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_pre.m <- with(impM, lm(andr_fat_mass_13_z ~ prenatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_pre.f <- with(impF, lm(andr_fat_mass_13_z ~ prenatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict android fat mass at 13 using postnatal stress only 
fit_fata_pos <- with(imp, lm(andr_fat_mass_13_z ~ postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_pos.m <- with(impM, lm(andr_fat_mass_13_z ~ postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_pos.f <- with(impF, lm(andr_fat_mass_13_z ~ postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict android fat mass at 13 using an additive prenatal and postnatal model
fit_fata_add  <- with(imp, lm(andr_fat_mass_13_z ~ prenatal_stress_z + postnatal_stress_z +
                      sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_add.m <- with(impM, lm(andr_fat_mass_13_z ~ prenatal_stress_z + postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_add.f <- with(impF, lm(andr_fat_mass_13_z ~ prenatal_stress_z + postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict android fat mass at 13 using prenatal by postnatal model
fit_fata_int <- with(imp, lm(andr_fat_mass_13_z ~ prenatal_stress_z * postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_int.m <- with(impM, lm(andr_fat_mass_13_z ~ prenatal_stress_z * postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
fit_fata_int.f <- with(impF, lm(andr_fat_mass_13_z ~ prenatal_stress_z * postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
fata_cov <- pool_fit(fit_fata_cov, "Android fat - base model")
fata_pre <- pool_fit(fit_fata_pre, "Androif fat - pren model")
fata_pos <- pool_fit(fit_fata_pos, "Android fat - post model")
fata_add <- pool_fit(fit_fata_add, "Android fat - pre+post model")
fata_int <- pool_fit(fit_fata_int, "Android fat - pre*post model")
fata_cov.m <- pool_fit(fit_fata_cov.m, "Male: android fat - base model")
fata_pre.m <- pool_fit(fit_fata_pre.m, "Male: androif fat - pren model")
fata_pos.m <- pool_fit(fit_fata_pos.m, "Male: android fat - post model")
fata_add.m <- pool_fit(fit_fata_add.m, "Male: android fat - pre+post model")
fata_int.m <- pool_fit(fit_fata_int.m, "Male: android fat - pre*post model")
fata_cov.f <- pool_fit(fit_fata_cov.f, "Female: android fat - base model")
fata_pre.f <- pool_fit(fit_fata_pre.f, "Female: androif fat - pren model")
fata_pos.f <- pool_fit(fit_fata_pos.f, "Female: android fat - post model")
fata_add.f <- pool_fit(fit_fata_add.f, "Female: android fat - pre+post model")
fata_int.f <- pool_fit(fit_fata_int.f, "Female: android fat - pre*post model")
# combine them into a single frame
andfat   <- rbind(fata_cov, fata_pre, fata_pos, fata_add, fata_int)
andfat.m <- rbind(fata_cov.m, fata_pre.m, fata_pos.m, fata_add.m, fata_int.m)
andfat.f <- rbind(fata_cov.f, fata_pre.f, fata_pos.f, fata_add.f, fata_int.f)
# ============================================================================ #

# ============================================================================ #
####################### SENSITIVITY (2): RESPONDENTS ONLY ######################
# ============================================================================ #

impR <- miceadds::subset_datlist(imp, subset = (!is.na(imp$data$intern_score_13_z) & !is.na(imp$data$tot_fat_percent_13_z)),  
                                 toclass = 'mids') 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict internalizing at 13 using covariates only 
fit_int_cov.r <- with(impR, lm(intern_score_13_z ~ 
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using prenatal stress only 
fit_int_pre.r <- with(impR, lm(intern_score_13_z ~ prenatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using postnatal stress only 
fit_int_pos.r <- with(impR, lm(intern_score_13_z ~ postnatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using an additive prenatal and postnatal model
fit_int_add.r <- with(impR, lm(intern_score_13_z ~ prenatal_stress_z + postnatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict internalizing at 13 using prenatal by postnatal model
fit_int_int.r <- with(impR, lm(intern_score_13_z ~ prenatal_stress_z * postnatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
int_cov.r <- pool_fit(fit_int_cov.r, "(Resp) Internalizing - base model")
int_pre.r <- pool_fit(fit_int_pre.r, "(Resp) Internalizing - pren model")
int_pos.r <- pool_fit(fit_int_pos.r, "(Resp) Internalizing - post model")
int_add.r <- pool_fit(fit_int_add.r, "(Resp) Internalizing - pren+post model")
int_int.r <- pool_fit(fit_int_int.r, "(Resp) Internalizing - pren*post model")
# combine them into a single frame
intern.resp <- rbind(int_cov.r,int_pre.r,int_pos.r,int_add.r,int_int.r)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict fat mass percentage at 13 using covariates only 
fit_fat_cov.r <- with(impR, lm(tot_fat_percent_13_z ~ 
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal stress only 
fit_fat_pre.r <- with(impR, lm(tot_fat_percent_13_z ~ prenatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using postnatal stress only 
fit_fat_pos.r <- with(impR, lm(tot_fat_percent_13_z ~ postnatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using an additive prenatal and postnatal model
fit_fat_add.r <- with(impR, lm(tot_fat_percent_13_z ~ prenatal_stress_z + postnatal_stress_z +
                               sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Predict fat mass percentage at 13 using prenatal by postnatal model
fit_fat_int.r <- with(impR, lm(tot_fat_percent_13_z ~ prenatal_stress_z * postnatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
# Pool the results
fat_cov.r <- pool_fit(fit_fat_cov.r, "(Resp) Fat mass % - base model")
fat_pre.r <- pool_fit(fit_fat_pre.r, "(Resp) Fat mass % - pren model")
fat_pos.r <- pool_fit(fit_fat_pos.r, "(Resp) Fat mass % - post model")
fat_add.r <- pool_fit(fit_fat_add.r, "(Resp) Fat mass % - pren+post model")
fat_int.r <- pool_fit(fit_fat_int.r, "(Resp) Fat mass % - pren*post model")
# combine them into a single frame
fatmas.resp <- rbind(fat_cov.r, fat_pre.r, fat_pos.r, fat_add.r, fat_int.r)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 covariats only 
fit_grps_cov.r <- with(impR,  nnet::multinom(risk_groups_perc_REC ~
                  sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.r <- with(impR,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z +
                  sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.r <- with(impR,  nnet::multinom(risk_groups_perc_REC ~ postnatal_stress_z +
                  sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.r <- with(impR,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z + postnatal_stress_z +
                  sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int.r <- with(impR,  nnet::multinom(risk_groups_perc_REC ~ prenatal_stress_z * postnatal_stress_z +
                  sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
# Pool the results
grps_cov.r <- pool_fit(fit_grps_cov.r, "(Resp) Comorbidity - base model", logm =T, rev =T)
grps_pre.r <- pool_fit(fit_grps_pre.r, "(Resp) Comorbidity - pren model", logm =T, rev =T)
grps_pos.r <- pool_fit(fit_grps_pos.r, "(Resp) Comorbidity - post model", logm =T, rev =T)
grps_add.r <- pool_fit(fit_grps_add.r, "(Resp) Comorbidity - pren+post model", logm =T, rev =T)
grps_int.r <- pool_fit(fit_grps_int.r, "(Resp) Comorbidity - pren*post model", logm =T, rev =T)
# combine them into a single frame
comorb.resp <- rbind(grps_cov.r, grps_pre.r, grps_pos.r, grps_add.r, grps_int.r)
# ============================================================================ #

# ============================================================================ #
###################### SENSITIVITY (3): ROBUST REGRESSION ######################
# ============================================================================ #

robfit_int <- with(imp, robustbase::lmrob(intern_score_13_z ~ prenatal_stress_z + postnatal_stress_z +
                  sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
robfit_fat <- with(imp, robustbase::lmrob(tot_fat_percent_13_z ~ prenatal_stress_z + postnatal_stress_z +
                   sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking))
intern.robust <- pool_fit(robfit_int, "(Robust) Internalizing - pren+post model", r_squared = F)
fatmas.robust <- pool_fit(robfit_fat, "(Robust) Fat mass % - pren+post model", r_squared = F)

############ original group healthy as reference 
# ============================================================================ #
# Save results into two excel files

modls <- list("intern" = intern, "fatmas" = fatmas, "comorb" = comorb, 
              "intern_m" = intern.m, "intern_m" = intern.f, 
              "fatmas_m" = fatmas.m, "fatmas_f" = fatmas.f, 
              "comorb_m" = comorb.m, "comorb_f" = comorb.f, 
              "intern_dm" = intern.dm, "fatmas_dm" = fatmas.dm, "comorb_dm" = comorb.dm,
              "intern_cma" = intern.cma, "fatmas_cma" = fatmas.cma, "other_cma" = mediation, 
              "android" = andfat, "android_m" = andfat.m, "android_m" = andfat.f, 
              "intern_select" = intern.resp, "fatmas_select" = fatmas.resp, "comorb_select" = comorb.resp, 
              "intern_robust_reg" = intern.robust, "fatmas_robust_reg" = fatmas.robust)

openxlsx::write.xlsx(modls, file = file.path(dirname, paste0(c("GenR_Results", Sys.Date(), ".xlsx"))), overwrite = T)
                     