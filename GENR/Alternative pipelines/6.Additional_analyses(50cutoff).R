# extend sample size by using 50% cutoff

# Load required packages
utilis <- c('nnet', 'CMAverse', 'mice', 'robustbase', 'openxlsx')
lapply(utilis, require, character.only = T);

#-------------------------------------------------------------------------------
# define a function that randomly shuffles internalizing and fat mass values and 
# returns the new size of the "randomly multimorbid" group.
permute <- function(df) { 
  # create empty dataset for permutation
  perm <- data.frame(1:nrow(df))
  
  perm$new_int <- sample(df$int)
  perm$new_fat <- sample(df$fat)
  perm$new_groups <- ifelse(perm$new_int == 0 & perm$new_fat == 0, 0, 
                            ifelse(perm$new_int == 1 & perm$new_fat == 0, 1, 
                                   ifelse(perm$new_int == 0 & perm$new_fat == 1, 2, 3)))
  new_n <- unname(summary(as.factor(perm$new_groups))[4])
  return(new_n)
}
#-------------------------------------------------------------------------------
# define a function that construncts the groups used as outcome for the third set 
# of models 
construct_grp <- function(int_var, fm_var, df, cutoff = 0.8, permute = T) {
  df$int = ifelse(df[, int_var] > quantile(df[, int_var], probs = cutoff, na.rm = T), 1, 0) 
  df$fat = ifelse(df[, fm_var]  > quantile(df[, fm_var],  probs = cutoff, na.rm = T), 1, 0) 
  
  df$risk_groups = rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    if ( is.na(df$int[i]) | is.na(df$fat[i]) )  { df$risk_groups[i] = NA
    } else if (df$int[i] == 0 & df$fat[i] == 0) { df$risk_groups[i] = 0   # Healthy
    } else if (df$int[i] == 1 & df$fat[i] == 0) { df$risk_groups[i] = 1   # High internalizing  only
    } else if (df$int[i] == 0 & df$fat[i] == 1) { df$risk_groups[i] = 2   # High fat mass only
    } else {                                      df$risk_groups[i] = 3 } # Multimorbid
  }
  # # Let's first factor that bad boy 
  df$risk_groups = factor(df$risk_groups, 
                          levels = c(0:3), 
                          labels = c("healthy", "internalizing_only", "cardiometabolic_only", "multimorbid"))
  message(paste("\nCombining:", int_var, "and", fm_var, "\n"))
  print(summary(df$risk_groups))
  corz = cor(df[, c(int_var, fm_var)], use = 'complete.obs')
  plot(df[, int_var], df[, fm_var], main = paste("Corr =", round(corz[1,2], 2)), xlab = int_var, ylab = fm_var, 
       col = c("darkgreen", "blue", "darkgoldenrod2", "red")[df$risk_groups])
  
  if (permute == T) {
    count <- 0 
    iterations <- 1000
    set.seed(310896)
    
    for (i in 1:iterations) {
      origN <- unname(summary(df$risk_groups)[4])
      randN <- permute(df)
      if (randN > origN) {
        count <- count + 1
      }
    }
    
    pval = format(round(count / iterations, 10), nsmall = 3)
    if (pval == 0.000) { pval = "< .001"}
    cat("Permutation p-value:", pval)
  }
  
  return(df$risk_groups)
}

# ============================================================================ #
# Load mids object
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_sample.rds' file 
  dirname <- dirname(imp_path)
  imp <- readRDS(imp_path)
}

# construct and add the risk group variables in  each imputed dataset. 
# to do so we will use the function construct_grp defined in 0-Functions.R
# We will be doing this only on the selected sample to save some time. The same code can 
# be adapted to use the pre-selection "imputation" variable. 

# convert imputations into a long dataframe
long <- mice::complete(imp, action = "long", include = TRUE) 
# initialize empty long dataframe for return
return_dat <- data.frame() 
#  apply compute_grp to each imputed dataset
for (n_imp in unique(long$.imp)) {
  dset = long[long$.imp == n_imp, ] # one iteration only
  message(paste("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Imputation nr: ", n_imp, 
                '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'))
  # Contruct the risk group varible for each fat measure
  dset$risk_groups_50 <- construct_grp('intern_score_13', 'tot_fat_percent_13', cutoff=0.5, dset)
  
  # Append results to long-type dataframe
  return_dat <- rbind.data.frame(return_dat, dset) 
}

# Pool descriptives 
grp <- return_dat[, c(".imp", "risk_groups_50")]
grp_summary <- with(grp, by(grp, .imp, function(x) summary(x[, 2]))) # exclude .imp
grp_df <- rbind.data.frame(grp_summary[2:31])
grp_pooled <- data.frame(rowMeans(grp_df))

imp <- mice::as.mids(return_dat) # reconvert long dataframe to mids object

rm(dset, long, return_dat, grp, grp_summary, grp_df)

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
# Predict group probabilities at 13 covariats only 
fit_grps_cov <- with(imp,  nnet::multinom(risk_groups_50 ~
                                            sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov <- pool_fit(fit_grps_cov, "Comorbidity - base model", logm =T, rev =T)
rm(fit_grps_cov)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre <- with(imp,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z +
                                            sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre <- pool_fit(fit_grps_pre, "Comorbidity - pren model", logm =T, rev =T)
rm(fit_grps_pre)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos <- with(imp,  nnet::multinom(risk_groups_50 ~ postnatal_stress_z +
                                            sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos <- pool_fit(fit_grps_pos, "Comorbidity - post model", logm =T, rev =T)
rm(fit_grps_pos)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add <- with(imp,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z + postnatal_stress_z +
                                            sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add <- pool_fit(fit_grps_add, "Comorbidity - pren+post model", logm =T, rev =T)
rm(fit_grps_add)
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int <- with(imp,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z * postnatal_stress_z +
                                            sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_int <- pool_fit(fit_grps_int, "Comorbidity - pren*post model", logm =T, rev =T)
rm(fit_grps_int)

# combine them into a single frame
comorb <- rbind(grps_cov,grps_pre,grps_pos,grps_add,grps_int)
rm(grps_cov,grps_pre,grps_pos,grps_add,grps_int)

# ============================================================================ #
# ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 
# Risk groups
fit_grp_dm <- with(imp, nnet::multinom(risk_groups_perc_REC ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                                         post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                         sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
dm_grp <- pool_fit(fit_grp_dm, "Comorbidity - domain model", logm = T, rev=T)
rm(fit_grp_dm)

# ============================================================================ #
####################### ADDITIONAL (1): STRATIFY BY SEX ########################
# ============================================================================ #

impM <- miceadds::subset_datlist(imp, subset = imp$data$sex == 1,  toclass = 'mids') 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 using covariates only 
fit_grps_cov.m <- with(impM,  nnet::multinom(risk_groups_50 ~ 
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov.m <- pool_fit(fit_grps_cov.m, "Male: comorbidity - base model", logm =T, rev =T)
rm(fit_grps_cov.m)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.m <- with(impM,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre.m <- pool_fit(fit_grps_pre.m, "Male: comorbidity - pren model", logm =T, rev =T)
rm(fit_grps_pre.m)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.m <- with(impM,  nnet::multinom(risk_groups_50 ~ postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos.m <- pool_fit(fit_grps_pos.m, "Male: comorbidity - post model", logm =T, rev =T)
rm(fit_grps_pos.m)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.m <- with(impM,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z + postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add.m <- pool_fit(fit_grps_add.m, "Male: comorbidity - pren+post model",logm =T, rev =T)
rm(fit_grps_add.m)
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int.m <- with(impM,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z * postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_int.m <- pool_fit(fit_grps_int.m, "Male: comorbidity - pren*post model",logm =T, rev =T)
rm(fit_grps_int.m)

# combine them into a single frame
comorb.m <- rbind(grps_cov.m, grps_pre.m, grps_pos.m, grps_add.m, grps_int.m)
rm(grps_cov.m, grps_pre.m, grps_pos.m, grps_add.m, grps_int.m)

# ============================================================================ #
# ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 
# Risk groups
fit_grp_dm.m <- with(impM, nnet::multinom(risk_groups_50 ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                                            post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                            age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
dm_grp.m <- pool_fit(fit_grp_dm.m, "Male: comorbidity - domain model", logm = T, rev=T)
rm(fit_grp_dm.m)

# =====================================
rm(impM)
# ================================= FEMALES ================================== #

impF <- miceadds::subset_datlist(imp, subset = imp$data$sex == 2,  toclass = 'mids') 

# Predict group probabilities at 13 using covariates only 
fit_grps_cov.f <- with(impF,  nnet::multinom(risk_groups_50 ~ 
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov.f <- pool_fit(fit_grps_cov.f, "Female: comorbidity - base model", logm =T, rev =T)
rm(fit_grps_cov.f)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.f <- with(impF,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre.f <- pool_fit(fit_grps_pre.f, "Female: comorbidity - pren model", logm =T, rev =T)
rm(fit_grps_pre.f)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.f <- with(impF,  nnet::multinom(risk_groups_50 ~ postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos.f <- pool_fit(fit_grps_pos.f, "Female: comorbidity - post model", logm =T, rev =T)
rm(fit_grps_pos.f)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.f <- with(impF,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z + postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add.f <- pool_fit(fit_grps_add.f, "Female: comorbidity - pren+post model",logm =T, rev =T)
rm(fit_grps_add.f)
# Predict group probabilities at 13 using prenatal by postnatal model
fit_grps_int.f <- with(impF,  nnet::multinom(risk_groups_50 ~ prenatal_stress_z * postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_int.f <- pool_fit(fit_grps_int.f, "Female: comorbidity - pren*post model",logm =T, rev =T)
rm(fit_grps_int.f)

# combine them into a single frame
comorb.f <- rbind(grps_cov.f,grps_pre.f,grps_pos.f,grps_add.f,grps_int.f)
rm(grps_cov.f,grps_pre.f,grps_pos.f,grps_add.f,grps_int.f)

# ============================================================================ #
# ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 

# Risk groups
fit_grp_dm.f <- with(impF, nnet::multinom(risk_groups_50 ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                                            post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                            age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
dm_grp.f <- pool_fit(fit_grp_dm.f, "Female: comorbidity - domain model", logm = T, rev=T)
rm(fit_grp_dm.f)

comorb.dm <- rbind(dm_grp, dm_grp.m, dm_grp.f)

rm(impF)

# ============================================================================ #
# ============================================================================ #
# Save results into two excel files

modls <- list("comorb" = comorb,"comorb_m" = comorb.m, "comorb_f" = comorb.f, 
              "comorb_dm" = comorb.dm, "comorb_summary" = grp_pooled)

openxlsx::write.xlsx(modls, file = file.path(dirname, paste0(Sys.Date(), "_ctff50_Results.xlsx")), overwrite = T)
