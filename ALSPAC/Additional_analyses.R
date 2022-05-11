# Load required packages
utilis <- c('nnet', 'CMAverse', 'mice', 'openxlsx')
lapply(utilis, require, character.only = T);

# ============================================================================ #
# Load mids object
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_sample.rds' file 
  dirname <- dirname(imp_path)
  imp <- readRDS(imp_path)
}

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

# Reverse the comparison group
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 covariats only 
fit_grps_cov <- with(imp,  nnet::multinom(risk_groups_perc ~
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov <- pool_fit(fit_grps_cov, "Comorbidity - base model", logm =T)
rm(fit_grps_cov)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre <- with(imp,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre <- pool_fit(fit_grps_pre, "Comorbidity - pren model", logm =T)
rm(fit_grps_pre)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos <- with(imp,  nnet::multinom(risk_groups_perc ~ postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos <- pool_fit(fit_grps_pos, "Comorbidity - post model", logm =T)
rm(fit_grps_pos)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add <- with(imp,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z + postnatal_stress_z +
                     sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add <- pool_fit(fit_grps_add, "Comorbidity - pren+post model", logm =T)
rm(fit_grps_add)

# combine them into a single frame
comorb <- rbind(grps_cov,grps_pre,grps_pos,grps_add)
rm(grps_cov,grps_pre,grps_pos,grps_add)

#  ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 
fit_grp_dm <- with(imp, nnet::multinom(risk_groups_perc ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                   post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                   sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
dm_grp <- pool_fit(fit_grp_dm, "Comorbidity - domain model", logm = T)
rm(fit_grp_dm)

####################### ADDITIONAL (1): STRATIFY BY SEX ########################
# ============================================================================ #

impM <- miceadds::subset_datlist(imp, subset = imp$data$sex == "Male",  toclass = 'mids') 
impF <- miceadds::subset_datlist(imp, subset = imp$data$sex == "Female",  toclass = 'mids')
# also create the response only dataset 
impR <- miceadds::subset_datlist(imp, subset = (!is.na(imp$data$intern_score_13_z) & !is.na(imp$data$tot_fat_percent_13_z)),  
                                 toclass = 'mids') 
rm(imp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 using covariates only 
fit_grps_cov.m <- with(impM,  nnet::multinom(risk_groups_perc ~ 
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov.m <- pool_fit(fit_grps_cov.m, "Male: comorbidity - base model", logm =T)
rm(fit_grps_cov.m)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.m <- with(impM,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre.m <- pool_fit(fit_grps_pre.m, "Male: comorbidity - pren model", logm =T)
rm(fit_grps_pre.m)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.m <- with(impM,  nnet::multinom(risk_groups_perc ~ postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos.m <- pool_fit(fit_grps_pos.m, "Male: comorbidity - post model", logm =T)
rm(fit_grps_pos.m)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.m <- with(impM,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z + postnatal_stress_z +
                                               age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add.m <- pool_fit(fit_grps_add.m, "Male: comorbidity - pren+post model",logm =T)
rm(fit_grps_add.m)

# combine them into a single frame
comorb.m <- rbind(grps_cov.m, grps_pre.m, grps_pos.m, grps_add.m)
rm(grps_cov.m, grps_pre.m, grps_pos.m, grps_add.m)

# ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 
fit_grp_dm.m <- with(impM, nnet::multinom(risk_groups_perc ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
dm_grp.m <- pool_fit(fit_grp_dm.m, "Male: comorbidity - domain model", logm = T)
rm(fit_grp_dm.m)

# =====================================
datM <- miceadds::mids2datlist(impM)
rm(impM)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 using covariates only 
fit_grps_cov.f <- with(impF,  nnet::multinom(risk_groups_perc ~ 
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov.f <- pool_fit(fit_grps_cov.f, "Female: comorbidity - base model", logm =T)
rm(fit_grps_cov.f)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.f <- with(impF,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre.f <- pool_fit(fit_grps_pre.f, "Female: comorbidity - pren model", logm =T)
rm(fit_grps_pre.f)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.f <- with(impF,  nnet::multinom(risk_groups_perc ~ postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos.f <- pool_fit(fit_grps_pos.f, "Female: comorbidity - post model", logm =T)
rm(fit_grps_pos.f)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.f <- with(impF,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z + postnatal_stress_z +
                       age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add.f <- pool_fit(fit_grps_add.f, "Female: comorbidity - pren+post model",logm =T)
rm(fit_grps_add.f)

# combine them into a single frame
comorb.f <- rbind(grps_cov.f,grps_pre.f,grps_pos.f,grps_add.f)
rm(grps_cov.f,grps_pre.f,grps_pos.f,grps_add.f)

# ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 
fit_grp_dm.f <- with(impF, nnet::multinom(risk_groups_perc ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                     post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                     age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
dm_grp.f <- pool_fit(fit_grp_dm.f, "Female: comorbidity - domain model", logm = T)
rm(fit_grp_dm.f)

# =====================================
datF <- miceadds::mids2datlist(impF)
rm(impF)

comorb.dm <- rbind(dm_grp, dm_grp.m, dm_grp.f)
rm(dm_grp, dm_grp.m, dm_grp.f)

# ============================================================================ #
####################### SENSITIVITY (2): RESPONDENTS ONLY ######################
# ============================================================================ #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict group probabilities at 13 covariats only 
fit_grps_cov.r <- with(impR,  nnet::multinom(risk_groups_perc ~
                                               sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_cov.r <- pool_fit(fit_grps_cov.r, "(Resp) Comorbidity - base model", logm =T)
rm(fit_grps_cov.r)
# Predict group probabilities at 13 using prenatal stress only 
fit_grps_pre.r <- with(impR,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z +
                                               sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pre.r <- pool_fit(fit_grps_pre.r, "(Resp) Comorbidity - pren model", logm =T)
rm(fit_grps_pre.r)
# Predict group probabilities at 13 using postnatal stress only 
fit_grps_pos.r <- with(impR,  nnet::multinom(risk_groups_perc ~ postnatal_stress_z +
                                               sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_pos.r <- pool_fit(fit_grps_pos.r, "(Resp) Comorbidity - post model", logm =T)
rm(fit_grps_pos.r)
# Predict group probabilities at 13 using an additive prenatal and postnatal model
fit_grps_add.r <- with(impR,  nnet::multinom(risk_groups_perc ~ prenatal_stress_z + postnatal_stress_z +
                                               sex + age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking, model = T))
grps_add.r <- pool_fit(fit_grps_add.r, "(Resp) Comorbidity - pren+post model", logm =T)
rm(fit_grps_add.r)

# combine them into a single frame
comorb.resp <- rbind(grps_cov.r, grps_pre.r, grps_pos.r, grps_add.r)
rm(grps_cov.r, grps_pre.r, grps_pos.r, grps_add.r)

# =================================
datR <- miceadds::mids2datlist(impR)
rm(impR)

# ============================================================================ #
################## ADDITIONAL (3): CAUSAL MEDIATION ANALYSIS ###################
# ============================================================================ #

# ------------------------------------------------------------------------------
cma_loop <- function(outcome, datalist, sex=F, allow_interaction = F, n_imputations = 30) {
  # initiate lists where model fits per imputation will be stored
  yregr_fit <- list(); mregr_fit <- list(); decom_fit <- list()
  drink_fit <- list(); smoke_fit <- list()
  if (sex==F) { baseC <- c("sex", "age_child", "ethnicity", "m_bmi_before_pregnancy")
  } else { baseC <- c("age_child", "ethnicity", "m_bmi_before_pregnancy") }
  # loop through imputed datasets
  for (i in 1:n_imputations) {
    message("\nImputation number: ", i) # update
    est <- CMAverse::cmest(data = datalist[[i]], full = T,
                           exposure = "prenatal_stress_z",  mediator = c("postnatal_stress_z"),
                           outcome = outcome,  yreg = 'linear', # outcome regression model (linear or multimod)
                           EMint = allow_interaction, # exposure-mediator interaction in yreg (i.e. prenatal*postnatal). 
                           model = "gformula", 
                           # exposure-outcome confounder(s), exposure-mediator confounder(s) and mediator-outcome confounder(s) not affected by the exposure
                           basec = baseC, 
                           # mediator-outcome confounder(s) affected by the exposure following temporal order
                           postc = c("m_drinking", "m_smoking"), # regression model for each variable in postc (gformula): 
                           postcreg = list(glm(m_drinking ~ prenatal_stress_z + ethnicity + m_bmi_before_pregnancy, family = gaussian(), data = datalist[[i]]), 
                                           glm(m_smoking ~ prenatal_stress_z + ethnicity + m_bmi_before_pregnancy, family = gaussian(), data = datalist[[i]])), 
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
pool_cma_decomp <- function(fit, dataset, n_parameters = 7) {
  m = length(fit@decomp); # number of imputations
  n = nrow(dataset[[1]]) # sample size
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
fit_int_cma.m  <- cma_loop('intern_score_13_z', datM, sex=T)
intern.cma.m   <- pool_cma_decomp(fit_int_cma.m, datM)
yreg_int_cma.m <- pool_cma_regs(fit_int_cma.m@yreg, "Inernalizing - yreg")
medreg_cma.m  <- pool_cma_regs(fit_int_cma.m@mreg, "Mediator regression")
drink_cma.m   <- pool_cma_regs(fit_int_cma.m@c_dr, "Drink regression")
smoke_cma.m   <- pool_cma_regs(fit_int_cma.m@c_sm, "Smoke regression")
rm(fit_int_cma.m)

fit_fat_cma.m  <- cma_loop('tot_fat_percent_13_z', datM, sex=T)
fatmas.cma.m   <- pool_cma_decomp(fit_fat_cma.m, datM)
yreg_fat_cma.m <- pool_cma_regs(fit_fat_cma.m@yreg, "Fat mass % - yreg")
rm(fit_fat_cma.m)

mediation.m <- rbind(yreg_int_cma.m, yreg_fat_cma.m, medreg_cma.m, drink_cma.m, smoke_cma.m)
rm(datM, yreg_int_cma.m, yreg_fat_cma.m, medreg_cma.m, drink_cma.m, smoke_cma.m)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit the CMA model to all imputed datasets
fit_int_cma.f  <- cma_loop('intern_score_13_z', datF, sex=T)
intern.cma.f   <- pool_cma_decomp(fit_int_cma.f, datF)
yreg_int_cma.f <- pool_cma_regs(fit_int_cma.f@yreg, "Inernalizing - yreg")
medreg_cma.f  <- pool_cma_regs(fit_int_cma.f@mreg, "Mediator regression")
drink_cma.f   <- pool_cma_regs(fit_int_cma.f@c_dr, "Drink regression")
smoke_cma.f   <- pool_cma_regs(fit_int_cma.f@c_sm, "Smoke regression")
rm(fit_int_cma.f)

fit_fat_cma.f  <- cma_loop('tot_fat_percent_13_z', datF, sex=T)
fatmas.cma.f   <- pool_cma_decomp(fit_fat_cma.f, datF)
yreg_fat_cma.f <- pool_cma_regs(fit_fat_cma.f@yreg, "Fat mass % - yreg")
rm(fit_fat_cma.f)

mediation.f <- rbind(yreg_int_cma.f, yreg_fat_cma.f, medreg_cma.f, drink_cma.f, smoke_cma.f)
rm(datF, yreg_int_cma.f, yreg_fat_cma.f, medreg_cma.f, drink_cma.f, smoke_cma.f)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit the CMA model to all imputed datasets (RESPONDENTS ONLY)
fit_int_cma.r  <- cma_loop('intern_score_13_z', datR)
intern.cma.r   <- pool_cma_decomp(fit_int_cma.r, datR)
yreg_int_cma.r <- pool_cma_regs(fit_int_cma.r@yreg, "Inernalizing - yreg")
medreg_cma.r  <- pool_cma_regs(fit_int_cma.r@mreg, "Mediator regression")
drink_cma.r   <- pool_cma_regs(fit_int_cma.r@c_dr, "Drink regression")
smoke_cma.r   <- pool_cma_regs(fit_int_cma.r@c_sm, "Smoke regression")
rm(fit_int_cma.r)

fit_fat_cma.r  <- cma_loop('tot_fat_percent_13_z', datR)
fatmas.cma.r   <- pool_cma_decomp(fit_fat_cma.r, datR)
yreg_fat_cma.r <- pool_cma_regs(fit_fat_cma.r@yreg, "Fat mass % - yreg")
rm(fit_fat_cma.r)

mediation.r <- rbind(yreg_int_cma.r, yreg_fat_cma.r, medreg_cma.r, drink_cma.r, smoke_cma.r)
rm(datR, yreg_int_cma.r, yreg_fat_cma.r, medreg_cma.r, drink_cma.r, smoke_cma.r)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save results
modls <- list("comorb" = comorb, "comorb_m" = comorb.m, "comorb_f" = comorb.f, "comorb_dm" = comorb.dm,"comorb_select" = comorb.resp, 
              "intern_cma_m" = intern.cma.m, "fatmas_cma_m" = fatmas.cma.m, "other_cma_m" = mediation.m,
              "intern_cma_f" = intern.cma.f, "fatmas_cma_f" = fatmas.cma.f, "other_cma_f" = mediation.f, 
              "intern_cma_r" = intern.cma.r, "fatmas_cma_r" = fatmas.cma.r, "other_cma_r" = mediation.r)

openxlsx::write.xlsx(modls, file = file.path(dirname, paste0(Sys.Date(), "_ALSPAC_addResults.xlsx")), overwrite = T)
