# Load required packages
utilis <- c('nnet', 'CMAverse', 'mice', 'miceadds','robustbase', 'openxlsx','dplyr')
lapply(utilis, require, character.only = T);

# ============================================================================ #
# Load mids object
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_sample.rds' file 
  dirname <- dirname(imp_path)
  imp <- readRDS(imp_path)
}

# ============================================================================ #
# Scale maternal BMI for additional followup (but not doing that anymore)
# datl <- mids2datlist(imp)
# datl <- miceadds::scale_datlist(datl, orig_var = 'm_bmi_before_pregnancy', trafo_var = 'm_bmi_before_pregnancy_z')
# imp <- datlist2mids(datl)

# Construct recoding for interaction analyses 
long.impdata <- complete(imp, 'long', include = TRUE) %>%
  mutate(ethnicity_rec = if_else(ethnicity == 1, 0, 1)) %>%
  mutate(sex_rec = if_else(sex == 1, 0, 1)) #  %>%
  # mutate(m_bmi_before_pregnancy_z = scale(m_bmi_before_pregnancy))
# Convert back to mids
imp <- as.mids(long.impdata)

# ============================================================================ #
# Fitting and pooling functions 
pool_mod <- function(outc, exp, data=imp, subs='') {
  
  # Define covariate set
  if (subs %in% c('(Male)','(Female)')) { add_sex = ''} else { add_sex = '+ sex_rec' }
  covs = paste(add_sex, '+ age_child + ethnicity + m_bmi_before_pregnancy + m_smoking + m_drinking')
  
  if (outc=='risk_groups_perc') {
    fit <- with(data, nnet::multinom(as.formula(paste(outc,'~',exp,covs))))
    p_fit <- mice::pool(fit) # pool results 
    mod <- summary(p_fit) #, 'all', conf.int = 0.95) # extract relevant information
    mod$OR  <- exp(mod$estimate)
    mod$lci <- exp(mod$estimate - 1.96*mod$std.error)
    mod$uci <- exp(mod$estimate + 1.96*mod$std.error)
    # add a column for AIC
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) 
    # Round everything
    mod[,-c(1,2)] <-round(mod[,-c(1,2)], 3)
    levels(mod$y.level) <- c("Internalizing", "Adiposity", "Comorbidity") # make group comparisons easier to read
  } else { 
    fit <- with(data, lm(as.formula(paste(outc,'~',exp,covs)))) 
    p_fit <- mice::pool(fit) # pool results 
    mod <- summary(p_fit) #, 'all', conf.int = 0.95) # extract relevant information
    mod$lci <- (mod$estimate - 1.96*mod$std.error)
    mod$uci <- (mod$estimate + 1.96*mod$std.error)
    # Add R squarred
    mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(pool.r.squared(fit, adjusted = T)[1], rep(NA, nrow(mod)-1)) # adjusted R2
    # Round everything
    mod[,-1] <-round(mod[,-1], 3)
  }
  
  # Add confidence intervals
  # names(mod)[which(names(mod) %in% c('2.5 %','97.5 %'))] <- c('lci','uci')
  
  exp_name <-  ifelse(exp=='','base',ifelse(exp=='prenatal_stress_z', 'pren',
                                     ifelse(exp=='postnatal_stress_z','post',
                                     ifelse(exp=='prenatal_stress_z + postnatal_stress_z','pren+post',
                                     ifelse(exp=='prenatal_stress_z * postnatal_stress_z','pren*post',
                                     ifelse(grep('pre_life_events',exp)==1,'domain',
                                     ifelse(grep('* sex',exp)==1,'sex.int',
                                     ifelse(grep('* ethnicity',exp)==1,'ethn.int','other'))))))))
  out_name <- ifelse(outc=='intern_score_13_z', 'Internalizing', 
              ifelse(outc=='tot_fat_percent_13_z', 'Adiposity', 
              ifelse(outc=='risk_groups_perc', 'Comorbidity', 
              ifelse(outc=='andr_fat_mass_13_z', 'Android fat', 'Other'))))
  
  mod_name <- paste(subs,out_name,'-',exp_name,'model')
  
  mod <- cbind(data.frame("model" = rep(mod_name, nrow(mod))), mod)
  # And one space in between models
  mod[nrow(mod)+1,] <- NA
  
  return(mod)
}
# ------------------------------------------------------------------------------
exp_mod <- function(outc, exps = c('','prenatal_stress_z','postnatal_stress_z',
                                      'prenatal_stress_z + postnatal_stress_z', 
                                      'prenatal_stress_z * postnatal_stress_z'), dat=imp, subs='') {
  mod <- data.frame()
  for (exp in exps) { 
    p <- pool_mod(outc, exp, dat, subs=subs)
    mod <- rbind(mod, p) 
  }
  # Calculate FDR per predictor
  mod[,'FDR'] <- NA
  terms <- as.vector(unique(mod[,'term'])) # get predictors
  terms <- terms[!is.na(terms)] # clean NAs
  for (t in terms) {
    fdrs <- p.adjust(mod[!is.na(mod[,'term']) & mod[,'term']==t, 'p.value'], method = 'fdr')
    # cat(t,'\t', length(fdrs),'\n')
    mod[!is.na(mod[,'term']) & mod[,'term']==t, 'FDR'] <- fdrs
  }
  mod[,'FDR'] <- round(mod[,'FDR'], 3)
  # add a column to highlight significant terms
  mod[,'sign_raw'] <- ifelse(mod[,'p.value'] < 0.05, '*', '') 
  mod[,'sign_fdr'] <- ifelse(mod[,'FDR'] < 0.05, '*', '')
  
  return(mod)
}
# ------------------------------------------------------------------------------
stack_set <- function(m1,m2,m3=NULL,m4=NULL,m5=NULL) {
  stack <- rbind(m1,m2,m3,m4,m5); rm(m1,m2,m3,m4,m5)
  return(stack)
}
  
# ============================================================================ #
# MAIN ANALYSIS: STEPWISE APPROACH
# ============================================================================ #
int <- exp_mod('intern_score_13_z')
fat <- exp_mod('tot_fat_percent_13_z')
main_s <- stack_set(int, fat)
main_c <- exp_mod('risk_groups_perc')

# ============================================================================ #
#  ADDITIONAL (2): STRESS DOMAIN CONTRIBUTION 
# ============================================================================ #
dm_predictors <- c('pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk + post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization')

int_dm <- exp_mod('intern_score_13_z', dm_predictors)
fat_dm <- exp_mod('tot_fat_percent_13_z', dm_predictors)
doms_s <- stack_set(int_dm, fat_dm)
doms_c <- exp_mod('risk_groups_perc', dm_predictors)

# ============================================================================ #
# SENSITIVITY (1): ANDROID FAT MASS 
# ============================================================================ #
fat_and <- exp_mod('andr_fat_mass_13_z')

# ============================================================================ #
#  ADDITIONAL (1): STRATIFY BY SEX
# ============================================================================ #
# test interaction
m_0 <- c('prenatal_stress_z * sex_rec','postnatal_stress_z * sex_rec') # normal coding male == 0

int_m <- exp_mod('intern_score_13_z', m_0)
fat_m <- exp_mod('tot_fat_percent_13_z', m_0)
sexi_s <- stack_set(int_m, fat_m)

sexi_c <- exp_mod('risk_groups_perc', m_0)

# further stratify the models ##################################################
impM <- miceadds::subset_datlist(imp, subset = imp$data$sex == 1,  toclass = 'mids') 

int_m <- exp_mod('intern_score_13_z', dat=impM, subs='(Male)')
fat_m <- exp_mod('tot_fat_percent_13_z', dat=impM, subs='(Male)')
idm_m <- exp_mod('intern_score_13_z', dm_predictors, dat=impM, subs='(Male)')
fdm_m <- exp_mod('tot_fat_percent_13_z', dm_predictors, dat=impM, subs='(Male)')
anf_m <- exp_mod('andr_fat_mass_13_z', dat=impM, subs='(Male)')
male_s <- stack_set(int_m, fat_m, idm_m, fdm_m, anf_m)
com_m <- exp_mod('risk_groups_perc', dat=impM, subs='(Male)')
cdm_m <- exp_mod('risk_groups_perc', dm_predictors, dat=impM, subs='(Male)')
male_c <- stack_set(com_m,cdm_m)
rm(impM)
################################################################################
impF <- miceadds::subset_datlist(imp, subset = imp$data$sex == 2,  toclass = 'mids') 

int_f <- exp_mod('intern_score_13_z', dat=impF, subs='(Female)')
fat_f <- exp_mod('tot_fat_percent_13_z', dat=impF, subs='(Female)')
idm_f <- exp_mod('intern_score_13_z', dm_predictors, dat=impF, subs='(Female)')
fdm_f <- exp_mod('tot_fat_percent_13_z', dm_predictors, dat=impF, subs='(Female)')
anf_f <- exp_mod('andr_fat_mass_13_z', dat=impF, subs='(Female)')
fema_s <- stack_set(int_f, fat_f, idm_f, fdm_f, anf_f)
com_f <- exp_mod('risk_groups_perc', dat=impF, subs='(Female)')
cdm_f <- exp_mod('risk_groups_perc', dm_predictors, dat=impF, subs='(Female)')
fema_c <- stack_set(com_f,cdm_f)
rm(impF)

# ============================================================================ #
#  SENSITIVITY (2): RESPONDENTS ONLY
# ============================================================================ #
impR <- miceadds::subset_datlist(imp, subset = (!is.na(imp$data$intern_score_13_z) & !is.na(imp$data$tot_fat_percent_13_z)),  
                                 toclass = 'mids') 

int_r <- exp_mod('intern_score_13_z', dat=impR, subs='(Resp)')
fat_r <- exp_mod('tot_fat_percent_13_z', dat=impR, subs='(Resp)')
resp_s <- stack_set(int_r, fat_r)
resp_c <- exp_mod('risk_groups_perc', dat=impR, subs='(Resp)')
rm(impR)

# ============================================================================ #
#  ADDITIONAL FOLLOWUP: ETHNICITY INTERACTION
# ============================================================================ #
e_0 <- c('prenatal_stress_z * ethnicity','postnatal_stress_z * ethnicity') # normal coding
e_1 <- c('prenatal_stress_z * ethnicity_rec','postnatal_stress_z * ethnicity_rec') # reversed coding

int_e1 <- exp_mod('intern_score_13_z', e_0)
int_er <- exp_mod('intern_score_13_z', e_1)

fat_e1 <- exp_mod('tot_fat_percent_13_z', e_0)
fat_er <- exp_mod('tot_fat_percent_13_z', e_1)

grp_e1 <- exp_mod('risk_groups_perc', e_0)
grp_er <- exp_mod('risk_groups_perc', e_1)

ethn_s <- stack_set(int_e1, fat_e1, int_er, fat_er)
ethn_c <- stack_set(grp_e1, grp_er)

# ============================================================================ #
##################### MAIN (2): CAUSAL MEDIATION ANALYSIS ######################
# ============================================================================ #

# Transform into a datalist
dat <- miceadds::mids2datlist(imp)
rm(imp)

# ------------------------------------------------------------------------------
cma_loop <- function(outcome, datalist, ymod = "linear", allow_interaction = F, n_imputations = 30) {
  # initiate lists where model fits per imputation will be stored
  yregr_fit <- list(); mregr_fit <- list(); decom_fit <- list()
  drink_fit <- list(); smoke_fit <- list()
  # loop through imputed datasets
  for (i in 1:n_imputations) {
    message("\nImputation number: ", i) # update
    est <- CMAverse::cmest(data = dat[[i]],
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
fit_int_cma  <- cma_loop('intern_score_13_z', datl)
fit_fat_cma  <- cma_loop('tot_fat_percent_13_z', datl)
intern.cma   <- pool_cma_decomp(fit_int_cma)
fatmas.cma   <- pool_cma_decomp(fit_fat_cma)
mediation <- stack_set(intern.cma, fatmas.cma)

yreg_int_cma <- pool_cma_regs(fit_int_cma@yreg, "Inernalizing - yreg")
yreg_fat_cma <- pool_cma_regs(fit_fat_cma@yreg, "Fat mass % - yreg")
rm(fit_fat_cma)
medreg_cma   <- pool_cma_regs(fit_int_cma@mreg, "Mediator regression")
drink_cma    <- pool_cma_regs(fit_int_cma@c_dr, "Drink regression")
smoke_cma    <- pool_cma_regs(fit_int_cma@c_sm, "Smoke regression")
rm(fit_int_cma)

med_other <- stack_set(yreg_int_cma, yreg_fat_cma, medreg_cma, drink_cma, smoke_cma)
# ============================================================================ #
# maternal BMI
# bmi1 <- exp_mod('intern_score_13_z', c('prenatal_stress_z*m_bmi_before_pregnancy_z','postnatal_stress_z*m_bmi_before_pregnancy_z'))
# bmi2 <- exp_mod('tot_fat_percent_13_z', c('prenatal_stress_z*m_bmi_before_pregnancy_z','postnatal_stress_z*m_bmi_before_pregnancy_z'))
# bmi3 <- exp_mod('risk_groups_perc', c('prenatal_stress_z*m_bmi_before_pregnancy_z','postnatal_stress_z*m_bmi_before_pregnancy_z'))
# mbmi <- rbind(bmi1, bmi2)
# ============================================================================ #

modls <- list("main_single" = main_s, "main_comorb" = main_c,
              "domains_single" = doms_s, "domains_comorb" = doms_c, 
              "fat_android" = fat_and, 
              "sex.int_single" = sexi_s, "sex.int_comorb" = sexi_c,
              "male_single" = male_s, "male_comorb" = male_c,
              "female_single" = fema_s, "female_comorb" = fema_c,
              "resp_single" = resp_s, "resp_comorb" = resp_c,
              "eth.int_single" = ethn_s, "eth.int_comorb" = ethn_c, 
              "mediation" = mediation, "med_other" = med_other)

openxlsx::write.xlsx(modls, file = file.path(dirname, paste0(Sys.Date(), "_GenR_Results.xlsx")), overwrite = T)
