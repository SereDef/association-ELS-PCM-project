
# Causal Mediation Analysis Approach

# Load required packages
utilis <- c('nnet', 'CMAverse', 'mice', 'openxlsx')
lapply(utilis, require, character.only = T);

# ============================================================================ #
# Load mids object
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'imputation_list_ELSPCM.rds' file 
  imp <- readRDS(imp_path)
}

# outlier exclusion
# imp <- miceadds::subset_datlist( imp, subset = ! imp$data$IDC %in% c(5105), toclass = 'mids')

# Transform into a datalist
dat <- miceadds::mids2datlist(imp)

# ============================================================================ #
# ============================================================================ #

 pool_fit <- function(outcome, datalist, ymod = "linear", n_imputations = 30) {
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
           EMint = T, # exposure-mediator interaction in yreg (i.e. prenatal*postnatal). 
           model = "gformula", 
           # exposure-outcome confounder(s), exposure-mediator confounder(s) and mediator-outcome confounder(s) not affected by the exposure
           basec = c("sex", "age_child", "ethnicity", "m_bmi_berore_pregnancy"), 
           # mediator-outcome confounder(s) affected by the exposure following temporal order
           postc = c("m_drinking", "m_smoking"), # regression model for each variable in postc (gformula): 
           postcreg = list(glm(m_drinking ~ prenatal_stress_z + ethnicity + m_bmi_berore_pregnancy, family = gaussian(), data = dat[[i]]), 
                           glm(m_smoking ~ prenatal_stress_z + ethnicity + m_bmi_berore_pregnancy, family = gaussian(), data = dat[[i]])), 
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
# ============================================================================ #

pool_regs <- function(fit) {
  p = mice::pool(mice::as.mira(fit)) 
  s = summary(p)
  s$sign <- ifelse(s$p.value < 0.05, '*', '') # add a column to highlight significant terms
  print(s)
  return(s)
}
# ============================================================================ #
# Rubin’s Rules # https://bookdown.org/mwheymans/bookmi/rubins-rules.html

pool_decomp <- function(fit, n_parameters = 15) {
  m = length(fit@decomp); # number of imputations
  n = nrow(dat[[1]]) # sample size
  p = n_parameters
  
  # Extract estimates and SEs for all imputations
  est <- sapply(fit@decomp, "[", , "Estimate") 
  se <- lapply(fit@decomp, "[", , "Std.error") 
  
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
# ============================================================================ #
# Define a customized function that pools regression results and builds a dataframe 
# that is easy to read and save
modeltable <- function(fit, logm = FALSE, rev = FALSE) {
  p_fit <- mice::pool(fit) # pool results 
  mod <- summary(p_fit) # extract relevant information
  mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
  if (logm == FALSE) {
    mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)) # adjusted R2
  } else {
    if (rev == T) { levels(mod$y.level) <- c("M:healthy", "M:intern", "M:fatmas") 
    } else { levels(mod$y.level) <- c("H:intern", "H:fatmas", "H:multim") } # make group comparisons easier to read
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
  }
  return(mod)
}
# ============================================================================ #

fit_is <- pool_fit('intern_score_z', dat)
fit_fm <- pool_fit('fat_mass_z', dat)
fit_gr <- pool_fit('risk_groups_rec', dat, ymod = "multinomial")

y_is <- pool_regs(fit_is@yreg)
y_fm <- pool_regs(fit_fm@yreg)
y_gr <- pool_regs(fit_gr@yreg)

medr <- pool_regs(fit_is@mreg)

drink <- pool_regs(fit_is@c_dr)
smoke <- pool_regs(fit_is@c_sm)

dec_is <- pool_decomp(fit_is)
dec_fm <- pool_decomp(fit_fm)
dec_gr <- pool_decomp(fit_gr, 17)

# ------------------------------------------------------------------------------

# Finally,
# Let's have a look at the contribution of individual domains of stress to the probability 
# of belonging to each group. Again fist a minumally adjusted and then a fully adjusted 
# model is run. 

# Fully adjusted model
fit_int_dm_full <- with(imp, nnet::multinom(intern_score_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                                              post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                              sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, model = T))
dm_int <- modeltable(fit_int_dm_full, logm = T)

# Fully adjusted model
fit_ftm_dm_full <- with(imp, nnet::multinom(fat_mass_z ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                                              post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                              sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, model = T))
dm_ftm <- modeltable(fit_ftm_dm_full, logm = T)

# Fully adjusted model
fit_grp_dm_full <- with(imp, nnet::multinom(risk_groups_rec ~ pre_life_events + pre_contextual_risk +  pre_parental_risk + pre_interpersonal_risk +
                                              post_life_events + post_contextual_risk + post_parental_risk + post_interpersonal_risk + post_direct_victimization +
                                              sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, model = T))
dm_grp <- modeltable(fit_grp_dm_full, logm = T)


# ============================================================================ #
modls <- list("1.intern" = y_is, "1.1.intern_dec" = dec_is,
              "2.fatmas" = y_fm, "2.1.fatmas_dec" = dec_fm, 
              "3.groups" = y_gr, "3.1.groups_dec" = dec_gr, 
              "med_reg" = medr, "drink_reg" = drink, "smoke_reg" = smoke, 
              "4.DM_int" = dm_int, "5.DM_fat" = dm_ftm, "6.DM_grp" = dm_grp) 

openxlsx::write.xlsx(modls, file = paste0(dirname(imp_path), "/Results_cma.xlsx"), overwrite = T)
