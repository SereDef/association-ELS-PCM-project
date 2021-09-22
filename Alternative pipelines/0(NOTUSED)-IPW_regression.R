# Define a function that created IP weights and runs regression
reg_long <- function(outcome, datalist, n_imputations = 30) {
  modform <- paste(outcome, '~ Stress * Period + (1 | IDC)')
  fitlist <- list() 
  # loop through imputed datasets
  for (i in 1:n_imputations) {
  # transform dataset 
  dl <- reshape(datalist[[i]], idvar = "IDC", timevar = "Period", direction="long", 
                varying = c('prenatal_stress_z', 'postnatal_stress_z'),
                v.names = 'Stress')
  
  
  fit <- lme4::lmer(modform, data = dl) 
  
  fitlist[[i]] <- fit
  }
  # Transform list into a mira object in order to use mice::pool() function
  #fitmira <- mice::as.mira(fitlist) 
  #pooled_fit <- mice::pool(fitmira) # pool results 
  
  return(fitlist)
}

is <- reg_long('intern_score_z', dat)




# Define a function that applies the IP weighted regression to each imputed dataset
# and pools estimates across them
pool_fit <- function(outcome, datalist, n_imputations = 30) {
  # initiate list where model fits will be stored
  fitlist <- list() 
  # loop through imputed datasets
  for (i in 1:n_imputations) {
    # Fit IP weighted regression model for outcome & model of choice
    wstb_pre <- ipw::ipwpoint( exposure =  prenatal_stress_z, # IP weights for prenatal stress
       family = "gaussian",
       numerator = ~ 1,
       denominator = ~ sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, 
       data = datalist[[i]])
    wstb_pos <- ipw::ipwpoint( exposure =  postnatal_stress_z, # IP weights for postnatal stress
       family = "gaussian",
       numerator = ~ prenatal_stress_z,
       denominator = ~ prenatal_stress_z + sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, 
       data = datalist[[i]])
    # Formulate the model
    modform <- paste(outcome, '~ prenatal_stress_z + postnatal_stress_z') 
    # Fit the regression
    if (outcome == 'risk_groups_rec') { 
      fit <- nnet::multinom(as.formula(modform), data = datalist[[i]], w = (wstb_pre$ipw.weights * wstb_pos$ipw.weights))
    } else { 
      fit <- lm(as.formula(modform), data = datalist[[i]], w = (wstb_pre$ipw.weights * wstb_pos$ipw.weights)) }
    # Store results to the list 
    fitlist[[i]] <- fit 
  }
  # Transform list into a mira object in order to use mice::pool() function
  fitmira <- mice::as.mira(fitlist) 
  pooled_fit <- mice::pool(fitmira) # pool results 
  
  return(pooled_fit)
}

# Define a function that pools additional info and adds them to the output
modeltable <- function(fit, logm = FALSE) {
  p_fit <- mice::pool(fit) # pool results 
  mod <- summary(p_fit) # extract relevant information
  mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
  if (logm == FALSE) {
    mod$rsq <- c(mice::pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(mice::pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)) # adjusted R2
  } else {
    # levels(mod$y.level) <- c("intern", "fatmas", "multim") # make group comparisons easier to read
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
  }
  return(mod)
}

# ============================================================================ #
# Load mids object
if (exists("imp") == F) { 
  imp_path <- file.choose() # choose the 'ELSPCM_imputed.rds' file 
  imp <- readRDS(imp_path)
}
# Transform into a datalist
dat <- miceadds::mids2datlist(imp)

# ============================================================================ #
#long trajectory regression 

# fit_is <- pool_fit('intern_score_z', dat)
# mod1 <- summary(fit_is); mod1$sign <- ifelse(mod1$p.value < 0.05, '*', '') # add a column to highlight significant terms
# mod1
# 
# fit_ft <- pool_fit('fat_mass_z', dat)
# mod2 <- summary(fit_ft); mod2$sign <- ifelse(mod2$p.value < 0.05, '*', '') # add a column to highlight significant terms
# mod2

# dl <- reshape(dat[[30]], idvar = "IDC", timevar = "period", direction = "long", 
#               varying = c('prenatal_stress_z', 'postnatal_stress_z'),
#               v.names = 'stress')
# dl$period <- as.factor(dl$period)
# ggplot2::qplot(x = stress, y = intern_score_z, data = dl, color = period) +
#   ggplot2::geom_smooth(method = "lm") 
# ggplot2::qplot(x = stress, y = fat_mass_z, data = dl, color = period) +
#   ggplot2::geom_smooth(method = "lm") 

# ============================================================================ #

# Fit IP weighted model for internalizing 
fit_is <- pool_fit('intern_score_z', dat)
mod1 <- summary(fit_is); mod1$sign <- ifelse(mod1$p.value < 0.05, '*', '') # add a column to highlight significant terms
mod1

# Fit IP weighted  model for fat mass
fit_ft <- pool_fit('fat_mass_z', dat)
mod2 <- summary(fit_ft); mod2$sign <- ifelse(mod2$p.value < 0.05, '*', '') # add a column to highlight significant terms
mod2

fit_gr <- pool_fit('risk_groups_rec', dat)
mod3 <- summary(fit_gr); mod3$sign <- ifelse(mod3$p.value < 0.05, '*', '') # add a column to highlight significant terms
mod3

# -----------------------------------------------
# Compare model with covariates for internalizing
fit_is_full <- with(imp, lm(intern_score_z ~ prenatal_stress_z + postnatal_stress_z + sex + age_child + 
                              ethnicity +  m_bmi_berore_pregnancy + m_smoking + m_drinking))
mod1.1 <- modeltable(fit_is_full)

# Compare model with covariates for fat mass
fit_fm_full <- with(imp, lm(fat_mass_z ~ prenatal_stress_z + postnatal_stress_z + sex + age_child + 
                              ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking))
mod2.1 <- modeltable(fit_fm_full)

fit_gr_full <- with(imp, { nnet::multinom(risk_groups ~ prenatal_stress_z + postnatal_stress_z +
                              sex + age_child + ethnicity + m_bmi_berore_pregnancy + m_smoking + m_drinking, model = T) });
mod3.1 <- modeltable(fit_gr_full, logm = T)

# ============================================================================ #
modls <- list("1.intern_ipw" = mod1, "1.intern_ful" = mod1.1,
              "2.fatmas_ipw" = mod2, "2.fatmas_ful" = mod2.1, 
              "3.groups_ipw" = mod3, "3.groups_ful" = mod3.1) 

openxlsx::write.xlsx(modls, file = paste0(dirname(imp_path), "/Results_ipw.xlsx"), overwrite = T)

# Try out longitudinal modeling (IGNORE: EQUIVALENT TO PREVIOUS FOR 2 TIMEPOINTS)
# fit_is_long <- pool_fit('intern_score_z', dat, long = T)
# mod <- summary(fit_is_long); mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
# mod
# fit_ft_long <- pool_fit('fat_mass_z', dat, long = T)
# mod <- summary(fit_ft_long); mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
# mod

# ============================================================================ #

est <- CMAverse::cmest(data = dat[[i]], model = "gformula", outcome = outcome, 
        exposure = "prenatal_stress_z",  mediator = c("postnatal_stress_z"), 
        basec = c("sex", "age_child", "ethnicity", "m_bmi_berore_pregnancy"), 
        postc = c('m_drinking', "m_smoking"), postcreg = list('linear', 'linear'),
        EMint = TRUE,
        mreg = list("linear"), yreg = "linear",
        astar = 0, a = 1, mval = list(0),
        estimation = "imputation", inference = "bootstrap", nboot = 1000)

s = summary(est1)

# Define a function that applies the IP weighted regression to each imputed dataset
# and pools estimates across them
pool_med <- function(outcome, datalist, n_imputations = 30, long = F) {
  fitlist <- list() # initiate list where model fits will be stored
  for (i in 1:n_imputations) { # loop through imputed datasets
    est <- CMAverse::cmest(data = datalist[[i]], model = "gformula", outcome = outcome, 
                           exposure = "prenatal_stress_z",  mediator = c("postnatal_stress_z"), 
                           basec = c("sex", "age_child", "ethnicity", "m_bmi_berore_pregnancy", "m_smoking"), 
                           postc = c('m_drinking'), postcreg = list('linear'),
                           EMint = TRUE,
                           mreg = list("linear"), yreg = "linear",
                           astar = 0, a = 1, mval = list(0),
                           estimation = "imputation", inference = "bootstrap", nboot = 20)
    # fit <- cbind( 'Estimate' = est$effect.pe, 
    #              'Std.error' = est$effect.se, 
    #                '95% CIL' = est$effect.ci.low, 
    #                '95% CIU' = est$effect.ci.high, 
    #                  'P.val' = est$effect.pval)
    fit <- est$effect.pe
    fitlist[[i]] <- fit # store results to the list 
  }
  pooled_fit <- sapply(fitlist, mean)
  # Transform list into a mira object in order to use mice::pool() function
  # fitmira <- mice::as.mira(fitlist) 
  # pooled_fit <- mice::pool(fitmira) # pool results 
  
  return(pooled_fit)
}

med <- pool_med("intern_score_z", dat)

# Define a function that pools additional info and adds them to the output
modeltable <- function(fit, logm = FALSE) {
  p_fit <- mice::pool(fit) # pool results 
  mod <- summary(p_fit) # extract relevant information
  mod$sign <- ifelse(mod$p.value < 0.05, '*', '') # add a column to highlight significant terms
  if (logm == FALSE) {
    mod$rsq <- c(pool.r.squared(fit)[1], rep(NA, nrow(mod)-1)) # add a column for R2
    mod$rsq_adj <- c(pool.r.squared(fit, adjusted = TRUE)[1], rep(NA, nrow(mod)-1)) # adjusted R2
  } else {
    levels(mod$y.level) <- c("intern", "fatmas", "multim") # make group comparisons easier to read
    mod$AIC <- c(mean(p_fit$glanced$AIC), rep(NA, nrow(mod)-1)) # add a column for AIC
  }
  return(mod)


# ============================================================================ #
# Try out Parametric G-formula for rep exposures with dataset 1
dl <- reshape(data = dat[[1]], idvar = "IDC", timevar = "time", direction="long", 
              varying = c('prenatal_stress_z', 'postnatal_stress_z'),
              v.names = 'stress')
dl$time <- dl$time - 1 # The time index must start at 0 (indexing baseline) and increase in increments of 1.

g <- gfoRmula::gformula(obs_data = dl, id = "IDC", time_name = "time", 
                        covnames = c("m_smoking","m_drinking"), 
                        covtypes = c("categorical", "categorical"), 
                        covparams = list(covmodels = c(m_smoking, m_drinking)),
                        outcome_name = "intern_score_z",
                        outcome_type = 'continuous_eof', 
                        ymodel = as.formula('intern_score_z ~ stress + sex + age_child + ethnicity +  m_bmi_berore_pregnancy + m_smoking + m_drinking'),
                        basecovs = c('sex', 'age_child','ethnicity','m_bmi_berore_pregnancy'), 
                        seed = 310896)


# ============================================================================ #
# ============================================================================ #
# summary(wstb_pre$ipw.weights)
#   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.01638   0.71993   0.84833   1.09814   1.03503 141.00299 
# summary(wstb_post$ipw.weights)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.4169  0.9605  0.9939  1.0009  1.0271  3.4468 

# Plot inverse probability weights
# graphics.off()
# ipw::ipwplot(weights = wstb_pre$ipw.weights, logscale = FALSE,
#              main = "Stabilized weights", xlim = c(0, 142))
# ipw::ipwplot(weights = wstb_post$ipw.weights, logscale = FALSE,
#             main = "Stabilized weights", xlim = c(0.4, 3.5))

# Examine numerator and denominator models.
# summary(wstb_pre$num.mod)
# summary(wstb_pre$den.mod)
# summary(wstb_post$num.mod)
# summary(wstb_post$den.mod)

# model3 <- nnet::multinom(risk_groups ~ prenatal_stress + postnatal_stress, data = imp, w = IPW) 
# summary(model3)
