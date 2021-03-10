# Hola,
# The following code runs a set of analyses for the project "Early-life stress as a risk factor 
# for poor mental and physical health in children: A population-based study" looking at the 
# association between ELS and psycho-cardio-metabolic multi-morbidity at age 10.  

# All it requires is just one file: the imputed dataset we built using the "Imputation.R"
# script, that you can find in this repository: https://github.com/SereDef/cumulative-ELS-score. 

# Ok, let's get started!

# check if the path to the dataset is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

#------------------------------------------------------------------------------#
# Load the (complete) dataset
datarisk <- readRDS(paste(pathtodata, 'ELS_PCM_imputed.rds', sep = ""))
# imputation <- readRDS(paste(pathtodata,'imputation_list.rds', sep = ""))

# Compute groups 
attach(datarisk)
datarisk$int = ifelse(intern_score_z > quantile(intern_score_z, probs = 0.8), 1, 0) # 590 risk, 2780 no risk
datarisk$fat = ifelse(fat_mass_z > quantile(fat_mass_z, probs = 0.8), 1, 0) # 674 risk, 2696 no risk

datarisk$risk_groups = rep(NA, 3371)
for (i in 1:3371) {
  if (datarisk$int[i] == 0 & datarisk$fat[i] == 0) { datarisk$risk_groups[i] = 0 # healthy
  } else if (datarisk$int[i] == 1 & datarisk$fat[i] == 0) { datarisk$risk_groups[i] = 1 # depression only
  } else if (datarisk$int[i] == 0 & datarisk$fat[i] == 1) { datarisk$risk_groups[i] = 2 # fat only
  } else { datarisk$risk_groups[i] = 3 } # multimorbid
}

datarisk$risk_groups = as.factor(datarisk$risk_groups)
summary(datarisk$risk_groups)
##    0    1    2    3 
## 2267  429  513  161 

datarisk$presum <- rowSums(datarisk[c("pre_life_events","pre_contextual_risk", 
                                      "pre_personal_stress", "pre_interpersonal_stress")], na.rm = F)
datarisk$postsum <- rowSums(datarisk[c("post_life_events","post_contextual_risk", 
                                       "post_parental_risk", "post_interpersonal_risk", 
                                       "post_direct_victimization")], na.rm = F)
datarisk$elssum <- rowSums(datarisk[c("presum","postsum")],  na.rm = F)

attach(datarisk)

lin = lm(risk_groups ~ presum + postsum + presum*postsum, data = datarisk)
summary(lin)

library(MASS)
ord = polr(risk_groups ~ presum + postsum + presum*postsum, data = datarisk, Hess=TRUE)
(ctable <- coef(summary(ord)))

## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = round(p, 4)))


saveRDS(datarisk, paste(pathtodata,'attempt.rds', sep = ""))

################################################################################
################################################################################
################################################################################
library(mlogit)

groups <- readRDS(paste(pathtodata, 'attempt.rds', sep = ""))

# Let's first factor that bad boy 
groups$risk_groups = factor(groups$risk_groups, 
                            levels = c(0:3), 
                            labels = c("healthy", "internalizing_only", "cardiometabolic_only", "multimorbid"))

# NOTE: the ratio of number of cases is important! 
# It is always easy to predict the larger category because probabilistically it is more likely
# rule of thumb is no higher that 1 to 5 ratio. 

#display groups 
attach(groups); plot(intern_score_z, fat_mass_z, 
                     col=c("cornflowerblue","darkgrey", "darkgoldenrod2","chartreuse4")[risk_groups]); detach(groups)

#display groups in relation to stress 
attach(groups); plot(presum, postsum, 
                     col=c("cornflowerblue","darkgrey", "darkgoldenrod2","chartreuse4")[risk_groups]); detach(groups)

# additivity 
correl = cor(groups[, 97:105]) # select just domains
symnum(correl) # this is really cool to quickly visualize collinearity (* or B for > 9)

# create the dataset for analysis
essential = groups[, c(97:105, 127, 128, 126)]
# we have to reshape the data to a format mlogit likes (re-stacking)
longdata = mlogit.data(essential, choice = "risk_groups", shape = "wide")
# not the reshape package because the data is actually already in the wide format
# but what is going to happen is that now groups are coded in a true -false manner
# for cor each participant the row is repeated  4 times 
# this is the dummy coding we are going to use 

# lets run the model! 
model = mlogit(risk_groups ~ 1 | presum + postsum + (presum*postsum), 
               data = longdata, reflevel = "healthy")

summary(model)
## Frequencies of alternatives:choice
## healthy cardiometabolic_only   internalizing_only          multimorbid
## 0.672700             0.152226             0.127300             0.047774
# percentages indicating proportion of each category
# this is problematic beause they are unbalanced 

# nr method
# 6 iterations, 0h:0m:0s
# g'(-H)^-1g = 1.04E-05
# successive function values within tolerance limits
# 
# Coefficients :
#                                  Estimate Std. Error  z-value  Pr(>|z|)
# (Intercept):cardiometabolic_only -2.53021    0.15556 -16.2651 < 2.2e-16 ***
# (Intercept):internalizing_only   -2.72761    0.16239 -16.7964 < 2.2e-16 ***
# (Intercept):multimorbid          -5.33946    0.32415 -16.4721 < 2.2e-16 ***
# presum:cardiometabolic_only       1.44587    0.26117   5.5362 3.091e-08 ***
# presum:internalizing_only         0.25419    0.28180   0.9020  0.367043
# presum:multimorbid                1.71288    0.41268   4.1506 3.316e-05 ***
# postsum:cardiometabolic_only      1.02465    0.20317   5.0434 4.573e-07 ***
# postsum:internalizing_only        1.32527    0.19630   6.7512 1.467e-11 ***
# postsum:multimorbid               2.56026    0.30535   8.3848 < 2.2e-16 ***
# presum:cardiometabolic_only      -0.88523    0.22636  -3.9106 9.205e-05 ***
# presum:internalizing_only        -0.16391    0.20160  -0.8130  0.416196
# presum:multimorbid               -0.91979    0.28082  -3.2754  0.001055 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# is prenatal and postantal stess even useful ? 
# Log-Likelihood: -3068.5
# McFadden R^2:  0.052456 # or pseudo-R2
# Likelihood ratio test : chisq = 339.74 (p.value = < 2.22e-16) # 12 dfs 
