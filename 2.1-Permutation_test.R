
# Here we test whether the multimorbidity group size is bigger than what you would 
# expect by chance alone. We do it with a permutation test. We randomly permute the 
# int and fat variables (redoting above 80th percentile internalizing or fat mass)
# and use the new random variables to re-compute the multimorbid groups size. 
# To obtain a p-value we than simply count how many times a group of equal of greater
# size occurs by chance and divide that count by the number of iterations. 

data <- readRDS('/Users/Serena/Desktop/Data/ELSPCM_imputed.rds')

data$int = ifelse(data$intern_score_z > quantile(data$intern_score_z, probs = 0.8), 1, 0) # 590 risk, 2780 no risk
data$fat = ifelse(data$fat_mass_z > quantile(data$fat_mass_z, probs = 0.8), 1, 0) # 674 risk, 2696 no risk

data$groups <- ifelse(data$int == 0 & data$fat == 0, 0, 
                      ifelse(data$int == 1 & data$fat == 0, 1, 
                             ifelse(data$int == 0 & data$fat == 1, 2, 3)))

summary(as.factor(data$groups))

# define a function that permutes (randomly shuffles )
permute <- function(df) { 
  # create empty dataset for permutation
  perm <- data.frame(1:dim(df)[1])
  
  perm$new_int <- sample(df$int)
  perm$new_fat <- sample(df$fat)
  perm$new_groups <- ifelse(perm$new_int == 0 & perm$new_fat == 0, 0, 
                            ifelse(perm$new_int == 1 & perm$new_fat == 0, 1, 
                                   ifelse(perm$new_int == 0 & perm$new_fat == 1, 2, 3)))
  summary(as.factor(perm$new_groups))[4]
  }

count <- 0 
itarations <- 1000
set.seed(3100896)

for (i in 1:itarations) {
  if (permute(data) > summary(as.factor(data$groups))[4]) {
    count <- count + 1
  }
}

pval = round(count / itarations, 3)

set.seed(3100896)
################################################################################
################################################################################
################################################################################
# Create two more groups to account for the effect of tress on low fat mass

# Compute groups 
data$int = ifelse(data$intern_score_z > quantile(data$intern_score_z, probs = 0.8), 1, 0) # 590 risk, 2780 no risk
data$fat = ifelse(data$fat_mass_z > quantile(data$fat_mass_z, probs = 0.8), 1, 0) # 674 risk, 2696 no risk
data$lowfat = ifelse(data$fat_mass_z < quantile(data$fat_mass_z, probs = 0.2), 1, 0)


data$groups = rep(NA, 3371)
for (i in 1:3371) {
  if (data$int[i] == 0 & data$fat[i] == 0 & data$lowfat[i] == 0) { data$groups[i] = 0 # healthy
  } else if (data$int[i] == 0 & data$fat[i] == 0 & data$lowfat[i] == 1) { data$groups[i] = 1 # low fat only 
  } else if (data$int[i] == 1 & data$fat[i] == 0 & data$lowfat[i] == 1) {  data$groups[i] = 2 # anorexia type (high internalizing low fat)
  } else if (data$int[i] == 1 & data$fat[i] == 0 & data$lowfat[i] == 0) { data$groups[i] = 3 # High intenalizing  only
  } else if (data$int[i] == 0 & data$fat[i] == 1 & data$lowfat[i] == 0) { data$groups[i] = 4 # High fat mass only
  } else { data$groups[i] = 5 } # multimorbid
}

data$groups = as.factor(data$groups)
summary(data$groups)

# # Let's first factor that bad boy 
data$groups = factor(data$groups,
                     levels = c(0:5),
                     labels = c("healthy", "low_fat", "anorexic_type", "internalizing_only", "cardiometabolic_only", "multimorbid"))

# display groups
attach(data); plot(intern_score_z, fat_mass_z,
                     col=c("cornflowerblue","black", "red", "darkgrey", "darkgoldenrod2","chartreuse4")[groups]); detach(data)

saveRDS(data, "/Users/Serena/Desktop/attempt.rds")
