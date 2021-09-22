require(openxlsx)

source('0-Functions.R') # where flowchart(), select_sibling() and the domains are defined

summdf <- function(object) {
  # take summary object, clean the strings and note them as row.names, return a data.frame
  m <- apply(object, 2, function(y) as.numeric(sub(".*:", "", y))) 
  m <- as.data.frame(m, row.names = c('Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.', 'NAs'))
  m[8,] <- apply(sample, 2, sd, na.rm = T)
  row.names(m)[8] <- 'SD'
  return(m[, -1])
}

# Load datasets
ELSPCM <- readRDS(paste0(pathtodata, 'ELSPCM_dataset.rds'))

################################################################################
# For the sake of time efficiency (and my mental health) let's select only those 
# variables that are needed for imputation and subsequent sample selection. 
# I will use the variable names defined in 0-Functions.R. 
# Once I am at it, I also order them by domain. This is important because mice is sensitive
# to the order of the variables in the set (even though this may be a version-specific issue)

ELSPCM_essentials <- ELSPCM[, c('IDC', 
                                # all variables for prenatal risk
                                pre_LE, pre_CR, pre_PR, pre_IR,
                                # all variables for postnatal risk
                                post_LE, post_CR, post_PR, post_IR, post_DV,
                                # all domain scores
                                'pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk', 
                                'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization',
                                # cumulative prenatal and postnatal stress exposure
                                'prenatal_stress', 'postnatal_stress',
                                # outcome variables and covariates + additional auxiliary variables for imputation
                                outcomes, covars, auxil, exclusion_criteria)]

sample <- flowchart(ELSPCM_essentials, return_selected_sample = T)
# ============================================================================ #
# Flowchart
fc <- capture.output(flowchart(ELSPCM_essentials))[1:(which(fc=='$final_sample')+2)]
fcm <- as.data.frame(t(matrix(unlist(fc), ncol = 13)[1:2, ])) 
fcm <- data.frame(fcm[,-1], row.names = fcm[,1])
names(fcm) = 'N'
fcm$N <- as.numeric(sub("\\[1]", "", fcm$N))

# Sample summary
s <- summdf(summary(sample))

# Group specific summary
bys <- by(sample, sample$risk_groups, summary)

ht <- summdf(bys[["healthy"]])
it <- summdf(bys[["internalizing_only"]])
ft <- summdf(bys[["cardiometabolic_only"]])
mm <- summdf(bys[["multimorbid"]])

# sex and ethnicity
bysex <- by(sample, sample$sex, summary)

boys <- summdf(bysex[[1]])
girl <- summdf(bysex[[2]])

# Correlation matrix
cors <- as.data.frame(cor(sample[, !names(sample) %in% c('IDC', 'twin', 'mother', 
                                                         'sex', 'ethnicity', 
                                                         'risk_groups', 'risk_groups_rec')], 
            use = 'pairwise.complete.obs'))

################################################################################

# Export the outputs of summary statistics into an xlsx file with one model per sheet

stats <- list("flowchart" = fcm, "summ_full" = s, "corr_mat" = cors,
              "summ_health" = ht, "summ_intern" = it, "summ_fatmas" = ft, "summ_multim" = mm, 
              "summ_boys" = boys, "summ_girls" = girl)

openxlsx::write.xlsx(stats, file = paste0(pathtodata, "Descriptives.xlsx"), 
                     row.names = T, overwrite = T)
