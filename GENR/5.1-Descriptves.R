
source('0-Functions.R') # where flowchart(), select_sibling() and the variables are defined

summdf <- function(object) {
  # take summary object, clean the strings and note them as row.names, return a data.frame
  m <- apply(object, 2, function(y) as.numeric(sub(".*:", "", y))) 
  m <- as.data.frame(m, row.names = c('Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.', 'NAs'))
  m[8,] <- apply(sample, 2, sd, na.rm = T)
  row.names(m)[8] <- 'SD'
  return(m[, -1])
}

pool_descriptives <- function(implist, column_names, categorical = c()) {
  num_pool <- lapply(implist, function(m) matrix(as.numeric(sapply(strsplit(m, ":"), "[[", 2)), nrow = dim(m)[1], ncol=dim(m)[2]))
  pool_mean <- Reduce("+",num_pool)/length(num_pool)
  colnames(pool_mean) <- column_names
  if (length(categorical)<1) {
    rownames(pool_mean) <- c( 'Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.','Max.')
  } else { rownames(pool_mean) <- categorical }
  
  return(data.frame(pool_mean))
}

# Load datasets
imp_samp <- readRDS(file.path(pathtoresults, 'imputation_list_sample.rds'))
imp_full <- readRDS(file.path(pathtoresults, 'imputation_list_full.rds'))

# Extract the original set (with NAs)
full <- complete(imp_full, 0) 
sample <- complete(imp_samp, 0) 
# Stack imputed datasets in long format, excluding the original data
impdat <- complete(imp_samp, action="long", include = FALSE)

################################################################################

# Flowchart capture the output of the flowchart function defined in 0-Functions
fc <- capture.output(flowchart(full, return_selected_sample = T))
# select on the relevant output are reshape it in a more readable format
fcm <- as.data.frame(t(matrix(unlist(fc[1:(which(fc =='$final_sample')+2)]), ncol = 9)[1:2, ])) 
fcm <- data.frame(fcm[,-1], row.names = fcm[,1])
names(fcm) = 'N'
fcm$N <- as.numeric(sub("\\[1]", "", fcm$N))

# ------------------------------------------------------------------------------
# Sample summary (before imputation!)
s <- summdf(summary(sample))
# for (n in summary(as.factor(sample$m_smoking))) { cat(paste0(n, " (", round((n/4268)*100), "%) \n"))}
# for (n in summary(as.factor(sample$m_drinking))) { cat(paste0(n, " (", round((n/4268)*100), "%) \n"))}

# ------------------------------------------------------------------------------
# Sample summary (after imputation!)
# split the sample in only continuous and other categorical variables 
cont <- impdat[, -c(which(colnames(impdat) %in% c('sex', "ethnicity", risk_grps)))]
cate <- impdat[, c(".imp", "ethnicity", risk_grps[1:3])]
drsm <- impdat[,  c(".imp", "m_drinking", "m_smoking")]
# compute mean and standard deviation in each imputed dataset, dividing variables into
# continuous and categorical and slitting categorical further according to the number of categories
cont_summary <- with(cont, by(cont, .imp, function(x) summary(x[, -c(1, 2)]))) # exclude .imp and .id cols
grps_summary <- with(cate, by(cate, .imp, function(x) summary(x[, -c(1, 2)]))) # exclude .imp and ethnicity
ethn_summary <- with(cate, by(cate, .imp, function(x) summary(x[, 2]))) # select ethnicity only
drnk_summary <- with(drsm, by(drsm, .imp, function(x) summary(as.factor(x[, 2]))))
smok_summary <- with(drsm, by(drsm, .imp, function(x) summary(as.factor(x[, 3]))))

# Pool descriptive s
cont_pooled <- pool_descriptives(cont_summary, colnames(cont[-c(1,2)]))
grps_pooled <- pool_descriptives(grps_summary, risk_grps[1:3], categorical = c( 'healthy', 'internalizing_only', 'cardiometabolic_only', 'multimorbid'))
ethn_pooled <- data.frame(Reduce("+",ethn_summary)/length(ethn_summary))
names(ethn_pooled) <- "ethnicity"
drnk_pooled <- data.frame(Reduce("+",drnk_summary)/length(drnk_summary))
smok_pooled <- data.frame(Reduce("+",smok_summary)/length(smok_summary))

# ------------------------------------------------------------------------------
notcat <- names(sample)[names(sample) %notin% c('IDC', 'twin', 'mother', 'sex', 'ethnicity', risk_grps)]
# Correlation matrix in the original set
cors <- as.data.frame(cor(sample[, notcat], use = 'pairwise.complete.obs'))
#  Correlation matrix in the imputed set
cors_imp <- miceadds::micombine.cor(mi.res = sample, variables = notcat) 

# ------------------------------------------------------------------------------
# Group specific summary (before imputation only)
bys <- by(sample, sample$risk_groups_perc, summary)

ht <- summdf(bys[["healthy"]])
it <- summdf(bys[["internalizing_only"]])
ft <- summdf(bys[["cardiometabolic_only"]])
mm <- summdf(bys[["multimorbid"]])

# sex 
bysex <- by(sample, sample$sex, summary)

boys <- summdf(bysex[[1]])
girl <- summdf(bysex[[2]])

################################################################################

# Export the outputs of summary statistics into an xlsx file with one model per sheet
stats <- list("flowchart" = fcm, "summ_orig" = s, "summ_imp" = cont_pooled, 
              "summ_imp_gr" = grps_pooled, "summ_imp_eth" = ethn_pooled, "corr_mat" = cors, "corr_imp" = cors_imp,
              "summ_health" = ht, "summ_intern" = it, "summ_fatmas" = ft, "summ_multim" = mm, 
              "summ_boys" = boys, "summ_girls" = girl)

openxlsx::write.xlsx(stats, file = file.path(pathtoresults, "Descriptives.xlsx"), 
                     row.names = T, overwrite = T)
