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
imp_samp <- readRDS(file.path(pathtoresults, 'imputation_list_sample.rds'))
imp_full <- readRDS(file.path(pathtoresults, 'imputation_list_full.rds'))

full <- complete(imp_full, 0) 
sample <- complete(imp_samp, 0) 

################################################################################

# Flowchart
fc <- capture.output(flowchart(full, return_selected_sample = T))
fc <- fc[1:(which(fc =='$final_sample')+2)]
fcm <- as.data.frame(t(matrix(unlist(fc), ncol = 9)[1:2, ])) 
fcm <- data.frame(fcm[,-1], row.names = fcm[,1])
names(fcm) = 'N'
fcm$N <- as.numeric(sub("\\[1]", "", fcm$N))

# Sample summary
s <- summdf(summary(sample))

# Stack imputed datasets in long format, exclude the original data
impdat <- complete(imp_samp, action="long", include = FALSE)

cont <- impdat[, -c(which(colnames(impdat) %in% c('sex', "ethnicity", risk_grps)))]
cate <- impdat[, c(".imp", "ethnicity", risk_grps[1:3])]

# compute mean and standard deviation in each imputed dataset
pool_mean <- with(cont, by(cont, .imp, function(x) summary(x)))
pool_numb <- with(cate, by(cate, .imp, function(x) summary(x)))

num_pool <- lapply(pool_mean, function(m) matrix(as.numeric(sapply(strsplit(m, ":"), "[[", 2)), nrow = dim(m)[1], ncol=dim(m)[2]))
pool_mean <- Reduce("+",num_pool)/length(num_pool)
colnames(pool_mean) <- colnames(cont)
rownames(pool_mean) <- c( 'Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.','Max.')

count_pool <- lapply(pool_numb, function(m) matrix(as.numeric(sapply(strsplit(m, ":"), "[[", 2)), nrow = dim(m)[1], ncol=dim(m)[2]))


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
