
# Load required packages
utilis <- c('foreign', 'stats', 'mice', 'miceadds')
lapply(utilis, require, character.only = T);

# Defining the path to the data
# check if the path to the data is already in memory, otherwise ask for it. 

if (exists("pathtodata") == F) { 
  message("ATTENTION! You will be prompted with a window, please navigate to the directory
          where your input files are stored. Choose any file in the directory to continue.
          Note: the code assumes that all (raw) data is stored in ONE folder.")
  f <- file.choose() 
  pathtodata <- dirname(f) }

while (exists("pathtoresults") == F) { 
  res <- readline(prompt = "Do you want to save results in the same directory? [y/n] ")
  if (!res %in% c('y', 'yes', 'n', 'no')) {
    message('You did not answer my question. RUDE! Let us start over, shall we')
    res <- readline(prompt = "Do you want to save results in the same directory? [y/n] ")
  } else if (res == 'y' | res == 'yes') { 
    pathtoresults <- file.path(pathtodata, 'Results')
    dir.create(pathtoresults, showWarnings = FALSE)
  } else if (res == 'n' | res == 'no') { 
    pathtoresults <- readline(prompt = 'Enter the full path where results should be saved: ') }
}

# This will come in handy for exclusion
'%notin%' <- Negate('%in%')

# Organize variable names into domains to specify them later more easily
pre_LE <- c('family_member_died','friend_relative_died', 'family_member_ill_pregnancy','admitted_to_hospital', 
            'health', 'unemployed', 'work_study_problems', 'moved_house', 'blood_loss', 'examination', 
            'baby_worried', 'pregnancy_worried', 'obstetric_care', 'pregnancy_planned', 'victim_robbery')
pre_CR <- c('financial_problems', 'trouble_pay_pregnancy', 'income_reduced', 'housing_defects', 'housing_adequacy', 
            'housing_basic_living', 'm_education_pregnancy', 'p_education_pregnancy')
pre_PR <- c('m_depression_pregnancy', 'm_anxiety_pregnancy', 'm_interp_sensitivity_pregnancy', 
            'p_depression_pregnancy', 'p_anxiety_pregnancy', 'p_interp_sensitivity_pregnancy', 'm_violence_people', 
            'm_violence_property', 'm_criminal_record', 'p_criminal_record') # without age, as the same variable is already in postnatal PR score
pre_IR <- c('difficulties_contacts','difficulties_partner','difficulties_family_friend','marital_status_pregnancy',
            'divorce_pregnancy','family_support','family_acceptance','family_affection','family_acception','family_trust',
            'family_painful_feelings','family_decisions','family_conflict','family_decisions_problems',
            'family_plans','family_talk_sadness', 'family_talk_worries', 'family_size_pregnancy')
post_LE <- c('sick_or_accident','family_member_ill','smbd_important_ill','parent_died','smbd_important_died','pet_died',
             'school_workload','repeated_grade','lost_smth_important','moved','changed_school','friend_moved','fire_or_burglary')
post_CR <- c('material_deprivation','financial_difficulties','neiborhood_problems','trouble_pay_childhood','income_once',
             'income_chronic','unemployed_once','unemployed_chronic', 'm_education','p_education')
post_PR <- c('tension_at_work','m_interpersonal_sensitivity','p_interpersonal_sensitivity','m_depression','p_depression',
             'm_anxiety','p_anxiety','m_age','p_age')
post_IR <- c('marital_problems','marital_status','family_size','m_fad_5yrs','m_fad_9yrs','p_fad_9yrs','conflict_family_member',
             'conflict_smbd_else','conflict_in_family','divorce_childhood','argument_friend')
post_DV <- c('m_harsh_parent','p_harsh_parent','bullying','physical_violence','physical_threats','sexual_harrasment',
             'sexual_behavior','rumors_or_gossip')
domains <- c('pre_life_events', 'pre_contextual_risk', 'pre_parental_risk', 'pre_interpersonal_risk', 
             'post_life_events', 'post_contextual_risk', 'post_parental_risk', 'post_interpersonal_risk', 'post_direct_victimization')
outcomes_13y <- c('intern_score_13', 'total_fat_13', 'andr_fat_mass_13', 'tot_fat_percent_13')
outcomes_09y <- c('intern_score_09', 'total_fat_09', 'andr_fat_mass_09', 'tot_fat_percent_09')
risk_grps <-c("risk_groups_tot", "risk_groups_andr", "risk_groups_perc", "risk_groups_tot_REC", "risk_groups_andr_REC", "risk_groups_perc_REC")
covars   <- c('sex', 'age_child', 'm_bmi_before_pregnancy', 'm_smoking', 'm_drinking')
auxil    <- c('m_bmi_pregnancy','m_dep_cont_pregnancy', 'p_dep_cont_pregnancy', # for postnatal only 
              'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs',               # for prenatal only 
              'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_age_cont')
exclusion_criteria <- c('pre_percent_missing', 'post_percent_missing', 'twin', 'mother')


#### ------------------------------ FUNCTIONS ----------------------------- ####

# Read in the data and fix all SPSS weird missing codes into NAs
readquick <- function(filename, rootdir = pathtodata, exclude_col = "") { # only works for SPSS files
  dat <- read.spss(file.path(rootdir, filename), 
                   use.value.labels = F, to.data.frame = T)
  # Get rid of all capital letters in column names (so you don't have to worry)
  names(dat) <- tolower(names(dat))
  # Replace values of 777, 888 or 999 with NAs unless they are IDCs or IDMs 
  # If you do not want this to happen for any other column use the exclude_col argument. 
  for (i in 1:length(dat)) {
    if (colnames(dat)[i] == "idm" | colnames(dat)[i] == "idc" | colnames(dat)[i] == exclude_col) {
      dat[,i] <- dat[,i]
    } else {
      dat[,i] <- ifelse(dat[,i] == 777 | dat[,i] == 888 | dat[,i] == 999, NA, dat[,i]) }
  } 
  return(dat)
}

#-------------------------------------------------------------------------------
## Select only one sibling (based on data availability or, if missing are the same, randomly).
select_sibling <- function(dt, column_selection = c()) {
  
  if (length(column_selection) > 0) {
    dt <- dt[, c('IDC', 'mother', column_selection)]
  } # if no selection is specified, missingness in the entire dataframe is used
  
  # Get rid of empty NA values for mother
  dt <- dt[!is.na(dt$mother),]
  # Determine a list of mothers that have more than one child in the set.
  m = dt$mother[duplicated(dt$mother)] # i.e.  which mother IDs recur more than once
  
  # Initialize list to fill with IDs of siblings to exclude
  worse_sibling = list()
  
  # Loop through the duplicated mother IDs and identify siblings 
  for (i in 1:length(m)) {
    siblings = dt[dt$mother == m[i], ] # identify the couples of siblings
    # Code does not account for more than 3 siblings
    if (nrow(siblings) > 3)  { message("ATTENTION: some mothers have more than 3 siblings, remove them manually") 
    
    # First, mothers with 3 siblings, let's select the "best" one and get rid of the other two
    } else if (nrow(siblings) == 3) {
      nmiss = c( sum(is.na(siblings[1,])), sum(is.na(siblings[2,])), sum(is.na(siblings[3,])) )
      if (which.min(nmiss) == 1) { worse_sibling = c(worse_sibling, siblings[2,'IDC'], siblings[3,'IDC'])
      } else if (which.min(nmiss) == 2) { worse_sibling = c(worse_sibling, siblings[1,'IDC'], siblings[3,'IDC'])
      } else {  worse_sibling = c(worse_sibling, siblings[1,'IDC'], siblings[2,'IDC']) }
    
    # Otherwise, select the "worse" sibling (with more missing) and add to it the the black list
    } else if ( sum(is.na(siblings[1,])) > sum(is.na(siblings[2,])) ) {
      worse_sibling = c(worse_sibling, siblings[1,'IDC'])
    } else if ( sum(is.na(siblings[1,])) == sum(is.na(siblings[2,])) ) {
      worse_sibling = c(worse_sibling, siblings[sample(1:2, 1),'IDC'])
    } else { worse_sibling = c(worse_sibling, siblings[2,'IDC']) }
  }
  
  return(worse_sibling)
}

#-------------------------------------------------------------------------------
flowchart <- function(df, return_selected_sample = F) {
  
  fc <- list(initial_sample = nrow(df))
  # 
  step1 <- df[df$pre_percent_missing < 50.0,]
  loss <- nrow(step1) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_pren = loss, after_pren_selection = nrow(step1))
  #
  step2 <- step1[step1$post_percent_missing < 50.0,]
  loss <- nrow(step2) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_post = loss, after_post_selection = nrow(step2))
  #
  # step3 <- step2[!is.na(step2$intern_score_z),] 
  # loss <- nrow(step3) - as.numeric(fc[length(fc)])
  # fc <- c(fc, no_inte = loss, after_inte_selection = nrow(step3))
  # #
  # step4 <- step3[!is.na(step3$fat_mass_z),]
  # loss <- nrow(step4) - as.numeric(fc[length(fc)])
  # fc <- c(fc, no_fatm = loss, after_fatm_selection = nrow(step4))
  #
  step5 <- step2[step2$twin == 0,]
  loss <- nrow(step5) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_twin = loss, after_twin_selection = nrow(step5))
  # 
  worse_sib_list <- select_sibling(step5, column_selection = c(pre_LE, pre_CR, pre_PR, pre_IR,
                                                               post_LE, post_CR, post_PR, post_IR, post_DV,
                                                               outcomes_13y, outcomes_09y, covars, auxil))
  finalsample <- step5[step5$IDC %notin% worse_sib_list, ]
  loss <- nrow(finalsample) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_siblings = loss, final_sample = nrow(finalsample))
  
  print(fc)
  
  if (return_selected_sample == T) { return(finalsample)    }
  if (return_selected_sample == F) { return(worse_sib_list) }
  
}


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

#-------------------------------------------------------------------------------
passive_imp_formula <- function(domain, add = "") {
  conc <- paste(domain, collapse = " + ")
  str <- paste0("~I( (", conc, ") / ", length(domain), ")")
  if (add != "") {
    str <- paste0("~I( (", conc, " + ", add, ") / ", length(domain)+1, ")")
  }
  return(str)
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------