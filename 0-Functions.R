# First, let's point to the necessary libraries
library(foreign)
library(stats)
library(mice)
library(miceadds)

# This will come in handy for exclusion
'%notin%' <- Negate('%in%')

# Defining the path to the data
# check if the path to the data is already in memory, otherwise ask for it. 
if (exists("pathtodata") == F) { pathtodata = readline(prompt="Enter path to data: ") }

# ATTENTION!!! If prompted with an "Enter path to data:" message -> Enter the location
# of your datafiles. The code assumes that all (raw) data is stored in ONE folder.
# Do not forget the final slash in your path, and, speaking of slashes, beware of 
# OS sensitive changes when you want to modify the structure of your dirs!


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
outcomes <- c('intern_score_z', 'fmi_z', 'risk_groups') # fat_mass_z
covars   <- c('sex', 'age_child', 'm_bmi_berore_pregnancy', 'm_smoking', 'm_drinking')
auxil    <- c('m_bmi_pregnancy','m_dep_cont_pregnancy', 'p_dep_cont_pregnancy', # for postnatal only 
              'm_bmi_5yrs', 'm_dep_cont_3yrs', 'p_dep_cont_3yrs',               # for prenatal only 
              'ethnicity', 'parity', 'gest_age_birth', 'gest_weight', 'm_age_cont')
exclusion_criteria <- c('pre_percent_missing', 'post_percent_missing', 'twin', 'mother')


#### ------------------------------ FUNCTIONS ----------------------------- ####

# Read in the data and fix all SPSS weird missing codes into NAs
readquick <- function(filename, rootdir = pathtodata, exclude_col = "") { # only works for SPSS files
  dat <- read.spss(paste(rootdir, filename, sep=""), 
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
  step3 <- step2[!is.na(step2$intern_score_z),] 
  loss <- nrow(step3) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_inte = loss, after_inte_selection = nrow(step3))
  #
  step4 <- step3[!is.na(step3$fmi_z),] # fat_mass_z
  loss <- nrow(step4) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_fatm = loss, after_fatm_selection = nrow(step4))
  #
  step5 <- step4[step4$twin == 0,]
  loss <- nrow(step5) - as.numeric(fc[length(fc)])
  fc <- c(fc, no_twin = loss, after_twin_selection = nrow(step5))
  # 
  worse_sib_list <- select_sibling(step5, column_selection = c(pre_LE, pre_CR, pre_PR, pre_IR,
                                                               post_LE, post_CR, post_PR, post_IR, post_DV,
                                                               outcomes, covars, auxil))
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