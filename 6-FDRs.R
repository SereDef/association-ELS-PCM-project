# FDR 
library("readxl")
library("fuzzySim")

multiplesheets <- function(fname) {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)

  # assigning names to data frames
  names(data_frame) <- sheets
  
  # return data frame
  return(data_frame)
}

# specifying the path name
path <- "/Users/Serena/Desktop/mockData/Results/2022-04-15_GENR_Results.xlsx"
path <- "/Users/Serena/Desktop/Paperdraft/ALSPAC-results/2022-04-15_ALSPAC_Results.xlsx"

f = multiplesheets(path)

checkFDR <- function(model) {
  for (s in 1:length(f)) {
    message("\n", names(f)[s])
    # remove NA rows 
    f[[s]] <- f[[s]][rowSums(is.na(f[[s]])) != ncol(f[[s]]),]
    # adapt names 
    if (names(f[[s]][2]) == 'y.level') {
      f[[s]]['TERM'] <- paste(f[[s]][,1], f[[s]][,2], f[[s]][,3], sep=" - ")
    } else { f[[s]]['TERM'] <- paste(f[[s]][,1], f[[s]][,2], sep=" - ") }
    
    if (grepl("cma", names(f)[s], fixed = TRUE) & !grepl("other", names(f)[s], fixed = TRUE) ) {
      rep = fuzzySim::FDR(pvalues = data.frame(var = f[[s]]$TERM, pval = f[[s]]$P.val));
      print(rep$exclude[rep$exclude$p.value < 0.05, ])
    } else {
      mod <- f[[s]][grepl(model, f[[s]]$model, fixed = TRUE), ]
      rep = fuzzySim::FDR(pvalues = data.frame(var = mod$TERM, pval = mod$p.value));
      print(rep$exclude[rep$exclude$p.value < 0.05, ])
    }
  }
}


checkFDR("pren model")
checkFDR(" post model")
checkFDR("pren+post model")
