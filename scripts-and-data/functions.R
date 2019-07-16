## Functions for the productivity analysis of ENHG derivational suffixes for my Master's thesis.

# General "make" functions ------------

make_spc <- function(dataset){
  ## This function takes in a dataset exported from ANNIS and reformatted and returns
  ## an spc object, ready for processing with zipfR.
  
  if(!all(is.na(dataset))){
  
    ## Create dataframes with the frequencies of each token (lemmatised).
    data.tokens <- as.data.frame(table(as.vector(dataset$lemma)))
    
    ## Create a frequency spectrum for current subset by taking the frequency of those frequencies.
    data.freq.spec <- as.data.frame(table(data.tokens$Freq))
    
    ## Rename the headers the way zipfR wants them.
    names(data.freq.spec) <- c("m", "Vm")
    
    ## Write dataframe to external .spc file.
    write.table(data.freq.spec, "r-analysis/data.spc", sep="\t", row.names=FALSE)
    
    ## Read data.spc back in with read.spc; now the spectrum is ready for use with zipfR.
    data.spc <- read.spc("r-analysis/data.spc")
    
    return(data.spc)
  
  }else{
    return(NA)
  }
}


make_freq_list <- function(concordance){
  ## Make a  frequency list from the concordance (ordered alphabetically by lemma).
  if(!all(is.na(concordance))){
    concordance.freqlist <- as.data.frame(table(as.vector(concordance$lemma)))
    names(concordance.freqlist) <- c("lemma", "frequency")
    ## This next line is here to convert the "lemma" column from factors to just a character vector.
    concordance.freqlist[,] <- lapply(concordance.freqlist, function(x) type.convert(as.character(x), as.is = TRUE))
    return(concordance.freqlist)
  }else{
    return(NA)  
  }
}


make_rank_freq_profile <- function(frequency.list){
  ## Make a rank/frequency profile (including lemmas) from the frequency list by 
  ## sorting that dataframe and adding a column for rank (which increments)
  
  if(!all(is.na(frequency.list))){
    concordance.rankprofile <- frequency.list[order(-frequency.list$frequency),]
    rank <- seq(1, length.out=nrow(concordance.rankprofile))
    concordance.rankprofile <- cbind(concordance.rankprofile, rank)
    return(concordance.rankprofile)
  }else{
    return(NA)
  }
}

make_empirical_growth_df <- function(concordance){
  
  if(!all(is.na(concordance))){
    conclength <- length(concordance$id)
    emp <- data.frame(matrix(ncol=2))
    names(emp) <- c("N", "V")
    discovered <- c()
    
    for(idx in 1:conclength){
      ## Check if the lemma has already been seen or not; if not, store it in
      ## the initially empty vector "discovered".
      
      lemma <- as.character(concordance$lemma)
      
      if(!(is.element(lemma[idx], discovered))){
        discovered <- c(discovered, lemma[idx])
      }
      
      ## Write the results to the "emp" dataframe; N = idx, V = length(discovered)
      emp[idx,] <- c(idx, length(discovered))
      
    }
    return(emp)
  }else{
    return(NA)
  }
}

# Potential productivity ------------

pp <- function(hapax.count, token.count){
  return(hapax.count/token.count)
}


strip_subsets_to_min_token <- function(subset.info.df, concs.list){
  ## Determines the smallest token size of all subsets in
  ## the given list, randomises the order of the concordances
  ## in each subset and takes only the same number as the min
  ## token size, resulting in equally sized subsets.
  
  ## Get smallest number of tokens.
  mintok <- min(subset.info.df$tokens)
  
  ## Randomise the concordances in each list element and then
  ## take the first mintok of them.
  ## (Otherwise it's only chronological, if you do it by period,
  ## so you always miss the end of the period)
  set.seed(42)
  gctoken.concs <- list()
  
  for(idx in 1:length(concs.list)){
    rando <- concs.list[[idx]][sample(1:nrow(concs.list[[idx]])), ]
    gctoken.concs[[idx]] <- head(rando, mintok)
  }
  
  return(gctoken.concs)
}


equate_subsets_to_n_tokens <- function(conc.info.df, concordance, num.subsets){
  ## Get total number of tokens in concordance.
  tokens <- conc.info.df$tokens
  
  ## Divide this number by desired number of subsets to get
  ## num tokens per subset (rounded down -- there might be a
  ## couple tokens left over at the end).
  tokens.per.subset <- floor(tokens/num.subsets)
  
  ## Split concordance into a list with num.subsets elements,
  ## each of which contains tokens.per.subset successive tokens.
  equated.concs <- list()
  for(sub in 1:num.subsets){
    start.idx <- tokens.per.subset*(sub-1)+1
    end.idx <- tokens.per.subset*sub
    equated.concs[[sub]] <- concordance[start.idx:end.idx, ]
    
  }
  
  ## Get the first and last year in each new subset and label
  ## list elements accordingly.
  labels <- c()
  for(idx in 1:length(equated.concs)){
    minyear <- min(equated.concs[[idx]]$date)
    maxyear <- max(equated.concs[[idx]]$date)
    labels <- c(labels, paste0(minyear, "-", maxyear))
  }
  names(equated.concs) <- labels
  return(equated.concs)
}

pp_of_equally_sized_conc_list <- function(conc.list){
  ## Only works for lists with the same number of tokens per
  ## list element.
  
  ## Save token count.
  token.count <- length(conc.list[[1]]$lemma)
  
  ## Get hapax count per subset.
  hapax.list <- c()
  for(idx in 1:length(conc.list)){
    freqtable <- as.data.frame(table(as.vector(conc.list[[idx]]$lemma)))
    spectrum <- as.data.frame(table(as.vector(freqtable$Freq)))
    hapax.list <- c(hapax.list, spectrum[1,2])
  }
  
  pp.list <- hapax.list/token.count
  
  pp.df <- cbind(names(conc.list), pp.list)
  colnames(pp.df) <- c("subset", "pp")
  
  return(pp.df)
}


# Concordance functions ------------

gen_fzm_options <- function(spectrum){
  models <- list()

  ## Make a new fZM model for each combination of parameters, save each one to a list.
  ## Then, I can go through this list and find the best gof.
  
  for(costfunc in c("chisq", "linear", "smooth.linear", "mse", "exact")){
    for(exactval in c(T, F)){
      for(methodalg in c("Custom", "NLM", "Nelder-Mead", "SANN")){
        
        model <- tryCatch(
          {
            ## The code to run if no error arises.
            lnre("fzm", spectrum, cost=costfunc, exact=exactval, method=methodalg)
          },
          error = function(cond){
            ## Return output if there is an error, e.g. parameter out of range or
            ## estimation fails.
            # message("Here's the original error message:")
            # message(cond)
            return(NA)
          }
        )
          
        modelname <- paste(costfunc, 
                           as.character(exactval),
                           methodalg,
                           sep="_")
        models[[ modelname ]] <- model
      }
    }
  }
  return(models)
}


make_fzm_info_df <- function(list.of.fzms){
  ## Create dataframe with information about each fZM model contained in input list.
  ## If model estimation didn't work, then the model information is 
  ## included as NA.
  
  fzm.info <- data.frame()
  
  for(idx in 1:length(list.of.fzms)){
    
    temp.info <- setNames(data.frame(matrix(ncol = 9, nrow = 1)), 
                          c("id", "cost", "exact", "method", "alpha", "S", "gof.X2", "gof.df", "gof.p"))
    
    curr.combo <- names(list.of.fzms)[idx]
    curr.combo.sep <- strsplit(curr.combo, "_")

    temp.info["id"] = curr.combo
    temp.info["cost"] = curr.combo.sep[[1]][1]
    temp.info["exact"] = curr.combo.sep[[1]][2]
    temp.info["method"] = curr.combo.sep[[1]][3]
        
    if(!all(is.na(list.of.fzms[[idx]]))){
      ## If list element not NA, i.e. fZM model estimation worked.
      temp.info["alpha"] = list.of.fzms[[idx]]$param$alpha
      temp.info["S"] = list.of.fzms[[idx]]$S
      temp.info["gof.X2"] = list.of.fzms[[idx]]$gof$X2
      temp.info["gof.df"] = list.of.fzms[[idx]]$gof$df
      temp.info["gof.p"] = list.of.fzms[[idx]]$gof$p
    }else{
      temp.info["alpha"] = NA
      temp.info["S"] = NA
      temp.info["gof.X2"] = NA
      temp.info["gof.df"] = NA
      temp.info["gof.p"] = NA
    }
    fzm.info <- rbind(fzm.info, temp.info)
  }
  return(fzm.info)
}



make_conc_basic_info_df <- function(spectrum){
  ## Makes a dataframe with overall type, token, and hapax counts, as well as
  ## data from the fZM model calculated for the entire concordance.
  
  #Initialize empty data frame.
  new.df <- data.frame(matrix(ncol = 3, nrow = 0))
  
  ## Add data together to make new info dataframe.
  info <- rbind(new.df, c(N(spectrum), 
                          V(spectrum), 
                          Vm(spectrum, 1)))
  names(info) <- c("tokens", "types", "hapaxes")
  
  return(info)
}


# Subperiod functions ------------

get_regular_timespans <- function(timestep){
  ## Returns a vector starting with 1482 and ending with 1915, and starting at
  ## 1500, the given time step proceeds forward in that interval.
  # c(1482, seq(1500,1914,timestep), 1915)
  c(seq(1482, 1914, timestep), 1915)
}


get_timespans_as_char_vec <- function(timevec){
  ## Based on the input vector of years, creates a list of these spans as strings,
  ## for use in dataframes containing info about the subsets and as legends of plots.
  spans <- c()
  
  for(idx in 1:(length(timevec)-1)){
    span <- paste0(as.character(timevec[idx]), "-", as.character(timevec[idx+1]-1))
    spans <- append(spans, span)
  }
  return(spans)
}


get_period_subsets <- function(timespans, dataset){
  ## This function takes in a vector containing the years in which each interval
  ## should begin and end, i.e. it begins with 1482 and ends with 1915 with however
  ## many divisions in between. Each time period is understood as beginning with and
  ## including the first date and extending to but not including the second date (hence 1915,
  ## not 1914, as 1914 is the last year with data).
  ## The function also takes the dataset in full detail, in dataframe form (not spc).
  ## The function returns a list of dataframes which are the provided data in "dataset"
  ## divided into the intervals provided in "timespans".
  
  period.subsets <- list()
  for(idx in 1:(length(timespans)-1)){
    
    current.subset <- subset(dataset, date >= timespans[idx] & date < timespans[idx+1])
    # print(paste(as.character(timespans[idx]), "-", as.character(timespans[idx+1]-1)))
    period.subsets[[idx]] <- current.subset
    
  }
  return(period.subsets)
}


get_document_subsets <- function(doctitles.list, concordance){

  doc.subsets <- list()
  for(idx in 1:(length(doctitles.list))){
    
    current.subset <- subset(concordance, concordance$title == doctitles.list[idx])

    if(dim(current.subset)[1] == 0){
      doc.subsets[[idx]] <- NA
    }else{
      doc.subsets[[idx]] <- current.subset
    }
  }
  return(doc.subsets)
}


make_spcs_for_subsets <- function(subsets.list, data.vector){
  ## This function creates frequency spectra for every subset of the concordance
  ## dataframe contained in the input list and returns a list of the same length
  ## with all of the frequency spectra contained inside it.
  
  subset.spectra <- list()

  for(idx in 1:length(subsets.list)){
    working.dataset <- subsets.list[[idx]]
    working.spc <- make_spc(working.dataset, data.vector)
    subset.spectra[[idx]] <- working.spc
  }

  return(subset.spectra)
}


make_sp_info_df <- function(subset.spectra, label.vec){
  ## Create dataframe with information about each subset.
  ## Type, token, hapax information available for all subperiods.
  
  subset.info <- data.frame()
  
  for(idx in 1:length(subset.spectra)){
    temp.info <- setNames(data.frame(matrix(ncol = 6, nrow = 1)), 
                          c("subset", "tokens", "types", "hapaxes", "pp", "ttr"))
    
    temp.info["subset"] = label.vec[idx]
    
    if(!all(is.na(subset.spectra[[idx]]))){
      temp.info["tokens"] = N(subset.spectra[[idx]])
      temp.info["types"] = V(subset.spectra[[idx]])
      temp.info["hapaxes"] = Vm(subset.spectra[[idx]], 1)
      temp.info["pp"] = Vm(subset.spectra[[idx]], 1) / N(subset.spectra[[idx]])
      temp.info["ttr"] = V(subset.spectra[[idx]]) / N(subset.spectra[[idx]])
    }else{
      temp.info["tokens"] = NA
      temp.info["types"] = NA
      temp.info["hapaxes"] = NA
      temp.info["pp"] = NA
      temp.info["ttr"] = NA
    }
    subset.info <- rbind(subset.info, temp.info)
    
  }
  return(subset.info)
}


make_spcs_into_bin_vgcs <- function(subset.spectra){
  ## This function takes the list of spectra for each subset and converts
  ## each of them into interpolated VGC dataframes, ready to be plotted.
  
  subset.bin.vgcs <- list()
  
  for(idx in 1:length(subset.spectra)){
    subset.bin.vgcs[[idx]] <- vgc.interp(subset.spectra[[idx]], 1:N(subset.spectra[[idx]]))
  }
  
  return(subset.bin.vgcs)
}


make_sp_fzms <- function(subset.spectra){
  ## Function to get the default lnre fZM models for each subset, with error handling
  ## (instead of saving a fZM model for a subperiod where parameters cannot
  ## be estimated, that list element instead just becomes NA.)
  
  subset.fzms <- list()
  
  for(idx in 1:length(subset.spectra)){
    subset.fzm <- tryCatch(
      {
        ## The 'try' part: this is the output if there is no error.
        lnre("fzm", subset.spectra[[idx]])  # exact=exact.boolean removed bc not default
      },
      error=function(cond){
        ## Return output if there is an error, e.g. parameter estimation fails.
        # message("Caught an error!")
        # message("Here's the original error message:")
        # message(cond)
        return(NA)
      }
    )
      
    subset.fzms[[idx]] <- subset.fzm
  }
  return(subset.fzms)
}

make_sp_fzms_exactf <- function(subset.spectra){
  ## Function to get the default lnre fZM models for each subset, with error handling
  ## (instead of saving a fZM model for a subperiod where parameters cannot
  ## be estimated, that list element instead just becomes NA.)
  
  subset.fzms <- list()
  
  for(idx in 1:length(subset.spectra)){
    subset.fzm <- tryCatch(
      {
        ## The 'try' part: this is the output if there is no error.
        lnre("fzm", subset.spectra[[idx]], exact=F)
      },
      error=function(cond){
        ## Return output if there is an error, e.g. parameter estimation fails.
        # message("Caught an error!")
        # message("Here's the original error message:")
        # message(cond)
        return(NA)
      }
    )
    
    subset.fzms[[idx]] <- subset.fzm
  }
  return(subset.fzms)
}

gen_sp_fzm_options <- function(subset.spectra){
  ## Function to get all lnre fZM models for each subset, with error handling
  ## (instead of saving a fZM model for a subperiod where parameters cannot
  ## be estimated, that list element instead just becomes NA.)
  
  subset.fzms <- list()
  
  for(idx in 1:length(subset.spectra)){
    for(costfunc in c("chisq", "linear", "smooth.linear", "mse", "exact")){
      for(exactval in c(T, F)){
        for(methodalg in c("Custom", "NLM", "Nelder-Mead", "SANN")){
          subset.fzm <- tryCatch(
            {
              ## Code to run if no error arises.
              lnre("fzm", subset.spectra[[idx]], cost=costfunc, exact=exactval, method=methodalg)
            },
            error=function(cond){
              ## Return output if there is an error, e.g. parameter estimation fails.
              # message("Caught an error!")
              # message("Here's the original error message:")
              # message(cond)
              return(NA)
            }
          )
          
          modelname <- paste(as.character(idx),
                             costfunc,
                             as.character(exactval),
                             methodalg,
                             sep="_")
          
          subset.fzms[[ modelname ]] <- subset.fzm
          
        }
      }
    }
  }
  return(subset.fzms)
}

make_sp_fzm_info_df <- function(list.of.fzms, times.chars){
  ## Create dataframe with information about each fZM model contained in input list
  ## (subperiod edition).
  ## If model estimation didn't work, then the model information is 
  ## included as NA.
  
  fzm.info <- data.frame()
  
  for(idx in 1:length(list.of.fzms)){
    
    temp.info <- setNames(data.frame(matrix(ncol = 10, nrow = 1)), 
                          c("id", "period", "cost", "exact", "method", "alpha", "S", "gof.X2", "gof.df", "gof.p"))
    
    curr.combo <- names(list.of.fzms)[idx]
    curr.combo.sep <- strsplit(curr.combo, "_")
    
    temp.info["id"] = curr.combo
    temp.info["period"] = times.chars[ as.numeric(curr.combo.sep[[1]][1]) ]
    temp.info["cost"] = curr.combo.sep[[1]][2]
    temp.info["exact"] = curr.combo.sep[[1]][3]
    temp.info["method"] = curr.combo.sep[[1]][4]  
    
    if(!all(is.na(list.of.fzms[[idx]]))){
      ## If list element not NA, i.e. fZM model estimation worked.
      temp.info["alpha"] = list.of.fzms[[idx]]$param$alpha
      temp.info["S"] = list.of.fzms[[idx]]$S
      temp.info["gof.X2"] = list.of.fzms[[idx]]$gof$X2
      temp.info["gof.df"] = list.of.fzms[[idx]]$gof$df
      temp.info["gof.p"] = list.of.fzms[[idx]]$gof$p
    }else{
      temp.info["alpha"] = NA
      temp.info["S"] = NA
      temp.info["gof.X2"] = NA
      temp.info["gof.df"] = NA
      temp.info["gof.p"] = NA
    }
    fzm.info <- rbind(fzm.info, temp.info)
  }
  return(fzm.info)
}



make_sp_fzm_vgcs <- function(subset.fzms, upperlimit){
  
  subset.fzm.vgcs <- list()
  
  for(idx in 1:length(subset.fzms)){
    
    if(!all(is.na(subset.fzms[[idx]]))){
      subset.fzm.vgcs[[idx]] <- lnre.vgc(subset.fzms[[idx]], 1:upperlimit, variances=TRUE)
    }else{
      subset.fzm.vgcs[[idx]] <- NA
    }
    
  }
  return(subset.fzm.vgcs)
}


get_subset_concs_from_period_vector <- function(period.vector, dataset){
  ## This function takes in the vector with years for each time period
  ## and returns the list containing frequency spectra for the requested
  ## periods (relying on helper functions defined above).
  
  f.subsets.conc <- get_period_subsets(period.vector, dataset)

  return(f.subsets.conc)
}


apply_func_to_subset_list <- function(func.to.apply, subset.list){
  subset.newlist <- list()
  for(idx in 1:length(subset.list)){
    subset.listitem <- func.to.apply(subset.list[[idx]])
    subset.newlist[[idx]] <- subset.listitem
  }
  return(subset.newlist)
}

# Other -----------------------

normalise <- function(corp.norm, corp.actual, count.actual){
  count.norm <- (count.actual * corp.norm) / corp.actual
  return(count.norm) 
}

round2 = function(x, n) {
  ## from https://stackoverflow.com/questions/12688717/round-up-from-5 
  ## x is the number to round, n is the number of digits to be rounded to
  
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

add_norm_to_doc_info <- function(info.df, tokens.df, corpus.norm){
  
  for(idx in 1:length(info.df$subset)){
    corpus.actual <- tokens.df$dipl.tokens[ info.df$subset[idx] == tokens.df$title ]
    
    info.df$corpus.actual[idx] <- corpus.actual
    info.df$corpus.norm[idx] <- corpus.norm
    
    info.df$tokens.norm[idx] <- round2( normalise(corpus.norm, corpus.actual, info.df$tokens[idx]) , 0)
    info.df$types.norm[idx] <- round2( normalise(corpus.norm, corpus.actual, info.df$types[idx]) , 0)
    info.df$hapaxes.norm[idx] <- round2( normalise(corpus.norm, corpus.actual, info.df$hapaxes[idx]) , 0)
  }
  return(info.df)
}

rate_of_addition <- function(rankprof.list){
  
  roa <- setNames(data.frame(matrix(ncol = 2, nrow = length(rankprof.list) )), 
                  c("subset", "new.types"))
  
  ## We take the first subset as the initial lexicon.
  init.lex <- rankprof.list[[1]]$lemma
  # print(init.lex)
  
  roa$subset <- names(rankprof.list)[1]
  roa$new.types <- length(init.lex)
  
  ## For all the rest of the subsets, loop over their elements and check if each
  ## is already contained in initial lexicon. If not, add it to new.types.
  ## Documents without tokens (NA in rank profile) are considered as adding no new types.
  
  for(idx in 2:length(rankprof.list)){
    
    if(!all(is.na(rankprof.list[idx]))){
      
      current.types <- rankprof.list[[idx]]$lemma
      new.types <- c()
      
      for(wd in current.types){
        ## Check which of the types in this subset are not in the initial lexicon.
        if(!is.element(wd, init.lex)){
          new.types <- c(new.types, wd)
        }
      }
      
    }else{
      ## If no tokens for that subperiod, i.e. that list entry is NA.
      new.types <- c()
    }
    
    roa$subset[idx] <- names(rankprof.list)[idx]
    roa$new.types[idx] <- length(new.types)
    
    # Save the initial lexicon size as an attribute of the dataframe.
    attr(roa, 'init.lex') <- length(init.lex)
    
  }
  return(roa)
}

# Save/export results -----------------------


save_mc_s_as_csv <- function(){
  write.csv(er_mc_approx_overview, file = "r-analysis/mc-analysis/er_mc_s_overview.csv", row.names=FALSE)
  write.csv(er_mc_approx_stats_no_outliers, file = "r-analysis/mc-analysis/er_mc_s_stats_no_outliers.csv", row.names=FALSE)
  write.csv(er_mc_approx_stats_w_outliers, file = "r-analysis/mc-analysis/er_mc_s_stats_w_outliers.csv", row.names=FALSE)
  
  write.csv(hkeit_mc_approx_overview, file = "r-analysis/mc-analysis/hkeit_mc_s_overview.csv", row.names=FALSE)
  write.csv(hkeit_mc_approx_stats_no_outliers, file = "r-analysis/mc-analysis/hkeit_mc_s_stats_no_outliers.csv", row.names=FALSE)
  write.csv(hkeit_mc_approx_stats_w_outliers, file = "r-analysis/mc-analysis/hkeit_mc_s_stats_w_outliers.csv", row.names=FALSE)
  
  write.csv(ung_mc_approx_overview, file = "r-analysis/mc-analysis/ung_mc_s_overview.csv", row.names=FALSE)
  write.csv(ung_mc_approx_stats_no_outliers, file = "r-analysis/mc-analysis/ung_mc_s_stats_no_outliers.csv", row.names=FALSE)
  write.csv(ung_mc_approx_stats_w_outliers, file = "r-analysis/mc-analysis/ung_mc_s_stats_w_outliers.csv", row.names=FALSE)
}

save_mc_pp_as_csv <- function(){
  write.csv(er_ppmc_stats, file = "r-analysis/mc-analysis/er_mc_pp_stats.csv", row.names=FALSE)
  write.csv(hkeit_ppmc_stats, file = "r-analysis/mc-analysis/hkeit_mc_pp_stats.csv", row.names=FALSE)
  write.csv(ung_ppmc_stats, file = "r-analysis/mc-analysis/ung_mc_pp_stats.csv", row.names=FALSE)
}

save_info_as_csv <- function(){
  write.csv(er_conc_info, file = "r-analysis/info-tables/er_conc_info.csv", row.names=FALSE)
  write.csv(er_sp_info, file = "r-analysis/info-tables/er_sp_info.csv", row.names=FALSE)
  write.csv(er_doc_info, file = "r-analysis/info-tables/er_doc_info.csv", row.names=FALSE)
  
  write.csv(hkeit_conc_info, file = "r-analysis/info-tables/hkeit_conc_info.csv", row.names=FALSE)
  write.csv(hkeit_sp_info, file = "r-analysis/info-tables/hkeit_sp_info.csv", row.names=FALSE)
  write.csv(hkeit_doc_info, file = "r-analysis/info-tables/hkeit_doc_info.csv", row.names=FALSE)
  
  write.csv(ung_conc_info, file = "r-analysis/info-tables/ung_conc_info.csv", row.names=FALSE)
  write.csv(ung_sp_info, file = "r-analysis/info-tables/ung_sp_info.csv", row.names=FALSE)
  write.csv(ung_doc_info, file = "r-analysis/info-tables/ung_doc_info.csv", row.names=FALSE)
}

save_spcs_as_csv <- function(){
  write.csv(er_conc_spc, file = "r-analysis/spcs/er_conc_spc.csv", row.names=FALSE)
  write.csv(er_sp_spc_tall, file = "r-analysis/spcs/er_sp_spc_tall.csv", row.names=FALSE)
  write.csv(er_doc_spc_tall, file = "r-analysis/spcs/er_doc_spc_tall.csv", row.names=FALSE)
  
  write.csv(hkeit_conc_spc, file = "r-analysis/spcs/hkeit_conc_spc.csv", row.names=FALSE)
  write.csv(hkeit_sp_spc_tall, file = "r-analysis/spcs/hkeit_sp_spc_tall.csv", row.names=FALSE)
  write.csv(hkeit_doc_spc_tall, file = "r-analysis/spcs/hkeit_doc_spc_tall.csv", row.names=FALSE)
  
  write.csv(ung_conc_spc, file = "r-analysis/spcs/ung_conc_spc.csv", row.names=FALSE)
  write.csv(ung_sp_spc_tall, file = "r-analysis/spcs/ung_sp_spc_tall.csv", row.names=FALSE)
  write.csv(ung_doc_spc_tall, file = "r-analysis/spcs/ung_doc_spc_tall.csv", row.names=FALSE)
}

save_fzm_vgcs_as_csv <- function(){
  ## For entire concordance.
  write.csv(ung_conc_fzm_approx_vgc, file = "r-analysis/vgcs/ung_conc_fzm_approx_vgc.csv", row.names=FALSE)
  write.csv(er_conc_fzm_approx_vgc, file = "r-analysis/vgcs/er_conc_fzm_approx_vgc.csv", row.names=FALSE)
  write.csv(hkeit_conc_fzm_approx_vgc, file = "r-analysis/vgcs/hkeit_conc_fzm_approx_vgc.csv", row.names=FALSE)
}

save_emp_vgcs_as_csv <- function(){
  write.csv(er_sp_emp_vgc_tall, file = "r-analysis/vgcs/er_sp_emp_vgc_tall.csv", row.names=FALSE)
  write.csv(er_conc_emp_vgc, file = "r-analysis/vgcs/er_conc_emp_vgc.csv", row.names=FALSE)
  write.csv(er_doc_emp_vgc_tall, file = "r-analysis/vgcs/er_doc_emp_vgc_tall.csv", row.names=FALSE)
  
  write.csv(hkeit_sp_emp_vgc_tall, file = "r-analysis/vgcs/hkeit_sp_emp_vgc_tall.csv", row.names=FALSE)
  write.csv(hkeit_conc_emp_vgc, file = "r-analysis/vgcs/hkeit_conc_emp_vgc.csv", row.names=FALSE)
  write.csv(hkeit_doc_emp_vgc_tall, file = "r-analysis/vgcs/hkeit_doc_emp_vgc_tall.csv", row.names=FALSE)
  
  write.csv(ung_sp_emp_vgc_tall, file = "r-analysis/vgcs/ung_sp_emp_vgc_tall.csv", row.names=FALSE)
  write.csv(ung_conc_emp_vgc, file = "r-analysis/vgcs/ung_conc_emp_vgc.csv", row.names=FALSE)
  write.csv(ung_doc_emp_vgc_tall, file = "r-analysis/vgcs/ung_doc_emp_vgc_tall.csv", row.names=FALSE)
}

save_bin_vgcs_as_csv <- function(){
  write.csv(er_conc_bin_vgc, file = "r-analysis/vgcs/er_conc_bin_vgc.csv", row.names=FALSE)
  write.csv(hkeit_conc_bin_vgc, file = "r-analysis/vgcs/hkeit_conc_bin_vgc.csv", row.names=FALSE)
  write.csv(ung_conc_bin_vgc, file = "r-analysis/vgcs/ung_conc_bin_vgc.csv", row.names=FALSE)
}

save_freqlists_as_csv <- function(){
  write.csv(er_conc_freqlist, file = "r-analysis/freqlists/er_conc_freqlist.csv", row.names=FALSE)
  write.csv(er_sp_freqlist_tall, file = "r-analysis/freqlists/er_sp_freqlist_tall.csv", row.names=FALSE)
  write.csv(er_doc_freqlist_tall, file = "r-analysis/freqlists/er_doc_freqlist_tall.csv", row.names=FALSE)
  
  write.csv(hkeit_conc_freqlist, file = "r-analysis/freqlists/hkeit_conc_freqlist.csv", row.names=FALSE)
  write.csv(hkeit_sp_freqlist_tall, file = "r-analysis/freqlists/hkeit_sp_freqlist_tall.csv", row.names=FALSE)
  write.csv(hkeit_doc_freqlist_tall, file = "r-analysis/freqlists/hkeit_doc_freqlist_tall.csv", row.names=FALSE)  
  
  write.csv(ung_conc_freqlist, file = "r-analysis/freqlists/ung_conc_freqlist.csv", row.names=FALSE)
  write.csv(ung_sp_freqlist_tall, file = "r-analysis/freqlists/ung_sp_freqlist_tall.csv", row.names=FALSE)
  write.csv(ung_doc_freqlist_tall, file = "r-analysis/freqlists/ung_doc_freqlist_tall.csv", row.names=FALSE)
}

save_rankprofiles_as_csv <- function(){
  write.csv(er_conc_rankprof, file = "r-analysis/rank-freq-profiles/er_conc_rankprof.csv", row.names=FALSE)
  write.csv(er_sp_rankprof_tall, file = "r-analysis/rank-freq-profiles/er_sp_rankprof_tall.csv", row.names=FALSE)
  write.csv(er_doc_rankprof_tall, file = "r-analysis/rank-freq-profiles/er_doc_rankprof_tall.csv", row.names=FALSE)
  
  write.csv(hkeit_conc_rankprof, file = "r-analysis/rank-freq-profiles/hkeit_conc_rankprof.csv", row.names=FALSE)
  write.csv(hkeit_sp_rankprof_tall, file = "r-analysis/rank-freq-profiles/hkeit_sp_rankprof_tall.csv", row.names=FALSE)
  write.csv(hkeit_doc_rankprof_tall, file = "r-analysis/rank-freq-profiles/hkeit_doc_rankprof_tall.csv", row.names=FALSE)  
  
  write.csv(ung_conc_rankprof, file = "r-analysis/rank-freq-profiles/ung_conc_rankprof.csv", row.names=FALSE)
  write.csv(ung_sp_rankprof_tall, file = "r-analysis/rank-freq-profiles/ung_sp_rankprof_tall.csv", row.names=FALSE)
  write.csv(ung_doc_rankprof_tall, file = "r-analysis/rank-freq-profiles/ung_doc_rankprof_tall.csv", row.names=FALSE)
}
