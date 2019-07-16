## General functions ----

mc_stats <- function(mc_df, tidied.mc.list){
  ## Takes a tidied list where all NAs (and optionally outliers) have already been
  ## removed, and returns a dataframe that contains the step number, the mean, and
  ## the 95% confidence interval for the MC distribution at each measurement point.
  means.vec <- montecarlo_s_means(tidied.mc.list)
  ci.df <- montecarlo_s_ci(tidied.mc.list)
  
  mc.stats.df <- data.frame(matrix(nrow=length(tidied.mc.list), ncol=5))
  names(mc.stats.df) <- c("step", "n.seen", "mean", "lower", "upper")
  mc.stats.df$step <- 1:length(tidied.mc.list)
  mc.stats.df$n.seen <- as.numeric(names(mc_df))
  mc.stats.df$mean <- means.vec
  mc.stats.df$lower <- ci.df$lower
  mc.stats.df$upper <- ci.df$upper
  
  return(mc.stats.df)
} 

mc_overview <- function(mc_df, mc_list_wo_outliers, mc_list_w_outliers, num.iterations){
  ## A function that looks at how many NAs and outliers there are for a given model
  ## and outputs a dataframe with an overview of these numbers.
  
  overview_df <- data.frame(matrix(nrow=length(mc_list_wo_outliers), ncol=6))
  names(overview_df) <- c("step", "n.seen", "successes", "failures", "outliers", "data.pts")
  
  for(idx in 1:length(mc_list_wo_outliers)){
    overview_df$step[idx] <- idx
    overview_df$n.seen[idx] <- as.numeric(names(mc_df)[1])*idx
    overview_df$successes[idx] <- length(mc_list_w_outliers[[idx]])
    overview_df$failures[idx] <- num.iterations-length(mc_list_w_outliers[[idx]])
    overview_df$outliers[idx] <- length(mc_list_w_outliers[[idx]]) - length(mc_list_wo_outliers[[idx]])
    overview_df$data.pts[idx] <- length(mc_list_wo_outliers[[idx]])
  }
  
  return(overview_df)
}

take_selection_of_mc_list <- function(overview_df, tidied_mc_list){
  elemts_to_take <- min(overview_df$data.pts)
  selection <- list()
  
  for(idx in 1:length(tidied_mc_list)){
    set.seed(idx)
    random.idxs <- sample.int(length(tidied_mc_list[[idx]]), size=elemts_to_take)
    selected_elements <- tidied_mc_list[[idx]][random.idxs]
    selection[[idx]] <- selected_elements
  }
  
  return(selection)
}

mc_analysis <- function(data_info, mc_df, exact.bool){
  if(exact.bool==T){
    exactness <- "exact"
  }else{
    exactness <- "approx"
  }
  
  ## Clean the dataframe, removing NAs, and creating one cleaned version that also has
  ## outliers removed.
  mc_list_no_outl <- montecarlo_s_clean(mc_df, T)
  mc_list_w_outl <- montecarlo_s_clean(mc_df, F)
  
  ## Create a dataframe that gets counts for successful parameter estimations, how many
  ## outliers, etc.
  mc_overview <- mc_overview(mc_df, mc_list_no_outl, mc_list_w_outl, 2000)
  
  ## Take a random selection of the data for each measurement point that is the same size as
  ## the minimum number of datapoints of all points (after outlier removal), so that the statistics are
  ## computed over the same number of data points for each measurement point.
  mc_list_no_outl_selection <- take_selection_of_mc_list(mc_overview, mc_list_no_outl)
  mc_list_w_outl_selection <- take_selection_of_mc_list(mc_overview, mc_list_w_outl)
  
  ## Create a dataframe that contains statistics about the simulation at each point (mean and bounds of
  ## 95% confidence intervals).
  mc_list_no_outl_selection_stats <- mc_stats(mc_df, mc_list_no_outl_selection)
  mc_list_w_outl_selection_stats <- mc_stats(mc_df, mc_list_w_outl_selection)
  
  ## Option to uncomment and also calculate stats for the full dataset, where each measurement
  ## point has a different number of data points (the results are basically the same!)
  # mc_exact_no_outl_full_stats <- mc_stats(mc_20x2000_exact_no_outl)
  # mc_exact_w_outl_full_stats <- mc_stats(mc_20x2000_exact_w_outl)
  
  ## Save important files: ## overview tables, full stats, selection stats
  saveRDS(mc_overview, paste0("r-analysis/rds/", data_info["suffix.plain"], 
                              "_mc_", exactness, "_overview.rds"))
  saveRDS(mc_list_no_outl_selection_stats, paste0("r-analysis/rds/", data_info["suffix.plain"], 
                                                  "_mc_", exactness, "_stats_no_outliers.rds"))
  saveRDS(mc_list_w_outl_selection_stats, paste0("r-analysis/rds/", data_info["suffix.plain"], 
                                                 "_mc_", exactness, "_stats_w_outliers.rds"))
}

run_mc_analysis2 <- function(data_list, mc_df, exact.bool){
  
  mc_analysis(data_list, mc_df, exact.bool)
  
  regex_patt <- paste0("^", data_list["suffix.plain"], "_mc.*rds$")

  files <- list.files(pattern = regex_patt, recursive=TRUE)

  for(idx in 1:length(files)){
    var_name <- substr(files[idx], 16, nchar(files[idx])-4)

    assign(var_name,
           readRDS(files[idx]),
           envir = parent.frame())
  }
}

## Functions for S ----

montecarlo_s_df <- function(data.info, conc, num.steps, num.iterations, exact.bool){
  ## Randomises a concordance num.iterations times, divides it into num.steps
  ## equal steps, and progresses through each randomised version, calculating an
  ## fZM model at each step.
  
  ## Get helpful values.
  tokens <- length(conc$id)
  tokens.per.step <- floor(tokens/num.steps)  
  n.vals <- c(1:num.steps)*tokens.per.step
  
  ## Set up dataframe that will hold the data (columns: each value of N;
  ## rows: each iteration).
  mc.df <- data.frame(matrix(ncol=num.steps, nrow=num.iterations))
  names(mc.df) <- n.vals
  
  ## The randomising and sampling procedure is repeated many times and 
  ## leads to a distribution of that
  ## constant at each measurement point.
  
  for(iteridx in 1:num.iterations){
    
    ## -> Randomly permute the order in which the words appear in a text.
    set.seed(iteridx)  ## For reproducibility.
    random.conc <- conc[sample(1:nrow(conc)), ]
    
    ## -> Following permutation, calculate the values of a given constant for a 
    ## pre-specified number K0 of text lengths.  
    
    ## Get a list of concordances, progressively adding new tokens.per.step 
    ## onto each.
    
    mc.conc.list <- list()
    for(idx in 1:num.steps){
      end.idx <- n.vals[idx]
      mc.conc.list[[idx]] <- random.conc[1:end.idx, ]
    }
    
    ## Make concs into spectra.
    mc.conc.spc <- apply_func_to_subset_list(make_spc, mc.conc.list)
    
    ## Calculate the fZM for each N in mc.conc.list.
    
    if(exact.bool==T){
      ## Exact S:
      mc.conc.fzm <- make_sp_fzms(mc.conc.spc)
    }else{
      ## Approximate S:
      mc.conc.fzm <- make_sp_fzms_exactf(mc.conc.spc)
    }
    
    ## Save the S values (or NA if model estimation fails) as a row to mc.df.
    S.vec <- c()
    
    for(idx in 1:length(mc.conc.fzm)){
      if(!all(is.na(mc.conc.fzm[[idx]]))){
        S.vec <- c(S.vec, mc.conc.fzm[[idx]]$S)
      }else{
        S.vec <- c(S.vec, NA)
      }
    }
    
    mc.df[iteridx, ] <- S.vec
    message(paste(data.info["suffix.plain"], ": Finished iteration", iteridx, "of", num.iterations, "at", Sys.time()))
  }
  
  return(mc.df)
}


remove_outliers <- function(x, na.rm = TRUE, ...) {
  ## Source: https://stackoverflow.com/questions/4787332/how-to-remove-outliers-from-a-dataset
  ## Info: https://www.purplemath.com/modules/boxwhisk3.htm
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


montecarlo_s_clean <- function(mc.df, rm.outliers.bool){
  ## Cleans up the dataframe, removing NAs from where the parameter estimation failed
  ## and optionally removing outliers too (if rm.outliers.bool==T, then outliers will
  ## be removed).
  
  mc.clean.list <- list()
  
  for(colidx in 1:length(mc.df)){
    ## For every column in the dataframe, i.e. every measurement point,
    ## remove NAs.
    col.no.na <- mc.df[ , colidx][!is.na(mc.df[ ,colidx])]
    
    if(rm.outliers.bool==T){
      ## If outliers should be removed, then do that too.
      ## remove_outliers() replaces them with NA, so we just get
      ## rid of all NAs again to leave the final column.
      col.no.outliers <- remove_outliers(col.no.na)
      col.no.outliers <- col.no.outliers[!is.na(col.no.outliers)]
      mc.clean.list[[colidx]] <- col.no.outliers
    }else{
      mc.clean.list[[colidx]] <- col.no.na
    }
  }
  return(mc.clean.list)
}


montecarlo_s_means <- function(mc.list){
  ## Calculates the mean for each measurement point (list element)
  ## and returns it as a vector.
  
  means <- c()
  for(idx in 1:length(mc.list)){
    means <- c(means, mean(mc.list[[idx]]))
  }
  return(means)
}

montecarlo_s_ci <- function(mc.list){
  ## Calculates the 95% CI for each measurement point (list element)
  ## and returns upper and lower bounds in a dataframe.
  ## (Where the 2.5% doesn't divide cleanly, I round the number of elements to exclude
  ## on either side of the interval to the nearest integer.)
  
  mc.cis <- data.frame(matrix(ncol=2, nrow=length(mc.list)))
  names(mc.cis) <- c("lower", "upper")
  
  for(idx in 1:length(mc.list)){
    elemts_outside_interval <- round2(length(mc.list[[idx]])*0.025, 0)
    lower.idx <- elemts_outside_interval
    upper.idx <- length(mc.list[[idx]]) - elemts_outside_interval
    
    ## Sort list in ascending order.
    estimates.sorted <- sort(mc.list[[idx]])
    
    lower.bound <- estimates.sorted[lower.idx]
    upper.bound <- estimates.sorted[upper.idx]
    
    # print(paste("Lower bound at index", lower.idx, "is", lower.bound))
    # print(paste("Upper bound at index", upper.idx, "is", upper.bound))
    # print(lower.bound < upper.bound)
    
    mc.cis$lower[idx] <- lower.bound
    mc.cis$upper[idx] <- upper.bound
  }

  return(mc.cis)
  
}


## Functions for PP ----

montecarlo_pp_df <- function(data.info, conc, num.steps, num.iterations){
  ## Randomises a concordance num.iterations times, divides it into num.steps
  ## equal steps, and progresses through each randomised version, calculating PP
  ## at each step.
  
  ## Get helpful values.
  tokens <- length(conc$id)
  tokens.per.step <- floor(tokens/num.steps)  
  n.vals <- c(1:num.steps)*tokens.per.step
  
  ## Set up dataframe that will hold the data (columns: each value of N;
  ## rows: each iteration).
  mc.df <- data.frame(matrix(ncol=num.steps, nrow=num.iterations))
  names(mc.df) <- n.vals
  
  ## The randomising and sampling procedure is repeated many times and 
  ## leads to a distribution of that value at each measurement point.
  
  for(iteridx in 1:num.iterations){
    
    ## -> Randomly permute the order in which the words appear in a text.
    set.seed(iteridx)  ## For reproducibility.
    random.conc <- conc[sample(1:nrow(conc)), ]
    
    ## -> Following permutation, calculate the values of a given constant for a 
    ## pre-specified number K0 of text lengths.  
    
    ## Get a list of concordances, progressively adding new tokens.per.step 
    ## onto each.
    
    mc.conc.list <- list()
    for(idx in 1:num.steps){
      end.idx <- n.vals[idx]
      mc.conc.list[[idx]] <- random.conc[1:end.idx, ]
    }
    
    ## Make concs into spectra.
    mc.conc.spc <- apply_func_to_subset_list(make_spc, mc.conc.list)
    
    ## Calculate PP for each N in mc.conc.list and save the values
    ## as a row to mc.df.
    PP.vec <- c()
    
    for(idx in 1:length(mc.conc.spc)){
      pp <- Vm(mc.conc.spc[[idx]], 1) / N(mc.conc.spc[[idx]]) ## hapaxes / tokens
      # print(paste("Hapaxes", Vm(mc.conc.spc[[idx]], 1), ". Tokens:", N(mc.conc.spc[[idx]])))
      # print(paste("PP:", pp))
      PP.vec <- c(PP.vec, pp)
    }
    
    mc.df[iteridx, ] <- PP.vec
    message(paste(data.info["suffix.plain"], ": Finished PP iteration", iteridx, "of", num.iterations, "at", Sys.time()))
  }
  
  return(mc.df)
}


montecarlo_pp_means <- function(mc.df){
  ## Calculates the mean for each measurement point (col of mc.df)
  ## and returns it as a vector.
  
  means <- c()
  for(idx in 1:ncol(j)){
    means <- c(means, mean(mc.df[,idx]))
  }
  return(means)
}

montecarlo_pp_ci <- function(mc.df){
  ## Calculates the 95% CI for each measurement point (col of mc.df)
  ## and returns upper and lower bounds in a dataframe.
  ## (Where the 2.5% doesn't divide cleanly, I round the number of elements to exclude
  ## on either side of the interval to the nearest integer.)
  
  mc.cis <- data.frame(matrix(ncol=2, nrow=nrow(mc.df)))
  names(mc.cis) <- c("lower", "upper")
  
  for(idx in 1:ncol(mc.df)){
    elemts_outside_interval <- round2(nrow(mc.df)*0.025, 0)
    lower.idx <- elemts_outside_interval
    upper.idx <- nrow(mc.df) - elemts_outside_interval
    
    ## Sort list in ascending order.
    estimates.sorted <- sort(mc.df[,idx])

    lower.bound <- estimates.sorted[lower.idx]
    upper.bound <- estimates.sorted[upper.idx]
    
    # print(paste("Lower bound at index", lower.idx, "is", lower.bound))
    # print(paste("Upper bound at index", upper.idx, "is", upper.bound))
    # print(lower.bound < upper.bound)
    
    mc.cis$lower[idx] <- lower.bound
    mc.cis$upper[idx] <- upper.bound
  }

  return(mc.cis)
  
}

mc_pp_stats <- function(mc_df, data_info){
  ## Returns a dataframe that contains the step number, the mean, and
  ## the 95% confidence interval for the MC distribution at each 
  ## measurement point.
  
  means.vec <- montecarlo_s_means(mc_df)
  ci.df <- montecarlo_s_ci(mc_df)
  
  mc.pstats.df <- data.frame(matrix(nrow=ncol(mc_df), ncol=5))
  names(mc.pstats.df) <- c("step", "n.seen", "mean", "lower", "upper")
  mc.pstats.df$step <- 1:ncol(mc_df)
  mc.pstats.df$n.seen <- as.numeric(names(mc_df))
  mc.pstats.df$mean <- means.vec
  mc.pstats.df$lower <- ci.df$lower
  mc.pstats.df$upper <- ci.df$upper
  
  saveRDS(mc.pstats.df, 
          paste0("r-analysis/rds/", data_info["suffix.plain"], "_ppmc_stats.rds"))
}

run_ppmc_analysis <- function(data_list, mc_df){
  
  mc_pp_stats(mc_df, data_list)
  
  regex_patt <- paste0("^", data_list["suffix.plain"], "_ppmc.*rds$")

  files <- list.files(pattern = regex_patt, recursive=TRUE)

  for(idx in 1:length(files)){
    var_name <- substr(files[idx], 16, nchar(files[idx])-4)

    assign(var_name,
           readRDS(files[idx]),
           envir = parent.frame())
  }
}
