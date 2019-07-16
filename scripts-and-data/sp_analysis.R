sp_analysis <- function(data.info, conc, times.vector, fzm.extent, norm.corpus){
  
  # ** Descriptive ---------------------------------
  
  ## Make string vector of subperiod year spans (for legends, dataframes etc.)
  times.chars <- get_timespans_as_char_vec(times.vector)
  
  ## Divide the concordance into subperiods using times.vector, save all in a list.
  sp.concs <- get_subset_concs_from_period_vector(times.vector, conc)
  
  ## Get frequency lists, rank/frequency profiles, and frequency spectra 
  ## for each subperiod.
  sp.freqlists <- apply_func_to_subset_list(make_freq_list, sp.concs)
  sp.rankprofiles <- apply_func_to_subset_list(make_rank_freq_profile, sp.freqlists)
  sp.spcs <- apply_func_to_subset_list(make_spc, sp.concs)
  sp.emp.vgc <- apply_func_to_subset_list(make_empirical_growth_df, sp.concs)
  
  ## Bind all of these into one single dataframe each for export. Setting the names of list 
  ## elements first means that these names will be used as the identifiers in the new
  ## "period" column of the tall dataframe.
  names(sp.concs) <- times.chars
  names(sp.freqlists) <- times.chars
  names(sp.rankprofiles) <- times.chars
  names(sp.spcs) <- times.chars
  names(sp.emp.vgc) <- times.chars
  
  sp.concs.tall <- bind_rows(sp.concs, .id="period")
  sp.spcs.tall <- bind_rows(sp.spcs, .id="period")
  sp.rankprofiles.tall <- bind_rows(sp.rankprofiles, .id="period")
  sp.freqlists.tall <- bind_rows(sp.freqlists, .id="period")
  sp.emp.vgc.tall <- bind_rows(sp.emp.vgc, .id="period")
  
  ## Save info about each subperiod (types, tokens, hapaxes)
  ## in one info dataframe for all subperiods (to be exported later).
  sp.info <- make_sp_info_df(sp.spcs, times.chars)

  # ** Comparable potential productivity ----------------------
  ## Method 1: Strip away all data beyond smallest token size.
  
  sp.mintoken.concs <- strip_subsets_to_min_token(sp.info, sp.concs)
  names(sp.mintoken.concs) <- times.chars
  sp.mintoken.pp <- pp_of_equally_sized_conc_list(sp.mintoken.concs)
  
  # ** Interpolated VGCs ---------------------------------
  
  ## Make interpolated VGC for each subperiod.
  sp.bin.vgcs <- make_spcs_into_bin_vgcs(sp.spcs)
  ## TO-DO: Make tall for plotting.
  
  # ** Default fZMs ---------------------------------
  
  ## Make fZM model for each subperiod. If fZM model estimation fails and returns an error
  ## (e.g. if zipfR cannot estimate any parameters for a particular subperiod), the 
  ## corresponding list element will be NA.
  sp.fzms.exact <- make_sp_fzms(sp.spcs)
  sp.fzms.approx <- make_sp_fzms_exactf(sp.spcs)
  
  ## Add S from each of these models to sp.info.

  for(idx in 1:length(sp.fzms.exact)){
    if(!all(is.na(sp.fzms.exact[[idx]]))){
      sp.info$S.exact[idx] <- sp.fzms.exact[[idx]]$S
    }else{
      sp.info$S.exact[idx] <- NA
    }
    if(!all(is.na(sp.fzms.approx[[idx]]))){
      sp.info$S.approx[idx] <- sp.fzms.approx[[idx]]$S
    }else{
      sp.info$S.approx[idx] <- NA
    }
  }
  
  ## Create VGCs for fZM models of subperiods.
  ## Option: Instead of having the upper limit of extrapolation as set above in vgc.extent,
  ## the second argument here could also be the maximum sample size among all subsets, 
  ## i.e. max(sp.info$tokens)
  sp.fzm.exact.vgcs <- make_sp_fzm_vgcs(sp.fzms.exact, fzm.extent)
  sp.fzm.approx.vgcs <- make_sp_fzm_vgcs(sp.fzms.approx, fzm.extent)
  
  
  ## TO-DO: Make tall for plotting.
  # hkeit_sp_fzm_vgc_nona <- hkeit_sp_fzm_vgc[!is.na(hkeit_sp_fzm_vgc)]
  # hkeit_sp_fzm_vgc_tall <- bind_rows(hkeit_sp_fzm_vgc, .id="period")
  # ung_sp_fzm_vgc_tall <- bind_rows(ung_sp_fzm_vgc, .id="period")
  
  # ** All fZMs ---------------------------------
  
  ## Generate all fZMs with all options for each subperiod.
  # sp.fzms.all <- gen_sp_fzm_options(sp.spcs)
  
  ## Save result to dataframe. (If parameter estimation fails, NA inserted)
  # sp.fzms.all.info <- make_sp_fzm_info_df(sp.fzms.all, times.chars)
  
  
  # and VGCs for each fzm?
  ## TO-DO
  
  
  # ** ROA by subperiod ---------------------------------
  
  ## The first row of this dataframe contains the initial lexicon.
  sp.roa <- rate_of_addition(sp.rankprofiles)
  
  # ** Save one by one ---------------------------------
  
  saveRDS(sp.concs, "r-analysis/rds/sp_conc.rds")
  saveRDS(sp.concs.tall, "r-analysis/rds/sp_conc_tall.rds")
  saveRDS(sp.emp.vgc, "r-analysis/rds/sp_emp_vgc.rds")
  saveRDS(sp.emp.vgc.tall, "r-analysis/rds/sp_emp_vgc_tall.rds")
  saveRDS(sp.freqlists, "r-analysis/rds/sp_freqlist.rds")
  saveRDS(sp.freqlists.tall, "r-analysis/rds/sp_freqlist_tall.rds")
  saveRDS(sp.info, "r-analysis/rds/sp_info.rds")
  saveRDS(sp.rankprofiles, "r-analysis/rds/sp_rankprof.rds")
  saveRDS(sp.rankprofiles.tall, "r-analysis/rds/sp_rankprof_tall.rds")
  saveRDS(sp.spcs, "r-analysis/rds/sp_spc.rds")
  saveRDS(sp.spcs.tall, "r-analysis/rds/sp_spc_tall.rds")
  saveRDS(sp.mintoken.concs, "r-analysis/rds/sp_mintoken_concs.rds")
  saveRDS(sp.mintoken.pp, "r-analysis/rds/sp_mintoken_pp.rds")
  saveRDS(sp.bin.vgcs, "r-analysis/rds/sp_bin_vgc.rds")
  saveRDS(sp.fzms.exact, "r-analysis/rds/sp_fzm_exact.rds")
  saveRDS(sp.fzms.approx, "r-analysis/rds/sp_fzm_approx.rds")
  saveRDS(sp.fzm.exact.vgcs, "r-analysis/rds/sp_fzm_exact_vgc.rds")
  saveRDS(sp.fzm.approx.vgcs, "r-analysis/rds/sp_fzm_approx_vgc.rds")
  saveRDS(sp.roa, "r-analysis/rds/sp_roa.rds")
}

run_sp_analysis <- function(data.info, conc, times.vector, fzm.extent, norm.corpus){
  
  sp_analysis(data.info, conc, times.vector, fzm.extent, norm.corpus)
  
  ## Find all rds files generated by the conc_analysis function.
  files <- list.files(pattern ="^sp_.*rds$", recursive = TRUE)
  
  ## Load and save them all with the prefix of the respective suffix.
  for(idx in 1:length(files)){
    assign(paste0(data.info["suffix.plain"], "_", 
                  substr(files[idx], 16, nchar(files[idx])-4)), 
           readRDS(files[idx]),
           envir = parent.frame())
  }
}