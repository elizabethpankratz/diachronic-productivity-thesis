## Elizabeth Pankratz, July 2019

## Primary script for the productivity analyses used for my M.A. thesis.
## The samples analysed by this script are found in child directory "r-analysis"
## and the csv output of this script is also saved there in various subdirectories.

## Initial set-up ----

library(ggplot2)
library(extrafont)
# loadfonts()
library(zipfR)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(RColorBrewer)

source("functions.R")
source("conc_analysis.R")
source("sp_analysis.R")
source("doc_analysis.R")
source("montecarlo.R")


## Set working dir in current dir if sourcing script:
# setwd(getSrcDirectory()[1])
setwd(dirname(parent.frame(2)$ofile))

## Set working dir in current dir if running lines of script:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## Set up and read in datasets.
## From RIDGES:
data.ung <- c(csv = "r-analysis/ung.csv", 
              suffix = "-ung", 
              suffix.plain = "ung")
data.hkeit <- c(csv = "r-analysis/hkeit.csv", 
                suffix = "-heit/-keit", 
                suffix.plain = "hkeit")
data.er <- c(csv = "r-analysis/er.csv", 
             suffix = "-er", 
             suffix.plain = "er")

conc.hkeit = read.csv2(data.hkeit["csv"],header=T,sep=",")
conc.er = read.csv2(data.er["csv"],header=T,sep=",")
conc.ung = read.csv2(data.ung["csv"],header=T,sep=",")

doc.tokencounts <- read.csv2("r-analysis/corpus-size-by-title.csv", header=T, sep=",")

## From DECOW16B:
data.hkeit.decow <- c(csv = "r-analysis/hkeit_decow.csv",
                      suffix = "-heit/-keit (DECOW)",
                      suffix.plain = "hkeitDECOW")
conc.hkeit.decow <- read.csv2(data.hkeit.decow["csv"],header=T,sep=",")


## Choose subperiods for the sp-based analysis ----

## This is Moser's (1926) ENHG-internal divisions, plus a final stage for earlier Modern German.
# times.vector <- c(1482, 1520, 1620, 1760, 1915)

## Two divisions: Late ENHG, earlier modern German.
# times.vector <- c(1482, 1650, 1915)

## Approximately 100-year divisions.
times.vector <- c(1482, 1550, 1650, 1750, 1850, 1915)

## Approximately 50-year divisions.
# times.vector <- c(1482, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1915)


## Set constants ---------------------------------

## To which value of N should the extrapolated model go?
vgc.extent <- 3000

## To which number of tokens should the type/token/hapax counts be normalised?
norm.corpus.size <- 10000


# General analyses ---------------------------------

## Executing each call to run_x_analysis will run all of 
## the analyses contained
## in the x_analysis() function defined in each document
## and saves the output
## to the workspace (global environment).

## Entire concordance:

run_conc_analysis(data.hkeit, conc.hkeit, vgc.extent, 4)
run_conc_analysis(data.er, conc.er, vgc.extent, 4)
run_conc_analysis(data.ung, conc.ung, vgc.extent, 4)

# Subperiods:

run_sp_analysis(data.hkeit, conc.hkeit, times.vector, vgc.extent, norm.corpus.size)
run_sp_analysis(data.er, conc.er, times.vector, vgc.extent, norm.corpus.size)
run_sp_analysis(data.ung, conc.ung, times.vector, vgc.extent, norm.corpus.size)

# Per document:

run_doc_analysis(data.hkeit, conc.hkeit, doc.tokencounts, norm.corpus.size)
run_doc_analysis(data.er, conc.er, doc.tokencounts, norm.corpus.size)
run_doc_analysis(data.ung, conc.ung, doc.tokencounts, norm.corpus.size)


# Monte Carlo S ----

## Running montecarlo_s_df() takes a few hours each time, so we can also just 
## read in the saved RDS file of that analysis (which is a 20 x 2000 dataframe 
## containing the MC distribution of S at the 20 measurement points).


## MC distributions for ung: 20 steps x 2000 iterations
### ung_mc_approx <- montecarlo_s_df(data.ung, conc.ung, 20, 2000, F)
### saveRDS(ung_mc_approx, "r-analysis/rds/ung_mc_approx")
ung_mc_approx <- readRDS("r-analysis/rds/ung_mc_approx.rds")

## MC distributions for hkeit: 20 steps x 2000 iterations
### hkeit_mc_approx <- montecarlo_s_df(data.hkeit, conc.hkeit, 20, 2000, F)
### saveRDS(hkeit_mc_approx, "r-analysis/rds/hkeit_mc_approx.rds")
hkeit_mc_approx <- readRDS("r-analysis/rds/hkeit_mc_approx.rds")

## MC distributions for er: 20 steps x 2000 iterations
### er_mc_approx <- montecarlo_s_df(data.er, conc.er, 20, 2000, F)
### saveRDS(er_mc_approx, "r-analysis/rds/er_mc_approx.rds")
er_mc_approx <- readRDS("r-analysis/rds/er_mc_approx.rds")

run_mc_analysis2(data.ung, ung_mc_approx, F)
run_mc_analysis2(data.er, er_mc_approx, F)
run_mc_analysis2(data.hkeit, hkeit_mc_approx, F)


## MC distribution for large hkeit sample from DECOW: 100 steps x 2000 iterations.

### hkeit_decow_mc_s_100 <- montecarlo_s_df(data.hkeit.decow, conc.hkeit.decow, 100, 2000, F)
### saveRDS(hkeit_decow_mc_s_100, "r-analysis/rds/hkeit_decow_mc_approx_100.rds")
hkeit_decow_mc_s_100 <- readRDS("r-analysis/rds/hkeit_decow_mc_approx_100.rds")
run_mc_analysis2(data.hkeit.decow, hkeit_decow_mc_s_100, F)


# Monte Carlo PP ----

## Same as above, montecarlo_pp_df() takes a while to run, so it's quicker to read in
## the existing RDSs.

### hkeit_ppmc_data <- montecarlo_pp_df(data.hkeit, conc.hkeit, 20, 2000)
### saveRDS(hkeit_ppmc_data, "r-analysis/rds/hkeit_ppmc_data.rds")
hkeit_ppmc_data <- readRDS("r-analysis/rds/hkeit_ppmc_data.rds")

### er_ppmc_data <- montecarlo_pp_df(data.er, conc.er, 20, 2000)
### saveRDS(er_ppmc_data, "r-analysis/rds/er_ppmc_data.rds")
er_ppmc_data <- readRDS("r-analysis/rds/er_ppmc_data.rds")

### ung_ppmc_data <- montecarlo_pp_df(data.ung, conc.ung, 20, 2000)
### saveRDS(ung_ppmc_data, "r-analysis/rds/ung_ppmc_data.rds")
ung_ppmc_data <- readRDS("r-analysis/rds/ung_ppmc_data.rds")

run_ppmc_analysis(data.er, er_ppmc_data)
run_ppmc_analysis(data.hkeit, hkeit_ppmc_data)
run_ppmc_analysis(data.ung, ung_ppmc_data)


# Save results -------------------

save_info_as_csv()
save_spcs_as_csv()
save_fzm_vgcs_as_csv()
save_mc_s_as_csv()
save_mc_pp_as_csv()
save_emp_vgcs_as_csv()
save_bin_vgcs_as_csv()
save_freqlists_as_csv()
save_rankprofiles_as_csv()
