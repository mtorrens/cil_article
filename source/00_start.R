################################################################################
# Script: 00_start.R
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
################################################################################

################################################################################
# Setup
################################################################################
#Silent mode (deactivated by default)
if (! 'silent' %in% ls()) { silent <- FALSE }  # Set TRUE to run silently
if (! 'runALL' %in% ls()) { runALL <- FALSE }  # Set TRUE to run ALL scripts

#Version control
vN <- 'v16'
vR <- 'xvi'

# Initialize project
if (silent == FALSE) {
  cat('** Initialising...
* Project:  Treatment Effect Estimation with Confounder Importance Learning
* Author:   Miquel Torrens (c)
* Version:  ', vR, ' (September 2021)
* Packages: mvtnorm, parallel, doMC, glmnet, mombf, hdm, BACprior, bacr,
            selectiveInference, statmod, compiler, regimes [NON-CRAN],
            pracma [NOT LOADED], plotly, ipumsr\n', sep = '')
}

#Paths
PATH <- paste('~/Desktop/year3/bma_teff/', vN, '/', sep = '')  #Change to LOCAL
SRCDIR <- paste(PATH, 'syntax/', sep = '')
DATDIR <- paste(PATH, 'data/', sep = '')
OUTDIR <- paste(PATH, 'output/', sep = '')
TMPDIR <- paste(PATH, 'temp/', sep = '')
FIGDIR <- paste(PATH, 'figures/', sep = '')
INPDIR <- paste(PATH, 'input/', sep = '')
FCNDIR <- paste(SRCDIR, 'functions/', sep = '')

#Create directories if they don't exist
if (! dir.exists(DATDIR)) {
  dir.create(DATDIR); cat('Created directory:', DATDIR, '\n')
}
if (! dir.exists(OUTDIR)) {
  dir.create(OUTDIR); cat('Created directory:', OUTDIR, '\n')
}
if (! dir.exists(TMPDIR)) {
  dir.create(TMPDIR); cat('Created directory:', TMPDIR, '\n')
}
if (! dir.exists(FIGDIR)) {
  dir.create(FIGDIR); cat('Created directory:', FIGDIR, '\n')
}
if (! dir.exists(INPDIR)) {
  dir.create(INPDIR); cat('Created directory:', INPDIR, '\n')
}

################################################################################
# Dependencies
################################################################################
req <- c()
req[01] <- require('mvtnorm')
req[02] <- require('parallel')
req[03] <- require('doMC')
req[04] <- require('mombf')  # BMA
req[05] <- require('glmnet')  # Lasso and PL generally
req[06] <- require('selectiveInference')  # Lasso with exact post-sel. inference
req[07] <- require('hdm')  # Double Lasso (Belloni et al., 2014)
req[08] <- require('bacr')  #BAC (original from Wang) by Wang et al. (2012)
req[09] <- require('BACprior')  # BAC based on Lefebvre et al. (2014)
req[10] <- require('statmod')  # Required by HDconfounding
req[11] <- require('regimes')  # ACPME by Wilson et al. (2018) [NON-CRAN]
req[12] <- require('compiler')  # Compile functions
req[13] <- require('plotly')  # Contour/level plots
req[14] <- require('ipumsr')  # CPS data loading
#req[15] <- require('pracma')  #For pseudo-inverse computation
if (any(! req == TRUE)) {
  warning('(!) NOT all packages could be loaded; review dependency list.')
  warning('Non-installed required packages need to be installed manually.')
}; rm(req)

#Where to find the NON-CRAN packages
# HDconfounding: available at https://github.com/jantonelli111/HDconfounding
#BayesPen: available at https://github.com/AnderWilson/BayesPen
#regimes: available at https://github.com/AnderWilson/regimes

#Specific functions
source(paste(FCNDIR, 'functions_butler.R', sep = ''))
source(paste(FCNDIR, 'functions_gradients.R', sep = ''))
source(paste(FCNDIR, 'functions_newmethod.R', sep = ''))
source(paste(FCNDIR, 'support_bac.R', sep = ''))
source(paste(FCNDIR, 'functions_simulator.R', sep = ''))
source(paste(FCNDIR, 'functions_plot.R', sep = ''))
source(paste(FCNDIR, 'functions_cps.R', sep = ''))

# Finish
if (silent == FALSE) {
  cat('* Project PATH:', PATH, '\n** Ready!\n')
}

################################################################################
# Scripts
################################################################################
# Run every script (for reproducibility)
if (runALL == TRUE) {
  cat('** RUNNING ALL scripts:\n')
  cat('** SCRIPT: 01_paper_simulations.R\n')
  source(paste(SRCDIR, '01_paper_simulations.R', sep = ''))
  cat('** SCRIPT: 02_paper_figures.R\n')
  source(paste(SRCDIR, '02_paper_figures.R', sep = ''))
  cat('** SCRIPT: 03_generic_illustrations.R\n')
  source(paste(SRCDIR, '03_generic_illustrations.R', sep = ''))
  cat('** SCRIPT: 04a_cps_format.R\n')
  source(paste(SRCDIR, '04a_cps_format.R', sep = ''))
  cat('** SCRIPT: 04b_cps_transform.R\n')
  source(paste(SRCDIR, '04b_cps_transform.R', sep = ''))
  cat('** SCRIPT: 05_cps_analysis.R\n')
  source(paste(SRCDIR, '05_cps_analysis.R', sep = ''))
  cat('** SCRIPT: 06_cps_figures.R\n')
  source(paste(SRCDIR, '06_cps_figures.R', sep = ''))
}
# END OF SCRIPT
