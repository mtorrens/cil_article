################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source('/no_backup/jferrer/mtorrens/stats/cil_article/source/00_start.R')
# source(paste(SRCDIR, '01_paper_simulations.R', sep = ''))
################################################################################

################################################################################
# Global parameters
force.sim <- FALSE  # Force simulation even if file already exists

# Methods to be tried out
do.bac <- TRUE  # Wang, Dominici & Parmigiani (2012)
do.bma <- TRUE  # Regular BMA
do.ssl <- FALSE # Antonelli et al (2019) -- WARNING: it needs library HDconfounding
do.pcr <- FALSE # Wilson & Reich (2014) -- WARNING: it needs library BayesPen
do.ndl <- TRUE  # Naive and Double Lasso
do.acm <- TRUE  # Wilson et al. (2018) -- ACPME for MT
do.nme <- TRUE  # New method (EP version)
do.nm2 <- TRUE  # New method (EB version)

# Number of cores used in the simulations (override manually)
ncor <- 6
#ncor <- parallel::detectCores() - 1

# Number of simulated datasets on each run
N <- 250
################################################################################

################################################################################
# SETUPS
# Status: (P) Pending; (U) Unfinished; (R) Running; (D) Done
################################################################################
run00 <- TRUE   # (D) 0: [ST] n =  50; p+T =  25; a = 1  ; p_y = p_d = 6
run01 <- TRUE   # (D) 1: [ST] n = 100; p+T =  50; a = 1  ; p_y = p_d = 6
run02 <- TRUE   # (D) 2: [ST] n = 100; p+T =  50; a = 1/3; p_y = p_d = 6
run03 <- TRUE   # (D) 3: [ST] n = 100; p+T =  50; a = 0  ; p_y = p_d = 6
run04 <- FALSE  # (D) 4: [ST] n = 100; p+T = 100; a = 1  ; p_y = p_d = 6
run05 <- FALSE  # (D) 5: [ST] n = 100; p+T = 200; a = 1  ; p_y = p_d = 6
run06 <- FALSE  # (D) 6: [ST] n = 100; p+T = 100; a = 1  ; p_y = p_d = 18
run07 <- FALSE  # (D) 7: [ST] n = 100; p+T = 100; a = 1  ; p_y = p_d = 12
run08 <- FALSE  # (D) 8: [MT] n = 100; p+T = 100; a = 1s ; nts = 2,3,4,5
################################################################################

################################################################################
# SETTING 0: [ST] n = 50; p+T = 25; a = 1; p_y = p_d = 6
################################################################################
term <- '_ST50x24.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run00 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  beta.prior <- 'nlp'
  bma.mprior <- 'bbin'
  only.bacInf <- FALSE
  acpme.onlyDEF <- FALSE

  # Setup
  cat('** RUNNING: SETUP 0\n')
  n <- 50; p <- 24; a <- 1; pstar <- 6
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - length(ones))),
      bd = c(rep(0L, k - 1), ones, rep(0L, p - k + 1 - length(ones))),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, acpme.onlyDEF = acpme.onlyDEF,
      do.bac = do.bac, do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl,
      do.ndl = do.ndl, do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 1: [ST] n = 100; p+T = 50; a = 1; p_y = p_d = 6
################################################################################
term <- '_ST100x49.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run01 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  beta.prior <- 'nlp'
  bma.mprior <- 'bbin'
  only.bacInf <- FALSE
  acpme.onlyDEF <- FALSE

  # Setup
  cat('** RUNNING: SETUP 1\n')
  n <- 100; p <- 49; a <- 1; pstar <- 6
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - length(ones))),
      bd = c(rep(0L, k - 1), ones, rep(0L, p - k + 1 - length(ones))),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, acpme.onlyDEF = acpme.onlyDEF,
      do.bac = do.bac, do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl,
      do.ndl = do.ndl, do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 2: [ST] n = 100; p+T = 50; a = 1/3; p_y = p_d = 6
################################################################################
term <- '_ST100x49a03.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run02 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  beta.prior <- 'nlp'
  bma.mprior <- 'bbin'
  only.bacInf <- FALSE
  acpme.onlyDEF <- FALSE

  # Setup
  cat('** RUNNING: SETUP 2\n')
  n <- 100; p <- 49; a <- 1/3; pstar <- 6
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - length(ones))),
      bd = c(rep(0L, k - 1), ones, rep(0L, p - k + 1 - length(ones))),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, acpme.onlyDEF = acpme.onlyDEF,
      do.bac = do.bac, do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl,
      do.ndl = do.ndl, do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 3: [ST] n = 100; p+T = 50; a = 0; p_y = p_d = 6
################################################################################
term <- '_ST100x49a0.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run03 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  beta.prior <- 'nlp'
  bma.mprior <- 'bbin'
  only.bacInf <- FALSE
  acpme.onlyDEF <- FALSE

  # Setup
  cat('** RUNNING: SETUP 3\n')
  n <- 100; p <- 49; a <- 0; pstar <- 6
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - length(ones))),
      bd = c(rep(0L, k - 1), ones, rep(0L, p - k + 1 - length(ones))),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, acpme.onlyDEF = acpme.onlyDEF,
      do.bac = do.bac, do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl,
      do.ndl = do.ndl, do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 4: [ST] n = 100; p+T = 100; a = 1; p_y = p_d = 6
################################################################################
term <- '_ST100x99.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run04 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  beta.prior <- 'nlp'
  bma.mprior <- 'bbin'
  only.bacInf <- TRUE

  # Setup
  cat('** RUNNING: SETUP 4\n')
  n <- 100; p <- 99; a <- 1; pstar <- 6
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - length(ones))),
      bd = c(rep(0L, k - 1), ones, rep(0L, p - k + 1 - length(ones))),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, do.bac = do.bac,
      do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl, do.ndl = do.ndl,
      do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 5: [ST] n = 100; p+T = 200; a = 1; p_y = p_d = 6
################################################################################
term <- '_ST100x199.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run05 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  beta.prior <- 'nlp'
  bma.mprior <- 'bbin'
  only.bacInf <- TRUE

  # Setup
  cat('** RUNNING: SETUP 5\n')
  n <- 100; p <- 199; a <- 1; pstar <- 6
  ones <- rep(1L, pstar)
  res <- NULL
  #for (k in 7) {
  for (k in 1:(pstar + 1)) {
    # Simulation
    fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
    if (! file.exists(fin) | force.sim == TRUE) {
      now <- paste('[', Sys.time(), ' CET]', sep = '')
      cat(now, ' Batch: ', k, ' of ', pstar + 1, '... ', sep = '')
      sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
        by = c(ones, rep(0L, p - length(ones))),
        bd = c(rep(0L, k - 1), ones, rep(0L, p - k + 1 - length(ones))),
        N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
        max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
        bac.pkg = bac.pkg, only.bacInf = only.bacInf, do.bac = do.bac,
        do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl, do.ndl = do.ndl,
        do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)
  
      # Save partial progress
      fin <- paste(DATDIR, 'psim_', pstar - k + 1, 'of', pstar, term, sep = '')
      save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)
    } else {
      sim <- get(load(file = fin)); cat('Loaded file:', fin, '\n')
    }

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 6: [ST] n = 100; p+T = 100; a = 1; p_y = p_d = 18
################################################################################
term <- '_ST100x99p18.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run06 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  bma.mprior <- 'bbin'
  beta.prior <- 'nlp'
  only.bacInf <- TRUE

  # Setup
  cat('** RUNNING: SETUP 6\n')
  n <- 100; p <- 99; a <- 1; pstar <- 18
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar / 3 + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar / 3 + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - pstar)),
      bd = c(rep(0L, 3 * (k - 1)), ones, rep(0L, p - 3 * (k - 1) - pstar)),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, do.bac = do.bac,
      do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl, do.ndl = do.ndl,
      do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    nc <- pstar - 3 * (k - 1)
    fin <- paste(DATDIR, 'psim_', nc, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bma.mprior = bma.mprior, bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 7: [ST] n = 100; p+T = 100; a = 1; p_y = p_d = 12
################################################################################
term <- '_ST100x99p12.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run07 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  # Specifications
  mod1 <- 'lasso_bic'
  bac.pkg <- 'bacr'
  th.prior <- 'unif'
  bma.mprior <- 'bbin'
  beta.prior <- 'nlp'
  only.bacInf <- TRUE

  # Setup
  cat('** RUNNING: SETUP 7\n')
  n <- 100; p <- 99; a <- 1; pstar <- 12
  ones <- rep(1L, pstar)
  res <- NULL
  for (k in 1:(pstar / 2 + 1)) {
    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: ', k, ' of ', pstar / 2 + 1, '... ', sep = '')
    sim <- sim.teff(a = a, n = n, p = p, S = diag(p), phid = 1, phiy = 1,
      by = c(ones, rep(0L, p - pstar)),
      bd = c(rep(0L, 2 * (k - 1)), ones, rep(0L, p - 2 * (k - 1) - pstar)),
      N = N, ncores = ncor, th.prior = th.prior, mod1 = mod1, R = 1e4,
      max.mod = Inf, bma.mprior = bma.mprior, beta.prior = beta.prior,
      bac.pkg = bac.pkg, only.bacInf = only.bacInf, do.bac = do.bac,
      do.acm = do.acm, do.bma = do.bma, do.ssl = do.ssl, do.ndl = do.ndl,
      do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    nc <- pstar - 2 * (k - 1)
    fin <- paste(DATDIR, 'psim_', nc, 'of', pstar, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, pstar = pstar,
    mod1 = mod1, th.prior = th.prior, beta.prior = beta.prior,
    bac.pkg = bac.pkg, only.bacInf = only.bacInf)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# SETTING 8: [MT] n = 100; p+T = 100; a = 1s; nts = 2:5
################################################################################
term <- '_MTR100x95.RData'
file.out <- paste(DATDIR, 'psim_T', term, sep = '')
if (run08 == TRUE & (! file.exists(file.out) | force.sim == TRUE)) {
  cat('** RUNNING: SETUP 08\n')
  res <- NULL  # Initialise output object
  nts <- 2:5
  for (nt in nts) {
    # Model parameters
    n <- 100
    p <- 95
    a <- rep(1, nt)
    S0 <- diag(p)

    # Simulation
    now <- paste('[', Sys.time(), ' CET]', sep = '')
    cat(now, ' Batch: (nt = ', nt, '); -- ', sep = '')
    sim <- sim.teff.mt(a = a, n = n, p = p, S = diag(p), phiy = 1,
      phid = rep(1, length(a)), N = N, ncores = 6, th.prior = 'unif',
      mod1 = 'lasso_bic', max.mod = Inf, R = 1e4, bma.mprior = 'bbin',
      beta.prior = 'nlp', do.bac = do.acm, do.bma = do.bma, do.ssl = do.ssl,
      do.ndl = do.ndl, do.pcr = do.pcr, do.nme = do.nme, do.nm2 = do.nm2)

    # Save partial progress
    fin <- paste(DATDIR, 'psim_nt', nt, term, sep = '')
    save(sim, file = fin); cat('Saved file:', fin, '\n'); rm(fin)

    # Build output matrix
    res <- rbind.data.frame(res, sim)
  }

  # Output object
  out <- list(res = as.matrix(res), n = n, p = p, a = a, mod1 = 'lasso_bic',
    th.prior = 'unif', beta.prior = 'nlp', bma.mprior = 'bbin', bac.pkg = NA)

  # Save data
  save(out, file = file.out); cat('Saved file:', file.out, '\n')
} else {
  cat('   NOT run.\n')
  #load(file = file.out); cat('Loaded file:', file.out, '\n')
}; rm(file.out)
# END OF SCRIPT

