################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(FCNDIR, 'functions_cps.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * chop
# * clean.names
# * conv.dummy
# * gen.fake.preds
# * cps.methods
# * sample.expabsh
# * plot.2010vs2019
# * min.twoway.boxplot
################################################################################

################################################################################
# Fold assignment function
chop <- function(x, k) {
################################################################################
  rank(x) %% k + 1
}

################################################################################
# Extract strange characters from column names
clean.names <- function(x) {
################################################################################
  x <- tolower(x)
  x <- gsub('-', '_', x)
  x <- gsub('\\(', '_', x)
  x <- gsub(')', '_', x)
  x <- gsub('>', '_', x)
  x <- gsub('<', '_', x)
  x <- gsub(',', '_', x)
  x <- gsub('/', '_', x)
  x <- gsub('__', '_', x)
  return(x)
}

################################################################################
# Convert categorical variable into matrix of dummy variables
conv.dummy <- function(data, colname, refcat = NULL, ret.df = FALSE) {
################################################################################
  if (is.null(ncol(data))) {
    data <- data.frame(cat = data)
    colname <- 'cat'
  }   
  cats <- unique(c(refcat, as.character(data[, colname])))
  ncats <- length(cats)
  dummat <- matrix(NA, nrow = nrow(data), ncol = ncats)
  for (k in 1:ncats) {
    dummat[, k] <- as.numeric(data[, colname] == cats[k])
  }
  colnames(dummat) <- gsub(' ', '', paste(colname, cats, sep = '_'))
  if (! is.null(refcat)) {
    dummat <- dummat[, -1]
  }
  if (ret.df == TRUE) {
    dummat <- as.data.frame(dummat)
  }
  return(dummat)
}

##############################################################################
# Generate fake predictors for the CPS dataset
gen.fake.preds <- function(nfr, centre, D, seed = 666) {
##############################################################################
  # Fake predictors
  #nfr <- 25  # No. of fake regressors per treatment
  noms <- c()
  for (i in 1:ncol(D)) {
    noms <- c(noms, paste('FR', sprintf('%02.0f', 1:nfr), 'T', i, sep = ''))
  }

  # Simulate fake regressors
  set.seed(seed)
  M <- matrix(nrow = nrow(D), ncol = nfr * ncol(D))
  colnames(M) <- noms
  for (j in 1:ncol(M)) {
    idt <- as.numeric(substr(noms[j], nchar(noms[j]), nchar(noms[j])))
    out <- D[, idt]
    M[which(out == 0), j] <- rnorm(sum(out == 0), mean = -centre, sd = 1)
    M[which(out == 1), j] <- rnorm(sum(out == 1), mean = +centre, sd = 1)
  }

  # End
  return(invisible(M))
}

################################################################################
# Run various compared methods on the CPS dataset
cps.methods <- function(y, D, X, I = NULL, minimal = FALSE, eta = 0.05,
  lpen = 'lambda.min', do.acm = TRUE, out.post = FALSE) {
################################################################################
  # Global parameters
  pr.eta <- c(eta / 2, 1 - eta / 2)
  noms <- paste(100 * pr.eta, '%', sep = '')
  if (is.null(colnames(X))) { colnames(X) <- paste('V', 1:ncol(X), sep = '') }

  # Number of treatments
  nt <- ncol(D)
  ni <- ifelse(is.null(I), 0, ncol(I))
  nc <- nt + ni

  # Adapt design matrix
  A <- cbind(D, X[, -1])
  o <- colnames(A)[unname(which(is.na(coef(lm(y ~ A))))) - 1]  # Collinearities
  # if (is.null(o)) { o <- NA }
  X <- X[, which(! colnames(X) %in% o)]  # Train
  Z <- cbind(D, I, X)  # Train
  n <- nrow(X)  # Sample size

  # Regular OLS ################################################################
  Zi <- cbind(D, I, X[, -1])
  ols.fit <- lm(y ~ Zi)
  ols.sum <- summary(ols.fit)

  # Estimates
  ah.ols <- unname(coef(ols.fit)[1 + 1:nc])
  c1.ols <- unname(confint(ols.fit)[1 + 1:nc, 1])
  c2.ols <- unname(confint(ols.fit)[1 + 1:nc, 2])
  il.ols <- c2.ols - c1.ols
  pd.ols <- unname(ols.sum[['coefficients']][1 + 1:nc, 4])
  nr.ols <- sum(ols.sum[['coefficients']][, 4] < eta)

  # Regular LASSO ##############################################################
  ah.nle <- c1.nle <- c2.nle <- il.nle <- pd.nle <- rep(NA, nc)
  nr.nle <- NA
  if (minimal == FALSE) {
    nl.fit <- cv.glmnet(Zi, y, alpha = 1, intercept = TRUE, standardize = FALSE)

    # Estimates
    ah.nle <- coef(nl.fit, s = nl.fit[[lpen]])[1 + 1:nc]
    c1.nle <- rep(NA, nc)
    c2.nle <- rep(NA, nc)
    il.nle <- rep(NA, nc)
    pd.nle <- ifelse(ah.nle != 0, 1L, 0L)
    nr.nle <- sum(coef(nl.fit, s = nl.fit[[lpen]]) != 0)
  }

  # Double Lasso ###############################################################
  dl.fit <- rlassoEffects(y = y, x = Z, index = 1:nc, post = TRUE,
    intercept = TRUE, method = 'double selection')

  # Estimates
  ah.dle <- unname(coef(dl.fit)[1:nc])
  c1.dle <- unname(confint(dl.fit))[, 1]
  c2.dle <- unname(confint(dl.fit))[, 2]
  il.dle <- c2.dle - c1.dle
  pd.dle <- unname(summary(dl.fit)[['coefficients']][1:nc, 4])
  nr.dle <- sum(dl.fit[['selection.matrix']][, 1], na.rm = TRUE)

  # ACPME ######################################################################
  ah.bac <- c1.bac <- c2.bac <- il.bac <- pd.bac <- rep(NA, nc)
  nr.bac <- NA
  if (do.acm == TRUE) {
  #if (do.acm == TRUE & minimal == FALSE) {
    old.names <- colnames(X)
    colnames(X) <- clean.names(colnames(X))
    tmp.file <- paste(TMPDIR, 'tmp_acpme_', Sys.time(), '.txt', sep = '')
    sink(file = tmp.file)
    acm.fit <- try(acpme(y = y, Z = D, C = X[, -1], niter = 1e4), silent = FALSE)
    #acm.fit <- try(acpme(y = y, Z = D, C = X[, -1], niter = 1500), silent = T)
    sink(); unlink(tmp.file)
    colnames(X) <- old.names

    # Estimates
    if (class(acm.fit) != 'try-error') {
      acm.coef <- acm.fit[['beta']]
      ah.bac <- unname(colMeans(acm.coef)[1:nc])
      c1.bac <- unname(apply(acm.coef[, 1:nc], 2, quantile, probs = pr.eta[1]))
      c2.bac <- unname(apply(acm.coef[, 1:nc], 2, quantile, probs = pr.eta[2]))
      il.bac <- c2.bac - c1.bac
      pd.bac <- ifelse(sign(c1.bac) == sign(c2.bac), 1, 0)
      nr.bac <- mean(rowSums(acm.fit[['alpha']]))
    }
  }

  # PCR ########################################################################
  ah.pcr <- c1.pcr <- c2.pcr <- il.pcr <- pd.pcr <- rep(NA, nc)
  nr.pcr <- NA
  if (minimal == FALSE) {
    tmp.file <- paste(TMPDIR, 'tmp_pcr_', Sys.time(), '.txt', sep = '')
    sink(file = tmp.file)
    pcr.fit <- BayesPen.lm.confounders(y, D, X)
    sink(); unlink(tmp.file)

    # Choose min BIC
    prpcr <- mean(y) + cbind(D, X) %*% t(pcr.fit[['coefs']])
    nC <- apply(pcr.fit[['coefs']], 1, function(x) sum(x != 0))
    nr <- nrow(X)
    #y2 <- matrix(NA, nrow = length(y), ncol = ncol(prpcr))
    #for (j in 1:ncol(y2)) { y2[, j] <- y }
    y2 <- y
    bic <- nr * log(colSums((y2 - prpcr)^2) / length(y)) +
           nr * (log(2 * pi) + 1) + log(nr) * nC
    pcr.coef <- pcr.fit[['coefs']][which.min(bic), ]

    # Estimates
    ah.pcr <- unname(pcr.coef[1:nc])
    c1.pcr <- rep(NA, nc)
    c2.pcr <- rep(NA, nc)
    il.pcr <- rep(NA, nc)
    pd.pcr <- ifelse(ah.pcr != 0, 1L, 0L)
    nr.pcr <- sum(pcr.coef != 0)
  }

  # BMA ########################################################################
  ms <- mombf::modelSelection
  igp <- mombf::igprior(alpha = 0.01, lambda = 0.01)  # IGam(0.01, 0.01)
  nlp <- mombf::momprior(tau = 0.348)  # MOM(0.348)
  bbp <- mombf::modelbbprior(alpha.p = 1, beta.p = 1)  # BetaBin(1, 1)

  # Fit
  bms.fit <- ms(y, Z, priorCoef = nlp, priorDelta = bbp, priorVar = igp,
    niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)
  gibbs <- mombf::rnlp(msfit = bms.fit, priorCoef = nlp, niter = 1e4)
  pprob <- mombf::postProb(bms.fit)

  # Estimates
  ah.bma <- unname(colMeans(gibbs[, 1 + 1:nc]))
  c1.bma <- unname(apply(gibbs[, 1 + 1:nc], 2, quantile, probs = pr.eta[1]))
  c2.bma <- unname(apply(gibbs[, 1 + 1:nc], 2, quantile, probs = pr.eta[2]))
  il.bma <- c2.bma - c1.bma
  pd.bma <- unname(bms.fit[['margpp']][1:nc])
  nr.bma <- sum(unlist(lapply(
    strsplit(as.character(pprob[, 1]), ','), length)) * pprob[, 3])
  ps.bma <- NA
  if (out.post == TRUE) { ps.bma <- bms.fit }

  # New method (EP) ############################################################
  nmP.fit <- try(new.method.mt(y, D, X, I, th.search = 'EP', mod1 = 'lasso_bic',
    th.prior = 'unif', beta.prior = 'nlp', th.range = seq(-9, 9, 1.5),
    rho.min = NULL, rho.max = NULL, R = 1e4, max.mod = 5e3), silent = TRUE)

  # Estimates
  if (class(nmP.fit) != 'try-error') {
    aux <- nmP.fit[['pprob']]
    ah.nmP <- unname(nmP.fit[['teff.nlp']][1:nc])
    c1.nmP <- unname(nmP.fit[['teff.d.nlp']][noms, ])[1, 1:nc]
    c2.nmP <- unname(nmP.fit[['teff.d.nlp']][noms, ])[2, 1:nc]
    il.nmP <- c2.nmP - c1.nmP
    pd.nmP <- unname(nmP.fit[['mpprob']][1:nc])
    nr.nmP <- sum(unlist(lapply(
      strsplit(as.character(aux[, 1]), ','), length)) * aux[, 2])
    m1.nmP <- nmP.fit[['mod1.coef']]
    ps.nmP <- NA
    if (out.post == TRUE) { ps.nmP <- nmP.fit }
  } else {
    ah.nmP <- c1.nmP <- c2.nmP <- il.nmP <- pd.nmP <- rep(NA, nc)
    nr.nmP <- m1.nmP <- NA
  }

  # New method (EB) ############################################################
  ah.nmB <- c1.nmB <- c2.nmB <- il.nmB <- pd.nmB <- rep(NA, nt)
  nr.nmB <- m1.nmB <- ps.nmB <- NA
  if (minimal == FALSE) {
    nmB.fit <- try(new.method.mt(y, D, X, I, th.search = 'EB',
      mod1 = 'lasso_bic', th.prior = 'unif', beta.prior = 'nlp', rho.min = NULL,
      rho.max = NULL, R = 1e4, max.mod = 5e3,
      bvs.fit0 = nmP.fit[['init.msfit']], th.EP = nmP.fit[['th.hat']]),
      silent = TRUE)

    # Estimates
    if (class(nmP.fit) != 'try-error') {
      aux <- nmB.fit[['pprob']]
      ah.nmB <- unname(nmB.fit[['teff.nlp']][1:nt])
      c1.nmB <- unname(nmB.fit[['teff.d.nlp']][noms, ])[1, 1:nt]
      c2.nmB <- unname(nmB.fit[['teff.d.nlp']][noms, ])[2, 1:nt]
      il.nmB <- c2.nmB - c1.nmB
      pd.nmB <- unname(nmB.fit[['mpprob']][1:nt])
      nr.nmB <- sum(unlist(lapply(
        strsplit(as.character(aux[, 1]), ','), length)) * aux[, 2])
      m1.nmB <- nmB.fit[['mod1.coef']]
      ps.nmB <- NA
      if (out.post == TRUE) { ps.nmB <- nmB.fit }
    }
  }

  # Finish
  return(list(ah.ols = ah.ols, c1.ols = c1.ols, c2.ols = c2.ols,
              il.ols = il.ols, pd.ols = pd.ols, nr.ols = nr.ols,
              ah.nle = ah.nle, c1.nle = c1.nle, c2.nle = c2.nle,
              il.nle = il.nle, pd.nle = pd.nle, nr.nle = nr.nle,
              ah.dle = ah.dle, c1.dle = c1.dle, c2.dle = c2.dle,
              il.dle = il.dle, pd.dle = pd.dle, nr.dle = nr.dle,
              ah.bac = ah.bac, c1.bac = c1.bac, c2.bac = c2.bac,
              il.bac = il.bac, pd.bac = pd.bac, nr.bac = nr.bac,
              ah.pcr = ah.pcr, c1.pcr = c1.pcr, c2.pcr = c2.pcr,
              il.pcr = il.pcr, pd.pcr = pd.pcr, nr.pcr = nr.pcr,
              ah.bma = ah.bma, c1.bma = c1.bma, c2.bma = c2.bma,
              il.bma = il.bma, pd.bma = pd.bma, nr.bma = nr.bma,
              ah.nmP = ah.nmP, c1.nmP = c1.nmP, c2.nmP = c2.nmP,
              il.nmP = il.nmP, pd.nmP = pd.nmP, nr.nmP = nr.nmP,
              ah.nmB = ah.nmB, c1.nmB = c1.nmB, c2.nmB = c2.nmB,
              il.nmB = il.nmB, pd.nmB = pd.nmB, nr.nmB = nr.nmB,
              ps.bma = ps.bma, ps.nmP = ps.nmP, ps.nmB = ps.nmB,
              m1.nmP = m1.nmP, m1.nmB = m1.nmB))
}

################################################################################
# Compute contribution to deviation from average salary (given alpha)
sample.expabsh <- function(D1, Dp, a, abs.val = TRUE) {
################################################################################
  obs <- sample(1:nrow(D1), 1)  # Draw a random observations of the observed
  Db <- Dp[obs, ]  # Predicted value of treatments
  Di <- D1[obs, ]  # Actual observed treatment values
  ev <- as.numeric(t(Di - Db) %*% a)  # Measure of interest
  if (abs.val == TRUE) { ev <- abs(ev) }
  return(exp(ev))
}

################################################################################
# Plot for FIGURE 2 and FIGURE S3: compare methods in 2010 and 2019
plot.2010vs2019 <- function(obj, treat, ymin, ymax, xlab.in = TRUE,
  year.tag = TRUE, year.pos = NA) {
################################################################################
  # Extract objects (must be four)
  J <- obj[[1]][c(1, 2), ]#[c(1, 10), ]
  K <- obj[[2]][c(1, 2), ]#[c(1, 10), ]
  L <- obj[[3]][c(1, 2), ]#[c(1, 10), ]

  # List of methods
  mths <- c('ols', 'dle', 'bma', 'nmP')

  # Pick treatment
  ntr <- switch(treat,
    '1' = 'gender: female',
    '2' = 'race: black',
    '3' = 'ethnicity: hispanic',
    '4' = 'birthplace: Latin America')

  # Plot parameters
  xmin <- 0; xmax <- 2
  xlabs <- c('OLS', 'DML', 'BMA', 'CIL')
  ylabs <- ntr
  #ylabs <- paste('Est. effect:', ntr)
  toplab1 <- c(expression(bold('2010')), expression(bold('2019')))
  posxl <- seq(0.05, 0.85, 0.8 / (length(xlabs) - 1))
  # Positions of labels
  cutx1 <- c(0.19, 0.46, 0.73)
  cutx2 <- c(1.27, 1.54, 1.81)
  posxl1 <- c(0.055, 0.325, 0.595, 0.865)
  posxl2 <- c(1.135, 1.405, 1.675, 1.945)
  # Positions of points / segments
  posp11 <- posxl1 - 0.06
  posp12 <- posxl1
  posp13 <- posxl1 + 0.06
  posp21 <- posxl2 - 0.06
  posp22 <- posxl2
  posp23 <- posxl2 + 0.06

  # Plot frame
  plot(NA,
    ylim = c(ymin, ymax), ylab = ylabs,
    xlim = c(xmin, xmax), xaxt = 'n', xlab = '')
  if (xlab.in == TRUE) {
    axis(1, posxl1, xlabs, tick = FALSE, cex.axis = 1)#, line = 0.2)
    axis(1, posxl2, xlabs, tick = FALSE, cex.axis = 1)#, line = 0.2)
  }
  for (k in 1:3) { abline(v = cutx1[k], lty = 2) }
  for (k in 1:3) { abline(v = cutx2[k], lty = 2) }
  for (l in seq(-0.5, 0.5, 0.05)) { abline(h = l, lty = 3, col = 'gray') }
  abline(h = 0, lty = 3)
  abline(v = xmax * 1/2)

  # Year tags
  if (year.tag == TRUE) {
    if (treat == 2) {
      rect(posxl1[1] - 0.134, ymax - 0.00675, xmax / 2, ymax + 0.00675, col = 'white')
      rect(xmax / 2, ymax - 0.00675, posxl2[4] + 0.135, ymax + 0.00675, col = 'white')
    } else {
      rect(posxl1[1] - 0.134, ymax - 0.0132, xmax / 2, ymax + 0.013, col = 'white')
      rect(xmax / 2, ymax - 0.0132, posxl2[4] + 0.135, ymax + 0.013, col = 'white')      
    }
    text(cutx1[2], ymax, expression(bold('2010')))
    text(cutx2[2], ymax, expression(bold('2019')))
  }

  # Fill in
  for (mth in mths) {
    qui <- which(mths == mth)
    xv1 <- c(posp11[qui], posp12[qui], posp13[qui])
    xv2 <- c(posp21[qui], posp22[qui], posp23[qui])

    # Plot parameters
    c1s <- paste('c1.', mth, treat, sep = '')
    c2s <- paste('c2.', mth, treat, sep = '')
    ahs <- paste('ah.', mth, treat, sep = '')
    pds <- paste('pd.', mth, treat, sep = '')
    nrs <- paste('nr.', mth, sep = '')
    if (mth == 'ols') { nrs <- 'p' }

    if (mth %in% c('ols', 'dle')) {
      pchs <- cbind(ifelse(J[, pds] < 0.05, 21, 22),
        ifelse(K[, pds] < 0.05, 21, 22),
        ifelse(L[, pds] < 0.05, 21, 22))
    } else {
      pchs <- cbind(ifelse(J[, pds] > 0.5, 21, 22),
        ifelse(K[, pds] > 0.5, 21, 22),
        ifelse(L[, pds] > 0.5, 21, 22))
    }

    # Point estimates and CIs
    segments(x0 = xv1[1], x1 = xv1[1],
      y0 = J[1, c1s], y1 = J[1, c2s], lwd = 1.5)
    points(xv1[1], J[1, ahs], bg = 'black', pch = pchs[1, 1])
    segments(x0 = xv1[2], x1 = xv1[2],
      y0 = K[1, c1s], y1 = K[1, c2s], col = 'darkgray', lwd = 1.5)
    points(xv1[2], K[1, ahs], bg = 'darkgray', pch = pchs[1, 2])
    segments(x0 = xv1[3], x1 = xv1[3],
      y0 = L[1, c1s], y1 = L[1, c2s], col = 'darkgray', lwd = 1.5)
    points(xv1[3], L[1, ahs], bg = 'darkgray', pch = pchs[1, 3])

    # Point estimates and CIs
    segments(x0 = xv2[1], x1 = xv2[1],
      y0 = J[2, c1s], y1 = J[2, c2s], lwd = 1.5)
    points(xv2[1], J[2, ahs], bg = 'black', pch = pchs[2, 1])
    segments(x0 = xv2[2], x1 = xv2[2],
      y0 = K[2, c1s], y1 = K[2, c2s], col = 'darkgray', lwd = 1.5)
    points(xv2[2], K[2, ahs], bg = 'darkgray', pch = pchs[2, 2])
    segments(x0 = xv2[3], x1 = xv2[3],
      y0 = L[2, c1s], y1 = L[2, c2s], col = 'darkgray', lwd = 1.5)
    points(xv2[3], L[2, ahs], bg = 'darkgray', pch = pchs[2, 3])
  }
}

################################################################################
# Plot for FIGURE 4 (LEFT PANEL): customized boxplot
min.twoway.boxplot <- function(x1, x2, dot = 'median', titol = '', ylab, xlab,
  cex.axis = 1) {
################################################################################
  pm.10 <- ifelse(dot == 'median', median(x1), mean(x1))
  pm.19 <- ifelse(dot == 'median', median(x2), mean(x2))
  tq.10 <- quantile(x1, c(0.25, 0.75))
  tq.19 <- quantile(x2, c(0.25, 0.75))
  ci.10 <- quantile(x1, c(0.05, 0.95))
  ci.19 <- quantile(x2, c(0.05, 0.95))

  ylims = c(0.99, 1.18)
  plot(NA, xlim = c(-0.1, +0.1), ylim = ylims, main = titol, xaxt = 'n',
    xlab = '', ylab = ylab, cex.axis = cex.axis)
  abline(h = 0.95, lty = 3, col = 'gray')
  for (k in seq(1.05, 1.25, 0.05)) { abline(h = k, lty = 3, col = 'gray') }
  abline(h = 1, lty = 2)
  axis(1, c(-0.05, +0.05), xlab, cex.axis = cex.axis)
  segments(x0 = -0.05, x1 = -0.05, y0 = ci.10[1], y1 = ci.10[2], lwd = 3)
  rect(-0.07, tq.10[1], -0.03, tq.10[2], col = 'gray')
  points(-0.05, pm.10, pch = 19, cex = 1.5)
  segments(x0 = +0.05, x1 = +0.05, y0 = ci.19[1], y1 = ci.19[2], lwd = 3)
  rect(+0.03, tq.19[1], +0.07, tq.19[2], col = 'gray')
  points(+0.05, pm.19, pch = 19, cex = 1.5)
}
# END OF SCRIPT
