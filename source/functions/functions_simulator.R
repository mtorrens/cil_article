################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(FCNDIR, 'functions_simulator.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * sim.teff
# * sim.teff.mt
################################################################################

################################################################################
# Simulator function of treatment effects simulation designs
#   It generates N datasets according to the specified characteristics and
#   runs the different methods, extracting a set of key statistics for every
#   method on each run.
sim.teff <- function(a, bd, by, n, p, S, phid, phiy, N, do.bac = FALSE,
  do.ndl = FALSE, do.bma = FALSE, do.pcr = FALSE, do.ssl = FALSE,
  do.nme = FALSE, do.nm2 = FALSE, do.acm = FALSE, bma.mprior = 'bbin',
  lpen = 'lambda.min', mod1 = 'ginv', th.prior = 'unif', beta.prior = 'nlp',
  rho.min = NULL, rho.max = NULL, eta = 0.05, R = 1e5, max.mod = Inf,
  bac.pkg = 'BACprior', only.bacInf = FALSE, acpme.onlyDEF = TRUE,
  ncores = NA) {
################################################################################
  # Function shortcuts and prior setup
  ms <- mombf::modelSelection
  igp <- mombf::igprior(alpha = 0.01, lambda = 0.01)
  zup <- mombf::zellnerprior(tau = n)
  nlp <- mombf::momprior(tau = 0.348)
  if (bma.mprior == 'bbin') {
    mprior <- mombf::modelbbprior(alpha.p = 1, beta.p = 1)
  } else if (bma.mprior == 'unif') {
    mprior <- mombf::modelunifprior()
  } else if (bma.mprior == 'binom') {
    ph <- (length(which(by != 0)) + as.numeric(a != 0)) / (p + 1)
    mprior <- mombf::modelbinomprior(p = ph)
  }

  # Amount of confounders
  conf.idx <- which(by == 1 & bd == 1)  # Indexes of confounder variables
  inst.idx <- which(by == 0 & bd == 1)  # Indexes of instrument variables
  #nc <- sum((by != 0) + (bd != 0) == 2)
  nc <- length(conf.idx)
  ni <- length(inst.idx)

  # Parallelise
  if (is.na(ncores)) { ncores <- parallel::detectCores() - 1 }
  doMC::registerDoMC(cores = ncores)
  cat('Running on', ncores, 'cores:\n')

  # Simulations
  sim <- foreach(r = 1:N, .errorhandling = 'pass', .combine = 'rbind') %dopar% {
    # Print progress
    cat('\r* Iteration:', r, 'of', N)
    cat(' (', sprintf('%01.1f', round(100 * r / N, 1)), '%)', sep = '')

    # Simulate data
    set.seed(100 * r)
    X <- mvtnorm::rmvnorm(n = n, mean = rep(0, nrow(S)), sigma = S)
    epsd <- rnorm(n, mean = 0, sd = sqrt(phid))
    epsy <- rnorm(n, mean = 0, sd = sqrt(phiy))
    d <- X %*% bd + epsd
    y <- a * d + X %*% by + epsy
    Z <- cbind(d, X)

    ############################################################################
    # Oracle model
    orc.model <- paste('x', which(by != 0), sep = '')
    if (a != 0) { orc.model <- c('d1', orc.model) }

    # Oracle OLS estimation
    orc.lm <- lm(y ~ d + X[, which(by != 0)] - 1)
    sum.lm <- summary(orc.lm)
    ah.orc <- coef(orc.lm)[['d']]
    nr.orc <- length(which(by != 0)) + 1 #as.numeric(a != 0)  # No. regressors
    fp.orc <- 0  # False positives
    fn.orc <- sum(summary(orc.lm)[['coefficients']][, 4] > eta)  # False negs.
    ci.orc <- confint(orc.lm, level = 1 - eta)['d', ]  # Confidence intervals
    il.orc <- diff(range(ci.orc))  # Interval length
    cg.orc <- inside(a, ci.orc)
    pd.orc <- !inside(0, ci.orc)
    pc.orc <- sum(sum.lm[['coefficients']][conf.idx + 1, 4] < eta) / nc
    pi.orc <- 0

    # Oracle estimation including whatever affects d
    inc <- unique(c(which(by != 0), which(bd != 0)))
    tinc <- which(by != 0) + 1; if (a != 0) { tinc <- c(1, tinc) }
    ore.lm <- lm(y ~ d + X[, inc] - 1)
    sum.lm <- summary(ore.lm)
    ah.ore <- coef(ore.lm)[['d']]
    nr.ore <- length(inc) + as.numeric(a != 0)
    fp.ore <- sum(! which(bd != 0) %in% which(by != 0))  # FPs
    fn.ore <- sum(summary(ore.lm)[['coefficients']][tinc, 4] > eta)  # FNs
    ci.ore <- confint(ore.lm, level = 1 - eta)['d', ]  # Confidence intervals
    il.ore <- diff(range(ci.ore))
    cg.ore <- inside(a, ci.ore)
    pd.ore <- !inside(0, ci.ore)
    pc.ore <- sum(sum.lm[['coefficients']][conf.idx + 1, 4] < eta) / nc
    pi.ore <- sum(sum.lm[['coefficients']][inst.idx + 1, 4] < eta) / ni

    ############################################################################
    # Lasso-based methods
    if (do.ndl == TRUE) {
      # Naive Lasso
      nl.fit <- cv.glmnet(Z, y, alpha = 1, intercept = FALSE, standardize = FALSE)
      ah.nle <- coef(nl.fit, s = nl.fit[[lpen]])['V1', ]
      aux <- which(coef(nl.fit, s = nl.fit[[lpen]]) != 0) - 2
      nle.model <- paste('x', aux[aux >= 0], sep = '')
      nle.model[which(nle.model == 'x0')] <- 'd1'
      nr.nle <- sum(coef(nl.fit, s = nl.fit[[lpen]])[-1] != 0)
      fp.nle <- sum(! nle.model %in% orc.model)
      fn.nle <- sum(! orc.model %in% nle.model)      
      il.nle <- NA
      cg.nle <- NA
      pd.nle <- length(which(nle.model == 'd1'))
      pc.nle <- sum(conf.idx %in% aux) / nc
      pi.nle <- sum(inst.idx %in% aux) / ni

      # Attempt inference on treatment using method in Lee et al. (2016)
      nl.beta <- as.vector(coef(nl.fit, s = nl.fit[[lpen]]))[-1]
      nl.sel <- nl.beta != 0
      nsigma <- sd((y - as.matrix(Z[, nl.sel]) %*% nl.beta[nl.sel != 0]))
      nsigma <- nsigma * sqrt(n / (n - sum(nl.sel) - 1))
      nl.inf <- selectiveInference::fixedLassoInf(x = Z, y = y,
        beta = coef(nl.fit, x = Z, y = y, s = lpen)[-1],
        lambda = n * nl.fit[[lpen]], sigma = nsigma, alpha = 0.05)
      tvar <- which(nl.inf[['vars']] == 1)
      if (length(tvar) > 0) { 
        ci.nle <- nl.inf[['ci']][tvar, ]
      } else {
        ci.nle <- c(0, 0)
      }
      il.nle <- diff(range(ci.nle))
      cg.nle <- inside(a, ci.nle)

      # Double Lasso
      dl.fit <- rlassoEffect(X, y, d, method = 'double selection',
        post = TRUE, intercept = FALSE)
      ah.dle <- unname(dl.fit[['alpha']])
      nr.dle <- length(dl.fit[['coefficients.reg']]) - 1
      aux <- names(dl.fit[['coefficients.reg']][2:(nr.dle + 1)])
      eux <- substr(aux[grep('^xx', aux)], 3, nchar(aux[grep('^xx', aux)]))
      dle.model <- substr(aux, 2, nchar(aux))
      fp.dle <- sum(! dle.model %in% orc.model)
      fn.dle <- sum(! orc.model %in% dle.model)
      ci.dle <- confint(dl.fit, level = 1 - eta)[1, ]
      il.dle <- diff(range(ci.dle))
      cg.dle <- inside(a, ci.dle)
      pd.dle <- !inside(0, ci.dle)
      pc.dle <- sum(conf.idx %in% as.numeric(eux)) / nc
      pi.dle <- sum(inst.idx %in% as.numeric(eux)) / ni
    } else {
      ah.nle <- nr.nle <- fp.nle <- fn.nle <- il.nle <- cg.nle <- pd.nle <- pc.nle <- pi.nle <- NA
      ah.dle <- nr.dle <- fp.dle <- fn.dle <- il.dle <- cg.dle <- pd.dle <- pc.dle <- pi.dle <- NA
    }

    ############################################################################
    # BAC
    if (do.bac == TRUE & bac.pkg == 'bacr') {
      # Reshape design (custom)
      W <- cbind.data.frame(y, d, X)
      colnames(W) <- c('y', 'd', paste('V', 1:ncol(X), sep = ''))

      # Problems with function when p + 1 >= n
      if (ncol(W) - 1 >= nrow(W)) {
        bacF <- bac2
      } else {
        bacF <- bacr::bac
      }

      # Run algorithm
      sink(file = paste(TMPDIR, 'bac_esborrar', r, '.txt', sep = ''))  # There is no silent option :_)
      bac.fitI <- try(bacF(data = W, exposure = 'd', outcome = 'y',
        confounders = paste('V', 1:ncol(X), sep = ''), interactors = NULL,
        familyX = 'gaussian', familyY = 'gaussian', omega = Inf,
        num_its = 5e3, burnM = 5e2, burnB = 5e2, thin = 1), silent = TRUE)
      if (class(bac.fitI) == 'try-error') {
        ah.bc1 <- nr.bc1 <- fp.bc1 <- fn.bc1 <- il.bc1 <- cg.bc1 <- pd.bc1 <- pc.bc1 <- pi.bc1 <- NA
        ah.bcT <- nr.bcT <- fp.bcT <- fn.bcT <- il.bcT <- cg.bcT <- pd.bcT <- pc.bcT <- pi.bcT <- NA
        ah.bcI <- nr.bcI <- fp.bcI <- fn.bcI <- il.bcI <- cg.bcI <- pd.bcI <- pc.bcI <- pi.bcI <- NA
        sink(); unlink(paste(TMPDIR, 'bac_esborrar', r, '.txt', sep = ''))
      } else {
        if (only.bacInf == TRUE) {
          bac.fit1 <- bac.fitT <- bac.fitI
        } else {
          bac.fit1 <- try(bacF(data = W, exposure = 'd', outcome = 'y',
            confounders = paste('V', 1:ncol(X), sep = ''), interactors = NULL,
            familyX = 'gaussian', familyY = 'gaussian', omega = 1, 
            num_its = 5e3, burnM = 5e2, burnB = 5e2, thin = 1), silent = TRUE)
          bac.fitT <- try(bacF(data = W, exposure = 'd', outcome = 'y',
            confounders = paste('V', 1:ncol(X), sep = ''), interactors = NULL,
            familyX = 'gaussian', familyY = 'gaussian', omega = 10, 
            num_its = 5e3, burnM = 5e2, burnB = 5e2, thin = 1), silent = TRUE)
        }

        # Point estimates for each omega
        bac1 <- summary(bac.fit1)
        bacI <- summary(bac.fitI)
        bacT <- summary(bac.fitT)
        sink(); unlink(paste(TMPDIR, 'bac_esborrar', r, '.txt', sep = ''))
        ci.bc1 <- bac1[['CI']]
        ci.bcI <- bacI[['CI']]
        ci.bcT <- bacT[['CI']]

        # Posterior average model size
        aux1a <- bac.fit1[['models']][[1]]
        aux2a <- bac.fit1[['models']][[2]] / sum(bac.fit1[['models']][[2]])
        aux3a <- bac.fit1[['models']][[1]][, conf.idx]
        aux4a <- bac.fit1[['models']][[1]][, inst.idx]
        nr.bc1 <- sum(rowSums(aux1a) * aux2a)
        pc.bc1 <- sum(rowSums(aux3a) * aux2a) / nc
        pi.bc1 <- sum(rowSums(aux4a) * aux2a) / ni
        aux1b <- bac.fitI[['models']][[1]]
        aux2b <- bac.fitI[['models']][[2]] / sum(bac.fitI[['models']][[2]])
        aux3b <- bac.fitI[['models']][[1]][, conf.idx]
        aux4b <- bac.fitI[['models']][[1]][, inst.idx]
        nr.bcI <- sum(rowSums(aux1b) * aux2b)
        pc.bcI <- sum(rowSums(aux3b) * aux2b) / nc
        pi.bcI <- sum(rowSums(aux4b) * aux2b) / ni
        aux1c <- bac.fitT[['models']][[1]]
        aux2c <- bac.fitT[['models']][[2]] / sum(bac.fitT[['models']][[2]])
        aux3c <- bac.fitT[['models']][[1]][, conf.idx]
        aux4c <- bac.fitT[['models']][[1]][, inst.idx]
        nr.bcT <- sum(rowSums(aux1c) * aux2c)
        pc.bcT <- sum(rowSums(aux3c) * aux2c) / nc
        pi.bcT <- sum(rowSums(aux4c) * aux2c) / ni

        # False positives and negatives
        ids1 <- apply(aux1a, 1, function(x) {
          mid <- paste('x', which(x[-length(x)] != 0), sep = '')
          if (x[length(x)] != 0) { mid <- c('d1', mid) }
          return(mid)
        })
        ids2 <- apply(aux1b, 1, function(x) {
          mid <- paste('x', which(x[-length(x)] != 0), sep = '')
          if (x[length(x)] != 0) { mid <- c('d1', mid) }
          return(mid)
        })
        ids3 <- apply(aux1c, 1, function(x) {
          mid <- paste('x', which(x[-length(x)] != 0), sep = '')
          if (x[length(x)] != 0) { mid <- c('d1', mid) }
          return(mid)
        })
        fps1 <- unname(unlist(lapply(ids1, function(x) {
          sum(! x %in% orc.model)
        } )))
        fps2 <- unname(unlist(lapply(ids2, function(x) {
          sum(! x %in% orc.model)
        } )))
        fps3 <- unname(unlist(lapply(ids3, function(x) {
          sum(! x %in% orc.model)
        } )))
        if (nrow(aux1a) == 1) {
          fn.bc1 <- sum(! orc.model %in% ids1)
          fp.bc1 <- sum(fps1)
        } else {
          fns1 <- unname(unlist(lapply(ids1, function(x) {
            sum(! orc.model %in% x)
          } )))
          fn.bc1 <- sum(fns1 * aux2a)
          fp.bc1 <- sum(fps1 * aux2a)
        }
        if (nrow(aux1b) == 1) {
          fn.bcI <- sum(! orc.model %in% ids2)
          fp.bcI <- sum(fps2)
        } else {
          fns2 <- unname(unlist(lapply(ids2, function(x) {
            sum(! orc.model %in% x)
          } )))
          fn.bcI <- sum(fns2 * aux2b)
          fp.bcI <- sum(fps2 * aux2b)
        }
        if (nrow(aux1c) == 1) {
          fn.bcT <- sum(! orc.model %in% ids3)
          fp.bcT <- sum(fps3)
        } else {
          fns3 <- unname(unlist(lapply(ids3, function(x) {
            sum(! orc.model %in% x)
          } )))
          fn.bcT <- sum(fns3 * aux2c)
          fp.bcT <- sum(fps3 * aux2c)
        }

        # Summary
        ah.bc1 <- bac1[['posterior.mean']]
        ah.bcI <- bacI[['posterior.mean']]
        ah.bcT <- bacT[['posterior.mean']]
        cg.bc1 <- inside(a, ci.bc1)
        cg.bcI <- inside(a, ci.bcI)
        cg.bcT <- inside(a, ci.bcT)
        il.bc1 <- diff(range(ci.bc1))
        il.bcI <- diff(range(ci.bcI))
        il.bcT <- diff(range(ci.bcT))
        pd.bc1 <- !inside(0, ci.bc1)
        pd.bcI <- !inside(0, ci.bcI)
        pd.bcT <- !inside(0, ci.bcT)

        # If running two BACs is too much
        if (only.bacInf == TRUE) {
          ah.bc1 <- nr.bc1 <- fp.bc1 <- fn.bc1 <- il.bc1 <- cg.bc1 <- pd.bc1 <- pc.bc1 <- pi.bc1 <- NA        
          ah.bcT <- nr.bcT <- fp.bcT <- fn.bcT <- il.bcT <- cg.bcT <- pd.bcT <- pc.bcT <- pi.bcT <- NA        
        }
      }  
    } else if (do.bac == TRUE & bac.pkg == 'BACprior') {
      os <- c(1, 10, Inf)
      bac.fit <- BACprior.lm(Y = as.numeric(y), X = as.numeric(d), U = X,
        omega = os, return.best = TRUE)
      bac.coef <- bac.fit[['results']]

      # Point estimates for each omega
      bac1 <- unname(bac.coef[which(bac.coef[, 1] == 1), ])
      bacT <- unname(bac.coef[which(bac.coef[, 1] == 10), ])
      bacI <- unname(bac.coef[which(bac.coef[, 1] == Inf), ])
      ci.bc1 <- bac1[2] + qnorm(1 - eta / 2) * c(-bac1[3], bac1[3])
      ci.bcT <- bacT[2] + qnorm(1 - eta / 2) * c(-bacT[3], bacT[3])
      ci.bcI <- bacI[2] + qnorm(1 - eta / 2) * c(-bacI[3], bacI[3])

      # Posterior average model size
      aux1 <- bac.fit[['best.models']]
      aux2 <- bac.fit[['posterior.prob']]
      nr.bc1 <- sum(unname(rowSums(aux1)) * aux2[, 1]) + 1
      nr.bcT <- sum(unname(rowSums(aux1)) * aux2[, 2]) + 1
      nr.bcI <- sum(unname(rowSums(aux1)) * aux2[, 3]) + 1

      # False positives and negatives
      ids <- apply(aux1, 1, function(x) {
        c('d1', paste('x', which(x != 0), sep = ''))
      })
      fps <- unname(unlist(lapply(ids, function(x) {
        sum(! x %in% orc.model)
      } )))
      fns <- unname(unlist(lapply(ids, function(x) {
        sum(! orc.model %in% x)
      } )))
      fp.bc1 <- sum(fps * aux2[, 1])
      fp.bcT <- sum(fps * aux2[, 2])
      fp.bcI <- sum(fps * aux2[, 3])
      fn.bc1 <- sum(fns * aux2[, 1])
      fn.bcT <- sum(fns * aux2[, 2])
      fn.bcI <- sum(fns * aux2[, 3])

      # Summary
      ah.bc1 <- bac1[2]
      ah.bcT <- bacT[2]
      ah.bcI <- bacI[2]
      cg.bc1 <- inside(a, ci.bc1)
      cg.bcT <- inside(a, ci.bcT)
      cg.bcI <- inside(a, ci.bcI)
      il.bc1 <- diff(range(ci.bc1))
      il.bcT <- diff(range(ci.bcT))
      il.bcI <- diff(range(ci.bcI))
      pd.bc1 <- !inside(0, ci.bc1)
      pd.bcT <- !inside(0, ci.bcT)
      pd.bcI <- !inside(0, ci.bcI)

      # LEFT TO DO IF NECESSARY
      pc.bc1 <- pc.bcT <- pc.bcI <- NA
      pi.bc1 <- pi.bcT <- pi.bcI <- NA     
    } else {
      ah.bc1 <- nr.bc1 <- fp.bc1 <- fn.bc1 <- il.bc1 <- cg.bc1 <- pd.bc1 <- pc.bc1 <- pi.bc1 <- NA
      ah.bcT <- nr.bcT <- fp.bcT <- fn.bcT <- il.bcT <- cg.bcT <- pd.bcT <- pc.bcT <- pi.bcT <- NA
      ah.bcI <- nr.bcI <- fp.bcI <- fn.bcI <- il.bcI <- cg.bcI <- pd.bcI <- pc.bcI <- pi.bcI <- NA
    }

    ############################################################################
    # Penalised Credible Regions (Wilson and Reich, 2014)
    if (do.pcr == TRUE) {
    #if (do.pcr == TRUE & nrow(Z) > ncol(Z)) {
      # Choose based on BIC
      sink(file = paste(TMPDIR, 'tmp_esborrar', r, '.txt', sep = ''))  # There is no silent option :'-(
      pcr.fit <- try(BayesPen::BayesPen.lm.confounders(y, d, X), silent = TRUE)
      if (class(pcr.fit) == 'try-error') {
        pcr.fit <- BayesPen.lm.confounders2(y, d, X)
      }
      sink(); unlink(paste(TMPDIR, 'tmp_esborrar', r, '.txt', sep = ''))

      # Choose min BIC
      prpcr <- mean(y) + cbind(d, X) %*% t(pcr.fit[['coefs']])
      nC <- apply(pcr.fit[['coefs']], 1, function(x) sum(x != 0))
      nr <- nrow(X)
      bic <- nr * log(colSums(apply(prpcr, 2, function(x) { x - y })^2) / length(y)) +
             nr * (log(2 * pi) + 1) + log(nr) * nC
      wr14.refit.coef <- pcr.fit[['coefs']][which.min(bic), ]
      names(wr14.refit.coef) <- c('d1', paste('x', 1:p, sep = ''))
      wr14.model <- names(wr14.refit.coef)[which(wr14.refit.coef != 0)]

      # Relevant values
      ah.pcr <- unname(wr14.refit.coef[1])
      nr.pcr <- sum(wr14.refit.coef != 0)
      fp.pcr <- sum(! wr14.model %in% orc.model)
      fn.pcr <- sum(! orc.model %in% wr14.model)
      cg.pcr <- NA
      il.pcr <- NA
      pd.pcr <- ifelse(ah.pcr == 0, 0, 1)
    } else {
      ah.pcr <- nr.pcr <- fp.pcr <- fn.pcr <- il.pcr <- cg.pcr <- pd.pcr <- NA
    }
    pc.pcr <- pi.pcr <- NA  # TO-DO if ever necessary

    ############################################################################
    # BMA on y
    if (do.bma == TRUE) {
      # BMS
      bms1 <- ms(y, Z, priorCoef = zup, priorDelta = mprior, priorVar = igp,
                 niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)
      bms2 <- ms(y, Z, priorCoef = nlp, priorDelta = mprior, priorVar = igp,
                 niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)

      # BMA
      gibbs1 <- mombf::rnlp(msfit = bms1, priorCoef = zup, niter = 1e4)
      gibbs2 <- mombf::rnlp(msfit = bms2, priorCoef = nlp, niter = 1e4)
      ah.bm1 <- mean(gibbs1[, 'beta1'])
      ah.bm2 <- mean(gibbs2[, 'beta1'])

      # Model size
      aux1 <- mombf::postProb(bms1)
      aux2 <- mombf::postProb(bms2)
      naux1 <- lapply(strsplit(as.character(aux1[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - 1, sep = '')
        if (ans[1] == 'x0') { ans[1] <- 'd1' }
        return(ans)
      })
      naux2 <- lapply(strsplit(as.character(aux2[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - 1, sep = '')
        if (ans[1] == 'x0') { ans[1] <- 'd1' }
        return(ans)
      })
      fps1 <- unlist(lapply(naux1, function(x) { sum(! x %in% orc.model) }))
      fps2 <- unlist(lapply(naux2, function(x) { sum(! x %in% orc.model) }))
      fns1 <- unlist(lapply(naux1, function(x) { sum(! orc.model %in% x) }))
      fns2 <- unlist(lapply(naux2, function(x) { sum(! orc.model %in% x) }))

      nr.bm1 <- sum(unlist(lapply(
        strsplit(as.character(aux1[, 1]), ','), length)) * aux1[, 'pp'])
      nr.bm2 <- sum(unlist(lapply(
        strsplit(as.character(aux2[, 1]), ','), length)) * aux2[, 'pp'])
      fp.bm1 <- sum(fps1 * aux1[, 'pp'])
      fp.bm2 <- sum(fps2 * aux2[, 'pp'])
      fn.bm1 <- sum(fns1 * aux1[, 'pp'])
      fn.bm2 <- sum(fns2 * aux2[, 'pp'])

      # Confounder and instrument detection
      naux1b <- lapply(naux1, function(x) {
        covs <- substr(x[grep('^x', x)], 2, nchar(x[grep('^x', x)]))
        ansc <- sum(conf.idx %in% as.numeric(covs)) / nc
        ansi <- sum(inst.idx %in% as.numeric(covs)) / ni
        return(c(ansc, ansi))
      })
      naux2b <- lapply(naux2, function(x) {
        covs <- substr(x[grep('^x', x)], 2, nchar(x[grep('^x', x)]))
        ansc <- sum(conf.idx %in% as.numeric(covs)) / nc
        ansi <- sum(inst.idx %in% as.numeric(covs)) / ni
        return(c(ansc, ansi))
      })
      pc.bm1 <- sum(unlist(lapply(naux1b, '[[', 1)) * aux1[, 'pp'])
      pc.bm2 <- sum(unlist(lapply(naux2b, '[[', 1)) * aux2[, 'pp'])
      pi.bm1 <- sum(unlist(lapply(naux1b, '[[', 2)) * aux1[, 'pp'])
      pi.bm2 <- sum(unlist(lapply(naux2b, '[[', 2)) * aux2[, 'pp'])

      # Inclusion indicators
      pr.eta <- c(eta / 2, 1 - eta / 2)
      ci.bm1 <- unname(quantile(gibbs1[, 'beta1'], probs = pr.eta))
      ci.bm2 <- unname(quantile(gibbs2[, 'beta1'], probs = pr.eta))
      il.bm1 <- diff(range(ci.bm1))
      il.bm2 <- diff(range(ci.bm2))
      cg.bm1 <- inside(a, ci.bm1)
      cg.bm2 <- inside(a, ci.bm2)
      pd.bm1 <- unname(bms1[['margpp']]['x1'])
      pd.bm2 <- unname(bms2[['margpp']]['x1'])
      #pd.bm1 <- !inside(0, ci.bm1)

      # MAP models
      map1 <- as.character(aux1[1, 1])
      map2 <- as.character(aux2[1, 1])
      t1.bm1 <- as.numeric(binary.id(map1, p = ncol(X))[1] == 1)
      t1.bm2 <- as.numeric(binary.id(map2, p = ncol(X))[1] == 1)
      t2.bm1 <- as.numeric(unname(bms1[['margpp']][1] >= 0.5))
      t2.bm2 <- as.numeric(unname(bms2[['margpp']][1] >= 0.5))
      t3.bm1 <- as.numeric(unname(bms1[['margpp']][1] >= 0.95))
      t3.bm2 <- as.numeric(unname(bms2[['margpp']][1] >= 0.95))
    } else {
      ah.bm1 <- nr.bm1 <- fp.bm1 <- fn.bm1 <- il.bm1 <- cg.bm1 <- pd.bm1 <- pc.bm1 <- pi.bm1 <- NA
      ah.bm2 <- nr.bm2 <- fp.bm2 <- fn.bm2 <- il.bm2 <- cg.bm2 <- pd.bm2 <- pc.bm2 <- pi.bm2 <- NA
      t1.bm1 <- t1.bm2 <- t2.bm1 <- t2.bm2 <- t3.bm1 <- t3.bm2 <- NA
    }

    ############################################################################
    # New SSL (2019)
    if (do.ssl == TRUE) {
      erase.file <- paste(TMPDIR, 'tmp_esborrar', r, '.txt', sep = '')
      sink(file = erase.file)  # There is no silent option :'-(
      ssl.fit <- HDconfounding::SSL(y = y, z = d, x = X, z_type = 'continuous',
        nScans = 8e3, burn = 1e3)
      sink(); unlink(erase.file)

      # Numerical summaries
      ah.ssl <- ssl.fit[['TreatEffect']]
      ci.ssl <- unname(ssl.fit[['TreatEffectCI']])
      il.ssl <- diff(range(ci.ssl))
      cg.ssl <- inside(a, ci.ssl)
      pd.ssl <- !inside(0, ci.ssl)

      # Model size
      aux <- ssl.fit[['betaPostCI']]
      nr.ssl <- sum(sign(aux[1, ]) == sign(aux[2, ])) + pd.ssl
      aux <- paste('x', which(sign(aux[1, ]) == sign(aux[2, ])), sep = '')
      if (pd.ssl == TRUE) { aux <- c('d1', aux) }
      fp.ssl <- sum(! aux %in% orc.model)
      fn.ssl <- sum(! orc.model %in% aux)
    } else {
      ah.ssl <- nr.ssl <- fp.ssl <- fn.ssl <- il.ssl <- cg.ssl <- pd.ssl <- NA
    }
    pc.ssl <- pi.ssl <- NA  # TO-DO if ever necessary

    ############################################################################
    # Our new method (using Expectation Propagation)
    if (do.nme == TRUE) {
      # Expectation Propagation
      nme.fit <- new.method(y, d, X, th.search = 'EP', mod1 = mod1,
        th.prior = th.prior, beta.prior = beta.prior, rho.min = rho.min,
        rho.max = rho.max, R = R, max.mod = max.mod)

      # Record estimated parameters
      th0.ep <- nme.fit[['th0.hat']]
      th1.ep <- nme.fit[['th1.hat']]

      # Model size
      aux <- nme.fit[['pprob']]
      naux <- lapply(strsplit(as.character(aux[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - 1, sep = '')
        if (ans[1] == 'x0') { ans[1] <- 'd1' }
        return(ans)
      })
      fps <- unlist(lapply(naux, function(x) { sum(! x %in% orc.model) }))
      fns <- unlist(lapply(naux, function(x) { sum(! orc.model %in% x) }))
      nr.nmA <- sum(unlist(lapply(
        strsplit(as.character(aux[, 1]), ','), length)) * aux[, 2])
      fp.nmA <- sum(fps * aux[, 2])
      fn.nmA <- sum(fns * aux[, 2])

      # Confounder and instrument detection
      naux1 <- lapply(naux, function(x) {
        covs <- substr(x[grep('^x', x)], 2, nchar(x[grep('^x', x)]))
        ansc <- sum(conf.idx %in% as.numeric(covs)) / nc
        ansi <- sum(inst.idx %in% as.numeric(covs)) / ni
        return(c(ansc, ansi))
      })
      pc.nmA <- sum(unlist(lapply(naux1, '[[', 1)) * aux[, 'pp'])
      pi.nmA <- sum(unlist(lapply(naux1, '[[', 2)) * aux[, 'pp'])

      # Numerical summaries (NLP)
      pr.eta <- c(eta / 2, 1 - eta / 2)
      noms <- paste(100 * pr.eta, '%', sep = '')
      ah.nmA <- nme.fit[['teff.nlp']]
      ci.nmA <- unname(nme.fit[['teff.d.nlp']][noms])
      il.nmA <- diff(range(ci.nmA))
      cg.nmA <- inside(a, ci.nmA)
      pd.nmA <- unname(nme.fit[['mpprob']][1])

      # Plan B (Gaussian prior)
      ah.nmB <- nme.fit[['teff.zup']]
      ci.nmB <- unname(nme.fit[['teff.d.zup']][noms])
      il.nmB <- diff(range(ci.nmB))
      cg.nmB <- inside(a, ci.nmB)
      nr.nmB <- nr.nmA
      fp.nmB <- fp.nmA
      fn.nmB <- fn.nmA
      pd.nmB <- pd.nmA
      pc.nmB <- pc.nmA
      pi.nmB <- pi.nmA

      # MAP model
      mapA <- as.character(nme.fit[['pprob']][1, 1])
      t1.nmA <- as.numeric(binary.id(mapA, p = ncol(X))[1] == 1)
      t2.nmA <- as.numeric(unname(nme.fit[['mpprob']][1] >= 0.5))
      t3.nmA <- as.numeric(unname(nme.fit[['mpprob']][1] >= 0.95))      
    } else {
      ah.nmA <- nr.nmA <- fp.nmA <- fn.nmA <- il.nmA <- cg.nmA <- pd.nmA <- pc.nmA <- pi.nmA <- NA
      ah.nmB <- nr.nmB <- fp.nmB <- fn.nmB <- il.nmB <- cg.nmB <- pd.nmB <- pc.nmB <- pi.nmB <- NA
      th0.ep <- th1.ep <- NA
      t1.nmA <- t2.nmA <- t3.nmA <- NA
    }

    ############################################################################
    # New method using Empirical Bayes
    if (do.nm2 == TRUE) {
      # Recycle initial MCMC if we've done it before
      if (do.nme == TRUE) {
        bvs.fit0 <- nme.fit[['init.msfit']]
      } else {
        bvs.fit0 <- NULL
      }

      # Empirical Bayes
      nme.fit2 <- new.method(y, d, X, th.search = 'EB', mod1 = mod1,
        th.prior = th.prior, beta.prior = beta.prior, rho.min = rho.min,
        rho.max = rho.max, R = R, max.mod = max.mod, bvs.fit0 = bvs.fit0)

      # Record estimated parameters
      th0.eb <- nme.fit2[['th0.hat']]
      th1.eb <- nme.fit2[['th1.hat']]

      # Model size
      aux <- nme.fit2[['pprob']]
      naux <- lapply(strsplit(as.character(aux[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - 1, sep = '')
        if (ans[1] == 'x0') { ans[1] <- 'd1' }
        return(ans)
      })
      fps <- unlist(lapply(naux, function(x) { sum(! x %in% orc.model) }))
      fns <- unlist(lapply(naux, function(x) { sum(! orc.model %in% x) }))
      nr.n2A <- sum(unlist(lapply(
        strsplit(as.character(aux[, 1]), ','), length)) * aux[, 2])
      fp.n2A <- sum(fps * aux[, 2])
      fn.n2A <- sum(fns * aux[, 2])

      # Confounder and instrument detection
      naux1 <- lapply(naux, function(x) {
        covs <- substr(x[grep('^x', x)], 2, nchar(x[grep('^x', x)]))
        ansc <- sum(conf.idx %in% as.numeric(covs)) / nc
        ansi <- sum(inst.idx %in% as.numeric(covs)) / ni
        return(c(ansc, ansi))
      })
      pc.n2A <- sum(unlist(lapply(naux1, '[[', 1)) * aux[, 'pp'])
      pi.n2A <- sum(unlist(lapply(naux1, '[[', 2)) * aux[, 'pp'])

      # Numerical summaries (NLP)
      noms <- paste(100 * c(eta / 2, 1 - eta / 2), '%', sep = '')
      ah.n2A <- nme.fit2[['teff.nlp']]
      ci.n2A <- unname(nme.fit2[['teff.d.nlp']][noms])
      il.n2A <- diff(range(ci.n2A))
      cg.n2A <- inside(a, ci.n2A)
      pd.n2A <- unname(nme.fit2[['mpprob']][1])

      # Plan B (Gaussian prior)
      ah.n2B <- nme.fit2[['teff.zup']]
      ci.n2B <- unname(nme.fit2[['teff.d.zup']][noms])
      il.n2B <- diff(range(ci.n2B))
      cg.n2B <- inside(a, ci.n2B)
      nr.n2B <- nr.n2A
      fp.n2B <- fp.n2A
      fn.n2B <- fn.n2A
      pd.n2B <- pd.n2A
      pc.n2B <- pc.n2A
      pi.n2B <- pi.n2A

      # MAP model
      mapB <- as.character(nme.fit2[['pprob']][1, 1])
      t1.nmB <- as.numeric(binary.id(mapB, p = ncol(X))[1] == 1)
      t2.nmB <- as.numeric(unname(nme.fit2[['mpprob']][1] >= 0.5))
      t3.nmB <- as.numeric(unname(nme.fit2[['mpprob']][1] >= 0.95))      
    } else {
      ah.n2A <- nr.n2A <- fp.n2A <- fn.n2A <- il.n2A <- cg.n2A <- pd.n2A <- pc.n2A <- pi.n2A <- NA
      ah.n2B <- nr.n2B <- fp.n2B <- fn.n2B <- il.n2B <- cg.n2B <- pd.n2B <- pc.n2A <- pi.n2A <- NA
      th0.eb <- th1.eb <- NA
      t1.nmB <- t2.nmB <- t3.nmB <- NA
    }

    ############################################################################
    # ACPME
    ah.acm <- nr.acm <- fp.acm <- fn.acm <- il.acm <- cg.acm <- pd.acm <- pc.acm <- pi.acm <- NA
    ah.ac1 <- nr.ac1 <- fp.ac1 <- fn.ac1 <- il.ac1 <- cg.ac1 <- pd.ac1 <- pc.ac1 <- pi.ac1 <- NA
    ah.ac2 <- nr.ac2 <- fp.ac2 <- fn.ac2 <- il.ac2 <- cg.ac2 <- pd.ac2 <- pc.ac2 <- pi.ac2 <- NA
    if (do.acm == TRUE) {
      # Run method
      sink(file = paste(TMPDIR, 'tmp_acpme.txt', sep = ''))  # No silent option
      acm.fit <- try(
        regimes::acpme(y = y, Z = d, C = X, niter = 1e4, pen.type = 'eigen'),
        silent = TRUE)
      #acm.fit <- regimes::acpme(y = y, Z = d, C = X, niter = 1e4)
      if (acpme.onlyDEF == FALSE) {
        ac1.fit <- try(
          regimes::acpme(y = y, Z = d, C = X, niter = 1e4, pen.type = 'correlation'),
          silent = TRUE)
        ac2.fit <- try(
          regimes::acpme(y = y, Z = d, C = X, niter = 1e4, pen.type = 'projection'),
          silent = TRUE)
      }
      sink()

      if (class(acm.fit) != 'try-error') {
        # BMA point estimate
        ah.acm <- colMeans(acm.fit[['beta']])

        # Model size
        orc.idx <- as.numeric(substr( orc.model[grep('x', orc.model)], 2,
          nchar(orc.model[grep('x', orc.model)])))
        nr.acm <- mean(rowSums(acm.fit[['alpha']]))
        fp.acm <- mean(apply(acm.fit[['alpha']], 1, function(x)
          sum(! unname(which(x == 1)) %in% orc.idx)
        ))
        fn.acm <- mean(apply(acm.fit[['alpha']], 1, function(x)
          sum(! orc.idx %in% unname(which(x == 1)))
        ))

        # Confounder and instrument detection
        pc.acm <- mean(rowSums(acm.fit[['alpha']][, conf.idx]) / nc)
        pi.acm <- mean(rowSums(acm.fit[['alpha']][, inst.idx]) / ni)

        # Inclusion indicators
        pr.eta <- c(eta / 2, 1 - eta / 2)
        ci.acm <- unname(apply(acm.fit[['beta']], 2, quantile, probs = pr.eta))
        il.acm <- apply(ci.acm, 2, function(x) { diff(range(x)) })
        cg.acm <- inside(a, ci.acm[, 1])
        pd.acm <- !inside(0, ci.acm[, 1])
      }
      if (acpme.onlyDEF == FALSE) {
        if (class(ac1.fit) != 'try-error') {
          # BMA point estimate
          ah.ac1 <- colMeans(ac1.fit[['beta']])

          # Model size
          orc.idx <- as.numeric(substr( orc.model[grep('x', orc.model)], 2,
            nchar(orc.model[grep('x', orc.model)])))
          nr.ac1 <- mean(rowSums(ac1.fit[['alpha']]))
          fp.ac1 <- mean(apply(ac1.fit[['alpha']], 1, function(x)
            sum(! unname(which(x == 1)) %in% orc.idx)
          ))
          fn.ac1 <- mean(apply(ac1.fit[['alpha']], 1, function(x)
            sum(! orc.idx %in% unname(which(x == 1)))
          ))

          # Confounder and instrument detection
          pc.ac1 <- mean(rowSums(ac1.fit[['alpha']][, conf.idx]) / nc)
          pi.ac1 <- mean(rowSums(ac1.fit[['alpha']][, inst.idx]) / ni)

          # Inclusion indicators
          pr.eta <- c(eta / 2, 1 - eta / 2)
          ci.ac1 <- unname(apply(ac1.fit[['beta']], 2, quantile, probs = pr.eta))
          il.ac1 <- apply(ci.ac1, 2, function(x) { diff(range(x)) })
          cg.ac1 <- inside(a, ci.ac1[, 1])
          pd.ac1 <- !inside(0, ci.ac1[, 1])
        }
        if (class(ac2.fit) != 'try-error') {
          # BMA point estimate
          ah.ac2 <- colMeans(ac2.fit[['beta']])

          # Model size
          orc.idx <- as.numeric(substr( orc.model[grep('x', orc.model)], 2,
            nchar(orc.model[grep('x', orc.model)])))
          nr.ac2 <- mean(rowSums(ac2.fit[['alpha']]))
          fp.ac2 <- mean(apply(ac2.fit[['alpha']], 1, function(x)
            sum(! unname(which(x == 1)) %in% orc.idx)
          ))
          fn.ac2 <- mean(apply(ac2.fit[['alpha']], 1, function(x)
            sum(! orc.idx %in% unname(which(x == 1)))
          ))

          # Confounder and instrument detection
          pc.ac2 <- mean(rowSums(ac2.fit[['alpha']][, conf.idx]) / nc)
          pi.ac2 <- mean(rowSums(ac2.fit[['alpha']][, inst.idx]) / ni)

          # Inclusion indicators
          pr.eta <- c(eta / 2, 1 - eta / 2)
          ci.ac2 <- unname(apply(ac2.fit[['beta']], 2, quantile, probs = pr.eta))
          il.ac2 <- apply(ci.ac2, 2, function(x) { diff(range(x)) })
          cg.ac2 <- inside(a, ci.ac2[, 1])
          pd.ac2 <- !inside(0, ci.ac2[, 1])
        }
      }
    }

    ############################################################################
    # Return results
    return(c(nc = nc,
      ah.orc = ah.orc, cg.orc = cg.orc, il.orc = il.orc,
      nr.orc = nr.orc, fp.orc = fp.orc, fn.orc = fn.orc,
      pd.orc = pd.orc, pc.orc = pc.orc, pi.orc = pi.orc,
      ah.ore = ah.ore, cg.ore = cg.ore, il.ore = il.ore,
      nr.ore = nr.ore, fp.ore = fp.ore, fn.ore = fn.ore,
      pd.ore = pd.ore, pc.ore = pc.ore, pi.ore = pi.ore,
      ah.nle = ah.nle, cg.nle = cg.nle, il.nle = il.nle,
      nr.nle = nr.nle, fp.nle = fp.nle, fn.nle = fn.nle,
      pd.nle = pd.nle, pc.nle = pc.nle, pi.nle = pi.nle,
      ah.dle = ah.dle, cg.dle = cg.dle, il.dle = il.dle,
      nr.dle = nr.dle, fp.dle = fp.dle, fn.dle = fn.dle,
      pd.dle = pd.dle, pc.dle = pc.dle, pi.dle = pi.dle,
      ah.ssl = ah.ssl, cg.ssl = cg.ssl, il.ssl = il.ssl,
      nr.ssl = nr.ssl, fp.ssl = fp.ssl, fn.ssl = fn.ssl,
      pd.ssl = pd.ssl, pc.ssl = pc.ssl, pi.ssl = pi.ssl,
      ah.pcr = ah.pcr, cg.pcr = cg.pcr, il.pcr = il.pcr,
      nr.pcr = nr.pcr, fp.pcr = fp.pcr, fn.pcr = fn.pcr,
      pd.pcr = pd.pcr, pc.pcr = pc.pcr, pi.pcr = pi.pcr,
      ah.bc1 = ah.bc1, cg.bc1 = cg.bc1, il.bc1 = il.bc1,
      nr.bc1 = nr.bc1, fp.bc1 = fp.bc1, fn.bc1 = fn.bc1,
      pd.bc1 = pd.bc1, pc.bc1 = pc.bc1, pi.bc1 = pi.bc1,
      ah.bcT = ah.bcT, cg.bcT = cg.bcT, il.bcT = il.bcT,
      nr.bcT = nr.bcT, fp.bcT = fp.bcT, fn.bcT = fn.bcT,
      pd.bcT = pd.bcT, pc.bcT = pc.bcT, pi.bcT = pi.bcT,
      ah.bcI = ah.bcI, cg.bcI = cg.bcI, il.bcI = il.bcI,
      nr.bcI = nr.bcI, fp.bcI = fp.bcI, fn.bcI = fn.bcI,
      pd.bcI = pd.bcI, pc.bcI = pc.bcI, pi.bcI = pi.bcI,
      ah.acm = ah.acm, cg.acm = cg.acm, il.acm = il.acm,
      nr.acm = nr.acm, fp.acm = fp.acm, fn.acm = fn.acm,
      pd.acm = pd.acm, pc.acm = pc.acm, pi.acm = pi.acm,
      ah.ac1 = ah.ac1, cg.ac1 = cg.ac1, il.ac1 = il.ac1,
      nr.ac1 = nr.ac1, fp.ac1 = fp.ac1, fn.ac1 = fn.ac1,
      pd.ac1 = pd.ac1, pc.ac1 = pc.ac1, pi.ac1 = pi.ac1,
      ah.ac2 = ah.ac2, cg.ac2 = cg.ac2, il.ac2 = il.ac2,
      nr.ac2 = nr.ac2, fp.ac2 = fp.ac2, fn.ac2 = fn.ac2,
      pd.ac2 = pd.ac2, pc.ac2 = pc.ac2, pi.ac2 = pi.ac2,
      ah.bm1 = ah.bm1, cg.bm1 = cg.bm1, il.bm1 = il.bm1,
      nr.bm1 = nr.bm1, fp.bm1 = fp.bm1, fn.bm1 = fn.bm1,
      pd.bm1 = pd.bm1, pc.bm1 = pc.bm1, pi.bm1 = pi.bm1,
      ah.bm2 = ah.bm2, cg.bm2 = cg.bm2, il.bm2 = il.bm2,
      nr.bm2 = nr.bm2, fp.bm2 = fp.bm2, fn.bm2 = fn.bm2,
      pd.bm2 = pd.bm2, pc.bm2 = pc.bm2, pi.bm2 = pi.bm2,
      ah.nmA = ah.nmA, cg.nmA = cg.nmA, il.nmA = il.nmA,
      nr.nmA = nr.nmA, fp.nmA = fp.nmA, fn.nmA = fn.nmA,
      pd.nmA = pd.nmA, pc.nmA = pc.nmA, pi.nmA = pi.nmA,
      ah.nmB = ah.nmB, cg.nmB = cg.nmB, il.nmB = il.nmB,
      nr.nmB = nr.nmB, fp.nmB = fp.nmB, fn.nmB = fn.nmB,
      pd.nmB = pd.nmB, pc.nmB = pc.nmB, pi.nmB = pi.nmB,
      ah.n2A = ah.n2A, cg.n2A = cg.n2A, il.n2A = il.n2A,
      nr.n2A = nr.n2A, fp.n2A = fp.n2A, fn.n2A = fn.n2A,
      pd.n2A = pd.n2A, pc.n2A = pc.n2A, pi.n2A = pi.n2A,
      ah.n2B = ah.n2B, cg.n2B = cg.n2B, il.n2B = il.n2B,
      nr.n2B = nr.n2B, fp.n2B = fp.n2B, fn.n2B = fn.n2B,
      pd.n2B = pd.n2B, pc.n2B = pc.n2B, pi.n2B = pi.n2B,
      th0.ep = th0.ep, th1.ep = th1.ep, th0.eb = th0.eb,
      th1.eb = th1.eb, t1.bm1 = t1.bm1, t2.bm1 = t2.bm1,
      t3.bm1 = t3.bm1, t1.bm2 = t1.bm2, t2.bm2 = t2.bm2,
      t3.bm2 = t3.bm2, t1.nmA = t1.nmA, t2.nmA = t2.nmA,
      t3.nmA = t3.nmA, t1.nmB = t1.nmB, t2.nmB = t2.nmB, t3.nmB = t3.nmB))
  }; cat(' | Done.\n')

  # Finish
  return(invisible(sim))
}; sim.teff <- compiler::cmpfun(sim.teff)

################################################################################
# Analogous of sim.teff (MULTIPLE TREATMENT version)
sim.teff.mt <- function(a, n, p, S, phiy, phid, N, do.bac = FALSE,
  do.ndl = FALSE, do.bma = FALSE, do.pcr = FALSE, do.ssl = FALSE,
  do.nme = FALSE, do.nm2 = FALSE, eta = 0.05, bma.mprior = 'bbin',
  lpen = 'lambda.min', mod1 = 'ginv', th.prior = 'unif', beta.prior = 'nlp',
  rho.min = NULL, rho.max = NULL, R = 1e5, max.mod = Inf, ncores = NA) {
################################################################################
  # Coefficients
  nt <- length(a)
  by <- c(rep(1, 4 * 5), rep(0, p - 4 * 5))
  bd <- matrix(0, nrow = p, ncol = nt)

  # Function shortcuts
  ms <- mombf::modelSelection
  igp <- mombf::igprior(alpha = 0.01, lambda = 0.01)
  zup <- mombf::zellnerprior(tau = n)
  nlp <- mombf::momprior(tau = 0.348)
  if (bma.mprior == 'bbin') {
    mprior <- mombf::modelbbprior(alpha.p = 1, beta.p = 1)
  } else if (bma.mprior == 'unif') {
    mprior <- mombf::modelunifprior()
  } else if (bma.mprior == 'binom') {
    ph <- (length(which(by != 0)) + length(which(a != 0))) / (p + nt)
    mprior <- mombf::modelbinomprior(p = ph)
  }

  # Parallelise
  if (is.na(ncores)) { ncores <- parallel::detectCores() - 1 }
  doMC::registerDoMC(cores = ncores)
  cat('Running on', ncores, 'cores:\n')

  # Simulations
  sim <- foreach(r = 1:N, .errorhandling = 'pass', .combine = 'rbind') %dopar% {
    # Print progress
    cat('\r* Iteration:', r, 'of', N)
    cat(' (', sprintf('%01.1f', round(100 * r / N, 1)), '%)', sep = '')

    # Simulate data
    set.seed(100 * r)
    X <- mvtnorm::rmvnorm(n = n, mean = rep(0, nrow(S)), sigma = S)
    D <- matrix(nrow = n, ncol = nt)
    for (t in 1:nt) {
      # Confounders correlating to treatment
      v1 <- t * 4
      bd[v1 - 0:3, t] <- 1

      # Other treatments
      bd[5 * 4 + 1:v1, t] <- 1

      # Simulate treatments
      epsd <- rnorm(n, mean = 0, sd = sqrt(phid[t]))
      D[, t] <- X %*% bd[, t] + epsd
    }
    # Simulate output
    epsy <- rnorm(n, mean = 0, sd = sqrt(phiy))
    y <- D %*% a + X %*% by + epsy
    Z <- cbind(D, X)

    # Oracle model
    orc.model <- c(paste('d', 1:nt, sep = ''),
                   paste('x', which(by != 0), sep = ''))
    if (sum(by != 0) == 0) { orc.model <- orc.model[-length(orc.model)] }

    ############################################################################
    # Oracle OLS estimation
    if (sum(by != 0) == 0) {
      orc.lm <- lm(y ~ D - 1)
    } else {
      orc.lm <- lm(y ~ D + X[, which(by != 0)] - 1)
    }
    ah.orc <- unname(coef(orc.lm)[1:nt])
    nr.orc <- nt + length(which(by != 0))  # No. regressors
    fp.orc <- 0  # False positives
    fn.orc <- sum(summary(orc.lm)[['coefficients']][-(1:nt), 4] > eta)  # False negatives
    ci.orc <- confint(orc.lm, level = 1 - eta)[1:nt, ]
    il.orc <- unname(apply(ci.orc, 1, function(x) { diff(range(x)) }))  # Interval length
    cg.orc <- sapply(1:nt, function(x) { inside(a[x], ci.orc[x, ]) })
    pd.orc <- sapply(1:nt, function(x) { !inside(0, ci.orc[x, ]) })

    # Oracle estimation including whatever affects D
    actT <- unique(as.vector(unlist(apply(bd, 2, function(x) { which(x != 0) }))))
    incX <- unique(c(which(by != 0), actT))
    tinc <- c(1:nt, which(by != 0) + nt)
    ore.lm <- lm(y ~ D + X[, incX] - 1)
    ah.ore <- unname(coef(ore.lm)[1:nt])
    nr.ore <- length(incX) + nt
    fp.ore <- sum(! actT %in% which(by != 0))  # FPs
    fn.ore <- sum(summary(ore.lm)[['coefficients']][tinc[-(1:nt)], 4] > eta)  # FNs
    ci.ore <- confint(ore.lm, level = 1 - eta)[1:nt, ]
    il.ore <- unname(apply(ci.ore, 1, function(x) { diff(range(x)) }))  # Interval length
    cg.ore <- sapply(1:nt, function(x) { inside(a[x], ci.ore[x, ]) })
    pd.ore <- sapply(1:nt, function(x) { !inside(0, ci.ore[x, ]) })

    ############################################################################
    # Lasso-based methods
    if (do.ndl == TRUE) {
      # Naive Lasso
      nl.fit <- cv.glmnet(x = Z, y = y, intercept = FALSE, standardize = FALSE)
      aux <- which(coef(nl.fit, s = nl.fit[[lpen]]) != 0) - nt - 1
      nle.model <- paste('x', aux[aux > 0], sep = '')
      nle.model <- c(paste('d', aux[aux <= 0] + 2, sep = ''), nle.model)
      ah.nle <- coef(nl.fit, s = nl.fit[[lpen]])[1 + 1:nt]
      nr.nle <- length(nle.model)
      fp.nle <- sum(! nle.model %in% orc.model)
      fn.nle <- sum(! orc.model %in% nle.model)
      il.nle <- NA
      cg.nle <- NA
      pd.nle <- coef(nl.fit, s = nl.fit[[lpen]])[1 + 1:nt] != 0

      # Attempt inference on treatment using method in Lee et al. (2016)
      nl.beta <- as.vector(coef(nl.fit, s = nl.fit[[lpen]]))[-1]
      nl.sel <- nl.beta != 0
      nsigma <- sd((y - as.matrix(Z[, nl.sel]) %*% nl.beta[nl.sel != 0]))
      nsigma <- nsigma * sqrt(n / (n - sum(nl.sel) - 1))
      nl.inf <- selectiveInference::fixedLassoInf(x = Z, y = y,
        beta = coef(nl.fit, x = Z, y = y, s = lpen)[-1],
        lambda = n * nl.fit[[lpen]], sigma = nsigma, alpha = 0.05)
      tvar <- which(nl.inf[['vars']] %in% 1:nt)
      ci.nle <- matrix(0, nrow = nt, ncol = 2)
      if (length(tvar) > 0) {
        ci.nle[tvar, ] <- nl.inf[['ci']][tvar, ]
      }
      il.nle <- unname(apply(ci.nle, 1, function(x) { diff(range(x)) }))
      cg.nle <- sapply(1:nt, function(x) { inside(a[x], ci.nle[x, ]) })

      # # Resulting OLS
      # dl.fit <- lm(y ~ D + X[, inc.dl] - 1)
      dl.fit <- try(rlassoEffects(y = y, x = Z, index = 1:nt, post = TRUE,
        intercept = FALSE, method = 'double selection'), silent = TRUE)
      if (class(dl.fit) != 'try-error') {
        ah.dle <- unname(coef(dl.fit)[1:nt])
        nr.dle <- length(coef(dl.fit))  # MIGHT NEED dl.fti$selection.matrix
        dle.model <- paste('d', 1:nt, sep = '')
        #dle.model <- c(paste('d', 1:nt, sep = ''), paste('x', inc.dl, sep = ''))
        fp.dle <- sum(! dle.model %in% orc.model)
        fn.dle <- sum(! orc.model %in% dle.model)
        ci.dle <- confint(dl.fit, level = 1 - eta)[1:nt, ]
        il.dle <- apply(ci.dle, 1, function(x) { diff(range(x)) })
        cg.dle <- sapply(1:nt, function(x) { inside(a[x], ci.dle[x, ]) })
        pd.dle <- sapply(1:nt, function(x) { !inside(0, ci.dle[x, ]) })
      } else {
        ah.dle <- nr.dle <- fp.dle <- fn.dle <- il.dle <- cg.dle <- pd.dle <- NA
      }
    } else {
      ah.nle <- nr.nle <- fp.nle <- fn.nle <- il.nle <- cg.nle <- pd.nle <- NA
      ah.dle <- nr.dle <- fp.dle <- fn.dle <- il.dle <- cg.dle <- pd.dle <- NA
    }

    ############################################################################
    # ACPME (goes on BAC output slots)
    ah.bc1 <- nr.bc1 <- fp.bc1 <- fn.bc1 <- il.bc1 <- cg.bc1 <- pd.bc1 <- NA
    ah.bcT <- nr.bcT <- fp.bcT <- fn.bcT <- il.bcT <- cg.bcT <- pd.bcT <- NA
    ah.bcI <- nr.bcI <- fp.bcI <- fn.bcI <- il.bcI <- cg.bcI <- pd.bcI <- NA
    if (do.bac == TRUE & (p + nt) <= n) {
      # Run method
      sink(file = paste(TMPDIR, 'tmp_acpme.txt', sep = ''))  # No silent option
      acm.fit <- try(acpme(y = y, Z = D, C = X, niter = 1e4), silent = TRUE)
      sink()

      if (class(acm.fit) != 'try-error') {
        # BMA point estimate
        ah.bcI <- colMeans(acm.fit[['beta']])    

        # Model size
        orc.idx <- as.numeric(substr( orc.model[grep('x', orc.model)], 2,
          nchar(orc.model[grep('x', orc.model)])))
        nr.bcI <- mean(rowSums(acm.fit[['alpha']]))
        fp.bcI <- mean(apply(acm.fit[['alpha']], 1, function(x)
          sum(! unname(which(x == 1)) %in% orc.idx)
        ))
        fn.bcI <- mean(apply(acm.fit[['alpha']], 1, function(x)
          sum(! orc.idx %in% unname(which(x == 1)))
        ))

        # Inclusion indicators
        pr.eta <- c(eta / 2, 1 - eta / 2)
        ci.bcI <- unname(apply(acm.fit[['beta']], 2, quantile, probs = pr.eta))
        il.bcI <- apply(ci.bcI, 2, function(x) { diff(range(x)) })
        cg.bcI <- sapply(1:nt, function(x) { inside(a[x], ci.bcI[, x]) })
        pd.bcI <- sapply(1:nt, function(x) { ! inside(0, ci.bcI[, x]) })
      }
    }

    ############################################################################
    # Penalised Credible Regions (Wilson and Reich, 2014)
    if (do.pcr == TRUE) {
    #if (do.pcr == TRUE & nrow(Z) > ncol(Z)) {
      # Choose based on BIC
      sink(file = paste(TMPDIR, 'tmp_esborrar', r, '.txt', sep = ''))  # No silent option :_)
      pcr.fit <- try(BayesPen::BayesPen.lm.confounders(y, D, X), silent = TRUE)
      if (class(pcr.fit) == 'try-error') {
        pcr.fit <- BayesPen.lm.confounders2(y, D, X)
      }
      sink(); unlink(paste(TMPDIR, 'tmp_esborrar', r, '.txt', sep = ''))
      #print(nrow(pcr.fit[['coefs']]))
      #print(pcr.fit[['joint.path']])

      # Choose minimum BIC
      prpcr <- mean(y) + cbind(D, X) %*% t(pcr.fit[['coefs']])
      nC <- apply(pcr.fit[['coefs']], 1, function(x) sum(x != 0))
      nr <- nrow(X)
      aux <- apply(prpcr, 2, function(x) { x - y })
      bic <- nr * log(colSums(aux ** 2) / length(y)) +
      #bic <- nr * log(colSums((y - prpcr)^2) / length(y)) +
             nr * (log(2 * pi) + 1) + log(nr) * nC
      wr14.refit.coef <- pcr.fit[['coefs']][which.min(bic), ]
      names(wr14.refit.coef) <- c(paste0('d', 1:nt), paste0('x', 1:p))
      wr14.model <- names(wr14.refit.coef)[which(wr14.refit.coef != 0)]

      # Relevant values
      ah.pcr <- unname(wr14.refit.coef[1:nt])
      nr.pcr <- sum(wr14.refit.coef != 0)
      fp.pcr <- sum(! wr14.model %in% orc.model)
      fn.pcr <- sum(! orc.model %in% wr14.model)
      cg.pcr <- NA
      il.pcr <- NA
      pd.pcr <- sapply(ah.pcr, function(x) { ifelse(x == 0, 0, 1) })
    } else {
      ah.pcr <- nr.pcr <- fp.pcr <- fn.pcr <- il.pcr <- cg.pcr <- pd.pcr <- NA
    }

    ############################################################################
    # BMA on y
    if (do.bma == TRUE) {
      # BMS
      bms1 <- ms(y, Z, priorCoef = zup, priorDelta = mprior, priorVar = igp,
                 niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)
      bms2 <- ms(y, Z, priorCoef = nlp, priorDelta = mprior, priorVar = igp,
                 niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)

      # BMA
      gibbs1 <- mombf::rnlp(msfit = bms1, priorCoef = zup, niter = 1e4)
      gibbs2 <- mombf::rnlp(msfit = bms2, priorCoef = nlp, niter = 1e4)
      ah.bm1 <- colMeans(gibbs1[, 1 + 1:nt])
      ah.bm2 <- colMeans(gibbs2[, 1 + 1:nt])

      # Model size
      aux1 <- mombf::postProb(bms1)
      aux2 <- mombf::postProb(bms2)
      naux1 <- lapply(strsplit(as.character(aux1[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - nt, sep = '')
        for (l in 1:nt) {
          ans[which(ans == paste0('x', l - nt))] <- paste0('d', l)
        }
        return(ans)
      })
      naux2 <- lapply(strsplit(as.character(aux2[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - nt, sep = '')
        for (l in 1:nt) {
          ans[which(ans == paste0('x', l - nt))] <- paste0('d', l)
        }
        return(ans)
      })
      fps1 <- unlist(lapply(naux1, function(x) { sum(! x %in% orc.model) }))
      fps2 <- unlist(lapply(naux2, function(x) { sum(! x %in% orc.model) }))
      fns1 <- unlist(lapply(naux1, function(x) { sum(! orc.model %in% x) }))
      fns2 <- unlist(lapply(naux2, function(x) { sum(! orc.model %in% x) }))
      nr.bm1 <- sum(unlist(lapply(
        strsplit(as.character(aux1[, 1]), ','), length)) * aux1[, 3])
      nr.bm2 <- sum(unlist(lapply(
        strsplit(as.character(aux2[, 1]), ','), length)) * aux2[, 3])
      fp.bm1 <- sum(fps1 * aux1[, 3])
      fp.bm2 <- sum(fps2 * aux2[, 3])
      fn.bm1 <- sum(fns1 * aux1[, 3])
      fn.bm2 <- sum(fns2 * aux2[, 3])

      # Inclusion indicators
      pr.eta <- c(eta / 2, 1 - eta / 2)
      ci.bm1 <- unname(apply(gibbs1[, 1 + 1:nt], 2, quantile, probs = pr.eta))
      ci.bm2 <- unname(apply(gibbs2[, 1 + 1:nt], 2, quantile, probs = pr.eta))
      il.bm1 <- apply(ci.bm1, 2, function(x) { diff(range(x)) })
      il.bm2 <- apply(ci.bm2, 2, function(x) { diff(range(x)) })
      cg.bm1 <- sapply(1:nt, function(x) { inside(a[x], ci.bm1[, x]) })
      cg.bm2 <- sapply(1:nt, function(x) { inside(a[x], ci.bm2[, x]) })
      pd.bm1 <- unname(bms1[['margpp']][1:nt])
      pd.bm2 <- unname(bms2[['margpp']][1:nt])
      #pd.nmA <- sapply(1:nt, function(x) { !inside(0, ci.nmA[, x]) })

      # MAP models
      map1 <- as.character(aux1[1, 1])
      map2 <- as.character(aux2[1, 1])
      t1.bm1 <- as.numeric(binary.id(map1, p = ncol(X))[1:nt] == 1)
      t1.bm2 <- as.numeric(binary.id(map2, p = ncol(X))[1:nt] == 1)
      t2.bm1 <- as.numeric(unname(bms1[['margpp']][1:nt] >= 0.5))
      t2.bm2 <- as.numeric(unname(bms2[['margpp']][1:nt] >= 0.5))
      t3.bm1 <- as.numeric(unname(bms1[['margpp']][1:nt] >= 0.95))
      t3.bm2 <- as.numeric(unname(bms2[['margpp']][1:nt] >= 0.95))
    } else {
      ah.bm1 <- nr.bm1 <- fp.bm1 <- fn.bm1 <- il.bm1 <- cg.bm1 <- pd.bm1 <- NA
      ah.bm2 <- nr.bm2 <- fp.bm2 <- fn.bm2 <- il.bm2 <- cg.bm2 <- pd.bm2 <- NA
      t1.bm1 <- t1.bm2 <- t2.bm1 <- t2.bm2 <- t3.bm1 <- t3.bm2 <- NA
    }

    ############################################################################
    # New SSL (2019)
    ah.ssl <- nr.ssl <- fp.ssl <- fn.ssl <- il.ssl <- cg.ssl <- pd.ssl <- NA

    ############################################################################
    # Our new method (using Expectation Propagation)
    if (do.nme == TRUE) {
      # Expectation Propagation
      nme.fit <- new.method.mt(y, D, X, th.search = 'EP', mod1 = mod1,
        th.prior = th.prior, beta.prior = beta.prior, rho.min = rho.min,
        rho.max = rho.max, R = R, max.mod = max.mod)

      # Record estimated parameters
      th.ep <- nme.fit[['th.hat']]

      # Model size
      aux <- nme.fit[['pprob']]
      naux <- lapply(strsplit(as.character(aux[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - nt, sep = '')
        for (l in 1:nt) {
          ans[which(ans == paste0('x', l - nt))] <- paste0('d', l)
        }
        return(ans)
      })
      fps <- unlist(lapply(naux, function(x) { sum(! x %in% orc.model) }))
      fns <- unlist(lapply(naux, function(x) { sum(! orc.model %in% x) }))
      nr.nmA <- sum(unlist(lapply(
        strsplit(as.character(aux[, 1]), ','), length)) * aux[, 2])
      fp.nmA <- sum(fps * aux[, 2])
      fn.nmA <- sum(fns * aux[, 2])

      # Numerical summaries (NLP)
      noms <- paste(100 * c(eta / 2, 1 - eta / 2), '%', sep = '')
      ah.nmA <- nme.fit[['teff.nlp']]
      ci.nmA <- unname(nme.fit[['teff.d.nlp']][noms, ])
      il.nmA <- apply(ci.nmA, 2, function(x) { diff(range(x)) })
      cg.nmA <- sapply(1:nt, function(x) { inside(a[x], ci.nmA[, x]) })
      pd.nmA <- sapply(1:nt, function(x) { !inside(0, ci.nmA[, x]) })

      # Plan B (Gaussian prior)
      ah.nmB <- NA#nme.fit[['teff.zup']]
      ci.nmB <- NA#unname(nme.fit[['teff.d.zup']][noms])
      il.nmB <- NA#diff(range(ci.nmB))
      cg.nmB <- NA#inside(a, ci.nmB)
      nr.nmB <- NA#nr.nmA
      fp.nmB <- NA#fp.nmA
      fn.nmB <- NA#fn.nmA
      pd.nmB <- NA#pd.nmA

      # MAP model
      mapA <- as.character(nme.fit[['pprob']][1, 1])
      t1.nmA <- as.numeric(binary.id(mapA, p = ncol(X))[1:nt] == 1)
      t2.nmA <- as.numeric(unname(nme.fit[['mpprob']][1:nt] >= 0.5))
      t3.nmA <- as.numeric(unname(nme.fit[['mpprob']][1:nt] >= 0.95))
    } else {
      ah.nmA <- nr.nmA <- fp.nmA <- fn.nmA <- il.nmA <- cg.nmA <- pd.nmA <- NA
      ah.nmB <- nr.nmB <- fp.nmB <- fn.nmB <- il.nmB <- cg.nmB <- pd.nmB <- NA
      th.ep <- NA
      t1.nmA <- t2.nmA <- t3.nmA <- NA
    }

    ############################################################################
    # Our new method (using Empirical Bayes)
    if (do.nm2 == TRUE) {
      # Recycle initial MCMC if we've done it before
      if (do.nme == TRUE) {
        bvs.fit0 <- nme.fit[['init.msfit']]
        th.EP <- nme.fit[['th.hat']]
        EP.is <- nme.fit[['EP.is']]
      } else {
        bvs.fit0 <- NULL
        th.EP <- NULL
        EP.is <- NULL
      }

      # Empirical Bayes
      nme.fit2 <- new.method.mt(y, D, X, th.search = 'EB', mod1 = mod1,
        th.prior = th.prior, beta.prior = beta.prior, rho.min = rho.min,
        rho.max = rho.max, R = R, max.mod = max.mod, bvs.fit0 = bvs.fit0,
        th.EP = th.EP, EP.is = EP.is)

      # Record estimated parameters
      th.eb <- nme.fit2[['th.hat']]

      # Model size
      aux <- nme.fit2[['pprob']]
      naux <- lapply(strsplit(as.character(aux[, 1]), ','), function(x) {
        ans <- paste('x', as.numeric(x) - nt, sep = '')
        for (l in 1:nt) {
          ans[which(ans == paste0('x', l - nt))] <- paste0('d', l)
        }
        return(ans)
      })
      fps <- unlist(lapply(naux, function(x) { sum(! x %in% orc.model) }))
      fns <- unlist(lapply(naux, function(x) { sum(! orc.model %in% x) }))
      nr.n2A <- sum(unlist(lapply(
        strsplit(as.character(aux[, 1]), ','), length)) * aux[, 2])
      fp.n2A <- sum(fps * aux[, 2])
      fn.n2A <- sum(fns * aux[, 2])

      # Numerical summaries (NLP)
      noms <- paste(100 * c(eta / 2, 1 - eta / 2), '%', sep = '')
      ah.n2A <- nme.fit2[['teff.nlp']]
      ci.n2A <- unname(nme.fit2[['teff.d.nlp']][noms, ])
      il.n2A <- apply(ci.n2A, 2, function(x) { diff(range(x)) })
      cg.n2A <- sapply(1:nt, function(x) { inside(a[x], ci.n2A[, x]) })
      pd.n2A <- sapply(1:nt, function(x) { !inside(0, ci.n2A[, x]) })

      # Plan B (Gaussian prior)
      ah.n2B <- NA#nme.fit2[['teff.zup']]
      ci.n2B <- NA#unname(nme.fit2[['teff.d.zup']][noms])
      il.n2B <- NA#diff(range(ci.n2B))
      cg.n2B <- NA#inside(a, ci.n2B)
      nr.n2B <- NA#nr.n2A
      fp.n2B <- NA#fp.n2A
      fn.n2B <- NA#fn.n2A
      pd.n2B <- NA#pd.n2A

      # MAP model
      mapB <- as.character(nme.fit2[['pprob']][1, 1])
      t1.nmB <- as.numeric(binary.id(mapB, p = ncol(X))[1:nt] == 1)
      t2.nmB <- as.numeric(unname(nme.fit2[['mpprob']][1:nt] >= 0.5))
      t3.nmB <- as.numeric(unname(nme.fit2[['mpprob']][1:nt] >= 0.95))
    } else {
      ah.n2A <- nr.n2A <- fp.n2A <- fn.n2A <- il.n2A <- cg.n2A <- pd.n2A <- NA
      ah.n2B <- nr.n2B <- fp.n2B <- fn.n2B <- il.n2B <- cg.n2B <- pd.n2B <- NA
      th.eb <- NA
      t1.nmB <- t2.nmB <- t3.nmB <- NA
    }

    # Return results
    return(list(nt = nt,#nc = nc,
      ah.orc = ah.orc, cg.orc = cg.orc, il.orc = il.orc, nr.orc = nr.orc,
      fp.orc = fp.orc, fn.orc = fn.orc, pd.orc = pd.orc,
      ah.ore = ah.ore, cg.ore = cg.ore, il.ore = il.ore, nr.ore = nr.ore,
      fp.ore = fp.ore, fn.ore = fn.ore, pd.ore = pd.ore,
      ah.nle = ah.nle, cg.nle = cg.nle, il.nle = il.nle, nr.nle = nr.nle,
      fp.nle = fp.nle, fn.nle = fn.nle, pd.nle = pd.nle,
      ah.dle = ah.dle, cg.dle = cg.dle, il.dle = il.dle, nr.dle = nr.dle,
      fp.dle = fp.dle, fn.dle = fn.dle, pd.dle = pd.dle,
      ah.ssl = ah.ssl, cg.ssl = cg.ssl, il.ssl = il.ssl, nr.ssl = nr.ssl,
      fp.ssl = fp.ssl, fn.ssl = fn.ssl, pd.ssl = pd.ssl,
      ah.pcr = ah.pcr, cg.pcr = cg.pcr, il.pcr = il.pcr, nr.pcr = nr.pcr,
      fp.pcr = fp.pcr, fn.pcr = fn.pcr, pd.pcr = pd.pcr,
      ah.bc1 = ah.bc1, cg.bc1 = cg.bc1, il.bc1 = il.bc1, nr.bc1 = nr.bc1,
      fp.bc1 = fp.bc1, fn.bc1 = fn.bc1, pd.bc1 = pd.bc1,
      ah.bcT = ah.bcT, cg.bcT = cg.bcT, il.bcT = il.bcT, nr.bcT = nr.bcT,
      fp.bcT = fp.bcT, fn.bcT = fn.bcT, pd.bcT = pd.bcT,
      ah.bcI = ah.bcI, cg.bcI = cg.bcI, il.bcI = il.bcI, nr.bcI = nr.bcI,
      fp.bcI = fp.bcI, fn.bcI = fn.bcI, pd.bcI = pd.bcI,
      ah.bm1 = ah.bm1, cg.bm1 = cg.bm1, il.bm1 = il.bm1, nr.bm1 = nr.bm1,
      fp.bm1 = fp.bm1, fn.bm1 = fn.bm1, pd.bm1 = pd.bm1,
      ah.bm2 = ah.bm2, cg.bm2 = cg.bm2, il.bm2 = il.bm2, nr.bm2 = nr.bm2,
      fp.bm2 = fp.bm2, fn.bm2 = fn.bm2, pd.bm2 = pd.bm2,
      ah.nmA = ah.nmA, cg.nmA = cg.nmA, il.nmA = il.nmA, nr.nmA = nr.nmA,
      fp.nmA = fp.nmA, fn.nmA = fn.nmA, pd.nmA = pd.nmA,
      ah.nmB = ah.nmB, cg.nmB = cg.nmB, il.nmB = il.nmB, nr.nmB = nr.nmB,
      fp.nmB = fp.nmB, fn.nmB = fn.nmB, pd.nmB = pd.nmB,
      ah.n2A = ah.n2A, cg.n2A = cg.n2A, il.n2A = il.n2A, nr.n2A = nr.n2A,
      fp.n2A = fp.n2A, fn.n2A = fn.n2A, pd.n2A = pd.n2A,
      ah.n2B = ah.n2B, cg.n2B = cg.n2B, il.n2B = il.n2B, nr.n2B = nr.n2B,
      fp.n2B = fp.n2B, fn.n2B = fn.n2B, pd.n2B = pd.n2B,
      th.ep = th.ep, th.eb = th.eb,
      t1.bm1 = t1.bm1, t2.bm1 = t2.bm1, t3.bm1 = t3.bm1,
      t1.bm2 = t1.bm2, t2.bm2 = t2.bm2, t3.bm2 = t3.bm2,
      t1.nmA = t1.nmA, t2.nmA = t2.nmA, t3.nmA = t3.nmA,
      t1.nmB = t1.nmB, t2.nmB = t2.nmB, t3.nmB = t3.nmB))
  }; cat(' | Done.\n')

  # Finish
  return(invisible(sim))
}; sim.teff.mt.rand <- compiler::cmpfun(sim.teff.mt)
# END OF SCRIPT
