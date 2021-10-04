################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(FCNDIR, 'functions_newmethod.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * binary.id
# * compact.id
# * pinc
# * pinc.mt
# * grid.opt
# * grid.opt.mt
# * set.s
# * pprob.nm
# * pprob.nm.mt
# * bma.teff
# * bma.teff.mt
# * new.method
# * new.method.mt
################################################################################

################################################################################
# Converting model indexes to binary strings
binary.id <- function(x, p, nt = 1, conc = FALSE) {
################################################################################
  out <- rep(0, p + nt)
  out[as.numeric(unlist(strsplit(x, ',')))] <- 1
  if (conc == TRUE) { out <- paste(out, collapse = '') }
  return(out)
}

################################################################################
# Converting binary strings to model indexes
compact.id <- function(x) {
################################################################################
  paste(which(unlist(strsplit(as.character(x), '')) == 1), collapse = ',')
}

################################################################################
# Prior inclusion probability function (conditional on coefficient size)
pinc <- function(x, th0, th1, rho.min = 0, rho.max = 1, squared = FALSE) {
################################################################################
  # If using abs(beta) or beta^2
  if (squared == TRUE) {
    probs <- (1 + exp(-(th0 + th1 * (x ** 2)))) ** (-1)
  } else {
    probs <- (1 + exp(-(th0 + th1 * abs(x)))) ** (-1)
  }
  
  # Limit upper and lower bounds
  probs <- rho.min + (rho.max - rho.min) * probs

  # End
  return(probs)
}

################################################################################
# Prior inclusion probability function (conditional on coefficient size)
pinc.mt <- function(x, th, rho.min = 0, rho.max = 1, squared = FALSE) {
################################################################################
  #  If there are interactions (CAREFUL: main effects + ORDERED interactions)
  th <- as.numeric(th)
  if (ncol(x) / (length(th) - 1) > 1) {
    th <- c(th, rep(th[-1], each = ncol(x) / (length(th) - 1) - 1))
  }

  # If using abs(beta) or beta^2
  if (squared == TRUE) {
    probs <- (1 + exp(-apply(x, 1, function(z) {
        sum(c(1, z^2) * th)
      })))^(-1)
  } else {
    probs <- (1 + exp(-apply(x, 1, function(z) {
        sum(c(1, abs(z)) * th)
      })))^(-1)
  }

  # Limit upper and lower bounds
  probs <- rho.min + (rho.max - rho.min) * probs

  # End
  return(probs)
}

################################################################################
# Grid optimizer for empirical Bayes theta value search
grid.opt <- function(G0, betad, pj1, method = 'EB', ws = NA, th.grid,
  th.prior = 'unif', rho.min = 0, rho.max = 1, V = diag(2),
  ret.plotly = FALSE) {
################################################################################
  opt.th <- th.grid
  if (th.prior == 'tunif') {
    opt.th <- opt.th[which(opt.th[, 1] >= -log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 2] >= -log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 1] <= log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 2] <= log(length(betad))), ]
  }

  #opt.th <- expand.grid(th0s, th1s)
  if (method == 'EB') {
    # Evaluate -log of objective function
    opt.th[, 3] <- exp(-apply(opt.th, 1, Of.EB,
      betad = betad, G0 = G0, ws = ws, th.prior = th.prior,
      rho.min = rho.min, rho.max = rho.max, V = V))
  } else if (method == 'EP') {
    opt.th[, 3] <- exp(-apply(opt.th, 1, Of.EP,
      betad = betad, pj1 = pj1, th.prior = th.prior,
      rho.min = rho.min, rho.max = rho.max, V = V))
  } else {
    stop('supported theta search methods are "EB" and "EP".')
  }

  # Optimal values
  th0.opt <- opt.th[which.max(opt.th[, 3]), 1]
  th1.opt <- opt.th[which.max(opt.th[, 3]), 2]

  if (ret.plotly == TRUE) {
    require(plotly)
    colnames(opt.th) <- c('theta0', 'theta1', 'postdens')
    titol <- paste('argmax: th0=', round(th0.opt, 1),
      '; th1=', round(th1.opt, 1), sep = '')
    fig <- plotly::plot_ly(opt.th, x = ~theta0, y = ~theta1, z = ~postdens,
      type = 'contour', colors = 'Greys',
      line = list(color = 'black')) %>% layout(title = titol, showlegend = FALSE)
  } else {
    fig <- NULL
  }

  # Output
  return(list(th0.opt = th0.opt, th1.opt = th1.opt, th.fig = fig))
}; grid.opt <- compiler::cmpfun(grid.opt)

################################################################################
# Grid optimizer for empirical Bayes theta search (MULTIPLE TREATMENT version)
grid.opt.mt <- function(G0, betad, pj1, method = 'EB', ws = NA, th.grid,
  th.prior = 'unif', rho.min = 0, rho.max = 1, V = diag(2),
  ret.plotly = FALSE) {
################################################################################
  opt.th <- th.grid
  if (th.prior == 'tunif') {
    opt.th <- opt.th[which(opt.th[, 1] >= -log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 2] >= -log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 1] <= log(length(betad))), ]
    opt.th <- opt.th[which(opt.th[, 2] <= log(length(betad))), ]
  }

  #opt.th <- expand.grid(th0s, th1s)
  if (method == 'EB') {
    # Evaluate -log of objective function
    opt.th[, ncol(opt.th) + 1] <- exp(-apply(opt.th, 1, Of.EB.mt,
      betad = betad, G0 = G0, ws = ws, th.prior = th.prior,
      rho.min = rho.min, rho.max = rho.max, V = V))
  } else if (method == 'EP') {
    opt.th[, ncol(opt.th) + 1] <- exp(-apply(opt.th, 1, Of.EP.mt,
      betad = betad, pj1 = pj1, th.prior = th.prior,
      rho.min = rho.min, rho.max = rho.max, V = V))
  } else {
    stop('supported theta search methods are "EB" and "EP".')
  }

  # Optimal values
  th.opt <- opt.th[which.max(opt.th[, ncol(opt.th)]), -ncol(opt.th)]
  th.opt <- as.numeric(th.opt)
  names(th.opt) <- paste('th', 1:length(th.opt) - 1, sep = '')

  # Output
  return(list(th.opt = th.opt))
}; grid.opt.mt <- compiler::cmpfun(grid.opt.mt)

################################################################################
# Find recommended variance value for a Gaussian prior (DEPRECATED - NOT USED)
set.s <- function(betad, cl = 0.99, eps = NA, rho01 = 0, bounds = c(0.01, 100),
  dens = 0.1, ptheta = 'gauss', tol = 1e-6, max.it = 100) {
################################################################################
  
  p <- length(betad) + 1
  mabj <- mean(abs(betad))
  done <- FALSE
  ss <- seq(bounds[1], bounds[2], dens)
  if (is.na(eps)) { eps <- 1/p }

  it <- 1
  while (done == FALSE) {
    tp <- rep(NA, length(ss))
    for (s in ss) {
      if (ptheta == 'gauss') {
        v0 <- -log(1/eps - 1) / sqrt(s * (1 + mabj^2 + 2 * mabj * rho01))
        tp[which(ss == s)] <- pnorm(-v0) - pnorm(v0)
      } else if (ptheta == 'cauchy') {
        v0 <- -log(1/eps - 1) / (s * (1 + mabj))
        tp[which(ss == s)] <- pcauchy(-v0) - pcauchy(v0)
      }
    }
    best <- which.min(abs(tp - cl))
    s.opt <- ss[best]
    t.opt <- tp[best]

    ss <- seq(s.opt - dens, s.opt + dens, dens / 10)
    dens <- dens / 10
    if (abs(t.opt - cl) < tol) { done <- TRUE }
    if (it == max.it) {
      s.opt <- NA
      done <- TRUE
    } else {
      it <- it + 1
    }
  }
  
  return(s.opt)
}; set.s <- compiler::cmpfun(set.s)

################################################################################
# Posterior model probabilities under CIL
pprob.nm <- function(y, d, X, mod1 = 'ginv', th.search = 'EB',
  th.prior = 'unif', beta.prior = 'nlp', rho.min = NULL, rho.max = NULL,
  tau = 0.348, max.mod = Inf, lpen = 'lambda.1se', eps = 1e-10, R = 1e5,
  bvs.fit0 = NULL, seed = NA) {
################################################################################
  # Renaming
  ms <- mombf::modelSelection
  igp <- mombf::igprior(0.01, 0.01)
  bbp <- mombf::modelbbprior(1, 1)
  nlp <- mombf::momprior(tau = tau)
  zup <- mombf::zellnerprior(tau = nrow(X))
  ufp <- mombf::modelunifprior()
  mlik <- mombf::nlpMarginal

  # Pre-computations
  p <- ncol(X)
  Z <- cbind(d, X)

  # Other fixed parameters
  if (all(is.null(c(rho.min, rho.max)))) {
    #rho.min <- 1 / ncol(Z)
    rho.min <- 1 / (ncol(Z)^2 + 1)
    rho.max <- 1 - rho.min
  }

  # Coefficient prior
  if (beta.prior == 'nlp') {
    cp <- nlp
  } else if (beta.prior == 'zup') {
    cp <- zup  
  } else {
    stop('currently supported "beta.prior" arg.: "nlp" (MOM) and "zup" (UIP)')
  }

  # Exposure family (support for Binomial and Gaussian ONLY)
  type.exp <- ifelse(length(unique(d)) == 2, 'binomial', 'gaussian')

  # MODULE 1: Estimate coefficients on exposure model ##########################
  if (mod1 == 'ginv') {
    betad0 <- as.numeric(pracma::pinv(t(X) %*% X) %*% t(X) %*% d)
  } else if (mod1 %in% c('bvs', 'bma', 'bms')) {
    bvs.fit <- ms(d, X, priorCoef = cp, priorVar = igp, priorDelta = bbp,
      niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)
    post.th <- colMeans(rnlp(msfit = bvs.fit, priorCoef = cp, niter = 1e4))
    betad0 <- unname(post.th[2:(length(post.th) - 1)])
    #phid <- unname(post.theta[p + 2])
  } else if (mod1 == 'lasso') {
    lasso.fit <- glmnet::cv.glmnet(X, d, intercept = FALSE)
    betad0 <- unname(coef(lasso.fit, s = lasso.fit[[lpen]])[-1, 1])
  } else if (mod1 == 'lasso_dml') {
    #rlasso.fit <- hdm::rlasso(d ~ X, post = FALSE, intercept = FALSE)
    rlasso.fit <- hdm::rlasso(d ~ X, post = TRUE, intercept = FALSE)
    betad0 <- rlasso.fit[['coefficients']]
  } else if (mod1 == 'lasso_bic') {
    lb.fit <- lasso.bic(d, X, intercept = FALSE, family = type.exp)
    betad0 <- unname(lb.fit[['coef']][-1])
  } else if (mod1 == 'ridge') {
    ridge.fit <- glmnet::cv.glmnet(X, d, alpha = 0, intercept = FALSE)
    betad0 <- unname(coef(ridge.fit, s = ridge.fit[['lambda.min']])[-1, ])
  } else {
    stop('currently unsupported "mod1" method.')
  }

  # Set final estimates
  betad <- betad0
  # # Scale betad (this is DEPRECATED and should be inactive *always*)
  # if (scale.delta == TRUE) {
  #   betad <- betad0 / max(abs(betad0))
  # } else {
  #   betad <- betad0
  # }

  # THETA SEARCH ###############################################################
  # Initial fit: th0 = 0; th1 = 0 (UNIFORM MODEL PRIOR)
  if (is.null(bvs.fit0)) {
    bvs.fit0 <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = ufp,
      niter = R, verbose = FALSE, center = FALSE, scale = FALSE)
  }

  # Posterior model probabilities
  pprobs <- postProb(bvs.fit0)[, c(1, 3)]

  # Marginal inclusion probabilities
  pj1 <- bvs.fit0[['margpp']]
  pj1[which(pj1 == 1)] <- 1 - eps  # In case there are too extreme values
  pj1[which(pj1 == 0)] <- eps

  # Set of visited models
  if (is.na(max.mod)) {  # Prevent a dictionary too big
    max.mod <- which(cumsum(pprobs[, ncol(pprobs)]) >= 0.9)[1]
  }
  G0 <- unname(sapply(
    as.character(pprobs[, 1]), binary.id, p = p, conc = TRUE))
  G0 <- G0[1:min(max.mod, length(G0))]

  # Find recommended value for prior theta variance
  if (th.prior == 'gauss') {
    # That leaves Pr(p_j(theta) \in (eps, 1-eps)) = 1-cl (correlation 0)
    s <- set.s(betad, eps = 0.05, cl = 0.1, rho01 = 0, dens = 0.01)
    Rm <- matrix(c(1, 0, 0, 1), ncol = 2)
    if (is.na(s)) {
      # Span of search
      s <- set.s(betad, eps = 0.05, cl = 0.1, rho01 = 0,
        bounds = c(100, 1e3), dens = 0.1)
    }
  } else if (th.prior %in% c('unif', 'tunif')) {
    s <- 1
    Rm <- diag(2)
  } else {
    stop('currently supported "th.prior" arg.: "unif" and "gauss"')
  }

  # Determine the values of hat{theta} #########################################
  # Set limits for optimisation (in principle: unbounded)
  if (th.prior == 'tunif') {
    ups <- c(log(ncol(Z) - 1), log(ncol(Z) - 1))  # p_j(theta, betad=0) < 1-1/p
    lws <- -ups  # p_j(theta, betad=0) > 1/p
    scale.prpr <- FALSE  # NO LONGER NEEDED
  } else {
    lws <- c(-Inf, -Inf)
    ups <- c(+Inf, +Inf)
  }

  # EP approximation (to set starting point)
  ths <- seq(-40, 40, 1)
  EP.is <- grid.opt(G0, betad, pj1, th.grid = expand.grid(ths, ths),
    method = 'EP', th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
    V = s * Rm)
  #orca(EP.is$th.fig, paste('Desktop/surrealism.pdf', sep = ''))

  # Load gradient functions and optimise
  st <- c(EP.is[['th0.opt']], EP.is[['th1.opt']])
  opt.EP <- nlminb(st, objective = Of.EP, gradient = Gf.EP, 
    lower = lws, upper = ups, pj1 = pj1, betad = betad,
    th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, V = s * Rm)
  # opt.EP <- optim(st, fn = Of.EP, gr = Gf.EP, method = 'L-BFGS-B',
  #   lower = lws, upper = ups, 
  #   betad = betad, pj1 = pj1, th.prior = th.prior, V = s * Rm)

  # Set values if there is convergence
  if (opt.EP[['convergence']] %in% 0:1) {
    th0.EP <- opt.EP[['par']][1]
    th1.EP <- opt.EP[['par']][2]
  } else {
    stop('theta search under "EP" did not converge.')
  }

  # Empirical Bayes search (starting at EP solution)
  if (th.search == 'EB') {
    # Model search at EP optimum
    auxprpr <- pinc(betad, th0 = th0.EP, th1 = th1.EP,
      rho.min = rho.min, rho.max = rho.max)
    auxprpr <- c(1/2, auxprpr)
    auxprpr[which(auxprpr == 1)] <- 1 - eps
    auxprpr[which(auxprpr == 0)] <- eps
    auxmp <- mombf::modelbinomprior(p = auxprpr)
    update.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = auxmp,
      niter = round(R / 3), verbose = FALSE, center = FALSE, scale = FALSE)

    # Append model set
    upprobs <- postProb(update.fit)[, c(1, 3)]
    G1 <- unname(sapply(
      as.character(upprobs[, 1]), binary.id, p = p, conc = TRUE))
    G0 <- unique(c(G0, G1[1:min(max.mod, length(G1))]))

    # Pre-compute marginal likelihoods
    ws <- unname(unlist(sapply(G0, function(m) {
      js <- which(unlist(strsplit(m, '')) == 1)
      return(mlik(js, y, Z, priorCoef = cp, logscale = TRUE))
    })))

    # Numerical problems with MLs that are too small
    if (all(ws < -600)) {  # We start having problems at exp(-750)
      ct <- abs(max(ws) + 350)  # Constant
      ws <- exp(ws + ct)  # Multiply by constant (NO EFFECT ON ARG MAX)
    } else {
      ws <- exp(ws)
    }

    # Optimise starting at the EP optimum
    st <- c(th0.EP, th1.EP)
    opt.EB <- try(nlminb(st, objective = Of.EB, gradient = Gf.EB,
      lower = lws, upper = ups, betad = betad, G0 = G0,
      ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
      V = s * Rm), silent = TRUE)
    # opt.EB <- optim(st, fn = Of.EB, gr = Gf.EB, method = 'L-BFGS-B',
    #   lower = lws, upper = ups, betad = betad, G0 = G0,
    #   ws = ws, th.prior = th.prior, V = s * Rm)
    if (class(opt.EB) == 'try-error') {
      st <- c(EP.is[['th0.opt']], EP.is[['th1.opt']])  # Change start
      opt.EB <- try(nlminb(st, objective = Of.EB, gradient = Gf.EB,
        lower = lws, upper = ups, betad = betad, G0 = G0,
        ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
        V = s * Rm))
      if (class(opt.EB) == 'try-error') {  # If EB cannot be evaluated
        stop('search for EB optima encountered numerical errors.')
      }
    }

    # Optimal values
    if (opt.EB[['convergence']] %in% 0:1) {
      th0.hat <- opt.EB[['par']][1]
      th1.hat <- opt.EB[['par']][2]
    } else {
      stop('theta search under "EB" did not converge.')
    }
  } else if (th.search == 'EP') {
    # Optimal values
    th0.hat <- th0.EP
    th1.hat <- th1.EP
  } else {
    stop('theta search method unsupported -- choose amongst "EB" and "EP".')
  }

  # MODEL SEARCH (under fixed theta hat) #######################################
  # Prior probabilities under theta
  pg1cth <- pinc(betad, th0 = th0.hat, th1 = th1.hat,
    rho.min = rho.min, rho.max = rho.max)
  pg1cth <- c(1/2, pg1cth)
  pg1cth[which(pg1cth == 1)] <- 1 - eps
  pg1cth[which(pg1cth == 0)] <- eps

  # Define new asymmetric binomial prior
  nbp <- mombf::modelbinomprior(p = pg1cth)

  # New model search
  new.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = nbp,
    niter = R, verbose = FALSE, center = FALSE, scale = FALSE)

  # Posterior model probabilities
  pprobs <- mombf::postProb(new.fit)[, c(1, 3)]
  mpprobs <- new.fit[['margpp']]

  # Finish
  return(list(pprob.y = pprobs, mpprob.y = mpprobs, th0.hat = th0.hat,
    th1.hat = th1.hat, msfit = new.fit, init.msfit = bvs.fit0))
}; pprob.nm <- compiler::cmpfun(pprob.nm)

################################################################################
# Posterior model probabilities under CIL (MULTIPLE TREATMENT version)
pprob.nm.mt <- function(y, D, X, I = NULL, mod1 = 'ginv', th.search = 'EB',
  th.prior = 'unif', beta.prior = 'nlp', rho.min = NULL, rho.max = NULL,
  th.range = NULL, tau = 0.348, max.mod = Inf, lpen = 'lambda.1se', eps = 1e-10,
  R = 1e5, cross.teff = FALSE, bvs.fit0 = NULL, th.EP = NULL, EP.is = NULL,
  seed = NA) {
################################################################################
  # Renaming
  ms <- mombf::modelSelection
  igp <- mombf::igprior(0.01, 0.01)
  bbp <- mombf::modelbbprior(1, 1)
  nlp <- mombf::momprior(tau = tau)
  zup <- mombf::zellnerprior(tau = nrow(X))
  ufp <- mombf::modelunifprior()
  mlik <- mombf::nlpMarginal
  ginv <- pracma::pinv

  # Pre-computations
  p <- ncol(X)
  Z <- cbind(D, I, X)
  Di <- cbind(D, I)

  # Other fixed parameters
  if (all(is.null(c(rho.min, rho.max)))) {
    rho.min <- 1 / ncol(Z)
    rho.max <- 1 - rho.min
  }

  # Coefficient prior
  if (beta.prior == 'nlp') {
    cp <- nlp
  } else if (beta.prior == 'zup') {
    cp <- zup  
  } else {
    stop('currently supported "beta.prior" arg.: "nlp" (MOM) and "zup" (UIP)')
  }

  # MODULE 1: Estimate coefficients on exposure model
  betad0 <- matrix(NA, nrow = ncol(X), ncol = ncol(Di))
  for (i in 1:ncol(Di)) {
    if (cross.teff == TRUE) {
      A <- Z[, -i]
      b <- Z[, i]
      exc <- 1:(nt - 1)
    } else {
      A <- X
      b <- Di[, i]
      exc <- 2 * p  # phantom column: not exclude anything
      if (length(unique(Di[, i])) <= 3 &  # This controls sum-to-zero constraint
          all(unique(Di[, i]) %in% c(-1, 0, 1))) {
        b[which(b == -1)] <- 0
      }
    }

    # Exposure family (support for Binomial and Gaussian ONLY)
    type.exp <- ifelse(length(unique(b)) == 2, 'binomial', 'gaussian')
    intcpt <- !any(apply(X, 2, var) == 0)

    if (length(table(b)) > 1 & min(table(b)) > 1) {
      if (mod1 == 'ginv') {
        betad0[, i] <- as.matrix(ginv(t(A) %*% A) %*% t(A) %*% b)[-exc]
      } else if (mod1 %in% c('bvs', 'bma', 'bms')) {
        m1.fit <- ms(b, A, priorCoef = cp, priorVar = igp, priorDelta = bbp,
          niter = 1e4, verbose = FALSE, center = FALSE, scale = FALSE)
        post.th <- colMeans(rnlp(msfit = m1.fit, priorCoef = cp, niter = 1e4))
        betad0[, i] <- unname(post.th[2:(length(post.th) - 1)])[-exc]
      } else if (mod1 == 'lasso') {
        m1.fit <- cv.glmnet(A, b, intercept = intcpt)
        betad0[, i] <- unname(coef(m1.fit, s = m1.fit[[lpen]])[-1, 1])[-exc]
      } else if (mod1 == 'lasso_dml') {
        m1.fit <- hdm::rlasso(b ~ A, post = TRUE, intercept = intcpt)
        betad0[, i] <- unname(m1.fit[['coefficients']])[-exc]
      } else if (mod1 == 'lasso_bic') {
        m1.fit <- lasso.bic(b, A, intercept = TRUE, family = type.exp)
        if (intcpt == FALSE) { m1.fit[['coef']][2] <- m1.fit[['coef']][1] }
        betad0[, i] <- unname(m1.fit[['coef']][-1][-exc])
      } else if (mod1 == 'ridge') {
        m1.fit <- glmnet::cv.glmnet(A, b, alpha = 0, intercept = intcpt)
        betad0[, i] <- unname(
          coef(m1.fit, s = m1.fit[['lambda.min']])[-1, ])[-exc]
      } else {
        stop('currently unsupported "mod1" method.')
      }
    } else {
      betad0[, i] <- rep(0, nrow(betad0))  # If treatment has no variability
    }
  }

  # # Scale betad (this is DEPRECATED and should be inactive *always*)
  # if (scale.delta == TRUE) {
  #   #betad <- betad0 / max(abs(betad0))
  #   betad <- apply(betad0, 2, function(x) { x / max(abs(x)) })
  # } else {
  betad <- betad0
  # }
  colnames(betad) <- colnames(Di)
  rownames(betad) <- colnames(X)

  # THETA SEARCH ###############################################################
  # Initial fit: th0 = 0; th1 = 0 (UNIFORM MODEL PRIOR)
  if (is.null(bvs.fit0)) {
    bvs.fit0 <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = ufp,
      niter = R, verbose = FALSE, center = FALSE, scale = FALSE)
  }
  # Posterior model probabilities
  pprobs <- postProb(bvs.fit0)[, c(1, 3)]

  # Marginal inclusion probabilities
  pj1 <- bvs.fit0[['margpp']]
  pj1[which(is.nan(pj1))] <- 0  # Potential problems in mombf with undefined
  pj1[which(pj1 == 1)] <- 1 - eps  # In case there are extreme values
  pj1[which(pj1 == 0)] <- eps

  # Set of visited models
  if (is.na(max.mod)) {  # Prevent a dictionary too big
    max.mod <- which(cumsum(pprobs[, ncol(pprobs)]) >= 0.9)[1]
  }
  G0 <- unname(sapply(
    as.character(pprobs[, 1]), binary.id, p = p, nt = ncol(Di), conc = TRUE))
  G0 <- G0[1:min(max.mod, length(G0))]  # Prevent a dictionary too big

  # Find recommended value for prior theta variance
  if (th.prior == 'gauss') {
    # That leaves Pr(p_j(theta) \in (eps, 1-eps)) = 1-cl (correlation 0)
    s <- set.s(betad, eps = 0.05, cl = 0.1, rho01 = 0, dens = 0.01)
    Rm <- matrix(c(1, 0, 0, 1), ncol = 2)
    if (is.na(s)) {
      # Span of search
      s <- set.s(betad, eps = 0.05, cl = 0.1, rho01 = 0,
        bounds = c(100, 1e3), dens = 0.1)
    }
  } else if (th.prior %in% c('unif', 'tunif')) {
    s <- 1
    Rm <- diag(2)
  } else {
    stop('currently supported "th.prior" arg.: "unif" and "gauss"')
  }

  # Determine the values of hat{theta} #########################################
  # Set limits for optimisation (in principle: unbounded)
  if (th.prior == 'tunif') {
    ups <- rep(log(ncol(Z) - 1), ncol(D) + 1)  # p_j(theta, betad=0) < 1-1/p
    lws <- -ups  # p_j(theta, betad=0) > 1/p
    scale.prpr <- FALSE  # NO LONGER NEEDED
  } else {
    lws <- rep(-Inf, ncol(D) + 1)
    ups <- rep(+Inf, ncol(D) + 1)
  }

  # Range for grid search
  if (is.null(th.range)) {
    if (ncol(D) == 1) {
      th.range <- seq(-40, 40, 1)
    } else if (ncol(D) == 2) {
      th.range <- seq(-40, 40, 2)
    } else if (ncol(D) <= 4) {
      th.range <- seq(-10, 10, 2)
    } else if (ncol(D) <= 7) {
      th.range <- seq(-7.5, 7.5, 2.5)
    } else if (ncol(D) <= 14) {
      th.range <- seq(-2, 2, 2)
      #th.range <- seq(-5, 5, 2.5)
    } else {
      th.range <- c(-2, 2)
    }
  }

  if (is.null(th.EP)) {
    # EP approximation (to set starting point)
    ths <- vector(ncol(D) + 1, mode = 'list')
    for (i in 1:length(ths)) { ths[[i]] <- th.range }
    EP.is <- grid.opt.mt(G0, betad, pj1, th.grid = expand.grid(ths),
      method = 'EP', th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
      V = s * Rm)

    # Load gradient functions and optimise
    st <- unlist(EP.is[[1]])
    opt.EP <- nlminb(st, objective = Of.EP.mt, gradient = Gf.EP.mt, 
      lower = lws, upper = ups, pj1 = pj1, betad = betad, 
      th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, V = s * Rm)
    # opt.EP <- optim(st, fn = Of.EP, gr = Gf.EP, method = 'L-BFGS-B',
    #   lower = lws, upper = ups, 
    #   betad = betad, pj1 = pj1, th.prior = th.prior, V = s * Rm)

    # Set values if there is convergence
    if (opt.EP[['convergence']] %in% 0:1) {
      th.EP <- opt.EP[['par']]
    } else {
      stop('theta search under "EP" did not converge.')
    }
  }

  # Empirical Bayes search (starting at EP solution)
  if (th.search == 'EB') {
    # Model search at EP optimum
    auxprpr <- pinc.mt(betad, th = th.EP, rho.min = rho.min, rho.max = rho.max)
    auxprpr <- c(rep(1/2, ncol(betad)), auxprpr)
    auxprpr[which(auxprpr == 1)] <- 1 - eps
    auxprpr[which(auxprpr == 0)] <- eps
    auxmp <- mombf::modelbinomprior(p = auxprpr)
    update.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = auxmp,
      niter = round(R / 3), verbose = FALSE, center = FALSE, scale = FALSE)

    # Append model set
    upprobs <- postProb(update.fit)[, c(1, 3)]
    G1 <- unname(sapply(
      as.character(upprobs[, 1]), binary.id, p = p, nt = ncol(D), conc = TRUE))
    G0 <- unique(c(G0, G1[1:min(max.mod, length(G1))]))

    # Pre-compute marginal likelihoods
    ws <- unname(unlist(sapply(G0, function(m) {
      js <- which(unlist(strsplit(m, '')) == 1)
      return(mlik(js, y, Z, priorCoef = cp, logscale = FALSE))
    })))

    # Numerical problems with MLs that are too small
    if (all(ws < -650)) {  # We start having problems at exp(-750)
      ct <- abs(max(ws) + 650)  # Constant
      ws <- exp(ws + ct)  # Multiply by constant (NO EFFECT ON ARG MAX)
    } else {
      ws <- exp(ws)
    }

    # Optimise starting at the EP optimum
    opt.EB <- try(nlminb(th.EP, objective = Of.EB.mt, gradient = Gf.EB.mt,
      lower = lws, upper = ups, betad = betad, G0 = G0,
      ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
      V = s * Rm), silent = TRUE)
    # opt.EB <- optim(st, fn = Of.EB, gr = Gf.EB, method = 'L-BFGS-B',
    #   lower = lws, upper = ups, betad = betad, G0 = G0,
    #   ws = ws, th.prior = th.prior, V = s * Rm)
    if (class(opt.EB) == 'try-error') {
      opt.EB <- try(nlminb(EP.is[[1]], objective = Of.EB.mt,
        gradient = Gf.EB.mt, lower = lws, upper = ups, betad = betad,
        G0 = G0, ws = ws, th.prior = th.prior, rho.min = rho.min,
        rho.max = rho.max, V = s * Rm))
      #opt.EB <- opt.EP
    }

    # Optimal values
    if (opt.EB[['convergence']] %in% 0:1) {
      th.hat <- opt.EB[['par']]
    } else {
      stop('theta search under "EB" did not converge.')
    }
  } else if (th.search == 'EP') {
    # Optimal values
    th.hat <- th.EP
  } else {
    stop('theta search method unsupported -- choose amongst "EB" and "EP".')
  }

  # MODEL SEARCH (under fixed theta hat) #######################################
  # Prior probabilities under theta
  pg1cth <- c(rep(1/2, ncol(betad)),
              pinc.mt(betad, th = th.hat, rho.min = rho.min, rho.max = rho.max))
  pg1cth[which(pg1cth == 1)] <- 1 - eps
  pg1cth[which(pg1cth == 0)] <- eps

  # Define new asymmetric binomial prior
  nbp <- mombf::modelbinomprior(p = pg1cth)

  # New model search
  new.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = nbp,
    niter = R, verbose = FALSE, center = FALSE, scale = FALSE)

  # Posterior model probabilities
  pprobs <- mombf::postProb(new.fit)[, c(1, 3)]
  mpprobs <- new.fit[['margpp']]

  # Finish
  return(list(pprob.y = pprobs, mpprob.y = mpprobs, th.hat = th.hat,
              msfit = new.fit, init.msfit = bvs.fit0, EP.is = EP.is,
              mod1.coef = betad))
}; pprob.nm.mt <- compiler::cmpfun(pprob.nm.mt)

################################################################################
# BMA estimate computation provided CIL model post. probs.
bma.teff <- function(msfit, ret.mcmc = TRUE) {
################################################################################
  # Posterior draws
  draws1 <- mombf::rnlp(msfit = msfit, niter = 1e4)

  # BMA Treatment Effect Estimates
  idx <- ifelse(grepl('ntercept', colnames(draws1)[1]), 2, 1)
  teff1 <- mean(draws1[, idx])

  # BMA distribution
  bma.distr1 <- quantile(draws1[, idx], seq(0, 1, 0.005))

  # Output
  out <- list(bma.nlp = teff1, bma.zup = NA,
    bma.d.nlp = bma.distr1, bma.d.zup = NA,
    mcmc.nlp = NA, mcmc.zup = NA)
  if (ret.mcmc == TRUE) {
    out[['mcmc.nlp']] <- draws1[, idx]
    out[['mcmc.zup']] <- NA
  }

  # End
  return(invisible(out))
}; bma.teff <- compiler::cmpfun(bma.teff)

################################################################################
# BMA est. computation given CIL model post. probs. (MULTIPLE TREATMENT version)
bma.teff.mt <- function(msfit, nt = 1, ret.mcmc = TRUE) {
################################################################################
  # Posterior draws
  draws1 <- mombf::rnlp(msfit = msfit, niter = 1e4)

  # BMA Treatment Effect Estimates
  teff1 <- colMeans(draws1[, 1 + 1:nt])

  # BMA distribution
  bma.distr1 <- apply(draws1[, 1 + 1:nt], 2, quantile, seq(0, 1, 0.005))

  # Output
  out <- list(bma.nlp = teff1, bma.zup = NA,
    bma.d.nlp = bma.distr1, bma.d.zup = NA,
    mcmc.nlp = NA, mcmc.zup = NA)
  if (ret.mcmc == TRUE) {
    out[['mcmc.nlp']] <- draws1[, 1 + 1:nt]
    out[['mcmc.zup']] <- NA
  }

  # End
  return(invisible(out))
}; bma.teff.mt <- compiler::cmpfun(bma.teff.mt)

################################################################################
# Wrapper function to compute post. probs. and BMA estimate
new.method <- function(y, d, X, priors = NULL, R = 1e5, ...) {
################################################################################
  # Posterior model probabilities
  pprobs <- pprob.nm(y, d, X, R = R, ...)

  # BMA computation
  teff <- bma.teff(pprobs[['msfit']])

  # Output object
  out <- list(teff.nlp = teff[['bma.nlp']], teff.zup = teff[['bma.zup']],
              pprob = pprobs[['pprob.y']], mpprob = pprobs[['mpprob.y']],
              th0.hat = pprobs[['th0.hat']], th1.hat = pprobs[['th1.hat']],
              teff.d.nlp = teff[['bma.d.nlp']],
              teff.d.zup = teff[['bma.d.zup']],
              init.msfit = pprobs[['init.msfit']],
              msfit = pprobs[['msfit']])

  # Exit
  return(invisible(out))
}

################################################################################
# Wrapper to compute post. probs. and BMA estimate (MULTIPLE TREATMENT version)
new.method.mt <- function(y, D, X, I = NULL, priors = NULL, R = 1e5, ...) {
################################################################################
  # Posterior model probabilities
  pprobs <- pprob.nm.mt(y, D, X, I, R = R, ...)

  # BMA computation
  teff <- bma.teff.mt(pprobs[['msfit']], nt = ncol(cbind(D, I)))

  # Output object
  out <- list(teff.nlp = teff[['bma.nlp']], teff.zup = teff[['bma.zup']],
              pprob = pprobs[['pprob.y']], mpprob = pprobs[['mpprob.y']],
              th.hat = pprobs[['th.hat']], EP.is = pprobs[['EP.is']],
              mod1.coef = pprobs[['mod1.coef']],
              teff.d.nlp = teff[['bma.d.nlp']],
              teff.d.zup = teff[['bma.d.zup']],
              init.msfit = pprobs[['init.msfit']],
              msfit = pprobs[['msfit']])

  # Exit
  return(invisible(out))
}
# END OF SCRIPT
