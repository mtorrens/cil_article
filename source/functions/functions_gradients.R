################################################################################
# source('~/Desktop/year3/bma_teff/v15/syntax/00_start.R')
# source(paste(FCNDIR, 'functions_gradients.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * Of.EB
# * Gf.EB
# * Of.EP
# * Gf.EP
# * Of.EB.mt
# * Gf.EB.mt
# * Of.EP.mt
# * Gf.EP.mt
################################################################################

################################################################################
# Objective function of Empirical Bayes original formulation
Of.EB <- function(x, betad, G0, ws, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # log p(y | theta)
  p1 <- pinc(betad, th0 = x[1], th1 = x[2], rho.min = rho.min, rho.max = rho.max)
  p1 <- c(1/2, p1)
  sc <- sum(sapply(G0, function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    pg <- prod((p1 ** gm) * ((1 - p1) ** (1 - gm)))
    wg <- ws[which(G0 == m)]
    return(wg * pg)
  }))
  pt1 <- log(sum(sc))  # log p(y | theta)

  # log p(theta)
  if (th.prior %in% c('unif', 'tunif')) {
    pt2 <- 0
  } else if (th.prior == 'gauss') {
    x <- matrix(c(x[1], x[2]), ncol = 1)
    pt2 <- as.numeric(-(1 / 2) * t(x) %*% solve(V) %*% x)
    #pt2 <- - (x[1] ** 2 + x[2] ** 2) / (2 * s)
  }

  return(-(pt1 + pt2))
}; Of.EB <- compiler::cmpfun(Of.EB)

################################################################################
# Gradient function of Empirical Bayes original formulation
Gf.EB <- function(x, betad, G0, ws, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # Compute (unnormalised) model inclusion probabilities
  p1 <- pinc(betad, th0 = x[1], th1 = x[2], rho.min = rho.min, rho.max = rho.max)
  p1 <- c(1/2, p1)
  sc <- sapply(G0, function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    wg <- ws[which(G0 == m)]
    pg <- (p1 ** gm) * (1 - p1) ** (1 - gm)
    return(wg * prod(pg))
  })

  # Marginal inclusion probabilities
  pprobs <- data.frame('model_id' = G0, 'pprob' = sc / sum(sc))
  pprobs[, 1] <- sapply(pprobs[, 1], compact.id)
  mpprobs <- apply(t(sapply(pprobs[, 1], binary.id, p = nchar(G0[1]) - 1)), 2,
    function(x) { sum(x * pprobs[, ncol(pprobs)]) })

  # grad log p(y | theta)
  # g0 <- sum(mpprobs - p1)
  # g1 <- sum((mpprobs - p1) * c(0, abs(betad)))
  g0 <- sum((mpprobs - p1)[-1])
  g1 <- sum((mpprobs - p1)[-1] * abs(betad))

  # grad log p(theta)
  if (th.prior %in% c('unif', 'tunif')) {
    t0 <- t1 <- 0
  } else if (th.prior == 'gauss') {
    x <- matrix(c(x[1], x[2]), ncol = 1)
    aux <- -t(x) %*% solve(V)
    t0 <- aux[1]
    t1 <- aux[2]
    # t0 <- -x[1] / s
    # t1 <- -x[2] / s
  }

  # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  return(c(-(g0 + t0), -(g1 + t1)))
}; Gf.EB <- compiler::cmpfun(Gf.EB)

################################################################################
# Objective function of Expectation Propagation formulation
Of.EP <- function(x, betad, pj1, th.prior = 'unif', rho.min = 0, rho.max = 1,
  V = diag(2)) {
################################################################################
  # log p(y | theta)
  p1 <- pinc(betad, th0 = x[1], th1 = x[2], rho.min = rho.min, rho.max = rho.max)
  p1 <- c(1/2, p1)
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)
  pt1 <- sum(log(fj))

  # log p(theta)
  if (th.prior %in% c('unif', 'tunif')) {
    pt2 <- 0
  } else if (th.prior == 'gauss') {
    x <- matrix(c(x[1], x[2]), ncol = 1)
    pt2 <- as.numeric(-(1 / 2) * t(x) %*% solve(V) %*% x)
    #pt2 <- - (x[1] ** 2 + x[2] ** 2) / (2 * s)
  }

  # - log p(theta | y) \propto -(log p(y | theta) + log p(theta))
  return(-(pt1 + pt2))
}; Of.EP <- compiler::cmpfun(Of.EP)

################################################################################
# Gradient function of Expectation Propagation formulation
Gf.EP <- function(x, betad, pj1, th.prior = 'unif', rho.min = 0, rho.max = 1,
  V = diag(2)) {
################################################################################
  # grad log p(y | theta)
  p1 <- pinc(betad, th0 = x[1], th1 = x[2], rho.min = rho.min, rho.max = rho.max)
  p1 <- c(1/2, p1)
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)
  g0 <- (((1 - p1) * p1) / fj) * (2 * pj1 - 1)
  g1 <- c(0, abs(betad)) * g0

  # grad log p(theta)
  if (th.prior %in% c('unif', 'tunif')) {
    t0 <- t1 <- 0
  } else if (th.prior == 'gauss') {
    x <- matrix(c(x[1], x[2]), ncol = 1)
    aux <- -t(x) %*% solve(V)
    t0 <- aux[1]
    t1 <- aux[2]
    # t0 <- -x[1] / s
    # t1 <- -x[2] / s
  }

  # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  return(c(-(sum(g0[-1]) + t0), -(sum(g1[-1]) + t1)))
}; Gf.EP <- compiler::cmpfun(Gf.EP)

################################################################################
# Objective function of EB original formulation (MULTIPLE TREATMENT version)
Of.EB.mt <- function(x, betad, G0, ws, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  p1 <- c(rep(1/2, nchar(G0[1]) - length(p1)), p1)
  sc <- sum(sapply(G0, function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    pg <- prod((p1 ** gm) * (1 - p1) ** (1 - gm))
    wg <- ws[which(G0 == m)]
    return(wg * pg)
  }))
  pt1 <- log(sum(sc))  # log p(y | theta)

  # log p(theta)
  # if (th.prior %in% c('unif', 'tunif')) {
    pt2 <- 0
  # } else if (th.prior == 'gauss') {
  #   x <- matrix(c(x[1], x[2]), ncol = 1)
  #   pt2 <- as.numeric(-(1 / 2) * t(x) %*% solve(V) %*% x)
  #   #pt2 <- - (x[1] ** 2 + x[2] ** 2) / (2 * s)
  # }

  return(-(pt1 + pt2))
}; Of.EB.mt <- compiler::cmpfun(Of.EB.mt)

################################################################################
# Gradient function of EB original formulation (MULTIPLE TREATMENT version)
Gf.EB.mt <- function(x, betad, G0, ws, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # Compute (unnormalised) model inclusion probabilities
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  sD <- nchar(G0[1]) - length(p1)
  p1 <- c(rep(1/2, sD), p1)
  sc <- sapply(G0, function(m) {
    gm <- as.numeric(unlist(strsplit(m, '')))
    wg <- ws[which(G0 == m)]
    pg <- (p1 ** gm) * (1 - p1) ** (1 - gm)
    return(wg * prod(pg))
  })

  # Marginal inclusion probabilities
  pprobs <- data.frame('model_id' = G0, 'pprob' = sc / sum(sc))
  pprobs[, 1] <- sapply(pprobs[, 1], compact.id)
  mpprobs <- apply(t(sapply(pprobs[, 1], binary.id,
    p = nchar(G0[1]) - ncol(betad), nt = ncol(betad))), 2,
    function(x) { sum(x * pprobs[, ncol(pprobs)]) })

  # grad log p(y | theta)
  gs <- rep(NA, length(x))
  gs[1] <- sum((mpprobs - p1)[-(1:sD)])
  for (k in 2:length(gs)) {
    aux0 <- (ncol(betad) / (length(x) - 1) - 1)
    idxs <- c(k - 1)
    if (aux0 > 0) { idxs <- c(idxs, (length(x) - 1) + aux0 * (k - 2) + 1:aux0) }
    aux1 <- mpprobs - p1
    aux2 <- rowSums(abs(as.matrix(betad[, idxs])))
    gs[k] <- sum(aux1[-(1:sD)] * aux2)
  }

  # # grad log p(theta)
  # if (th.prior %in% c('unif', 'tunif')) {
    t0 <- t1 <- 0
  # } else if (th.prior == 'gauss') {
  #   x <- matrix(c(x[1], x[2]), ncol = 1)
  #   aux <- -t(x) %*% solve(V)
  #   t0 <- aux[1]
  #   t1 <- aux[2]
  #   # t0 <- -x[1] / s
  #   # t1 <- -x[2] / s
  # }

  # # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  # return(c(-(g0 + t0), -(g1 + t1)))
  return(-(gs + t0))
}; Gf.EB.mt <- compiler::cmpfun(Gf.EB.mt)

################################################################################
# Objective function of Exp. Prop. formulation (MULTIPLE TREATMENT version)
Of.EP.mt <- function(x, betad, pj1, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  p1 <- c(rep(1/2, length(pj1) - length(p1)), p1)
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)
  pt1 <- sum(log(fj))

  # log p(theta)
  if (th.prior %in% c('unif', 'tunif')) {
    pt2 <- 0
  }
  # } else if (th.prior == 'gauss') {
  #   x <- matrix(c(x[1], x[2]), ncol = 1)
  #   pt2 <- as.numeric(-(1 / 2) * t(x) %*% solve(V) %*% x)
  #   #pt2 <- - (x[1] ** 2 + x[2] ** 2) / (2 * s)
  # }

  # - log p(theta | y) \propto -(log p(y | theta) + log p(theta))
  return(-(pt1 + pt2))
}; Of.EP.mt <- compiler::cmpfun(Of.EP.mt)

################################################################################
# Gradient function of Exp. Prop. formulation (MULTIPLE TREATMENT version)
Gf.EP.mt <- function(x, betad, nt, pj1, th.prior = 'unif', rho.min = 0,
  rho.max = 1, V = diag(2)) {
################################################################################
  # grad log p(y | theta)
  p1 <- pinc.mt(betad, th = x, rho.min = rho.min, rho.max = rho.max)
  sD <- length(pj1) - length(p1)
  p1 <- c(rep(1/2, sD), p1)
  fj <- pj1 * p1 + (1 - pj1) * (1 - p1)

  gs <- rep(NA, length(x))
  gs[1] <- sum(((((1 - p1) * p1) / fj) * (2 * pj1 - 1))[-(1:sD)])
  for (k in 2:length(gs)) {
    aux0 <- (ncol(betad) / (length(x) - 1) - 1)
    idxs <- c(k - 1)
    if (aux0 > 0) { idxs <- c(idxs, (length(x) - 1) + aux0 * (k - 2) + 1:aux0) }
    aux1 <- (((1 - p1) * p1) / fj) * (2 * pj1 - 1)
    aux2 <- rowSums(abs(as.matrix(betad[, idxs])))
    gs[k] <- sum(aux1[-(1:sD)] * aux2)
  }

  t0 <- 0
  # # grad log p(theta)
  # if (th.prior %in% c('unif', 'tunif')) {
  #   t0 <- t1 <- 0
  # }
  # # } else if (th.prior == 'gauss') {
  # #   x <- matrix(c(x[1], x[2]), ncol = 1)
  # #   aux <- -t(x) %*% solve(V)
  # #   t0 <- aux[1]
  # #   t1 <- aux[2]
  # #   # t0 <- -x[1] / s
  # #   # t1 <- -x[2] / s
  # # }

  # # - grad log p(theta | y) = -(grad log p(y | theta) + grad log p(theta))
  # return(c(-(sum(g0[-1]) + t0), -(sum(g1[-1]) + t1)))
  return(-(gs + t0))
}; Gf.EP.mt <- compiler::cmpfun(Gf.EP.mt)
# END OF SCRIPT
