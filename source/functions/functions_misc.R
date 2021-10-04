################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(FCNDIR, 'functions_butler.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * inside
# * unlist.matrix
# * logit
# * lasso.bic
################################################################################

################################################################################
# Test if value "x" is inside interval ab = (a, b)
inside <- function(x, ab, include.ab = TRUE) {
################################################################################
  ifelse(include.ab == TRUE,
    (x >= ab[1] && x <= ab[2]), (x > ab[1] && x < ab[2]))
}

################################################################################
# Unlist matrix whose columns can be elements of a list
unlist.matrix <- function(X, kill.names = TRUE) {
################################################################################
  for (j in 1:ncol(X)) { X[, j] <- unlist(X[, j]) }
  if (kill.names == TRUE) {
    rownames(X) <- NULL
    #colnames(X) <- NULL
  }
  return(X)
}

################################################################################
# Logit function
logit <- function(p) {
################################################################################
  -log(1 / p - 1)
}

################################################################################
# LASSO estimation with BIC criterion for tunning parameter
lasso.bic <- function(y, x, intercept = TRUE, standardize = TRUE,
  family = 'gaussian') {
# (c) David Rossell
################################################################################
  require(glmnet)
  fit <- glmnet(x = x, y = y, family = family, alpha = 1,
    intercept = intercept, standardize = standardize)
  if (intercept == TRUE) {
    pred <- cbind(1, x) %*% rbind(fit[['a0']], fit[['beta']])
  } else {
    pred <- x %*% fit[['beta']]
  }
  
  transf <- function(x, family) {  
    if (family == 'binomial') { x <- exp(pred) / (1 + exp(pred)) }
    return(x)
  }
  pred <- transf(pred, family = family)

  n <- length(y)
  p <- colSums(fit[['beta']] != 0) + ifelse(intercept == TRUE, 1, 0)
  bic <- n * log(colSums((y - pred)^2) / length(y)) +
         n * (log(2 * pi) + 1) + log(n) * p
  sel <- which.min(bic)
  beta <- c(fit[['a0']][sel], fit[['beta']][, sel])
  ypred <- pred[, sel]
  ans <- list(coef = beta, ypred = ypred, lambda.opt = fit[['lambda']][sel],
    lambda = data.frame(lambda = fit[['lambda']], bic = bic, nvars = p))
  return(ans)
}
# END OF SCRIPT
