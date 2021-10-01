################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '03_generic_illustrations.R', sep = ''))
################################################################################

################################################################################
# FIGURE 3: Representation of CIL prior for values of theta
################################################################################
# Parameters
J <- 99
f <- function(f1, rho, th) {
  rho + (1 - 2 * rho) * (1 + exp(-th[1] - th[2] * abs(f1)))^(-1)
}
rho <- (J^2 + 1)^(-1)
th1 <- c(-1, -2)
th2 <- c(-1, 2)
th3 <- c(-1, 4)
f1s <- seq(0, 3, 0.01)

# Plot
ffig <- paste(FIGDIR, 'cil_shape.pdf', sep = '')
pdf(ffig, height = 4, width = 6)
xl <- expression(f['j,1'])
yl <- expression(pi[j](theta))
par(mfrow = c(1, 1), mar = c(4, 4.5, 0.5, 0.5), bg = 'white')
plot(f1s, f(f1s, rho, th1), type = 'n', ylim = c(0, 1), xlab = xl, ylab = yl)
grid(); abline(v = 0, col = 'darkgray');
abline(h = 0, col = 'darkgray'); abline(h = 1, col = 'darkgray')
abline(h = rho,lty = 3);
abline(h = 1 - rho,lty = 3); abline(h = f(f1s, rho, th1)[1],lty = 3)
lines(f1s, f(f1s, rho, th1), lty = 6)
points(f1s[1], f(f1s, rho, th1)[1], pch = 16)
lines(f1s, f(f1s, rho, th2), lty = 1)
lines(f1s, f(f1s, rho, th3), lty = 5)
legend('bottomright', c(expression(theta[1] == +2), expression(theta[1] == +4),
  expression(theta[1] == -2)), lty = c(1, 5:6), bg = 'white')
dev.off(); cat('Printed figure:', ffig, '\n'); rm(ffig)

################################################################################
# FIGURE S1: Non-local prior illustration
################################################################################
xs <- seq(-3,3,0.01)
xl <- expression(alpha[t])
yl <- 'density'
ffig <- paste(FIGDIR, 'nlpmom.pdf', sep = '')
pdf(ffig, height = 3, width = 5)
par(mfrow = c(1, 1), mar = c(4, 4, 0.5, 0.5), bg = 'white')
plot(xs, mombf::dmom(x = xs, tau = 0.348), type = 'n',
  ylab = yl, xlab = xl, ylim = c(0, 0.53))
abline(h = 0, col = 'gray')
abline(v = 0, col = 'gray', lty = 2)
lines(xs, mombf::dmom(x = xs, tau = 0.348))
dev.off(); cat('Printed figure:', ffig, '\n'); rm(ffig)

################################################################################
# FIGURE S3: Representation of differences between EP and EB
################################################################################
# Design paramaters
r <- 25
n <- 100; nc <- 3; p <- 49; S <- diag(p)
bY <- c(1, 1, 1, 1, 1, 1, rep(0, p - 6))
bd <- c(rep(0, 6 - nc), rep(1, 6), rep(0, p - 12 + nc))
a <- 1; phiy <- phid <- 1

# Simulation environment
R <- 1e4
burnin <- round(0.1 * R)
eta <- 0.05
tau <- 0.348
mod1 <- 'lasso_bic'

# Dataset generation
set.seed(100 * r)
X <- mvtnorm::rmvnorm(n = n, mean = rep(0, nrow(S)), sigma = S)
epsd <- rnorm(n, mean = 0, sd = sqrt(phid))
epsy <- rnorm(n, mean = 0, sd = sqrt(phiy))
d <- X %*% bd + epsd
y <- a * d + X %*% bY + epsy
Z <- cbind(d, X)

# Shortcuts for mombf
ms <- mombf::modelSelection
mlik <- mombf::nlpMarginal
igp <- mombf::igprior(0.01, 0.01)
bbp <- mombf::modelbbprior(1, 1)
nlp <- mombf::momprior(tau = tau)
zup <- mombf::zellnerprior(tau = nrow(X))
priors <- list(nlp, zup, igp)

# General parameters of the function
th.search <- 'EB'
th.prior <- 'unif'
beta.prior <- 'nlp'
rho.min <- NULL
rho.max <- NULL
max.mod <- Inf
lpen <- 'lambda.1se'
eps <- 1e-10
bvs.fit0 <- NULL
seed <- NA
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
  rho.min <- 1 / (ncol(Z)^2 + 1)
  rho.max <- 1 - rho.min
}
cp <- nlp

# Exposure family (support for Binomial and Gaussian ONLY)
type.exp <- ifelse(length(unique(d)) == 2, 'binomial', 'gaussian')
lb.fit <- lasso.bic(d, X, intercept = FALSE, family = type.exp)
betad <- betad0 <- unname(lb.fit[['coef']][-1])
bvs.fit0 <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = ufp,
  niter = R, verbose = TRUE, center = FALSE, scale = FALSE)
pprobs <- postProb(bvs.fit0)[, c(1, 3)]

# Marginal inclusion probabilities
pj1 <- bvs.fit0[['margpp']]
pj1[which(pj1 == 1)] <- 1 - eps  # In case there are too extreme values
pj1[which(pj1 == 0)] <- eps
G0 <- unname(sapply(
  as.character(pprobs[, 1]), binary.id, p = p, conc = TRUE))
G0 <- G0[1:min(max.mod, length(G0))]

# Find recommended value for prior theta variance
s <- 1
Rm <- diag(2)
lws <- c(-Inf, -Inf)
ups <- c(+Inf, +Inf)
ths <- seq(-7, 7, 1)

# Here onwards we just run a code chunk extracted from the CIL function
ths1 <- seq(-6, 1, 1)
ths2 <- seq(-1, 8, 1)
EP.is <- grid.opt(G0, betad, pj1, th.grid = expand.grid(ths1, ths2),
  method = 'EP', th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
  V = s * Rm, ret.plotly = TRUE)

ffig <- paste(FIGDIR, 'thEP.pdf', sep = '')
plotly::orca(EP.is[['th.fig']], 'thEP.pdf')
file.copy('thEP.pdf', ffig); unlink('thEP.pdf')
cat('Printed figure:', ffig, '\n'); rm(ffig)

st <- c(EP.is[['th0.opt']], EP.is[['th1.opt']])
opt.EP <- nlminb(st, objective = Of.EP, gradient = Gf.EP, 
  lower = lws, upper = ups, pj1 = pj1, betad = betad,
  th.prior = th.prior, rho.min = rho.min, rho.max = rho.max, V = s * Rm)

# Set values if there is convergence
th0.EP <- opt.EP[['par']][1]  # -2.340109
th1.EP <- opt.EP[['par']][2]  # 3.07873

auxprpr <- pinc(betad, th0 = th0.EP, th1 = th1.EP,
    rho.min = rho.min, rho.max = rho.max)
auxprpr <- c(1/2, auxprpr)
auxprpr[which(auxprpr == 1)] <- 1 - eps
auxprpr[which(auxprpr == 0)] <- eps
auxmp <- mombf::modelbinomprior(p = auxprpr)
update.fit <- ms(y, Z, priorCoef = cp, priorVar = igp, priorDelta = auxmp,
  niter = R, verbose = TRUE, center = FALSE, scale = FALSE)
upprobs <- postProb(update.fit)[, c(1, 3)]
G1 <- unname(sapply(
  as.character(upprobs[, 1]), binary.id, p = p, conc = TRUE))
G0 <- unique(c(G0, G1[1:min(max.mod, length(G1))]))

# Pre-compute marginal likelihoods
ws <- unname(unlist(sapply(G0, function(m) {
  js <- which(unlist(strsplit(m, '')) == 1)
  return(mlik(js, y, Z, priorCoef = cp, logscale = TRUE))
})))
ws <- exp(ws)

ffig <- paste(FIGDIR, 'thEB.pdf', sep = '')
EB.is <- grid.opt(G0, betad, pj1, ws = ws, th.grid = expand.grid(ths1, ths2), 
  method = 'EB', th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
  V = s * Rm, ret.plotly = TRUE)
plotly::orca(EB.is[['th.fig']], 'thEB.pdf')
file.copy('thEB.pdf', ffig); unlink('thEB.pdf')
cat('Printed figure:', ffig, '\n'); rm(ffig)

st <- c(th0.EP, th1.EP)
opt.EB <- try(nlminb(st, objective = Of.EB, gradient = Gf.EB,
  lower = lws, upper = ups, betad = betad, G0 = G0,
  ws = ws, th.prior = th.prior, rho.min = rho.min, rho.max = rho.max,
  V = s * Rm), silent = TRUE)
th0.hat <- opt.EB[['par']][1]  # -2.426168
th1.hat <- opt.EB[['par']][2]  # 3.193549
# END OF SCRIPT
