################################################################################
# source('~/Desktop/stats/cil_article/source/00_start.R')
# source('/no_backup/jferrer/mtorrens/stats/cil_article/source/00_start.R')
# source(paste(SRCDIR, '07_abortion_analysis.R', sep = ''))
################################################################################

################################################################################
# Belloni's abortion data
# LevittExampleREStudTable.m:
#   created the lagged responses, treatments and controls
################################################################################

################################################################################
# Setup
################################################################################
# Parameterts
plot.vars <- FALSE  # TRUE if an exploratory plot of the vars is desired

# Global seed
set.seed(1)  # Reproducibility

################################################################################
# 1. VIOLENT CRIME
cat('* VIOLENT CRIME: ******************************************************\n')
################################################################################
# Read data
file1 <- paste(INPDIR, 'abortion/viol.txt', sep = '')
file2 <- paste(INPDIR, 'abortion/vnames.csv', sep = '')
data <- read.table(file1, sep = '\t', header = FALSE)
vnames <- read.csv(file2, sep = ',', header = FALSE)
names(data) <- gsub(' ', '', vnames[, 3])

# Determine inputs and outputs
y <- data[, 'violent']
x <- as.matrix(data[, -1])
x <- x[, colnames(x) != 'year01']  # Remove reference category
x <- x[, colnames(x) != 'xV0^2']  # Remove fully collinear covariate
#x <- x[, colnames(x) != 'xM0^2']  # Remove fully collinear covariate

# Visual
if (plot.vars == TRUE) { plot.xy(xvar = x[, 'abortion'], yvar = y) }
cor(x[, 'abortion'], y)
# [1] -0.2757895
summary(cor(x[, 1], x[, -1])[1, ])  # Are there strong correlations?
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.33820 -0.05313  0.01662  0.06459  0.13540  0.57874 

# Build modelling matrices
df <- data.frame(y, x)
sel <- grep('year', colnames(df))
for (j in 1:length(sel)) { df[, sel[j]] <- factor(df[, sel[j]]) }
D <- matrix(df[, 'abortion'], ncol = 1)
X <- df[, c(-1:-2), ]
f <- formula(paste('y ~ ', paste(names(df)[-1], collapse = ' + ')))

# Variables that will be forced into the model (year lags)
includevars <- c(TRUE, FALSE,
  rep(TRUE, length(sel)),
  rep(FALSE, ncol(df) - length(sel) - 2))

# OLS ##########################################################################
ls <- lm(f, data = df)
b.ls <- summary(ls)[['coefficients']]
b.ls <- cbind(b.ls[, 1], confint(ls), b.ls[, 4])
round(head(b.ls), 3)
#                     2.5 % 97.5 %      
# (Intercept)  0.059  0.032  0.086 0.000
# abortion    -0.035 -1.374  1.303 0.959
# year021     -0.080 -0.117 -0.044 0.000
# year031     -0.019 -0.068  0.030 0.439
# year041     -0.020 -0.100  0.060 0.620
# year051      0.034 -0.076  0.143 0.544

# BMA (Normal prior) ###########################################################
ms <- modelSelection(f, data = df, includevars = includevars, niter = 1e4,
  priorCoef = normalidprior(), priorDelta = modelbbprior())
b.bma <- coef(ms)
round(head(b.bma), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.051 -0.002  0.103  1.000
# abortion      -0.066 -0.288  0.000  0.212
# year021       -0.070 -0.087 -0.052  1.000
# year031       -0.009 -0.026  0.008  1.000
# year041       -0.009 -0.027  0.008  1.000
# year051        0.044  0.025  0.063  1.000

# BMA (non-local MOM prior)
ms2 <- modelSelection(f, data = df, includevars = includevars, niter = 1e4,
  priorCoef = momprior(), priorDelta = modelbbprior())
b.bma2 <- coef(ms2)
round(head(b.bma2), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.058  0.036  0.073  1.000
# abortion      -0.005 -0.141  0.000  0.034
# year021       -0.078 -0.101 -0.052  1.000
# year031       -0.025 -0.047  0.014  1.000
# year041       -0.026 -0.048  0.011  1.000
# year051        0.036  0.016  0.061  1.000

# CIL (EP + Normal) ############################################################
includevars <- includevars[1:(ncol(X) + 2)]
cilfit <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EP', priorCoef = normalidprior(), includevars = includevars)
b <- coef(cilfit[['msfit']])
round(head(b), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.054 -0.001  0.109   1.00
# D             -0.168 -0.321  0.000   0.98
# year021       -0.070 -0.087 -0.052   1.00
# year031       -0.008 -0.025  0.009   1.00
# year041       -0.006 -0.022  0.011   1.00
# year051        0.050  0.033  0.067   1.00

# CIL (EP + non-local MOM)
cilfit2 <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EP', priorCoef = momprior(), includevars = includevars)
b2 <- coef(cilfit2[['msfit']])
round(head(b2), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.068  0.052  0.084  1.000
# D             -0.110 -0.236  0.000  0.658
# year021       -0.084 -0.107 -0.061  1.000
# year031       -0.030 -0.051 -0.011  1.000
# year041       -0.029 -0.049 -0.011  1.000
# year051        0.036  0.017  0.058  1.000

# CIL (EB + Normal)
cilfit.eb <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EB', priorCoef = normalidprior(), includevars = includevars)
b.eb <- coef(cilfit.eb[['msfit']])
round(head(b.eb), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.054  0.000  0.110  1.000
# D             -0.168 -0.320  0.000  0.979
# year021       -0.070 -0.087 -0.052  1.000
# year031       -0.008 -0.025  0.009  1.000
# year041       -0.006 -0.022  0.011  1.000
# year051        0.050  0.033  0.068  1.000

# CIL (EB + non-local MOM)
cilfit2.eb <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EB', priorCoef = momprior(), includevars = includevars)
b2.eb <- coef(cilfit2.eb[['msfit']])
round(head(b2.eb), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.068  0.050  0.085  1.000
# D             -0.112 -0.237  0.000  0.658
# year021       -0.084 -0.107 -0.062  1.000
# year031       -0.031 -0.051 -0.011  1.000
# year041       -0.030 -0.050 -0.011  1.000
# year051        0.035  0.016  0.057  1.000

# Save results #################################################################
save(ms, ms2, cilfit, cilfit2, cilfit.eb, cilfit2.eb,
  b.ls, b.bma, b.bma2, b, b2, b.eb, b2.eb,
  file = paste(DATDIR, 'abortion/violent_output.RData', sep = ''), compress = TRUE)

# Post-processing analysis
xtable(b.bma[b.bma[, 4] > 0.1, ])
xtable(b.bma2[b.bma2[, 4] > 0.1, ])
xtable(b.eb[b.eb[, 4] > 0.1, ])
xtable(b2.eb[b2.eb[, 4] > 0.1, ])

head(postProb(ms))
head(postProb(ms2))
head(postProb(cilfit.eb))
head(postProb(cilfit2.eb))

# DML analysis #################################################################
dml <- rlassoEffect(x = x[, -1], y = y, d = D, method = 'double selection',
  I3 = includevars[-1:-2])

round(summary(dml)[[1]], 3)  # Inference for treatment effect
#    Estimate. Std. Error t value Pr(>|t|)
# d1    -0.209      0.129   -1.62    0.105
round(confint(dml), 3)
#     2.5 % 97.5 %
# d1 -0.462  0.044

# Results
sel <- names(dml[['coefficients.reg']])[-1:-2]
sel <- sub('x', '', sel)
lmfit <- lm(y ~ D + x[, sel])  # Final DML outcome model
tab <- summary(lmfit)[['coefficients']][, -3]
xtable(tab, digits = c(0, 2, 2, 3))

# BAC ##########################################################################
W <- cbind.data.frame(y, D, x[, -1])
colnames(W) <- c('y', 'd', paste('V', 1:ncol(X), sep = ''))
bac.fitI <- bac(data = W, exposure = 'd', outcome = 'y',
  confounders = paste('V', 1:ncol(X), sep = ''), interactors = NULL,
  familyX = 'gaussian', familyY = 'gaussian', num_its = 5e3, burnM = 5e2,
  burnB = 5e2, thin = 1)#, omega=Inf)
s <- summary(bac.fitI)
ci.bac <- c(s[[1]], s[[2]])
round(ci.bac, 3)
#         2.5%  97.5% 
# 0.700 -0.167  1.576 

# Results
margpp <- colMeans(bac.fitI[['models']][['models']])
#margpp <- colMeans(bac.fitI$models$models)
sum(margpp > 0.5)
# [1] 84

# ACPME ########################################################################
# Fit
x2 <- x; colnames(x2) <- paste('feat', 1:ncol(x2), sep = '')  # don't like names
acm.fit <- regimes::acpme(y = y, Z = D, C = x2[, -1], niter = 1e4)

# Result summary
colMeans(acm.fit[['beta']])  # Treatment effect estimate
# [1] -0.4230742
t(apply(acm.fit[['beta']], 2, quantile, probs = c(0.025, 0.975)))
#            2.5%      97.5%
# [1,] -0.5427336 -0.3085127
mean(rowSums(acm.fit[['alpha']]))
# [1] 51.8878

# Save
fout <- paste(DATDIR, 'abortion/violent_output_dml:bac:acpme.RData', sep = '')
save(dml, bac.fitI, acm.fit, file = fout, compress = TRUE)

################################################################################
# 2. PROPERTY CRIME
cat('* PROPERTY CRIME: *****************************************************\n')
################################################################################
# Read data
file1 <- paste(INPDIR, 'abortion/prop.txt', sep = '')
file2 <- paste(INPDIR, 'abortion/vnames.csv', sep = '')
data <- read.table(file1, sep = '\t', header = FALSE)
vnames <- read.csv(file2, sep = ',', header = FALSE)
names(data) <- gsub(' ', '', vnames[, 2])

# Determine inputs and outputs
y <- data[, 'property']
x <- as.matrix(data[, -1])
x <- x[, colnames(x) != 'year01']  # Remove reference category
x <- x[, colnames(x) != 'xV0^2']  # Remove fully collinear covariate
#x <- x[, colnames(x) != 'xM0^2']  # Remove fully collinear covariate

# Visual
if (plot.vars == TRUE) { plot.xy(xvar = x[, 'abortion'], yvar = y) }
cor(x[, 'abortion'], y)
# [1] -0.2739313
summary(cor(x[, 1], x[, -1])[1, ])  # Are there strong correlations?
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.375831 -0.040974  0.008925  0.066857  0.149867  0.634557 
# Variables 285–296 are very highly correlated
# This actually deflates somewhat post. probs. of CIL

# Build modelling matrices
df <- data.frame(y, x)
sel <- grep('year', colnames(df))
for (j in 1:length(sel)) { df[, sel[j]] <- factor(df[, sel[j]]) }
D <- matrix(df[, 'abortion'], ncol = 1)
X <- df[, c(-1:-2), ]
f <- formula(paste('y ~ ', paste(names(df)[-1], collapse = ' + ')))

# Variables that will be forced into the model (year lags)
includevars <- c(TRUE, FALSE,
  rep(TRUE, length(sel)),
  rep(FALSE, ncol(df) - length(sel) - 2))

# OLS ##########################################################################
ls <- lm(f, data = df)
b.ls <- summary(ls)[['coefficients']]
#b.ls <- cbind.data.frame(b.ls[, 1], confint(ls), b.ls[, 4])
kill <- which(! rownames(confint(ls)) %in% attr(b.ls, 'dimnames')[[1]])  # One NA
b.ls <- cbind.data.frame(b.ls[, 1], confint(ls)[-kill, ], b.ls[, 4])
round(head(b.ls), 3)
#             b.ls[, 1]  2.5 % 97.5 % b.ls[, 4]
# (Intercept)     0.045  0.024  0.066     0.000
# abortion       -0.188 -0.565  0.189     0.327
# year021        -0.008 -0.031  0.016     0.525
# year031        -0.024 -0.054  0.007     0.125
# year041        -0.008 -0.051  0.036     0.733
# year051        -0.001 -0.054  0.052     0.972

# BMA (Normal prior) ###########################################################
ms <- modelSelection(f, data = df, includevars = includevars, niter = 1e4,
  priorCoef = normalidprior(), priorDelta = modelbbprior())
b.bma <- coef(ms)
round(head(b.bma), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.029 -0.006  0.063  1.000
# abortion      -0.001  0.000  0.000  0.063
# year021       -0.007 -0.019  0.004  1.000
# year031       -0.028 -0.039 -0.016  1.000
# year041       -0.019 -0.030 -0.008  1.000
# year051       -0.017 -0.029 -0.006  1.000

# BMA (non-local MOM prior)
ms2 <- modelSelection(f, data = df, includevars = includevars, niter = 1e4,
  priorCoef = momprior(), priorDelta = modelbbprior())
b.bma2 <- coef(ms2)
round(head(b.bma2), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.034  0.023  0.044  1.000
# abortion       0.000  0.000  0.000  0.001
# year021       -0.017 -0.031  0.003  1.000
# year031       -0.034 -0.049 -0.019  1.000
# year041       -0.027 -0.041 -0.012  1.000
# year051       -0.025 -0.040 -0.012  1.000

# CIL (EP + Normal) ############################################################
includevars <- includevars[1:(ncol(X) + 2)]
cilfit <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EP', priorCoef = normalidprior(), includevars = includevars)
b <- coef(cilfit[['msfit']])
round(head(b), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.032 -0.005  0.069  1.000
# D             -0.048 -0.148  0.000  0.578
# year021       -0.006 -0.018  0.005  1.000
# year031       -0.025 -0.037 -0.013  1.000
# year041       -0.015 -0.029 -0.002  1.000
# year051       -0.012 -0.027  0.002  1.000

# CIL (EP + non-local MOM)
cilfit2 <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EP', priorCoef = momprior(), includevars = includevars)
b2 <- coef(cilfit2[['msfit']])
round(head(b2), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.036  0.028  0.045  1.000
# D             -0.001  0.000  0.000  0.014
# year021       -0.019 -0.032 -0.007  1.000
# year031       -0.036 -0.050 -0.022  1.000
# year041       -0.028 -0.043 -0.015  1.000
# year051       -0.027 -0.041 -0.014  1.000

# CIL (EB + Normal)
cilfit.eb <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EB', priorCoef = normalidprior(), includevars = includevars)
b.eb <- coef(cilfit.eb[['msfit']])
round(head(b.eb), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.032 -0.003  0.068  1.000
# D             -0.048 -0.146  0.000  0.569
# year021       -0.006 -0.018  0.005  1.000
# year031       -0.026 -0.037 -0.014  1.000
# year041       -0.015 -0.028 -0.002  1.000
# year051       -0.012 -0.027  0.002  1.000

# CIL (EB + non-local MOM)
cilfit2.eb <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EB', priorCoef = momprior(), includevars = includevars)
b2.eb <- coef(cilfit2.eb[['msfit']])
round(head(b2.eb), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.036  0.027  0.045  1.000
# D             -0.001  0.000  0.000  0.015
# year021       -0.019 -0.032 -0.007  1.000
# year031       -0.036 -0.049 -0.022  1.000
# year041       -0.028 -0.042 -0.015  1.000
# year051       -0.027 -0.041 -0.013  1.000

# Save results #################################################################
save(ms, ms2, cilfit, cilfit2, cilfit.eb, cilfit2.eb,
  b.ls, b.bma, b.bma2, b, b2, b.eb, b2.eb,
  file = paste(DATDIR, 'abortion/property_output.RData', sep = ''), compress = TRUE)

# Post-processing analysis
xtable(b.bma[b.bma[, 4] > 0.1, ])
xtable(b.bma2[b.bma2[, 4] > 0.1, ])
xtable(b.eb[b.eb[, 4] > 0.1, ])
xtable(b2.eb[b2.eb[, 4] > 0.1, ])

head(postProb(ms))
head(postProb(ms2))
head(postProb(cilfit.eb))
head(postProb(cilfit2.eb))

# DML analysis #################################################################
dml <- rlassoEffect(x = x[, -1], y = y, d = D, method = 'double selection',
  I3 = includevars[-1:-2])

round(summary(dml)[[1]], 3)  # Inference for treatment effect
#    Estimate. Std. Error t value Pr(>|t|)
# d1    -0.037      0.042  -0.899    0.369
round(confint(dml), 3)
#     2.5 % 97.5 %
# d1 -0.119  0.044

# Results
sel <- names(dml[['coefficients.reg']])[-1:-2]
sel <- sub('x', '', sel)
lmfit <- lm(y ~ D + x[, sel])  # Final DML outcome model
tab <- summary(lmfit)[['coefficients']][, -3]
xtable(tab, digits = c(0, 2, 2, 3))

# BAC ##########################################################################
W <- cbind.data.frame(y, D, x[, -1])
colnames(W) <- c('y', 'd', paste('V', 1:ncol(X), sep = ''))
bac.fitI <- bac(data = W, exposure = 'd', outcome = 'y',
  confounders = paste('V', 1:ncol(X), sep = ''), interactors = NULL,
  familyX = 'gaussian', familyY = 'gaussian', num_its = 5e3, burnM = 5e2,
  burnB = 5e2, thin = 1)#, omega=Inf)
s <- summary(bac.fitI)
ci.bac <- c(s[[1]], s[[2]])
round(ci.bac, 3)
#          2.5%  97.5% 
# -0.218 -0.489  0.042 

# Results
margpp <- colMeans(bac.fitI[['models']][['models']])
#margpp <- colMeans(bac.fitI$models$models)
sum(margpp > 0.5)
# [1] 84

# ACPME ########################################################################
# Fit
x2 <- x; colnames(x2) <- paste('feat', 1:ncol(x2), sep = '')  # don't like names
acm.fit <- regimes::acpme(y = y, Z = D, C = x2[, -1], niter = 1e4)

# Result summary
colMeans(acm.fit[['beta']])  # Treatment effect estimate
# [1] -0.1434458
t(apply(acm.fit[['beta']], 2, quantile, probs = c(0.025, 0.975)))
#            2.5%       97.5%
# [1,] -0.2163965 -0.07103259
mean(rowSums(acm.fit[['alpha']]))
# [1] 54.0566

# Save
fout <- paste(DATDIR, 'abortion/property_output_dml:bac:acpme.RData', sep = '')
save(dml, bac.fitI, acm.fit, file = fout, compress = TRUE)

################################################################################
# 3. MURDER 
cat('* MURDER: *************************************************************\n')
################################################################################
# Read data
file1 <- paste(INPDIR, 'abortion/murd.txt', sep = '')
file2 <- paste(INPDIR, 'abortion/vnames.csv', sep = '')
data <- read.table(file1, sep = '\t', header = FALSE)
vnames <- read.csv(file2, sep = ',', header = FALSE)
names(data) <- gsub(' ', '', vnames[, 1])

# Determine inputs and outputs
y <- data[, 'murder']
x <- as.matrix(data[, -1])
x <- x[, colnames(x) != 'year01']  # Remove reference category
x <- x[, colnames(x) != 'xV0^2']  # Remove fully collinear covariate
#x <- x[, colnames(x) != 'xM0^2']  # Remove fully collinear covariate

# Visual
if (plot.vars == TRUE) { plot.xy(xvar = x[, 'abortion'], yvar = y) }
cor(x[, 'abortion'], y)
# [1] -0.1071187
summary(cor(x[, 1], x[, -1])[1, ])  # Are there strong correlations?
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.29494 -0.05059  0.01344  0.05757  0.11353  0.52052 

# Build modelling matrices
df <- data.frame(y, x)
sel <- grep('year', colnames(df))
for (j in 1:length(sel)) { df[, sel[j]] <- factor(df[, sel[j]]) }
D <- matrix(df[, 'abortion'], ncol = 1)
X <- df[, c(-1:-2), ]
f <- formula(paste('y ~ ', paste(names(df)[-1], collapse = ' + ')))

# Variables that will be forced into the model (year lags)
includevars <- c(TRUE, FALSE,
  rep(TRUE, length(sel)),
  rep(FALSE, ncol(df) - length(sel) - 2))

# OLS ##########################################################################
ls <- lm(f, data = df)
b.ls <- summary(ls)[['coefficients']]
kill <- which(! rownames(confint(ls)) %in% attr(b.ls, 'dimnames')[[1]])  # One NA
b.ls <- cbind.data.frame(b.ls[, 1], confint(ls)[-kill, ], b.ls[, 4])
#b.ls <- cbind.data.frame(b.ls[, 1], confint(ls), b.ls[, 4])
round(head(b.ls), 3)
#             b.ls[, 1]  2.5 % 97.5 % b.ls[, 4]
# (Intercept)     0.060 -0.021  0.142     0.148
# abortion        1.725 -3.696  7.146     0.531
# year021        -0.111 -0.230  0.008     0.068
# year031        -0.073 -0.204  0.057     0.270
# year041        -0.128 -0.314  0.058     0.178
# year051        -0.101 -0.389  0.186     0.488

# BMA (Normal prior) ###########################################################
ms <- modelSelection(f, data = df, includevars = includevars, niter = 1e4,
  priorCoef = normalidprior(), priorDelta = modelbbprior())
b.bma <- coef(ms)
round(head(b.bma), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.050 -0.158  0.263  1.000
# abortion       0.000  0.000  0.000  0.009
# year021       -0.082 -0.149 -0.013  1.000
# year031       -0.037 -0.105  0.032  1.000
# year041       -0.063 -0.132  0.006  1.000
# year051        0.001 -0.068  0.071  1.000

# BMA (non-local MOM prior)
ms2 <- modelSelection(f, data = df, includevars = includevars, niter = 1e4,
  priorCoef = momprior(), priorDelta = modelbbprior())
b.bma2 <- coef(ms2)
round(head(b.bma2), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.147  0.102  0.191      1
# abortion       0.000  0.000  0.000      0
# year021       -0.185 -0.264 -0.107      1
# year031       -0.145 -0.221 -0.071      1
# year041       -0.168 -0.244 -0.092      1
# year051       -0.114 -0.188 -0.045      1

# CIL (EP + Normal) ############################################################
includevars <- includevars[1:(ncol(X) + 2)]
cilfit <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EP', priorCoef = normalidprior(), includevars = includevars)
b <- coef(cilfit[['msfit']])
round(head(b), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.050 -0.167  0.263  1.000
# D             -0.031 -0.587  0.091  0.136
# year021       -0.082 -0.153 -0.012  1.000
# year031       -0.037 -0.106  0.033  1.000
# year041       -0.063 -0.133  0.007  1.000
# year051        0.002 -0.066  0.070  1.000

# CIL (EP + non-local MOM)
cilfit2 <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EP', priorCoef = momprior(), includevars = includevars)
b2 <- coef(cilfit2[['msfit']])
round(head(b2), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.148  0.104  0.192  1.000
# D             -0.001  0.000  0.000  0.003
# year021       -0.186 -0.264 -0.109  1.000
# year031       -0.147 -0.223 -0.072  1.000
# year041       -0.170 -0.247 -0.094  1.000
# year051       -0.117 -0.189 -0.050  1.000

# CIL (EB + Normal)
cilfit.eb <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EB', priorCoef = normalidprior(), includevars = includevars)
b.eb <- coef(cilfit.eb[['msfit']])
round(head(b.eb), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.051 -0.163  0.264  1.000
# D             -0.035 -0.599  0.084  0.136
# year021       -0.082 -0.151 -0.012  1.000
# year031       -0.037 -0.105  0.031  1.000
# year041       -0.063 -0.132  0.008  1.000
# year051        0.001 -0.070  0.072  1.000

# CIL (EB + non-local MOM)
cilfit2.eb <- cil(y = matrix(df[, 'y'], ncol = 1), D = D, X = X, R = 1e4,
  th.search = 'EB', priorCoef = momprior(), includevars = includevars)
b2.eb <- coef(cilfit2.eb[['msfit']])
round(head(b2.eb), 3)
#             estimate   2.5%  97.5% margpp
# (Intercept)    0.149  0.102  0.193  1.000
# D              0.000  0.000  0.000  0.003
# year021       -0.186 -0.264 -0.110  1.000
# year031       -0.147 -0.223 -0.075  1.000
# year041       -0.170 -0.248 -0.093  1.000
# year051       -0.116 -0.190 -0.046  1.000

# Save results #################################################################
save(ms, ms2, cilfit, cilfit2, cilfit.eb, cilfit2.eb,
  b.ls, b.bma, b.bma2, b, b2, b.eb, b2.eb,
  file = paste(DATDIR, 'abortion/murder_output.RData', sep = ''), compress = TRUE)

# Post-processing analysis
xtable(b.bma[b.bma[, 4] > 0.1, ])
xtable(b.bma2[b.bma2[, 4] > 0.1, ])
xtable(b.eb[b.eb[, 4] > 0.1, ])
xtable(b2.eb[b2.eb[, 4] > 0.1, ])

head(postProb(ms))
head(postProb(ms2))
head(postProb(cilfit.eb))
head(postProb(cilfit2.eb))

# DML analysis #################################################################
dml <- rlassoEffect(x = x[, -1], y = y, d = D, method = 'double selection',
  I3 = includevars[-1:-2])

round(summary(dml)[[1]], 3)  # Inference for treatment effect
#    Estimate. Std. Error t value Pr(>|t|)
# d1    -0.116      0.424  -0.273    0.785
round(confint(dml), 3)
#     2.5 % 97.5 %
# d1 -0.948  0.716

# Results
sel <- names(dml[['coefficients.reg']])[-1:-2]
sel <- sub('x', '', sel)
lmfit <- lm(y ~ D + x[, sel])  # Final DML outcome model
tab <- summary(lmfit)[['coefficients']][, -3]
xtable(tab, digits = c(0, 2, 2, 3))

# BAC ##########################################################################
W <- cbind.data.frame(y, D, x[, -1])
colnames(W) <- c('y', 'd', paste('V', 1:ncol(X), sep = ''))
bac.fitI <- bac(data = W, exposure = 'd', outcome = 'y',
  confounders = paste('V', 1:ncol(X), sep = ''), interactors = NULL,
  familyX = 'gaussian', familyY = 'gaussian', num_its = 5e3, burnM = 5e2,
  burnB = 5e2, thin = 1)#, omega=Inf)
s <- summary(bac.fitI)
ci.bac <- c(s[[1]], s[[2]])
round(ci.bac, 3)
#         2.5%  97.5% 
# -0.775 -4.407  2.714 

# Results
margpp <- colMeans(bac.fitI[['models']][['models']])
#margpp <- colMeans(bac.fitI$models$models)
sum(margpp > 0.5)
# [1] 98

# ACPME ########################################################################
# Fit
x2 <- x; colnames(x2) <- paste('feat', 1:ncol(x2), sep = '')  # don't like names
acm.fit <- regimes::acpme(y = y, Z = D, C = x2[, -c(1, 292)], niter = 1e4)

# Result summary
colMeans(acm.fit[['beta']])  # Treatment effect estimate
# [1] -0.5083816
t(apply(acm.fit[['beta']], 2, quantile, probs = c(0.025, 0.975)))
#            2.5%      97.5%
# [1,] -0.9092138 -0.1066207
mean(rowSums(acm.fit[['alpha']]))
# [1] 41.4338

# Save
fout <- paste(DATDIR, 'abortion/murder_output_dml:bac:acpme.RData', sep = '')
save(dml, bac.fitI, acm.fit, file = fout, compress = TRUE)

# # End computaion
# stop('controlled error.')

################################################################################
# PLOTS
################################################################################
#rm(list = ls()[ls() != 'PATH'])
fv1 <- paste(DATDIR, 'abortion/violent_output.RData', sep = '')
fv2 <- paste(DATDIR, 'abortion/violent_output_dml:bac:acpme.RData', sep = '')
fp1 <- paste(DATDIR, 'abortion/property_output.RData', sep = '')
fp2 <- paste(DATDIR, 'abortion/property_output_dml:bac:acpme.RData', sep = '')
fm1 <- paste(DATDIR, 'abortion/murder_output.RData', sep = '')
fm2 <- paste(DATDIR, 'abortion/murder_output_dml:bac:acpme.RData', sep = '')

# Violent data
load(fv1); cat('Loaded file:', fv1, '\n')
load(fv2); cat('Loaded file:', fv2, '\n')

b.viol <- do.call(rbind, list(b.dnl = c(-0.13, -0.18, -0.08, 0),
  b.ls = b.ls[2, ],
  b.dml = c(summary(dml)[[1]][1], confint(dml), summary(dml)[[1]][2]),
  b.bac = c(summary(bac.fitI)[[1]], summary(bac.fitI)[[2]], 1),
  b.acm = c(colMeans(acm.fit[['beta']]), t(apply(acm.fit[['beta']], 2, quantile, probs = c(0.025, 0.975))), 1),
  b.bma = b.bma[2, ], b.bma2 = b.bma2[2, ],
  b = b[2, ], b2 = b2[2, ], b.eb = b.eb[2, ], b2.eb = b2.eb[2, ]))

# Property data
load(fp1); cat('Loaded file:', fp1, '\n')
load(fp2); cat('Loaded file:', fp2, '\n')

b.prop <- do.call(rbind, list(b.dnl = c(-0.90, -0.94, -0.87, 0),
  b.ls = b.ls[2, ],
  b.dml = c(summary(dml)[[1]][1], confint(dml), summary(dml)[[1]][2]),
  b.bac = c(summary(bac.fitI)[[1]], summary(bac.fitI)[[2]], 1),
  b.acm = c(colMeans(acm.fit[['beta']]), t(apply(acm.fit[['beta']], 2, quantile, probs = c(0.025, 0.975))), 1),
  b.bma = b.bma[2, ], b.bma2 = b.bma2[2, ],
  b = b[2, ], b2 = b2[2, ], b.eb = b.eb[2, ], b2.eb = b2.eb[2, ]))

# Murder data
load(fm1); cat('Loaded file:', fm1, '\n')
load(fm2); cat('Loaded file:', fm2, '\n')

b.murd <- do.call(rbind, list(b.dnl = c(-0.12, -0.21, -0.03, 0.01),
  b.ls = b.ls[2, ],
  b.dml = c(summary(dml)[[1]][1], confint(dml), summary(dml)[[1]][2]),
  b.bac = c(summary(bac.fitI)[[1]], summary(bac.fitI)[[2]], 1),
  b.acm = c(colMeans(acm.fit[['beta']]), t(apply(acm.fit[['beta']], 2, quantile, probs = c(0.025, 0.975))), 1),
  b.bma = b.bma[2, ], b.bma2 = b.bma2[2, ],
  b = b[2, ], b2 = b2[2, ], b.eb = b.eb[2, ], b2.eb = b2.eb[2, ]))

#stop('CONTROLLED ERROR.')
###

# VIOLENT CRIME
aux <- round(b.viol[1:9, 4], 3)
aux[nchar(aux) == 1] <- paste(aux[nchar(aux) == 1], '000', sep = '.')
for (k in 2:4) {
  txt <- paste(rep(0, 5 - k), sep = '', collapse = '')
  aux[nchar(aux) == k] <- paste(aux[nchar(aux) == k], txt, sep = '')
}
aux[which(aux == '0.000')] <- '<0.001'

fplot <- paste(FIGDIR, 'inference_viol.pdf', sep = '')
pdf(fplot, height = 5, width = 8)
par(mfrow = c(1, 1), mar = c(4, 6, 0.5, 4))
plot(NA, xlim = c(-1.4, 1.5), ylim = c(9.5, 0.5), yaxt = 'n',
  ylab = '', xlab = 'Treatment effect (violent crime)')
abline(v = 0, lty = 2)
for (j in seq(-2, 2, 0.5)) { abline(v = j, col = 'lightgray', lty = 3) }
for (i in 1:9) {
  abline(h = i, col = 'lightgray', lty = 3)
  points(b.viol[i, 1], i, pch = 16, cex = 1.25)
  segments(x0 = b.viol[i, 2], x1 = b.viol[i, 3], y0 = i, y1 = i, lwd = 1.25)
}
axis(2, at = 1:9, las = 2,
  c('Don. & Lev.', 'OLS (all)', 'DML', 'BAC', 'ACPME', 'BMA (Norm.)',
    'BMA (MOM)', 'CIL (Norm.)', 'CIL (MOM)'))
#axis(4, at = 1:8, las = 2, round(b.viol[8:1, 4], 3))
axis(4, at = 1:9, las = 2, aux, tick = FALSE)
dev.off(); cat('Plotted file:', fplot, '\n')

# PROPERTY CRIME
aux <- round(b.prop[1:9, 4], 3)
aux[nchar(aux) == 1] <- paste(aux[nchar(aux) == 1], '000', sep = '.')
for (k in 2:4) {
  txt <- paste(rep(0, 5 - k), sep = '', collapse = '')
  aux[nchar(aux) == k] <- paste(aux[nchar(aux) == k], txt, sep = '')
}
aux[which(aux == '0.000')] <- '<0.001'

fplot <- paste(FIGDIR, 'inference_prop.pdf', sep = '')
pdf(fplot, height = 5, width = 8)
par(mfrow = c(1, 1), mar = c(4, 6, 0.5, 4))
plot(NA, xlim = c(-1.4, 1.5), ylim = c(9.5, 0.5), yaxt = 'n',
  ylab = '', xlab = 'Treatment effect (violent crime)')
abline(v = 0, lty = 2)
for (j in seq(-2, 2, 0.5)) { abline(v = j, col = 'lightgray', lty = 3) }
for (i in 1:9) {
  abline(h = i, col = 'lightgray', lty = 3)
  points(b.prop[i, 1], i, pch = 16, cex = 1.25)
  segments(x0 = b.prop[i, 2], x1 = b.prop[i, 3], y0 = i, y1 = i, lwd = 1.25)
}
axis(2, at = 1:9, las = 2,
  c('Don. & Lev.', 'OLS (all)', 'DML', 'BAC', 'ACPME', 'BMA (Norm.)',
    'BMA (MOM)', 'CIL (Norm.)', 'CIL (MOM)'))
#axis(4, at = 1:8, las = 2, round(b.prop[8:1, 4], 3))
axis(4, at = 1:9, las = 2, aux, tick = FALSE)
dev.off(); cat('Plotted file:', fplot, '\n')

# MURDER
aux <- round(b.murd[1:9, 4], 3)
aux[nchar(aux) == 1] <- paste(aux[nchar(aux) == 1], '000', sep = '.')
for (k in 2:4) {
  txt <- paste(rep(0, 5 - k), sep = '', collapse = '')
  aux[nchar(aux) == k] <- paste(aux[nchar(aux) == k], txt, sep = '')
}
aux[which(aux == '0.000')] <- '<0.001'

fplot <- paste(FIGDIR, 'inference_murd.pdf', sep = '')
pdf(fplot, height = 5, width = 8)
par(mfrow = c(1, 1), mar = c(4, 6, 0.5, 4))
plot(NA, xlim = c(-1.4, 1.5), ylim = c(9.5, 0.5), yaxt = 'n',
  ylab = '', xlab = 'Treatment effect (violent crime)')
abline(v = 0, lty = 2)
for (j in seq(-2, 2, 0.5)) { abline(v = j, col = 'lightgray', lty = 3) }
for (i in 1:9) {
  abline(h = i, col = 'lightgray', lty = 3)
  points(b.murd[i, 1], i, pch = 16, cex = 1.25)
  segments(x0 = b.murd[i, 2], x1 = b.murd[i, 3], y0 = i, y1 = i, lwd = 1.25)
}
axis(2, at = 1:9, las = 2,
  c('Don. & Lev.', 'OLS (all)', 'DML', 'BAC', 'ACPME', 'BMA (Norm.)',
    'BMA (MOM)', 'CIL (Norm.)', 'CIL (MOM)'))
#axis(4, at = 1:8, las = 2, round(b.murd[8:1, 4], 3))
axis(4, at = 1:9, las = 2, aux, tick = FALSE)
dev.off(); cat('Plotted file:', fplot, '\n')

