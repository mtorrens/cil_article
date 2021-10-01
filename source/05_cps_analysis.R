################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '05_cps_analysis.R', sep = ''))
################################################################################

################################################################################
# Preparation
################################################################################
# Global parameters
eta <- 0.05  # Significance level
lpen <- 'lambda.min'  # LASSO type penalty
force <- FALSE  # Force execution even if file exists

# Load data
file.in <- paste(DATDIR, 'cps/cps_model_data_year.RData', sep = '')
raw <- load(file = file.in); cat('Loaded file:', file.in, '\n'); rm(file.in)
y <- get(raw[1])  # 615366 x 1
D <- get(raw[2])  # 615366 x 1
X <- get(raw[3])  # 615366 x 240
state <- get(raw[4])  # 615366
year <- get(raw[5])  # 615366

# Parameters
years <- as.character(unique(year))
states <- as.character(unique(state))

# Supress some persistent redundant columns
cols.out <- c('movedhouse_Abroad', 'famrel_notmember', 'relhh_grandchild',
  'relhh_fosterchild', 'isrel2hh', 'propoi_incsurv1', 'propoi_incsurv2',
  'propoi_incdisa2', 'is2019', 'region_northeast', 'region_midwest',
  'region_south')
X <- X[, which(! colnames(X) %in% cols.out)]  # 615366 x 232

# Initially we had only one treatment: gender (female)
# We consider multiple treatments now: race(black), hispanic, bpl(latam)
new.treatments <- c('race_Black', 'hispanic', 'bpl_latam')
D <- cbind(D, X[, new.treatments])  # 615366 x 4
X <- X[, which(! colnames(X) %in% new.treatments)]  # 615366 x 228
colnames(D) <- c('female', new.treatments)
Xor <- X  # Save original X for later

# Columns to standardize later
std.cols <- c('logschllunch', 'logspmlunch', 'atelunch', 'logctccrd',
  'logfoodstamp', 'nchild', 'nsibs', 'stampno', 'logspmchxpns',
  'logspmnchild', 'logspmcapxpns', 'logfica', 'logspmchsup', 'logspmmedxpns',
  'hiunpers', 'grpwho1', 'logspmnadults', 'logfedtaxac', 'logstataxac',
  'stampmo', 'logotherpersincome', 'hhrespln', 'logotherfamincome', 'pernum',
  'numprec', 'lineno', 'famsize', 'eldch', 'yngch', 'wksunemly',
  'numweekswkd', 'hoursworked', 'uhrsworkly', 'age')

################################################################################
# Dummies + interactions + Sum to zero constraint
################################################################################
st.dum <- conv.dummy(state)
colnames(st.dum) <- gsub('cat_', 'st', colnames(st.dum))
for (i in 1:ncol(D)) {
  ints <- D[, i] * st.dum
  colnames(ints) <- paste(
    substr(colnames(st.dum), 3, nchar(colnames(st.dum))), 'T', i, sep = '')
  ints[which(ints[, paste('CAT', i, sep = '')] == 1), ] <- -1  # STZ
  ints <- ints[, -which(colnames(ints) == paste('CAT', i, sep = ''))]
  if (i == 1) {
    int.dum <- ints
  } else {
    int.dum <- cbind(int.dum, ints)
  }
}
st.dum <- st.dum[, -which(colnames(st.dum) == 'stCA')]  # Ref. CA
################################################################################

################################################################################
# FULL DATA - Y2Y - STATE INTERACTIONS - POSTERIOR SAMPLES
################################################################################
# Simulations
fout <- paste(DATDIR, 'cps_stateints_y2y_post_FULL.RData', sep = '')
if (! file.exists(fout) | force == TRUE) {
  #registerDoMC(cores = detectCores() - 1)
  registerDoMC(cores = 2)
  y2y.bmkP <- foreach(v = 1:length(years), .errorhandling = 'pass') %dopar% {
    # Choose year
    ye <- years[v]
    cat('\r*** Year:', ye)

    # Adapt design matrix
    n <- sum(year == ye)
    so <- which(year == ye)
    vo <- names(which(apply(X[so, ], 2, var) == 0))[-1]
    Xs <- cbind(X[so, -which(colnames(X) %in% vo)], st.dum[so, ])
    Ds <- D[so, ]
    Is <- int.dum[so, ]
    ys <- y[so, ]
    Zs <- cbind(Ds, Is, Xs)
    Ws <- cbind(ys, Zs)

    # In an earlier version we computed every year, now we only compare 2010-19
    if (ye %in% c('2010', '2019')) {
      # Compute main statistics
      stats <- cps.methods(ys, Ds, Xs, Is, minimal = TRUE, eta, lpen,
        out.post = TRUE)

      # Recover intermediate data
      mod1.coef <- stats[length(stats) - 1:0]
      post.bma <- stats[length(stats) - 4][[1]]
      post.nmP <- stats[length(stats) - 3][[1]][['msfit']]
      #post.nmB <- stats[length(stats) - 2][[1]][['msfit']]
      stats <- stats[-(length(stats) - 4:0)]

      # Outcome
      out <- list(c(year = ye, n = n, p = ncol(Zs),
        nact = unname(c(colSums(Ds),
          apply(int.dum[so, ], 2, function(x) { length(which(x == 1)) }))),
        unlist(stats)), mod1.coef[[1]])
      post <- list(post.bma = post.bma, post.nmP = post.nmP)#, post.nmB)
    } else {  # No computations for other years
      # Outcome (results are empty)
      out <- list(c(year = ye, n = n, p = ncol(Zs),
        nact = unname(c(colSums(Ds),
          apply(int.dum[so, ], 2, function(x) { length(which(x == 1)) })))))
      post <- list(post.bma = NA, post.nmP = NA)#, post.nmB)
    }

    # Save partial results
    filew <- paste(DATDIR, 'cps_stateints_y2y_post_', ye, '.RData', sep = '')
    save(out, post, file = filew); cat('\nSaved file:', filew, '\n'); rm(filew)

    # Return relevant input and output
    return(list(out = out, post = post))
  }; cat('** Done.\n')

  # Save total results
  save(y2y.bmkP, file = fout); cat('\nSaved file:', fout, '\n')
} else {
  y2y.bmkP <- get(load(file = fout)); cat('Loaded file:', fout, '\n')  
}; rm(fout)

##############################################################################
# OTHER FUNCTIONALS USING POSTERIOR MODEL SAMPLES
##############################################################################
# Define data
ye <- years[c(1, 10)]  # 2010 and 2019
so <- which(year %in% ye)
vo <- names(which(apply(X[so, ], 2, var) == 0))[-1]
Xs <- cbind(X[so, -which(colnames(X) %in% vo)], st.dum[so, ])
Ts <- cbind(D[so, ], int.dum[so, ])
ys <- y[so, ]

# Prediction: E(d_{n+1} | x_{n+1}) for each treatment
good <- colnames(Xs)[! grepl('bpl', colnames(Xs)) &
                     ! grepl('race', colnames(Xs)) &
                     ! grepl('nativity', colnames(Xs)) &
                     ! grepl('citizen',colnames(Xs))]
Dp <- matrix(NA, nrow = length(so), ncol = ncol(D))
for (tr in 1:ncol(D)) {
  auX <- cbind.data.frame(Ts[, tr], Xs[, good][, -1])
  colnames(auX) <- c('d', paste0('x', 1:(ncol(auX) - 1)))
  trM <- glm(d ~ ., data = auX, family = binomial(link = 'logit'))
  Dp[, tr] <- trM[['fitted.values']]
}

# Analysis per state: predictions conditional on each of them
good <- good[which(! (grepl('^st', good) & nchar(good) == 4))]
DpST <- vector(length = length(states), mode = 'list')
names(DpST) <- states
st0 <- state[so]
for (st in states) {
  cat('\r* State:', st, '-')
  qui <- which(st0 == st)
  sti <- which(states == st)
  DpST[[sti]] <- matrix(NA, nrow = sum(st0 == states[sti]), ncol = ncol(D))
  for (tr in 1:ncol(D)) {
    cat('', tr)
    auX <- cbind.data.frame(Ts[qui, tr], Xs[qui, good][, -1])
    auX <- auX[, which(apply(auX, 2, var) != 0)]
    colnames(auX) <- c('d', paste0('x', 1:(ncol(auX) - 1)))
    trM <- glm(d ~ ., data = auX, family = binomial(link = 'logit'))
    auX <- auX[, c('d', names(which(! is.na(coef(trM))))[-1])]  # some NAs
    trM <- glm(d ~ ., data = auX, family = binomial(link = 'logit'))
    DpST[[sti]][, tr] <- trM[['fitted.values']]
  }
  cat(' - Done\n')
}

# Computation of related functionals
fout <- paste(DATDIR, 'cps_stateints_postsamples_FULL.RData', sep = '')
if (! file.exists(fout) | force == TRUE) {
  registerDoMC(cores = 2)
  ps.bmk <- foreach(v = c(1, 10), .errorhandling = 'pass') %dopar% {
    # Posterior samples
    post <- y2y.bmkP[[v]][['post']]

    # We run rnlp for both methods and compare measures
    set.seed(666)
    nlp <- mombf::momprior(tau = 0.348)
    gs.bma <- mombf::rnlp(msfit = post[[1]], priorCoef = nlp, niter = 1e4)
    gs.nmP <- mombf::rnlp(msfit = post[[2]], priorCoef = nlp, niter = 1e4)

    # Define columns in the model(s)
    colsX <- colnames(gs.bma[, -1])[colnames(gs.bma[, -1]) %in% colnames(Xs)]
    colsD <- colnames(gs.bma[, -1])[colnames(gs.bma[, -1]) %in% colnames(Ts)]

    # Select data
    Xy <- Xs[which(year[so] %in% years[v]), colsX]
    Ty <- Ts[which(year[so] %in% years[v]), colsD]
    XyC <- apply(Xy, 2, function(x) { x - mean(x) })
    TyC <- apply(Ty, 2, function(x) { x - mean(x) })

    # Deviation from average salary: exp{|h(x1, beta1)|}
    colsDm <- colsD[1:4]
    Tsm <- Ts[, colsDm]
    m5.bma <- apply(gs.bma[, -1], 1, function(x) {
      ans <- sample.expabsh(D1 = Tsm, Dp = Dp, a = x[colsDm])
      return(ans)
    })
    m5.nmP <- apply(gs.nmP[, -1], 1, function(x) {
      ans <- sample.expabsh(D1 = Tsm, Dp = Dp, a = x[colsDm])
      return(ans)
    })

    # State by state analysis
    colsDm <- colsD[1:4]
    doMC::registerDoMC(cores = 1)
    m5.nmP2 <- foreach(st = states, .errorhandling = 'pass') %dopar% {
      Dp2 <- DpST[[which(states == st)]]
      Tsm2 <- Ts[which(st0 == st), colsDm]
      res <- apply(gs.nmP[, -1], 1, function(x) {
        inters <- x[grep(paste('^', st, 'T', sep = ''), names(x))]
        if (st == 'CA') { inters <- rep(0, length(colsDm)) }  # Ref. category
        newx <- x[colsDm] + inters
        ans <- sample.expabsh(D1 = Tsm2, Dp = Dp2, a = newx)
        return(ans)
      })
      return(res)
    }

    # Save results
    res <- list(gs.bma = gs.bma, gs.nmP = gs.nmP, m5.bma = m5.bma,
      m5.nmP = m5.nmP, m5.nmP2 = m5.nmP2)
    fout <- paste(DATDIR, 'cps_stateints_postsamples_', years[v],
      '.RData', sep = '')
    save(res, file = fout); cat('Saved file:', fout, '\n'); rm(fout)
    return(res)
  }; cat('** Done.\n')

  # Save total results
  save(ps.bmk, file = fout); cat('Saved file:', fout, '\n')
} else {
  ps.bmk <- get(load(file = fout)); cat('Loaded file:', fout, '\n')  
}; rm(fout)

##############################################################################
# FULL DATA - Y2Y - STATE INTERACTIONS - 100 FAKE PREDS - POSTERIOR SAMPLES
##############################################################################
# Simulations
fout <- paste(DATDIR, 'cps_stateintsFP_y2y_FULL.RData', sep = '')
if (! file.exists(fout) | force == TRUE) {
  # Generate 100 fake predictors (regular correlation)
  M <- gen.fake.preds(nfr = round(100 / ncol(D)), centre = 1.5, D = D)

  # Add them to the design matrix
  X <- cbind(Xor, M)  # 123265 x 326

  registerDoMC(cores = 2)
  y2y.1hRP <- foreach(v = c(1, 10), .errorhandling = 'pass') %dopar% {
    # Choose year
    ye <- years[v]
    cat('\r*** Year:', ye)

    # Adapt design matrix
    n <- sum(year == ye)
    so <- which(year == ye)
    vo <- names(which(apply(X[so, ], 2, var) == 0))[-1]
    Xs <- cbind(X[so, -which(colnames(X) %in% vo)], st.dum[so, ])
    Ds <- D[so, ]
    Is <- int.dum[so, ]
    ys <- y[so, ]
    Zs <- cbind(Ds, Is, Xs)
    Ws <- cbind(ys, Zs)

    # Compute main statistics
    stats <- cps.methods(ys, Ds, Xs, Is, minimal = TRUE, eta, lpen,
      out.post = TRUE)

    # Recover intermediate data
    mod1.coef <- stats[length(stats) - 1:0]
    post.bma <- stats[length(stats) - 4][[1]]
    post.nmP <- stats[length(stats) - 3][[1]][['msfit']]
    #post.nmB <- stats[length(stats) - 2][[1]][['msfit']]
    stats <- stats[-(length(stats) - 4:0)]

    # Outcome
    out <- list(c(year = ye, n = n, p = ncol(Zs),
      nact = unname(c(colSums(Ds),
        apply(int.dum[so, ], 2, function(x) { length(which(x == 1)) }))),
      unlist(stats)), mod1.coef[[1]])
    post <- list(post.bma = post.bma, post.nmP = post.nmP)#, post.nmB)

    # Save partial results
    filew <- paste(DATDIR, 'cps_stateintsFP_y2y_post_', ye, '.RData', sep = '')
    save(out, post, file = filew); cat('\nSaved file:', filew, '\n'); rm(filew)

    # Return relevant input and output
    return(list(out = out, post = post))
  }; cat('** Done.\n')

  # Save total results
  save(y2y.1hRP, file = fout); cat('\nSaved file:', fout, '\n')
} else {
  y2y.1hRP <- get(load(file = fout)); cat('Loaded file:', fout, '\n')  
}; rm(fout)

##############################################################################
# FULL DATA - Y2Y - STATE INTERACTIONS - 200 FAKE PREDS - POSTERIOR SAMPLES
##############################################################################
# Simulations
fout <- paste(DATDIR, 'cps_stateintsFP200_y2y_FULL.RData', sep = '')
if (! file.exists(fout) | force == TRUE) {
  # Generate 200 fake predictors (regular correlation)
  M1 <- gen.fake.preds(nfr = round(100 / ncol(D)), centre = 1.5, D = D)
  M2 <- gen.fake.preds(nfr = round(100 / ncol(D)), centre = 1.5, D = D,
    seed = 999)
  M <- cbind(M1, M2)
  # M <- gen.fake.preds(nfr = round(200 / ncol(D)), centre = 1.5, D = D)

  # Add them to the design matrix
  X <- cbind(Xor, M)  # 123265 x 326

  registerDoMC(cores = 2)
  y2y.2hRP <- foreach(v = c(1, 10), .errorhandling = 'pass') %dopar% {
    # Choose year
    ye <- years[v]
    cat('\r*** Year:', ye)

    # Adapt design matrix
    n <- sum(year == ye)
    so <- which(year == ye)
    vo <- names(which(apply(X[so, ], 2, var) == 0))[-1]
    Xs <- cbind(X[so, -which(colnames(X) %in% vo)], st.dum[so, ])
    Ds <- D[so, ]
    Is <- int.dum[so, ]
    ys <- y[so, ]
    Zs <- cbind(Ds, Is, Xs)
    Ws <- cbind(ys, Zs)

    # Compute main statistics
    stats <- cps.methods(ys, Ds, Xs, Is, minimal = TRUE, eta, lpen,
      out.post = TRUE)

    # Recover intermediate data
    mod1.coef <- stats[length(stats) - 1:0]
    post.bma <- stats[length(stats) - 4][[1]]
    post.nmP <- stats[length(stats) - 3][[1]][['msfit']]
    #post.nmB <- stats[length(stats) - 2][[1]][['msfit']]
    stats <- stats[-(length(stats) - 4:0)]

    # Outcome
    out <- list(c(year = ye, n = n, p = ncol(Zs),
      nact = unname(c(colSums(Ds),
        apply(int.dum[so, ], 2, function(x) { length(which(x == 1)) }))),
      unlist(stats)), mod1.coef[[1]])
    post <- list(post.bma = post.bma, post.nmP = post.nmP)#, post.nmB)

    # Save partial results
    filew <- paste(DATDIR, 'cps_stateintsFP200_y2y_post_', ye, '.RData', sep = '')
    save(out, post, file = filew); cat('\nSaved file:', filew, '\n'); rm(filew)

    # Return relevant input and output
    return(list(out = out, post = post))
  }; cat('** Done.\n')

  # Save total results
  save(y2y.2hRP, file = fout); cat('\nSaved file:', fout, '\n')
} else {
  y2y.2hRP <- get(load(file = fout)); cat('Loaded file:', fout, '\n')
}; rm(fout)
# END OF SCRIPT
