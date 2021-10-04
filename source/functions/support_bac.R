################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(FCNDIR, 'support_bac.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * GetACE
# * calc.ACE
# * get.predictorsY
# * calclogpost
# * update_alpha
# * find.models
# * bac2
################################################################################

GetACE <- function (nsample, model, population, exposure, outcome, predictorsY, 
    interactions, confounders, familyY, thin, burnB, data) 
{
    startvalue = NA
    formu = paste(outcome, " ~ ", paste(predictorsY[model == 
        1], collapse = "+"), sep = "")
    kk = coef(glm(formu, data = data, family = familyY))
    model[(which(model == 1))[which(is.na(kk[(names(kk) != "(Intercept)") & 
        (names(kk) != exposure)]))]] = 0
    formu = paste(outcome, " ~ ", paste(predictorsY[model == 
        1], collapse = "+"), sep = "")
    b0 = 0
    B0 = 0
    if (familyY == "gaussian") {
        para = data.matrix(MCMCregress(formu, data = data, mcmc = nsample * 
            thin, b0 = b0, B0 = B0, burnin = burnB, thin = thin)[, 
            1:(sum(model == 1) + 1)])
        if (nsample > 1) {
            para = t(para)
        }
    }
    if (familyY == "binomial") {
        para = t(data.matrix(MCMClogit(formu, data = data, mcmc = nsample * 
            thin, b0 = b0, B0 = B0, burnin = burnB, thin = thin, 
            beta.start = startvalue, tune = 0.4)))
    }
    if (familyY == "poisson") {
        para = t(data.matrix(MCMCpoisson(formu, data = data, 
            mcmc = nsample * thin, b0 = b0, B0 = B0, burnin = burnB, 
            thin = thin, beta.start = startvalue, tune = 0.5)))
    }
    data1 = data.matrix(data)
    n = nrow(data)
    m = length(confounders)
    mm = m + length(interactions)
    if (!is.null(interactions)) {
        data_c = cbind(rep(1, n), data1[, confounders[model[1:m] == 
            1]], rep(0, n), rep(0, n) * data1[, confounders[interactions[model[(m + 
            1):mm] == 1]]])
        data_t = cbind(rep(1, n), data1[, confounders[model[1:m] == 
            1]], rep(1, n), rep(1, n) * data1[, confounders[interactions[model[(m + 
            1):mm] == 1]]])
    } else {
        data_c = cbind(rep(1, n), data1[, confounders[model[1:m] == 
            1]], rep(0, n))
        data_t = cbind(rep(1, n), data1[, confounders[model[1:m] == 
            1]], rep(1, n))
    }
    data_t = data_t[population, ]
    data_c = data_c[population, ]
    predict1 = data_t %*% para
    predict0 = data_c %*% para
    if (familyY == "binomial") {
        predict1 = 1/(1 + exp((-1) * predict1))
        predict0 = 1/(1 + exp((-1) * predict0))
    }
    if (familyY == "poisson") {
        predict1 = exp(predict1)
        predict0 = exp(predict0)
    }
    weights = t(rdirichlet(nsample, rep((-1) + 2, sum(population)))) * 
        sum(population)
    ACE.boot = colMeans(weights * (predict1 - predict0))
    return(list(ACE.boot = ACE.boot, para = para))
}

calc.ACE <- function (population, exposure, outcome, predictorsY, interactions, 
    confounders, MY, familyY, thin, burnB, data) 
{
    temp = find.models(MY)
    ACE.boot = NULL
    for (h in 1:length(temp$count)) {
        model = temp$models[h, ]
        vv = GetACE(temp$counts[h], model, population, exposure, 
            outcome, predictorsY, interactions, confounders, 
            familyY, thin, burnB, data)
        ACE.boot = c(ACE.boot, vv$ACE.boot)
    }
    return(list(ACE.boot = ACE.boot))
}

get.predictorsY <- function (exposure, confounders, interactions) 
{
    if (!is.null(interactions)) {
        predictorsY = c(confounders, paste(exposure, confounders[interactions], 
            sep = ":"), exposure)
    }
    else {
        predictorsY = c(confounders, exposure)
    }
    return(predictorsY)
}

calclogpost <- function (Yvar, Xvar, family, data) 
{
    if (length(Xvar) > 0) {
        formu = paste(Yvar, " ~ ", paste(Xvar, collapse = "+"), 
            sep = "")
    }
    else {
        formu = paste(Yvar, " ~ 1", sep = "")
    }
    fit = glm(formu, data = data, family = family)
    return((-1) * BIC(fit)/2)
}

update_alpha <- function (Yvar, Xvar, M0, M1, a0, tuning, family, data) 
{
    a1 = calclogpost(Yvar, Xvar[M1 == 1], family, data)
    BF = tuning * exp(a1 - a0)
    flag = 1
    if (BF < 1) {
        if (BF < runif(1)) {
            flag = 0
        }
    }
    if (flag == 1) {
        return(list(M = M1, a0 = a1))
    }
    else {
        return(list(M = M0, a0 = a0))
    }
}

find.models <- function (MY) 
{
    bb = apply(MY, 1, paste, collapse = " ")
    count = table(bb)
    models = matrix(as.numeric((sapply(strsplit(names(count), 
        split = " "), unlist))), byrow = T, nrow = length(count))
    return(list(models = models, counts = count))
}

bac2 <- function (data, exposure, outcome, confounders, interactors, 
    familyX, familyY, omega = Inf, num_its, burnM, burnB, thin, 
    population = NULL) 
{
    interactions = match(interactors, confounders)
    if (length(interactions) == 0) {
        interactions = NULL
    }
    nits = num_its + burnM
    predictorsY = get.predictorsY(exposure, confounders, interactions)
    predictorsX = confounders
    n = nrow(data)
    m = length(confounders)
    mm = m + length(interactions)
    if (length(population) == 0) {
        population = rep(T, n)
    }
    MX = matrix(nrow = nits, ncol = m)
    MY = matrix(nrow = nits, ncol = mm + 1)
    # M0X = rep(1, m)  # ORIGINAL CODE
    # M0Y = rep(1, mm + 1)  # ORIGINAL CODE
    M0X = rep(0, m)  # ADDED CODE
    M0Y = c(rep(0, mm), 1)  # ADDED CODE
    # IDEA: INSTEAD OF STARTING AT FULL MODEL START AT NULL AND DON'T GET STUCK
    MX[1, ] = M0X
    MY[1, ] = M0Y
    a0X = calclogpost(exposure, predictorsX[M0X == 1], familyX, 
        data)
    a0Y = calclogpost(outcome, predictorsY[M0Y == 1], familyY, 
        data)
    ptime = proc.time()[3]
    for (k in 2:nits) {
        M0X = MX[k - 1, ]
        M1X = M0X
        if (omega == Inf) {
            candidate = (1:m)[MY[k - 1, 1:m] == 1]
            if (length(candidate) > 0) {
                if (length(candidate) > 1) {
                  j = sample(candidate, 1)
                }
                else {
                  j = candidate
                }
                M1X[j] = 1 - M1X[j]
                tuning = 1
                temp = update_alpha(exposure, predictorsX, M0X, 
                  M1X, a0X, tuning, familyX, data)
                M0X = temp$M
                a0X = temp$a0
            }
        }
        else {
            j = sample(1:m, 1)
            M1X[j] = 1 - M1X[j]
            if (MY[k - 1, j] == 1) {
                tuning = 1
            }
            else {
                tuning = omega^(M0X[j] - M1X[j])
            }
            temp = update_alpha(exposure, predictorsX, M0X, M1X, 
                a0X, tuning, familyX, data)
            M0X = temp$M
            a0X = temp$a0
        }
        MX[k, ] = M0X
        M0Y = MY[k - 1, ]
        M1Y = M0Y
        ind = rep(0, m)
        ind[interactions] = M0Y[-c(1:m, mm + 1)]
        if (omega == Inf) {
            candidate = c((1:m)[(ind == 0) & (MX[k, ] == 0)], 
                (1:mm)[-(1:m)][M0Y[interactions] == 1])
            if (length(candidate) > 0) {
                if (length(candidate) > 1) {
                  j = sample(candidate, 1)
                }
                else {
                  j = candidate
                }
                M1Y[j] = 1 - M1Y[j]
                tuning = 1
                temp = update_alpha(outcome, predictorsY, M0Y, 
                  M1Y, a0Y, tuning, familyY, data)
                M0Y = temp$M
                a0Y = temp$a0
            }
        }
        else {
            candidate = c((1:m)[(ind == 0)], (1:mm)[-(1:m)])
            if (length(candidate) > 0) {
                if (length(candidate) > 1) {
                  j = sample(candidate, 1)
                }
                else {
                  j = candidate
                }
                M1Y[j] = 1 - M1Y[j]
                tuning = 1
                if (j <= m) {
                  if (MX[k, j] == 1) {
                    tuning = omega^(M1Y[j] - M0Y[j])
                  }
                }
                temp = update_alpha(outcome, predictorsY, M0Y, 
                  M1Y, a0Y, tuning, familyY, data)
                M0Y = temp$M
                a0Y = temp$a0
            }
        }
        MY[k, ] = M0Y
    }
    MX = MX[-(1:burnM), ]
    MY = MY[-(1:burnM), ]
    include = which(MY[, ncol(MY)] == 1)  # ADDED CODE
    MX = MX[include, ]  # ADDED CODE
    MY = MY[include, ]  # ADDED CODE
    # IDEA: IN FULL MODEL d == 1 ALWAYS; NOW MAYBE NOT, BUT BAC REQUIRES IT -- KILL REST
    include = which(rowSums(MY) <= nrow(data))  # ADDED CODE
    MX = MX[include, ]  # ADDED CODE
    MY = MY[include, ]  # ADDED CODE
    # IDEA: JUST IN CASE MODEL IS TOO BIG (MCMC DOESN'T ALLOW p > n)
    tmp = calc.ACE(population, exposure, outcome, predictorsY, 
        interactions, confounders, MY, familyY, thin, burnB, 
        data)
    ACE.boot = tmp$ACE.boot
    ptime = round(proc.time()[3] - ptime, digits = 2)
    print(paste("It took you ", ptime, " seconds to run BAC", 
        sep = ""))
    para = list(familyX = familyX, familyY = familyY, thin = thin, 
        omega = omega, num_its = num_its, burnM = burnM, burnB = burnB, 
        population = population)
    models = find.models(MY)
    result = list(data = data, MX = MX, MY = MY, models = models, 
        exposure = exposure, outcome = outcome, confounders = confounders, 
        interactions = interactions, predictorsY = predictorsY, 
        ACE = ACE.boot, para = para)
    class(result) = "bacr"
    return(result)
}
# END SCRIPT
