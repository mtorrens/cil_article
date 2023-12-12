################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(FCNDIR, 'functions_plot.R', sep = ''))
################################################################################

################################################################################
# List of functions
################################################################################
# * decouple.var.mt
# * adapt.output.mt
# * plot.compare.st
# * plot.mtrand
################################################################################

################################################################################
# For multiple treatment outputs, it decouples object into T<=5 columns
decouple.var.mt <- function(res, col) {
################################################################################
  out1 <- out2 <- out3 <- out4 <- out5 <- rep(NA, nrow(res))
  for (i in 1:nrow(res)) {
    out1[i] <- as.character(res[i, col][[1]][1])
    if (length(res[i, col][[1]]) > 1) {
      out2[i] <- as.character(res[i, col][[1]][2])
    }
    if (length(res[i, col][[1]]) >= 3) {
      out3[i] <- as.character(res[i, col][[1]][3])
    }
    if (length(res[i, col][[1]]) >= 4) {
      out4[i] <- as.character(res[i, col][[1]][4])
    }
    if (length(res[i, col][[1]]) >= 5) {
      out5[i] <- as.character(res[i, col][[1]][5])
    }
  }
  out1 <- as.numeric(out1)
  out2 <- as.numeric(out2)
  out3 <- as.numeric(out3)
  out4 <- as.numeric(out4)
  out5 <- as.numeric(out5)
  return(list(out1 = out1, out2 = out2, out3 = out3, out4 = out4, out5 = out5))
}

################################################################################
# For MT, it adapts simulator output object into T<=5 matrices with results
adapt.output.mt <- function(res) {
################################################################################
  aux01 <- decouple.var.mt(res, 'ah.orc')
  aux02 <- decouple.var.mt(res, 'ah.ore')
  aux03 <- decouple.var.mt(res, 'ah.nmA')
  aux04 <- decouple.var.mt(res, 'ah.n2A')
  aux05 <- decouple.var.mt(res, 'ah.bm1')
  aux06 <- decouple.var.mt(res, 'ah.bm2')
  aux07 <- decouple.var.mt(res, 'ah.nle')
  aux08 <- decouple.var.mt(res, 'ah.dle')
  aux09 <- decouple.var.mt(res, 'ah.pcr')
  aux10 <- decouple.var.mt(res, 'ah.bcI')
  aux11 <- decouple.var.mt(res, 'ah.ssl')

  cua <- c('orc', 'ore', 'nle', 'dle', 'pcr', 'bcI',
           'ssl', 'bm1', 'bm2', 'nmA', 'n2A')
  resD1 <- as.data.frame(matrix(nrow = nrow(res), ncol = 12))
  colnames(resD1) <- c(colnames(res)[1], paste('ah', cua, sep = '.'))
  #colnames(resD1) <- c('nc', paste('ah', cua, sep = '.'))
  #resD1[, 1] <- unname(unlist(res[, 'nc']))
  resD1[, 1] <- unname(unlist(res[, 1]))
  resD1[, 2] <- aux01[[1]]
  resD1[, 3] <- aux02[[1]]
  resD1[, 4] <- aux07[[1]]
  resD1[, 5] <- aux08[[1]]
  resD1[, 6] <- aux09[[1]]
  resD1[, 7] <- aux10[[1]]
  resD1[, 8] <- aux11[[1]]
  resD1[, 9] <- aux05[[1]]
  resD1[, 10] <- aux06[[1]]
  resD1[, 11] <- aux03[[1]]
  resD1[, 12] <- aux04[[1]]

  resD2 <- as.data.frame(matrix(nrow = nrow(res), ncol = 12))
  colnames(resD2) <- colnames(resD1)
  #resD2[, 1] <- unname(unlist(res[, 'nc']))
  resD2[, 1] <- unname(unlist(res[, 1]))
  resD2[, 2] <- aux01[[2]]
  resD2[, 3] <- aux02[[2]]
  resD2[, 4] <- aux07[[2]]
  resD2[, 5] <- aux08[[2]]
  resD2[, 6] <- aux09[[2]]
  resD2[, 7] <- aux10[[2]]
  resD2[, 8] <- aux11[[2]]
  resD2[, 9] <- aux05[[2]]
  resD2[, 10] <- aux06[[2]]
  resD2[, 11] <- aux03[[2]]
  resD2[, 12] <- aux04[[2]]

  resD3 <- as.data.frame(matrix(nrow = nrow(res), ncol = 12))
  colnames(resD3) <- colnames(resD1)
  if (! is.na(aux01[[3]][1])) {
    resD3[, 1] <- unname(unlist(res[, 1]))
    resD3[, 2] <- aux01[[3]]
    resD3[, 3] <- aux02[[3]]
    resD3[, 4] <- aux07[[3]]
    resD3[, 5] <- aux08[[3]]
    resD3[, 6] <- aux09[[3]]
    resD3[, 7] <- aux10[[3]]
    resD3[, 8] <- aux11[[3]]
    resD3[, 9] <- aux05[[3]]
    resD3[, 10] <- aux06[[3]]
    resD3[, 11] <- aux03[[3]]
    resD3[, 12] <- aux04[[3]]
  }

  resD4 <- as.data.frame(matrix(nrow = nrow(res), ncol = 12))
  colnames(resD4) <- colnames(resD1)
  if (! is.na(aux01[[4]][1])) {
    resD4[, 1] <- unname(unlist(res[, 1]))
    resD4[, 2] <- aux01[[4]]
    resD4[, 3] <- aux02[[4]]
    resD4[, 4] <- aux07[[4]]
    resD4[, 5] <- aux08[[4]]
    resD4[, 6] <- aux09[[4]]
    resD4[, 7] <- aux10[[4]]
    resD4[, 8] <- aux11[[4]]
    resD4[, 9] <- aux05[[4]]
    resD4[, 10] <- aux06[[4]]
    resD4[, 11] <- aux03[[4]]
    resD4[, 12] <- aux04[[4]]
  }

  resD5 <- as.data.frame(matrix(nrow = nrow(res), ncol = 12))
  colnames(resD5) <- colnames(resD1)
  if (! is.na(aux01[[5]][1])) {
    resD5[, 1] <- unname(unlist(res[, 1]))
    resD5[, 2] <- aux01[[5]]
    resD5[, 3] <- aux02[[5]]
    resD5[, 4] <- aux07[[5]]
    resD5[, 5] <- aux08[[5]]
    resD5[, 6] <- aux09[[5]]
    resD5[, 7] <- aux10[[5]]
    resD5[, 8] <- aux11[[5]]
    resD5[, 9] <- aux05[[5]]
    resD5[, 10] <- aux06[[5]]
    resD5[, 11] <- aux03[[5]]
    resD5[, 12] <- aux04[[5]]
  }

  return(list(resD1, resD2, resD3, resD4, resD5))
}

################################################################################
# Plot for FIGURE 1, FIGURE S4, FIGURE S5 and FIGURE S6: compare simulations
plot.compare.st <- function(sims, alpha, qty = 'rmse', xtags = 0:6, titol = '',
  cex.main = 1, yrange = NULL, legpos = NULL, add.legend = TRUE,
  var.bac = FALSE, mt = FALSE, show.ssl = FALSE, show.pcr = FALSE,
  show.lasso = TRUE, show.EB = FALSE, legend.EB = FALSE) {
################################################################################
  # Columns of interest and general plotting parameter defaults
  mth <- c('orc', 'ore', 'nle', 'dle', 'pcr', 'bcI', 'acm', 'bm2', 'nmA', 'n2A', 'ssl')
  leg.pos <- 'topleft'
  xtags <- sort(unique(sims[, 'nc']))

  # Tweak according to what is being plotted
  if (qty == 'rmse') {
    qty.f <- function(x, a = alpha) { sqrt(mean((x - a) ** 2, na.rm = TRUE)) }
    cols <- paste('ah', mth, sep = '.')
    ylims <- c(0, 5)
  } else if (qty == 'ncovs') {
    qty.f <- function(x) { mean(x, na.rm = TRUE) }
    cols <- paste('nr', mth, sep = '.')
    ylims <- c(0, 20)
  } else if (qty == 'pincl_treat') {
    qty.f <- function(x) { mean(x, na.rm = TRUE) }
    cols <- paste('pd', mth, sep = '.')
    bayes <- which(mth %in% c('bm2', 'nmA', 'n2A'))
    cols[bayes] <- paste('t2', c(mth[bayes[1:2]], 'nmB'), sep = '.')
    ylims <- c(0.5, 1)
    leg.pos <- 'bottomleft'
  }

  # Override choice of y-axis limits
  if (! is.null(yrange)) { ylims <- yrange }
  if (! is.null(legpos)) { leg.pos <- legpos }

  # Compute quantity
  qty.orc <- tapply(sims[, cols[1]], sims[, 'nc'], qty.f)
  if (qty == 'pincl_treat') { qty.orc <- 1 }
  rqty.ore <- tapply(sims[, cols[2]], sims[, 'nc'], qty.f) / qty.orc
  rqty.nle <- tapply(sims[, cols[3]], sims[, 'nc'], qty.f) / qty.orc
  rqty.dle <- tapply(sims[, cols[4]], sims[, 'nc'], qty.f) / qty.orc
  rqty.pcr <- tapply(sims[, cols[5]], sims[, 'nc'], qty.f) / qty.orc
  rqty.bcI <- tapply(sims[, cols[6]], sims[, 'nc'], qty.f) / qty.orc
  rqty.acm <- tapply(sims[, cols[7]], sims[, 'nc'], qty.f) / qty.orc
  rqty.bm2 <- tapply(sims[, cols[8]], sims[, 'nc'], qty.f) / qty.orc
  rqty.nmA <- tapply(sims[, cols[9]], sims[, 'nc'], qty.f) / qty.orc
  rqty.n2A <- tapply(sims[, cols[10]], sims[, 'nc'], qty.f) / qty.orc
  rqty.ssl <- tapply(sims[, cols[11]], sims[, 'nc'], qty.f) / qty.orc
  rqty.orc <- tapply(sims[, cols[1]], sims[, 'nc'], qty.f) / qty.orc
  if (var.bac == TRUE) { rqty.bcI <- jitter(rqty.bcI, factor = 2.5) }

  # Plot preparing
  pchs <- c(3, 4, 5, 2, 6, 1, 4, 8)#, 1)
  ltys <- c(2, 1, 1, 1, 2, 1, 1, 1)#, 1)
  labs <- c('LASSO', 'DL', 'PCR', 'BAC', 'BMA',
            'CIL (EP)', 'CIL (EB)', 'EM-SSL')#, 'Ext. Oracle')
  # labs <- c('LASSO (PSI)', 'DML (LASSO)', 'PCR (BIC)', 'BAC (Inf)', 'BMA (NLP)',
  #           'CIL (EP)', 'CIL (EB)', 'EM-SSL')#, 'Ext. Oracle')
  colors <- c(rep('azure4', 5), rep('black', 2), 'azure4')#, 'azure4')
  if (qty == 'rmse') { labs[1] <- 'LASSO' }
  if (mt == TRUE) { labs[4] <- 'ACPME' }

  # Plot
  plot(NA, type = 'l', xaxt = 'n', col = 'white', xlab = '', ylab = '',
    main = titol, xlim = c(1, 7), ylim = ylims, cex.main = cex.main)
  axis(1, 1:7, xtags); grid()
  if (qty == 'pincl_treat') { 
    lines(rqty.orc, col = 'gray')
  } else {
    abline(h = 1, col = 'gray')
  }
  if (show.lasso == TRUE) {
    lines(rqty.nle, col = colors[1], lty = ltys[1], type = 'b', pch = pchs[1])
  }
  lines(rqty.dle, col = colors[2], lty = ltys[2], type = 'b', pch = pchs[2])
  lines(rqty.bcI, col = colors[4], lty = ltys[4], type = 'b', pch = pchs[4])
  lines(rqty.bm2, col = colors[5], lty = ltys[5], type = 'b', pch = pchs[5])
  lines(rqty.nmA, col = colors[6], lty = ltys[6], type = 'b', pch = pchs[6])
  if (show.EB == TRUE) {
    lines(rqty.n2A, col = colors[7], lty = ltys[7], type = 'b', pch = pchs[7])
  }
  if (show.pcr == TRUE) {
    lines(rqty.pcr, col = colors[3], lty = ltys[3], type = 'b', pch = pchs[3])
  }
  if (show.ssl == TRUE) {
    lines(rqty.ssl, col = colors[8], lty = ltys[8], type = 'b', pch = pchs[8])
  }
  #lines(rqty.ore, col = colors[9], lty = ltys[8], type = 'b', pch = pchs[9])

  # Add ACPME for SINGLE TREATMENT
  if (mt == FALSE) {
    lines(rqty.acm, col = 'azure4', lty = 1, type = 'b', pch = 0)
  }

  # Legend
  if (add.legend == TRUE) {
    ordre <- c(1:5, 8, 6:7)
    # if (mt == TRUE) {
    if (mt == TRUE | show.ssl == FALSE) {
      #kill <- c(4, 8)
      kill <- 8
      labs <- labs[-kill]
      pchs <- pchs[-kill]
      ltys <- ltys[-kill]
      colors <- colors[-kill]
      ordre <- 1:length(labs)
    }
    # if (show.ssl == FALSE) {
    #   kill <- 8
    #   labs <- labs[-kill]
    #   pchs <- pchs[-kill]
    #   ltys <- ltys[-kill]
    #   colors <- colors[-kill]
    #   ordre <- 1:length(labs)
    # }
    kill <- NULL
    if (show.lasso == FALSE) { kill <- c(kill, 1) }
    if (show.pcr == FALSE) { kill <- c(kill, 3) }
    if (show.EB == FALSE & legend.EB == FALSE) { kill <- c(kill, 7) }
    labs <- labs[-kill]
    pchs <- pchs[-kill]
    ltys <- ltys[-kill]
    colors <- colors[-kill]
    ordre <- 1:length(labs)

    # Drop the EP part of the label if there is only one CIL version
    if (show.EB == FALSE & legend.EB == FALSE) {
      labs[which(labs == 'CIL (EP)')] <- 'CIL'
    }

    # Add ACPME layer for SINGLE TREATMENT after all
    if (mt == TRUE) {
      new.labs <- labs[ordre]
      new.colors <- colors[ordre]
      new.ltys <- ltys[ordre]
      new.pchs <- pchs[ordre]
    } else {
      idx <- length(labs) # CIL index
      new.labs <- c(labs[ordre][-idx], 'ACPME', labs[ordre][idx])
      new.colors <- c(colors[ordre][-idx], 'azure4', colors[ordre][idx])
      new.ltys <- c(ltys[ordre][-idx], 1, ltys[ordre][idx])
      new.pchs <- c(pchs[ordre][-idx], 0, pchs[ordre][idx])
    }

    legend(leg.pos, new.labs, col = new.colors, lty = new.ltys,
      pch = new.pchs, bg = 'white', ncol = 2, merge = FALSE, cex = 0.8)
  }
}

################################################################################
# Plot for FIGURE 5: compare simulations (MULTIPLE TREATMENT version)
plot.mtrand <- function(x1, x2, x3, x4, alpha, qty = 'rmse', yrange = NULL,
  add.legend = TRUE, legpos = 'topleft', cex.main = 1, titol = '', xnom = '',
  ynom = '', show.EB = FALSE, show.pcr = FALSE, which.bma = 1,
  var.bac = FALSE) {
################################################################################
  if (is.null(x4)) { x4 <- as.data.frame(matrix(NA, nrow = 0, ncol = ncol(x1))) }
  nts <- c(rep(2, nrow(x1)), rep(3, nrow(x2)), rep(4, nrow(x3)), rep(5, nrow(x4)))
  res <- cbind(nts, rbind(x1, x2, x3, x4))

  # Columns of interest and general plotting parameter defaults
  if (which.bma == 2) {
    mth <- c('orc', 'ore', 'nle', 'dle', 'pcr', 'bcI', 'bm2', 'nmA', 'n2A')
  } else if (which.bma == 1) {
    mth <- c('orc', 'ore', 'nle', 'dle', 'pcr', 'bcI', 'bm1', 'nmA', 'n2A')
  } else {
    stop('arg. "which.bma" must be either "1" or "2" (numeric).')
  }
  leg.pos <- 'topleft'
  xtags <- as.character(1 + 1:length(unique(nts)))

  # Tweak according to what is being plotted
  if (qty == 'rmse') {
    qty.f <- function(x, a = alpha) { sqrt(mean((x - a) ** 2, na.rm = TRUE)) }
    cols <- paste('ah', mth, sep = '.')
    ylims <- c(0.5, 2.5)
  # } else if (qty == 'ncovs') {
  #   qty.f <- function(x) { mean(x, na.rm = TRUE) }
  #   cols <- paste('nr', mth, sep = '.')
  #   ylims <- c(0, 20)
  # } else if (qty == 'pincl_treat') {
  #   qty.f <- function(x) { mean(x, na.rm = TRUE) }
  #   cols <- paste('pd', mth, sep = '.')
  #   ylims <- c(0.5, 1)
  #   leg.pos <- 'bottomleft'
  }

  # Override choice of y-axis limits
  if (! is.null(yrange)) { ylims <- yrange }
  if (! is.null(legpos)) { leg.pos <- legpos }

  # Compute quantity
  qty.orc <- tapply(res[, cols[1]], res[, 'nts'], qty.f)
  if (qty == 'pincl_treat') { qty.orc <- 1 }
  rqty.ore <- tapply(res[, cols[2]], res[, 'nts'], qty.f) / qty.orc
  rqty.nle <- tapply(res[, cols[3]], res[, 'nts'], qty.f) / qty.orc
  rqty.dle <- tapply(res[, cols[4]], res[, 'nts'], qty.f) / qty.orc
  rqty.pcr <- tapply(res[, cols[5]], res[, 'nts'], qty.f) / qty.orc
  rqty.bcI <- tapply(res[, cols[6]], res[, 'nts'], qty.f) / qty.orc
  rqty.bm2 <- tapply(res[, cols[7]], res[, 'nts'], qty.f) / qty.orc
  rqty.nmA <- tapply(res[, cols[8]], res[, 'nts'], qty.f) / qty.orc
  rqty.n2A <- tapply(res[, cols[9]], res[, 'nts'], qty.f) / qty.orc
  rqty.orc <- tapply(res[, cols[1]], res[, 'nts'], qty.f) / qty.orc
  if (var.bac == TRUE) { rqty.bcI <- jitter(rqty.bcI, factor = 2.5) }

  # Plot preparing
  pchs <- c(3, 4, 5, 2, 6, 1, 4, 8)#, 1)
  ltys <- c(2, 1, 1, 1, 2, 1, 1, 1)#, 1)
  labs <- c('LASSO', 'DML', 'PCR', 'ACPME', 'BMA',
            'CIL', 'CIL (EB)', 'EM-SSL')#, 'Ext. Oracle')
  colors <- c(rep('azure4', 5), rep('black', 2), 'azure4')#, 'azure4')

  # Plot
  plot(NA, type = 'l', xaxt = 'n', col = 'white', xlab = xnom, ylab = ynom,
    xlim = c(0.75, length(xtags) + 0.25), ylim = ylims, main = titol,
    cex.main = cex.main)
  axis(1, 1:length(xtags), xtags); grid()
  if (qty == 'pincl_treat') { 
    lines(rqty.orc, col = 'gray')
  } else {
    abline(h = 1, col = 'gray')
  }

  # Lines
  lines(rqty.nle, col = colors[1], lty = ltys[1], type = 'b', pch = pchs[1])
  lines(rqty.dle, col = colors[2], lty = ltys[2], type = 'b', pch = pchs[2])
  lines(rqty.bcI, col = colors[4], lty = ltys[4], type = 'b', pch = pchs[4])
  lines(rqty.bm2, col = colors[5], lty = ltys[5], type = 'b', pch = pchs[5])
  lines(rqty.nmA, col = colors[6], lty = ltys[6], type = 'b', pch = pchs[6])
  if (show.EB == TRUE) {
    lines(rqty.n2A, col = colors[7], lty = ltys[7], type = 'b', pch = pchs[7])
  }
  if (show.pcr == TRUE) {
    lines(rqty.pcr, col = colors[3], lty = ltys[3], type = 'b', pch = pchs[3])
  }

  # Legend
  if (add.legend == TRUE) {
    ordre <- c(1:5, 8, 6:7)
    kill <- 8
    if (show.pcr == FALSE) { kill <- c(kill, 3) }
    if (show.EB == FALSE) { kill <- c(kill, 7) }
    labs <- labs[-kill]
    pchs <- pchs[-kill]
    ltys <- ltys[-kill]
    colors <- colors[-kill]
    ordre <- 1:length(labs)

    legend(leg.pos, labs[ordre], col = colors[ordre], lty = ltys[ordre],
      pch = pchs[ordre], bg = 'white', ncol = 2, merge = FALSE, cex = 0.8)
  }
}
# END OF SCRIPT
