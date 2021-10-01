################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '06_cps_figures.R', sep = ''))
################################################################################

################################################################################
# Required data files with computed results
################################################################################
# Y2Y - State ints - NO fake predictors - Model
f01a <- paste(DATDIR, 'cps_stateints_y2y_post_FULL.RData', sep = '')
# Y2Y - State ints - NO fake predictors - Post. samples
f01b <- paste(DATDIR, 'cps_stateints_postsamples_FULL.RData', sep = '')
# Y2Y - State ints - 100 fake predictors - Model
f02a <- paste(DATDIR, 'cps_stateintsFP_y2y_FULL.RData', sep = '')
# Y2Y - State ints - 200 fake predictors - Model
f03a <- paste(DATDIR, 'cps_stateintsFP200_y2y_FULL.RData', sep = '')
################################################################################

################################################################################
# Load all data
################################################################################
y2y.bmkP <- get(load(file = f01a)); cat('Loaded file:', f01a, '\n')
y2y.1hRP <- get(load(file = f02a)); cat('Loaded file:', f02a, '\n')
y2y.2hRP <- get(load(file = f03a)); cat('Loaded file:', f03a, '\n')
ps.bmk <- get(load(file = f01b)); cat('Loaded file:', f01b, '\n')

# Extract relevant objects control
y2y.bmkP <- lapply(y2y.bmkP, '[[', 'out')
y2y.1hRP <- lapply(y2y.1hRP, '[[', 'out')
y2y.2hRP <- lapply(y2y.2hRP, '[[', 'out')

# State labels
f00 <- paste(DATDIR, 'cps/cps_model_data_year.RData', sep = '')
raw <- load(file = f00); cat('Loaded file:', f00, '\n')
states <- as.character(unique(get(raw[4])))
years <- as.numeric(as.character(unique(get(raw[5]))))
rm(f00, f01a, f01b, f02a, f03a)
################################################################################

################################################################################
# Rearrange results
################################################################################
pe.bmk <- pe.1hR <- pe.2hR <- list()
for (v in 1:length(years)) {
  w <- ifelse(v == 1, 1, 2)  # Augmented data only in 2010 and 2019 (diff. vers.)
  pe.bmk[[w]] <- y2y.bmkP[[v]][[1]]
  v <- ifelse(v == 1, 1, 2)
  pe.1hR[[w]] <- y2y.1hRP[[v]][[1]]
  pe.2hR[[w]] <- y2y.2hRP[[v]][[1]]
}
pe.bmk <- as.data.frame(do.call(rbind, pe.bmk))
pe.1hR <- as.data.frame(do.call(rbind, pe.1hR))
pe.2hR <- as.data.frame(do.call(rbind, pe.2hR))
pe.bmk <- apply(pe.bmk, 2, as.numeric)
pe.1hR <- apply(pe.1hR, 2, as.numeric)
pe.2hR <- apply(pe.2hR, 2, as.numeric)

# Define relevant posterior measures
m5.nmP10 <- ps.bmk[[1]][['m5.nmP']]
m5.nmP19 <- ps.bmk[[2]][['m5.nmP']]
m7.nmP10 <- ps.bmk[[1]][['m5.nmP2']]  # State by state
m7.nmP19 <- ps.bmk[[2]][['m5.nmP2']]

# Object containing inference for all expansions
obj <- list(pe.bmk, pe.1hR, pe.2hR)
################################################################################

################################################################################
# Plots
################################################################################
# Directory
if (! dir.exists(paste(FIGDIR, 'cps/', sep = ''))) {
  dir.create(paste(FIGDIR, 'cps/', sep = ''))
}

################################################################################
# FIGURE 2
################################################################################
# 1(A). Paper panels
figT <- paste(FIGDIR, 'cps/figApp1top.pdf', sep = '')
figB <- paste(FIGDIR, 'cps/figApp1bot.pdf', sep = '')
pdf(figT, height = 4.5, width = 5.5)  # TOP PANEL
par(mfrow = c(1, 1), bg = 'white')
plot.2010vs2019(obj, treat = 1, ymin = -0.270, ymax = +0.065,
  xlab.in = TRUE, year.tag = TRUE)
dev.off(); cat('Printed figure:', figT, '\n')
pdf(figB, height = 4.5, width = 5.5)  # BOTTOM PANEL
par(mfrow = c(1, 1), bg = 'white')
plot.2010vs2019(obj, treat = 2, ymin = -0.095, ymax = +0.075,
  xlab.in = FALSE, year.tag = FALSE)
dev.off(); cat('Printed figure:', figB, '\n')

################################################################################
# FIGURE S3
################################################################################
# 1(B). Supplement panels
figT <- paste(FIGDIR, 'cps/figApp1Btop.pdf', sep = '')
figB <- paste(FIGDIR, 'cps/figApp1Bbot.pdf', sep = '')
pdf(figT, height = 4.5, width = 5.5)  # TOP PANEL
par(mfrow = c(1, 1), bg = 'white')
plot.2010vs2019(obj, treat = 3, ymin = -0.130, ymax = +0.130,
  xlab.in = TRUE, year.tag = TRUE)
dev.off(); cat('Printed figure:', figT, '\n')
pdf(figB, height = 4.5, width = 5.5)  # BOTTOM PANEL
par(mfrow = c(1, 1), bg = 'white')
plot.2010vs2019(obj, treat = 4, ymin = -0.183, ymax = +0.183,
  xlab.in = FALSE, year.tag = FALSE)
dev.off(); cat('Printed figure:', figB, '\n')

################################################################################
# FIGURE 4
################################################################################
# [LEFT PANEL] 2(A). Posterior distribution of exp{|h|}
figL <- paste(FIGDIR, 'cps/figApp2.pdf', sep = '')
pdf(figL, height = 6, width = 5)  # LEFT PANEL
par(mfrow = c(1, 1), mar = c(2.5, 2.5, 0.5, 0.5), bg = 'white')
min.twoway.boxplot(m5.nmP10, m5.nmP19, cex.axis = 1.5,
  ylab = '', xlab = c('2010', '2019'))
dev.off(); cat('Printed figure:', figL, '\n')

# 2(A). Values present in the text
round(summary(m5.nmP10), 3)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.009   1.032   1.053   1.079   1.300 
round(quantile(m5.nmP10, probs = c(0.025, 0.05, 0.95, 0.975)), 3)
#  2.5%    5%   95% 97.5% 
# 1.001 1.001 1.178 1.204 
round(summary(m5.nmP19), 3)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.003   1.009   1.015   1.022   1.079 
round(quantile(m5.nmP19, probs = c(0.025, 0.05, 0.95, 0.975)), 3)
#  2.5%    5%   95% 97.5% 
# 1.000 1.000 1.048 1.055 

################################################################################
# [RIGHT PANEL] 2(B). State-by-state figure
a10 <- unlist(lapply(m7.nmP10, mean))
a19 <- unlist(lapply(m7.nmP19, mean))
m10 <- unlist(lapply(m7.nmP10, median))
m19 <- unlist(lapply(m7.nmP19, median))
q10 <- do.call(rbind, lapply(m7.nmP10, quantile, probs = c(0.25, 0.75)))
q19 <- do.call(rbind, lapply(m7.nmP19, quantile, probs = c(0.25, 0.75)))
names(a10) <- names(a19) <- names(m10) <- names(m19) <- states
rownames(q10) <- rownames(q19) <- states

# Plot (manual)
o <- order(m19, decreasing = TRUE)
s1 <- seq(1, 51, 2)
s2 <- seq(2, 51, 2)
figR <- paste(FIGDIR, 'cps/figApp2B.pdf', sep = '')
pdf(figR, height = 6, width = 5)
par(mfrow = c(1, 1), mar = c(2, 2.5, 0.5, 2.5), las = 1, bg = 'white')
plot(rep(NA, 51), 1:51, type = 'l', yaxt = 'n', cex.axis = 1.1,
  xlim = c(1, 1.09), ylim = c(0, 51), ylab = '', xlab = '')
for (i in seq(1.01, 1.1, 0.01)) { abline(v = i, col = 'lightgray', lty = 3) }
abline(v = 1, lty = 2)
for (i in 1:51) {
  abline(h = i, col = 'lightgray', lty = 3)
  segments(x0 = q10[o, ][i, 1], x1 = q10[o, ][i, 2],
    y0 = i + 0.1, y1 = i + 0.1, col = 'darkgray')
  segments(x0 = q19[o, ][i, 1], x1 = q19[o, ][i, 2],
    y0 = i - 0.1, y1 = i - 0.1)
}
points(m10[o], 1:51 + 0.1, pch = 21, cex = 1, bg = 'darkgray')
points(m19[o], 1:51 - 0.1, pch = 16, cex = 1)
axis(2, s1, names(m10[o])[s1], cex.axis = 1, family = 'mono')
axis(4, s2, names(m10[o])[s2], cex.axis = 1, family = 'mono')
text(1.025, -0.5, '2019', cex = 1.1)
text(1.062, -0.5, '2010', cex = 1.1)
dev.off(); cat('Printed figure:', figR, '\n')

# 2(B). Values present in the text
head(round(sort(m10 - m19), 4))
#     ME     WY     SD     ND     KY     MT 
# 0.0019 0.0021 0.0026 0.0044 0.0056 0.0063 
tail(round(sort(m10 - m19), 4))
#     MD     CA     NY     TX     OH     IN 
# 0.0250 0.0260 0.0263 0.0326 0.0402 0.0410 
head(round(sort((m19 - m10) / m10), 4))
#      IN      OH      TX      NY      CA      MD 
# -0.0389 -0.0382 -0.0313 -0.0254 -0.0251 -0.0242 
tail(round(sort((m19 - m10) / m10), 4))
#      MT      KY      ND      SD      WY      ME 
# -0.0063 -0.0056 -0.0044 -0.0025 -0.0021 -0.0019 
# END OF SCRIPT
