################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '02_paper_figures.R', sep = ''))
################################################################################

################################################################################
# SETUPS
# (P: Pending; F: Finished; U: Unfinished; R: Running)
################################################################################
# (F) 00: [ST]  n =  50; p+T =  25; a = 1; p_y = p_d = 6
file00 <- paste(DATDIR, 'psim_T_ST50x24.RData', sep = '')
# (F) 01: [ST]  n = 100; p+T =  50; a = 1; p_y = p_d = 6
file01 <- paste(DATDIR, 'psim_T_ST100x49.RData', sep = '')
# (F) 02: [ST]  n = 100; p+T =  50; a = 1/3; p_y = p_d = 6
file02 <- paste(DATDIR, 'psim_T_ST100x49a03.RData', sep = '')
# (F) 03: [ST]  n = 100; p+T =  50; a = 0; p_y = p_d = 6
file03 <- paste(DATDIR, 'psim_T_ST100x49a0.RData', sep = '')
# (F) 04: [ST]  n = 100; p+T = 100; a = 1; p_y = p_d = 6
file04 <- paste(DATDIR, 'psim_T_ST100x99.RData', sep = '')
# (F) 05: [ST]  n = 100; p+T = 200; a = 1; p_y = p_d = 6
file05 <- paste(DATDIR, 'psim_T_ST100x199.RData', sep = '')
# (F) 06: [ST]  n = 100; p+T = 100; a = 1; p_y = p_d = 18
file06 <- paste(DATDIR, 'psim_T_ST100x99p18.RData', sep = '')
# (F) 07: [ST]  n = 100; p+T = 100; a = 1; p_y = p_d = 12
file07 <- paste(DATDIR, 'psim_T_ST100x99p12.RData', sep = '')
# (F) 08: [MT]  n = 100; p+T = 100; a = 1,1; p_y = p_d = 18
file08 <- paste(DATDIR, 'psim_T_MT100x98.RData', sep = '')
# (F) 09: [MT]  n = 100; p+T = 100; a = 1,1; p_y = p_d = 12
file09 <- paste(DATDIR, 'psim_T_MT100x98p12.RData', sep = '')
# (F) 10: [MT]  n = 100; p+T = 100; a = 1/3,0; p_y = p_d = 18
file10 <- paste(DATDIR, 'psim_T_MT100x98smtr.RData', sep = '')
# (F) 11: [MT]  n = 100; p+T = 100; a = 1/3,0; p_y = p_d = 12
file11 <- paste(DATDIR, 'psim_T_MT100x98p12smtr.RData', sep = '')
# (F) 12: [MT]  n = 100; p+T =  50; a = 1,1; p_y = p_d = 6
file12 <- paste(DATDIR, 'psim_T_MT100x48p6.RData', sep = '')
# (F) 13: [MTR] n = 100; p+T = 100; a = b1; nts = 2:5
file13 <- paste(DATDIR, 'psim_T_MTR100x95.RData', sep = '')
################################################################################

################################################################################
# Load all data
res00 <- get(load(file = file00))[[1]]; cat('Loaded file:', file00, '\n')
res01 <- get(load(file = file01))[[1]]; cat('Loaded file:', file01, '\n')
res02 <- get(load(file = file02))[[1]]; cat('Loaded file:', file02, '\n')
res03 <- get(load(file = file03))[[1]]; cat('Loaded file:', file03, '\n')
res04 <- get(load(file = file04))[[1]]; cat('Loaded file:', file04, '\n')
res05 <- get(load(file = file05))[[1]]; cat('Loaded file:', file05, '\n')
res06 <- get(load(file = file06))[[1]]; cat('Loaded file:', file06, '\n')
res07 <- get(load(file = file07))[[1]]; cat('Loaded file:', file07, '\n')
res08 <- get(load(file = file08))[[1]]; cat('Loaded file:', file08, '\n')
res09 <- get(load(file = file09))[[1]]; cat('Loaded file:', file09, '\n')
res10 <- get(load(file = file10))[[1]]; cat('Loaded file:', file10, '\n')
res11 <- get(load(file = file11))[[1]]; cat('Loaded file:', file11, '\n')
res12 <- get(load(file = file12))[[1]]; cat('Loaded file:', file12, '\n')
res13 <- get(load(file = file13))[[1]]; cat('Loaded file:', file13, '\n')
rm(list = paste('file', sprintf('%02.0f', 0:13), sep = ''))
################################################################################

################################################################################
# FIGURE 1: ST alpha strong/weak/null
################################################################################
# Figure file names
fig1Ap1 <- paste(FIGDIR, 'fig1Ap1.pdf', sep = '')
fig1Ap2 <- paste(FIGDIR, 'fig1Ap2.pdf', sep = '')
fig1Ap3 <- paste(FIGDIR, 'fig1Ap3.pdf', sep = '')
fig1Bp1 <- paste(FIGDIR, 'fig1Bp1.pdf', sep = '')
fig1Bp2 <- paste(FIGDIR, 'fig1Bp2.pdf', sep = '')
fig1Bp3 <- paste(FIGDIR, 'fig1Bp3.pdf', sep = '')
fig1Bp4 <- paste(FIGDIR, 'fig1Bp4.pdf', sep = '')
fig1Bp5 <- paste(FIGDIR, 'fig1Bp5.pdf', sep = '')
fig1Bp6 <- paste(FIGDIR, 'fig1Bp6.pdf', sep = '')

# First Part (individual panels)
pdf(fig1Ap1, height = 4, width = 3) # First panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res01, alpha = 1, qty = 'rmse', yrange = c(0, 5),
  cex.main = 1.5, add.legend = TRUE)
dev.off(); cat('Printed figure:', fig1Ap1, '\n')
pdf(fig1Ap2, height = 4, width = 3) # Second panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res02, alpha = 1/3, qty = 'rmse', yrange = c(0, 5),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Ap2, '\n')
pdf(fig1Ap3, height = 4, width = 3) # Third panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res03, alpha = 0, qty = 'rmse', yrange = c(0, 5),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Ap3, '\n')

# Second Part (individual panels)
pdf(fig1Bp1, height = 3.5, width = 3.5)  # First top panel
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res01, alpha = 1, qty = 'ncovs', yrange = c(0.75, 4.25),
  cex.main = 1.5, add.legend = TRUE)
dev.off(); cat('Printed figure:', fig1Bp1, '\n')
pdf(fig1Bp2, height = 3.5, width = 3.5)  # Second top panel
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res02, alpha = 1/3, qty = 'ncovs', yrange = c(0.75, 4.25),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Bp2, '\n')
pdf(fig1Bp3, height = 3.5, width = 3.5)  # Third top panel
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res03, alpha = 0, qty = 'ncovs', yrange = c(0.75, 4.25),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Bp3, '\n')
pdf(fig1Bp4, height = 3.5, width = 3.5)  # First bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res01, alpha = 1, qty = 'pincl_treat', yrange = c(0.6, 1),
  add.legend = FALSE, show.lasso = FALSE)
dev.off(); cat('Printed figure:', fig1Bp4, '\n')
pdf(fig1Bp5, height = 3.5, width = 3.5)  # Second bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res02, alpha = 1/3, qty = 'pincl_treat', yrange = c(0.6, 1),
  add.legend = FALSE, var.bac = TRUE, show.lasso = FALSE)
dev.off(); cat('Printed figure:', fig1Bp5, '\n')
pdf(fig1Bp6, height = 3.5, width = 3.5)  # Third bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res03, alpha = 0, qty = 'pincl_treat', yrange = c(0, 0.25),
  add.legend = FALSE, var.bac = TRUE, show.lasso = FALSE)
dev.off(); cat('Printed figure:', fig1Bp6, '\n')

# Erase file traces
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 2: ST changing p (growing design)
################################################################################
# Figure file names
fig2Ap1 <- paste(FIGDIR, 'fig2Ap1.pdf', sep = '')
fig2Ap2 <- paste(FIGDIR, 'fig2Ap2.pdf', sep = '')
fig2Ap3 <- paste(FIGDIR, 'fig2Ap3.pdf', sep = '')

# Plot per panels
pdf(fig2Ap1, height = 4, width = 3) # First panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res00, alpha = 1, qty = 'rmse', yrange = c(0.5, 5.5),
  cex.main = 1.35, add.legend = TRUE, legpos = 'topright', legend.EB = TRUE)
dev.off(); cat('Printed figure:', fig2Ap1, '\n')
pdf(fig2Ap2, height = 4, width = 3) # Second panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res04, alpha = 1, qty = 'rmse', yrange = c(0.5, 5.5),
  cex.main = 1.35, add.legend = FALSE)
dev.off(); cat('Printed figure:', fig2Ap2, '\n')
pdf(fig2Ap3, height = 4, width = 3) # Third panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res05, alpha = 1, qty = 'rmse', yrange = c(0.5, 5.5),
  cex.main = 1.35, show.EB = TRUE, add.legend = FALSE)
dev.off(); cat('Printed figure:', fig2Ap3, '\n')

# Erase file traces
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 3: ST p* -> Inf (increase in confounding)
################################################################################
# Figure file names
fig3Ap1 <- paste(FIGDIR, 'fig3Ap1.pdf', sep = '')
fig3Ap2 <- paste(FIGDIR, 'fig3Ap2.pdf', sep = '')
fig3Ap3 <- paste(FIGDIR, 'fig3Ap3.pdf', sep = '')

# Plot per panels
pcr.here <- FALSE  # Include PCR in the figure
pdf(fig3Ap1, height = 4, width = 3) # First panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res04, alpha = 1, qty = 'rmse', yrange = c(0, 12),
  cex.main = 1.5, add.legend = TRUE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig3Ap1, '\n')
pdf(fig3Ap2, height = 4, width = 3) # Second panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res07, alpha = 1, qty = 'rmse', yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig3Ap2, '\n')
pdf(fig3Ap3, height = 4, width = 3) # Third panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res06, alpha = 1, qty = 'rmse', yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig3Ap3, '\n')

# Erase file traces
rm(pcr.here)
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 4: MT p* growing (strong alpha, fixed p)
################################################################################
# Figure file names
fig4A <- paste(FIGDIR, 'fig4A.pdf', sep = '')
fig4Ap1 <- paste(FIGDIR, 'fig4Ap1.pdf', sep = '')
fig4Ap2 <- paste(FIGDIR, 'fig4Ap2.pdf', sep = '')
fig4Ap3 <- paste(FIGDIR, 'fig4Ap3.pdf', sep = '')
fig4Ap4 <- paste(FIGDIR, 'fig4Ap4.pdf', sep = '')

# Decouple outputs
aux08 <- adapt.output.mt(res08)
aux09 <- adapt.output.mt(res09)
res08D1 <- aux08[[1]]
res08D2 <- aux08[[2]]
res09D1 <- aux09[[1]]
res09D2 <- aux09[[2]]

# Plot per panels
pcr.here <- FALSE
pdf(fig4Ap1, height = 4, width = 4)  # First top panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res09D1, alpha = 1, qty = 'rmse', mt = TRUE, yrange = c(0, 12),
  cex.main = 1.5, add.legend = TRUE, legpos = 'topleft', show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig4Ap1, '\n')
pdf(fig4Ap2, height = 4, width = 4)  # Second top panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res08D1, alpha = 1, qty = 'rmse', mt = TRUE, yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig4Ap2, '\n')
pdf(fig4Ap3, height = 4, width = 4)  # First bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res09D2, alpha = 1, qty = 'rmse', mt = TRUE, yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig4Ap3, '\n')
pdf(fig4Ap4, height = 4, width = 4)  # Second bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res08D2, alpha = 1, qty = 'rmse', mt = TRUE, yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig4Ap4, '\n')

# Erase file traces
rm(pcr.here)
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 4 (B): MT p* = 6 (small design, strong alpha)
################################################################################
# Figure file names
fig4bA <- paste(FIGDIR, 'fig4bA.pdf', sep = '')
fig4bAp1 <- paste(FIGDIR, 'fig4bAp1.pdf', sep = '')
fig4bAp2 <- paste(FIGDIR, 'fig4bAp2.pdf', sep = '')

# Decouple outputs
aux12 <- adapt.output.mt(res12)
res12D1 <- aux12[[1]]
res12D2 <- aux12[[2]]

# Plot per panels
pdf(fig4bAp1, height = 4, width = 4)  # Left panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res12D1, alpha = 1, qty = 'rmse', mt = TRUE,
  yrange = c(0.5, 5.5), cex.main = 1.5, legpos = 'topleft')
dev.off(); cat('Printed figure:', fig4bAp1, '\n')
pdf(fig4bAp2, height = 4, width = 4)  # Right panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res12D2, alpha = 1, qty = 'rmse', mt = TRUE,
  yrange = c(0.5, 5.5), cex.main = 1.5, add.legend = FALSE)
dev.off(); cat('Printed figure:', fig4bAp2, '\n')

# Erase file traces
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 5: MT p* growing (weak alpha, fixed p)
################################################################################
# Figure file names
fig5A <- paste(FIGDIR, 'fig5A.pdf', sep = '')
fig5Ap1 <- paste(FIGDIR, 'fig5Ap1.pdf', sep = '')
fig5Ap2 <- paste(FIGDIR, 'fig5Ap2.pdf', sep = '')
fig5Ap3 <- paste(FIGDIR, 'fig5Ap3.pdf', sep = '')
fig5Ap4 <- paste(FIGDIR, 'fig5Ap4.pdf', sep = '')

# Decouple outputs
aux10 <- adapt.output.mt(res10)
aux11 <- adapt.output.mt(res11)
res10D1 <- aux10[[1]]
res10D2 <- aux10[[2]]
res11D1 <- aux11[[1]]
res11D2 <- aux11[[2]]

# Plot per panels
pcr.here <- FALSE
pdf(fig5Ap1, height = 4, width = 4)  # First top panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res11D1, alpha = 1/3, qty = 'rmse', mt = TRUE,
  yrange = c(0, 12), cex.main = 1.5, add.legend = TRUE, legpos = 'topleft',
  show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig5Ap1, '\n')
pdf(fig5Ap2, height = 4, width = 4)  # Second top panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res10D1, alpha = 1/3, qty = 'rmse', mt = TRUE,
  yrange = c(0, 12), cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig5Ap2, '\n')
pdf(fig5Ap3, height = 4, width = 4)  # First bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res11D2, alpha = 0, qty = 'rmse', mt = TRUE, yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig5Ap3, '\n')
pdf(fig5Ap4, height = 4, width = 4)  # Second bottom panel
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res10D2, alpha = 0, qty = 'rmse', mt = TRUE, yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig5Ap4, '\n')

# Erase file traces
rm(pcr.here)
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 6: MT with increasing number of treatments
################################################################################
figX <- paste(FIGDIR, 'figX.pdf', sep = '')
figXp1 <- paste(FIGDIR, 'figXp1.pdf', sep = '')
figXp2 <- paste(FIGDIR, 'figXp2.pdf', sep = '')

# Load data
term <- '_MTR100x95.RData'
f2 <- paste(DATDIR, 'psim_nt2', term, sep = '')
f3 <- paste(DATDIR, 'psim_nt3', term, sep = '')
f4 <- paste(DATDIR, 'psim_nt4', term, sep = '')
f5 <- paste(DATDIR, 'psim_nt5', term, sep = '')
res2 <- get(load(file = f2)); cat('Loaded file:', f2, '\n')
res3 <- get(load(file = f3)); cat('Loaded file:', f3, '\n')
res4 <- get(load(file = f4)); cat('Loaded file:', f4, '\n')
res5 <- get(load(file = f5)); cat('Loaded file:', f5, '\n')
ares2 <- adapt.output.mt(res2)
ares3 <- adapt.output.mt(res3)
ares4 <- adapt.output.mt(res4)
ares5 <- adapt.output.mt(res5)

# Assign results
res2T1 <- ares2[[1]]; res3T1 <- ares3[[1]]
res4T1 <- ares4[[1]]; res5T1 <- ares5[[1]]
res2T2 <- ares2[[2]]; res3T2 <- ares3[[2]]
res4T2 <- ares4[[2]]; res5T2 <- ares5[[2]]
res2T3 <- ares2[[3]]; res3T3 <- ares3[[3]]
res4T3 <- ares4[[3]]; res5T3 <- ares5[[3]]
res2T4 <- ares2[[4]]; res3T4 <- ares3[[4]]
res4T4 <- ares4[[4]]; res5T4 <- ares5[[4]]
res2T5 <- ares2[[5]]; res3T5 <- ares3[[5]]
res4T5 <- ares4[[5]]; res5T5 <- ares5[[5]]

x1 <- res2T1
x2 <- res3T1
x3 <- res4T1
x4 <- res5T1
y1 <- rbind(res2T1, res2T2)
y2 <- rbind(res3T1, res3T2, res3T3)
y3 <- rbind(res4T1, res4T2, res4T3, res4T4)
y4 <- rbind(res5T1, res5T2, res5T3, res5T4, res5T5)

# Plots
pcr.here <- FALSE
pdf(figXp1, height = 4, width = 4)
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.mtrand(x1, x2, x3, x4, alpha = 1, yrange = c(0.5, 7.5),
  cex.main = 1.5, add.legend = TRUE, legpos = 'topleft', show.pcr = pcr.here)
dev.off(); cat('Printed figure:', figXp1, '\n')

pdf(figXp2, height = 4, width = 4.25)
par(mfrow = c(1, 1), mar = c(4, 2, 0.5, 0.5), bg = 'white')
plot.mtrand(y1, y2, y3, y4, alpha = 1, yrange = c(0.5, 12),
  cex.main = 1.5, add.legend = TRUE, legpos = 'topleft', show.pcr = pcr.here,
  xnom = 'Number of Treatments')
dev.off(); cat('Printed figure:', figXp2, '\n')

# Erase file traces
rm(pcr.here)
rm(list = ls()[grep('fig', ls())])
# END OF SCRIPT
