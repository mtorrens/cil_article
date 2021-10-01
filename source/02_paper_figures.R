################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '02_paper_figures.R', sep = ''))
################################################################################

################################################################################
# SETUPS
# (P: Pending; D: Done; U: Unfinished; R: Running)
################################################################################
# (D) 00: [ST]  n =  50; p+T =  25; a = 1; p_y = p_d = 6
file00 <- paste(DATDIR, 'psim_T_ST50x24.RData', sep = '')
# (D) 01: [ST]  n = 100; p+T =  50; a = 1; p_y = p_d = 6
file01 <- paste(DATDIR, 'psim_T_ST100x49.RData', sep = '')
# (D) 02: [ST]  n = 100; p+T =  50; a = 1/3; p_y = p_d = 6
file02 <- paste(DATDIR, 'psim_T_ST100x49a03.RData', sep = '')
# (D) 03: [ST]  n = 100; p+T =  50; a = 0; p_y = p_d = 6
file03 <- paste(DATDIR, 'psim_T_ST100x49a0.RData', sep = '')
# (D) 04: [ST]  n = 100; p+T = 100; a = 1; p_y = p_d = 6
file04 <- paste(DATDIR, 'psim_T_ST100x99.RData', sep = '')
# (D) 05: [ST]  n = 100; p+T = 200; a = 1; p_y = p_d = 6
file05 <- paste(DATDIR, 'psim_T_ST100x199.RData', sep = '')
# (D) 06: [ST]  n = 100; p+T = 100; a = 1; p_y = p_d = 18
file06 <- paste(DATDIR, 'psim_T_ST100x99p18.RData', sep = '')
# (D) 07: [ST]  n = 100; p+T = 100; a = 1; p_y = p_d = 12
file07 <- paste(DATDIR, 'psim_T_ST100x99p12.RData', sep = '')
# (D) 08: [MT] n = 100; p+T = 100; a = 1s; nts = 2:5
file08a <- paste(DATDIR, 'psim_nt2_MTR100x95.RData', sep = '')
file08b <- paste(DATDIR, 'psim_nt3_MTR100x95.RData', sep = '')
file08c <- paste(DATDIR, 'psim_nt4_MTR100x95.RData', sep = '')
file08d <- paste(DATDIR, 'psim_nt5_MTR100x95.RData', sep = '')
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
res08a <- get(load(file = file08a)); cat('Loaded file:', file08a, '\n')
res08b <- get(load(file = file08b)); cat('Loaded file:', file08b, '\n')
res08c <- get(load(file = file08c)); cat('Loaded file:', file08c, '\n')
res08d <- get(load(file = file08d)); cat('Loaded file:', file08d, '\n')
rm(list = ls()[grep('^file', ls())])
################################################################################

################################################################################
# FIGURE 1 and S4: ST alpha strong/weak/null
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

# Figure 1 (individual panels)
pdf(fig1Ap1, height = 4, width = 3) # LEFT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res01, alpha = 1, qty = 'rmse', yrange = c(0, 5),
  cex.main = 1.5, add.legend = TRUE)
dev.off(); cat('Printed figure:', fig1Ap1, '\n')
pdf(fig1Ap2, height = 4, width = 3) # MIDDLE PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res02, alpha = 1/3, qty = 'rmse', yrange = c(0, 5),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Ap2, '\n')
pdf(fig1Ap3, height = 4, width = 3) # RIGHT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res03, alpha = 0, qty = 'rmse', yrange = c(0, 5),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Ap3, '\n')

# Figure S4 (individual panels)
pdf(fig1Bp1, height = 3.5, width = 3.5)  # TOP LEFT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res01, alpha = 1, qty = 'ncovs', yrange = c(0.75, 4.25),
  cex.main = 1.5, add.legend = TRUE)
dev.off(); cat('Printed figure:', fig1Bp1, '\n')
pdf(fig1Bp2, height = 3.5, width = 3.5)  # MIDDLE TOP PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res02, alpha = 1/3, qty = 'ncovs', yrange = c(0.75, 4.25),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Bp2, '\n')
pdf(fig1Bp3, height = 3.5, width = 3.5)  # RIGHT TOP PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res03, alpha = 0, qty = 'ncovs', yrange = c(0.75, 4.25),
  cex.main = 1.5, add.legend = FALSE, var.bac = TRUE)
dev.off(); cat('Printed figure:', fig1Bp3, '\n')
pdf(fig1Bp4, height = 3.5, width = 3.5)  # BOTTOM LEFT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res01, alpha = 1, qty = 'pincl_treat', yrange = c(0.6, 1),
  add.legend = FALSE, show.lasso = FALSE)
dev.off(); cat('Printed figure:', fig1Bp4, '\n')
pdf(fig1Bp5, height = 3.5, width = 3.5)  # BOTTOM MIDDLE PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res02, alpha = 1/3, qty = 'pincl_treat', yrange = c(0.6, 1),
  add.legend = FALSE, var.bac = TRUE, show.lasso = FALSE)
dev.off(); cat('Printed figure:', fig1Bp5, '\n')
pdf(fig1Bp6, height = 3.5, width = 3.5)  # BOTTOM RIGHT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 1, 0.5), bg = 'white')
plot.compare.st(res03, alpha = 0, qty = 'pincl_treat', yrange = c(0, 0.25),
  add.legend = FALSE, var.bac = TRUE, show.lasso = FALSE)
dev.off(); cat('Printed figure:', fig1Bp6, '\n')

# Erase file traces
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE S5: ST changing p (growing design)
################################################################################
# Figure file names
fig2Ap1 <- paste(FIGDIR, 'fig2Ap1.pdf', sep = '')
fig2Ap2 <- paste(FIGDIR, 'fig2Ap2.pdf', sep = '')
fig2Ap3 <- paste(FIGDIR, 'fig2Ap3.pdf', sep = '')

# Plot per panels
pdf(fig2Ap1, height = 4, width = 3) # LEFT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res00, alpha = 1, qty = 'rmse', yrange = c(0.5, 5.5),
  cex.main = 1.35, add.legend = TRUE, legpos = 'topright', legend.EB = TRUE)
dev.off(); cat('Printed figure:', fig2Ap1, '\n')
pdf(fig2Ap2, height = 4, width = 3) # MIDDLE PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res04, alpha = 1, qty = 'rmse', yrange = c(0.5, 5.5),
  cex.main = 1.35, add.legend = FALSE)
dev.off(); cat('Printed figure:', fig2Ap2, '\n')
pdf(fig2Ap3, height = 4, width = 3) # RIGHT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res05, alpha = 1, qty = 'rmse', yrange = c(0.5, 5.5),
  cex.main = 1.35, show.EB = TRUE, add.legend = FALSE)
dev.off(); cat('Printed figure:', fig2Ap3, '\n')

# Erase file traces
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE S6: ST p* -> Inf (increase in confounding)
################################################################################
# Figure file names
fig3Ap1 <- paste(FIGDIR, 'fig3Ap1.pdf', sep = '')
fig3Ap2 <- paste(FIGDIR, 'fig3Ap2.pdf', sep = '')
fig3Ap3 <- paste(FIGDIR, 'fig3Ap3.pdf', sep = '')

# Plot per panels
pcr.here <- FALSE  # Include PCR in the figure
pdf(fig3Ap1, height = 4, width = 3) # LEFT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res04, alpha = 1, qty = 'rmse', yrange = c(0, 12),
  cex.main = 1.5, add.legend = TRUE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig3Ap1, '\n')
pdf(fig3Ap2, height = 4, width = 3) # MIDDLE PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res07, alpha = 1, qty = 'rmse', yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig3Ap2, '\n')
pdf(fig3Ap3, height = 4, width = 3) # RIGHT PANEL
par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
plot.compare.st(res06, alpha = 1, qty = 'rmse', yrange = c(0, 12),
  cex.main = 1.5, add.legend = FALSE, show.pcr = pcr.here)
dev.off(); cat('Printed figure:', fig3Ap3, '\n')

# Erase file traces
rm(pcr.here)
rm(list = ls()[grep('fig', ls())])

################################################################################
# FIGURE 5: MT with increasing number of treatments
################################################################################
# Plot files
figX <- paste(FIGDIR, 'figX.pdf', sep = '')
figXp1 <- paste(FIGDIR, 'figXp1.pdf', sep = '')
figXp2 <- paste(FIGDIR, 'figXp2.pdf', sep = '')

# Transform input data
dres08a <- adapt.output.mt(res08a)
dres08b <- adapt.output.mt(res08b)
dres08c <- adapt.output.mt(res08c)
dres08d <- adapt.output.mt(res08d)

# Assign results
res2T1 <- dres08a[[1]]; res3T1 <- dres08b[[1]]
res4T1 <- dres08c[[1]]; res5T1 <- dres08d[[1]]
res2T2 <- dres08a[[2]]; res3T2 <- dres08b[[2]]
res4T2 <- dres08c[[2]]; res5T2 <- dres08d[[2]]
res2T3 <- dres08a[[3]]; res3T3 <- dres08b[[3]]
res4T3 <- dres08c[[3]]; res5T3 <- dres08d[[3]]
res2T4 <- dres08a[[4]]; res3T4 <- dres08b[[4]]
res4T4 <- dres08c[[4]]; res5T4 <- dres08d[[4]]
res2T5 <- dres08a[[5]]; res3T5 <- dres08b[[5]]
res4T5 <- dres08c[[5]]; res5T5 <- dres08d[[5]]

# We plot averages across present treatments
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
# pdf(figXp1, height = 4, width = 4)
# par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5), bg = 'white')
# plot.mtrand(x1, x2, x3, x4, alpha = 1, yrange = c(0.5, 7.5),
#   cex.main = 1.5, add.legend = TRUE, legpos = 'topleft', show.pcr = pcr.here)
# dev.off(); cat('Printed figure:', figXp1, '\n')

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
