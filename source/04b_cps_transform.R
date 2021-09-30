################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '04b_cps_transform.R', sep = ''))
################################################################################

################################################################################
# Load data and prune NIU observations and covariates
################################################################################
# Directory
new.dir <- paste(DATDIR, 'cps/', sep = '')
if (! dir.exists(new.dir)) {
  dir.create(new.dir); cat('Created directory:', new.dir, '\n')
}

# Raw data
file <- paste(INPDIR, 'cps2year/salary_2010to2019.RData', sep = '')
raw <- get(load(file = file)); cat('Loaded file:', file, '\n')
cps <- as.data.frame(raw)  # 723163(144853) x 162

# Adapt classes
rownames(cps) <- NULL
for (col in 1:ncol(cps)) {
  if (length(class(cps[, col])) > 1) {
    cps[, col] <- as.vector(cps[, col])
  }
}

# Small changes
cps[which(cps[, 'schllunch'] == 99999), 'schllunch'] <- 0

# Choose people who work between 10 and 60 hours per week
cps <- cps[which(cps[, 'hoursworked'] >= 20 & cps[, 'hoursworked'] <= 60), ]
cps <- cps[which(cps[, 'uhrsworkly'] >= 20 & cps[, 'uhrsworkly'] <= 60), ]
# Cuts n = 144853 to n = 131524

# Choose people who work have some income wage
cps <- cps[which(cps[, 'incomewage'] > 1000), ]
# Cuts n = 131524 to n = 125559

# Choose people that are paid
cps <- cps[which(cps[, 'classworker'] != 'Unpaid family worker'), ]
# Cuts n = 125559 to n = 125539

# Choose families with positive income
cps <- cps[which(cps[, 'famincome'] > 0), ]
cps <- cps[which(cps[, 'personalincome'] > 0), ]
cps <- cps[which(cps[, 'personalincome'] - cps[, 'incomewage'] >= 0), ]
cps <- cps[which(cps[, 'famincome'] - cps[, 'incomewage'] >= 0), ]
# Cuts n = 125539 to n = 124857

# Delete under age observations
cps <- cps[which(cps[, 'age'] >= 18), ]
# Cuts n = 124857 to n = 124458

# Small errors or weird observarions
cps <- cps[which(cps[, 'hhrespln'] != 99), ]
cps <- cps[which(cps[, 'nfams'] <= 3), ]
cps <- cps[which(cps[, 'eldch'] <= 45), ]
cps <- cps[which(cps[, 'yngch'] <= 45), ]
cps <- cps[which(cps[, 'foodstamp'] <= 15e3), ]
cps <- cps[which(cps[, 'schllunch'] <= 5e3), ]
cps <- cps[which(cps[, 'spmchxpns'] <= 1e5), ]
cps <- cps[which(cps[, 'spmchsup'] <= 6e4), ]
# Cuts n = 124458 to n = 123656

# Negative incomes
inc.cols <- c('personalincome', 'incomewage', 'incbus', 'incfarm', 'incss',
  'incwelfr', 'incretir', 'incssi', 'incint', 'incunemp', 'incwkcom', 'incvet',
  'incsurv', 'incdisab', 'incdivid', 'incrent', 'inceduc', 'incchild',
  'incasist', 'incother', 'inclongj', 'oincbus', 'oincfarm', 'oincwage',
  'incsurv1', 'incsurv2', 'incdisa1', 'incdisa2', 'famincome')
cps <- cps[which(rowSums(cps[, inc.cols] < 0) == 0), ]
# Cuts n = 123656 to n = 123265

# What to do with oversamples?
if (FALSE) {
  cps <- cps[which(cps[, 'oversample'] == 0), ]
  # Cuts n = 123265 to n = 81786
  cps[, 'oversample'] <- NULL
}

################################################################################
# Create skeletons
y <- matrix(nrow = nrow(cps), ncol = 1)
X <- matrix(nrow = nrow(cps), ncol = 1)
D <- matrix(nrow = nrow(cps), ncol = 1)

################################################################################
# Outcome
################################################################################
colnames(y) <- 'loghourwage'
y[, 1] <- log((cps[, 'incomewage'] * cps[, 'cpi99']) / cps[, 'hoursworked'])
#y[,1] <- scale(y[, 1], center = TRUE, scale = TRUE)

################################################################################
# Treatments
################################################################################
colnames(D) <- 'female'
D[, 1] <- as.numeric(as.character(cps[, 'female']))
#m0 <- lm(y ~ cbind(1, D) - 1)

################################################################################
# Covariates
################################################################################
# Intercept
X[, 1] <- 1
colnames(X)[1] <- 'intercept'

# Oversampling observation (ARTIFICIAL VARIABLE)
# (over-represents HISPANIC -> younger age)
if ('oversample' %in% colnames(cps)) {
  X <- cbind(X, as.integer(cps[, 'oversample']))
  colnames(X)[ncol(X)] <- 'oversample'
}

# Number of person records following (ARTIFICIAL VARIABLE)
X <- cbind(X, as.integer(cps[, 'numprec']))
colnames(X)[ncol(X)] <- 'numprec'

# Line number of household respondent (ARTIFICIAL VARIABLE)
X <- cbind(X, as.integer(cps[, 'hhrespln']))
colnames(X)[ncol(X)] <- 'hhrespln'

# Person number in sample unit (ARTIFICIAL VARIABLE)
X <- cbind(X, as.integer(cps[, 'pernum']))
colnames(X)[ncol(X)] <- 'pernum'

# Line number on original form (ARTIFICIAL VARIABLE)
X <- cbind(X, as.integer(cps[, 'lineno']))
colnames(X)[ncol(X)] <- 'lineno'

# Year
X <- cbind(X, ifelse(cps[, 'year'] == 2019, 1L, 0L))
colnames(X)[ncol(X)] <- 'is2019'

# Region
X <- cbind(X, conv.dummy(cps, 'region', 'west'))

# State
NULL  # CAREFUL! We ignore the U.S. STATE for now

# Metropolitan area
X <- cbind(X, conv.dummy(cps, 'metro', 'no'))

# Size of the area
X <- cbind(X, conv.dummy(cps, 'areasize', '0 (<100k)'))

# Age
X <- cbind(X, cps[, 'age'])
colnames(X)[ncol(X)] <- 'age'

# Race
X <- cbind(X, conv.dummy(cps, 'race', 'White'))

# Marital status
X <- cbind(X, conv.dummy(cps, 'marital', 'Married'))

# Education level
X <- cbind(X, conv.dummy(cps, 'edu', 'HSD'))

# School attendance
inschool <- c('High school', 'College or university')
X <- cbind(X, ifelse(cps[, 'schlcoll'] %in% inschool, 1L, 0L))
colnames(X)[ncol(X)] <- 'attendsedu'

# Hispanic
X <- cbind(X, as.numeric(as.character(cps[, 'hispanic'])))
colnames(X)[ncol(X)] <- 'hispanic'

# Birthplace
#X <- cbind(X, conv.dummy(cps, 'bpl', 'HSD'))
northam <- c('US mainland', 'US overseas', 'North America')
latam <- c('Caribbean', 'Central America', 'South America')
X <- cbind(X, ifelse(cps[, 'bpl'] %in% latam, 1L, 0L))
colnames(X)[ncol(X)] <- 'bpl_latam'
X <- cbind(X, ifelse(cps[, 'bpl'] == 'Africa', 1L, 0L))
colnames(X)[ncol(X)] <- 'bpl_africa'
X <- cbind(X, ifelse(cps[, 'bpl'] == 'Asia', 1L, 0L))
colnames(X)[ncol(X)] <- 'bpl_asia'
X <- cbind(X, ifelse(cps[, 'bpl'] == 'Europe', 1L, 0L))
colnames(X)[ncol(X)] <- 'bpl_europe'
X <- cbind(X, ifelse(cps[, 'bpl'] == 'Oceania', 1L, 0L))
colnames(X)[ncol(X)] <- 'bpl_oceania'

# Time residing in the USA
X <- cbind(X, conv.dummy(cps, 'timeus', '0 (native)'))

# Citizenship
X <- cbind(X, conv.dummy(cps, 'citizen', 'Born US'))

# Mother's birthplace
X <- cbind(X, ifelse(cps[, 'mbpl'] %in% latam, 1L, 0L))
colnames(X)[ncol(X)] <- 'mbpl_latam'
X <- cbind(X, ifelse(cps[, 'mbpl'] == 'Africa', 1L, 0L))
colnames(X)[ncol(X)] <- 'mbpl_africa'
X <- cbind(X, ifelse(cps[, 'mbpl'] == 'Asia', 1L, 0L))
colnames(X)[ncol(X)] <- 'mbpl_asia'
X <- cbind(X, ifelse(cps[, 'mbpl'] == 'Europe', 1L, 0L))
colnames(X)[ncol(X)] <- 'mbpl_europe'
X <- cbind(X, ifelse(cps[, 'mbpl'] == 'Oceania', 1L, 0L))
colnames(X)[ncol(X)] <- 'mbpl_oceania'

# Father's birthplace
X <- cbind(X, ifelse(cps[, 'fbpl'] %in% latam, 1L, 0L))
colnames(X)[ncol(X)] <- 'fbpl_latam'
X <- cbind(X, ifelse(cps[, 'fbpl'] == 'Africa', 1L, 0L))
colnames(X)[ncol(X)] <- 'fbpl_africa'
X <- cbind(X, ifelse(cps[, 'fbpl'] == 'Asia', 1L, 0L))
colnames(X)[ncol(X)] <- 'fbpl_asia'
X <- cbind(X, ifelse(cps[, 'fbpl'] == 'Europe', 1L, 0L))
colnames(X)[ncol(X)] <- 'fbpl_europe'
X <- cbind(X, ifelse(cps[, 'fbpl'] == 'Oceania', 1L, 0L))
colnames(X)[ncol(X)] <- 'fbpl_oceania'

# Nativity
X <- cbind(X, conv.dummy(cps, 'nativity', 'Both US'))

# Moved state
X <- cbind(X, ifelse(cps[, 'movedstate'] == 'No', 0L, 1L))
colnames(X)[ncol(X)] <- 'movedstate'

# Moved house
X <- cbind(X, conv.dummy(cps, 'movedhouse', 'Same house'))

# Class of worker
X <- cbind(X, conv.dummy(cps, 'classworker', 'Wage/salary'))

# Change of class from previous year
X <- cbind(X, as.integer(cps[, 'chclasswkr']))
colnames(X)[ncol(X)] <- 'chclasswkr'

# Presence in labor force
NULL  # Only people "At work"

# Labor force information collected by self or proxy response
X <- cbind(X, as.integer(cps[, 'lfproxy']))
colnames(X)[ncol(X)] <- 'lfproxy'

# Occupation
X <- cbind(X, conv.dummy(cps, 'occ', 'office'))

# Changed occupation since last year
X <- cbind(X, as.integer(cps[, 'choccly']))
colnames(X)[ncol(X)] <- 'choccly'

# Industry
X <- cbind(X, conv.dummy(cps, 'ind', 'terciary'))

# Changed industry since last year
X <- cbind(X, as.integer(cps[, 'chindly']))
colnames(X)[ncol(X)] <- 'chindly'

# Firm size
X <- cbind(X, conv.dummy(cps, 'firmsize', '5 (>1000)'))

# Hours worked per week
X <- cbind(X, as.numeric(cps[, 'hoursworked']))
colnames(X)[ncol(X)] <- 'hoursworked'

# Work status
X <- cbind(X, ifelse(cps[, 'wkstat'] == 'Part-time', 1L, 0L))
colnames(X)[ncol(X)] <- 'wkstat_pt'

# Has more than one job
X <- cbind(X, ifelse(cps[, 'pcthrsmainjob'] == 1, 0L, 1L))
colnames(X)[ncol(X)] <- 'hasmorejobs'

# Number of weeks worked per year
X <- cbind(X, cps[, 'wkswork1'])
colnames(X)[ncol(X)] <- 'numweekswkd'

# Usually works full-time
X <- cbind(X, as.integer(cps[, 'usftptlw']))
colnames(X)[ncol(X)] <- 'usftptlw'

# Worked last year
NULL  # Everybody said yes
#X <- cbind(X, as.integer(cps[, 'workly']))
#colnames(X)[ncol(X)] <- 'workly'

# Usual hours worked per week (last year)
X <- cbind(X, cps[, 'uhrsworkly'])
colnames(X)[ncol(X)] <- 'uhrsworkly'

# Num. weeks unemployed last year
X <- cbind(X, cps[, 'wksunemly'])
colnames(X)[ncol(X)] <- 'wksunemly'

# Did they look for work last year
#X[, 'nwlookwkly'] <- NULL  # No relevant info

# A number of binary variables coming up
bin.vars <- c('paidhour', 'union', 'ssocinc', 'pension', 'veteran', 'medicare',
  'medicaid', 'milhealth', 'phinsur', 'phiown', 'out', 'grpdeply', 'grpownly',
  'grpoutly', 'dpdeply', 'dpownly', 'dpoutly', 'diffany', 'diffmob', 'diffeye',
  'diffrem', 'diffhear', 'diffphys', 'diffcare', 'disabwrk', 'quitsick')
for (nomcol in bin.vars) {
  X <- cbind(X, as.integer(cps[, nomcol]))
  colnames(X)[ncol(X)] <- nomcol
}

# Related numerical variables
num.vars <- c('numempsly', 'hiunpers', 'grpwho1', 'dpwho1')
for (nomcol in num.vars) {
  X <- cbind(X, as.numeric(cps[, nomcol]))
  colnames(X)[ncol(X)] <- nomcol
}

# Employer paid for group health plan
X <- cbind(X, conv.dummy(cps, 'emphplan', 'no'))

# Type of employment-based coverage last year
X <- cbind(X, conv.dummy(cps, 'grptyply', 'none'))

# Type of direct-purchase insurance plan, previous year
X <- cbind(X, conv.dummy(cps, 'dptyply', 'none'))

# Health status
X <- cbind(X, conv.dummy(cps, 'health', 'excellent'))

# Tax filer status
X <- cbind(X, conv.dummy(cps, 'taxstatus', 'joint'))

# Child Tax Credit (log)
X <- cbind(X, log(as.numeric(cps[, 'ctccrd'] * cps[, 'cpi99']) + 1))
colnames(X)[ncol(X)] <- 'logctccrd'

# Earned Income Tax Credit
X <- cbind(X, as.numeric(cps[, 'eitcred'] > 0))  # Make it binary
colnames(X)[ncol(X)] <- 'haseitcred'

# Federal income tax liability, after all credits (log) [ ?]
X <- cbind(X, log(cps[, 'fedtaxac'] * cps[, 'cpi99'] + 1e4))
#X <- cbind(X, log(cps[, 'fedtaxac'] - min(cps[, 'fedtaxac']) + 1))
colnames(X)[ncol(X)] <- 'logfedtaxac'

# Social security retirement payroll deduction (log)
X <- cbind(X, log(as.numeric(cps[, 'fica'] * cps[, 'cpi99']) + 1))
colnames(X)[ncol(X)] <- 'logfica'

# State income tax liability, after all credits (log) [ ?]
X <- cbind(X, log(cps[, 'stataxac'] * cps[, 'cpi99'] + 1e4))
colnames(X)[ncol(X)] <- 'logstataxac'

# Federal income marginal tax rate
X <- cbind(X, as.numeric(cps[, 'margtax']))
colnames(X)[ncol(X)] <- 'margtax'
 
# Proportion of salary taxable
log(cps[, 'taxinc'] / cps[, 'personalincome'] + 1)
colnames(X)[ncol(X)] <- 'logproptaxincwageinc'

# Type of tenure
X <- cbind(X, conv.dummy(cps, 'hhtenure', 'own'))

# Type of bulding
X <- cbind(X, cps[, 'buildtype'])
colnames(X)[ncol(X)] <- 'buildtype'

# Active mortgage
X <- cbind(X, cps[, 'mortgage'])
colnames(X)[ncol(X)] <- 'mortgage'

# Phone availability
X <- cbind(X, cps[, 'phone'])
colnames(X)[ncol(X)] <- 'phone'

# Variables of family member counts
X <- cbind(X,
  as.numeric(cps[, 'nfams']),  # Number of families in household (HH)
  as.numeric(cps[, 'ncouples']),  # Number of couples in HH
  as.numeric(cps[, 'nmothers']),  # Number of mothers in HH
  as.numeric(cps[, 'nfathers']),  # Number of fathers in HH
  as.numeric(cps[, 'famsize']),  # Number of family members
  as.numeric(cps[, 'nchild']),  # Number of children in HH
  as.numeric(cps[, 'nchlt5']),  # Number of children under 5 in HH
  as.numeric(cps[, 'eldch']),  # Age of eldest child
  as.numeric(cps[, 'yngch']),  # Age of youngest child
  as.numeric(cps[, 'famunit']),  # Size of family unit
  as.numeric(cps[, 'nsibs']))  # Number of siblings in HH
noms.fam <- c('nfams', 'ncouples', 'nmothers', 'nfathers', 'famsize', 'nchild',
  'nchlt5', 'eldch', 'yngch', 'famunit', 'nsibs')
colnames(X)[ncol(X) - ((length(noms.fam) - 1):0)] <- noms.fam

# Family type
X <- cbind(X, conv.dummy(cps, 'ftype', 'primary'))

# Kind of family
X <- cbind(X, conv.dummy(cps, 'fkind', 'couple'))

# Person status in family
X <- cbind(X, conv.dummy(cps, 'famrel', 'reference'))

# Relation to head of household
X <- cbind(X, conv.dummy(cps, 'relhh', 'head'))

# Is the respondent related to the HH
X <- cbind(X, ifelse(cps[, 'isrelhh'] == 'related', 1L, 0L))
colnames(X)[ncol(X)] <- 'isrel2hh'

# Another set of binary variables (subsidies / public housing)
bin.vars <- c('pubhous', 'rentsub', 'heatsub', 'lunchsub', 'foodstmp')
for (nomcol in bin.vars) {
  X <- cbind(X, as.integer(cps[, nomcol]))
  colnames(X)[ncol(X)] <- nomcol
}

# Related numerical variables
num.vars <- c('atelunch', 'frelunch', 'stampno', 'stampmo')
for (nomcol in num.vars) {
  X <- cbind(X, as.numeric(cps[, nomcol]))
  colnames(X)[ncol(X)] <- nomcol
}

# Family market value of food stamps
X <- cbind(X, log(cps[, 'foodstamp'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logfoodstamp'

# Family market value of school lunch
X <- cbind(X, log(cps[, 'schllunch'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logschllunch'

# SPM unit's school lunch value 
X <- cbind(X, log(cps[, 'spmlunch'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logspmlunch'

# SPM unit's medical out-of-pocket and Medicare B subsidy
X <- cbind(X, log(cps[, 'spmmedxpns'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logspmmedxpns'

# SPM unit's child care expenses - not capped
X <- cbind(X, log(cps[, 'spmchxpns'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logspmchxpns'

# SPM unit's capped work and child care expenses
X <- cbind(X, log(cps[, 'spmcapxpns'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logspmcapxpns'

# SPM unit's child support paid
X <- cbind(X, log(cps[, 'spmchsup'] * cps[, 'cpi99'] + 1))
colnames(X)[ncol(X)] <- 'logspmchsup'

# SPM unit's number of adults
X <- cbind(X, as.numeric(cps[, 'spmnadults']))
colnames(X)[ncol(X)] <- 'logspmnadults'

# SPM unit's number of children
X <- cbind(X, as.numeric(cps[, 'spmnchild']))
colnames(X)[ncol(X)] <- 'logspmnchild'

# SPM unit's number of persons (spmnpers)
NULL  # It is sum of previous two

# Most relevant income columns
# Family income
X <- cbind(X,
  log(1 + (cps[, 'personalincome'] - cps[, 'incomewage']) * cps[, 'cpi99']))
colnames(X)[ncol(X)] <- 'logotherpersincome'
X <- cbind(X,
  log(1 + (cps[, 'famincome'] - cps[, 'incomewage']) * cps[, 'cpi99']))
colnames(X)[ncol(X)] <- 'logotherfamincome'

# Different income columns
inc.cols <- c('incbus', 'incfarm', 'incss', 'incwelfr', 'incretir', 'incssi',
  'incint', 'incunemp', 'incwkcom', 'incvet', 'incsurv', 'incdisab', 'incdivid',
  'incrent', 'inceduc', 'incchild', 'incasist', 'incother', #'inclongj',
  'oincbus', 'oincfarm', 'incsurv1', 'incsurv2', 'incdisa1', #'oincwage',
  'incdisa2')
totothinc <- rowSums(cps[, inc.cols]) + 1
aux <- cps[, inc.cols] / totothinc
colnames(aux) <- paste('propoi', inc.cols, sep = '_')

# We add what proportion of other income is due to each source
X <- cbind(X, as.matrix(aux))  # 123265 x 239

# Keep state for state-level analysis later
state <- cps[, 'state']
#stateD <- conv.dummy(cps, 'state')

# Keep year for year-level analysis later
year <- cps[, 'year']

# Print results
cat('* Design matrix dimension:', nrow(X), 'x', ncol(X), '\n')

################################################################################
# Save data
################################################################################
file.out <- paste(DATDIR, 'cps/cps_model_data_year.RData', sep = '')
save(y, D, X, state, year, file = file.out); cat('Saved file:', file.out, '\n')
# END OF SCRIPT
