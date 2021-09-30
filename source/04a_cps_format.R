################################################################################
# source('~/Desktop/year3/bma_teff/v16/syntax/00_start.R')
# source(paste(SRCDIR, '04a_cps_format.R', sep = ''))
################################################################################

################################################################################
# Obtaining the data
################################################################################
# ALTERNATIVE OPTION TO MANUAL DOWNLOAD:
# Download https://www2.census.gov/programs-surveys/cps/datasets/2019/march/asecpub19csv.zip
# Download https://www2.census.gov/programs-surveys/cps/datasets/2010/march/asec2010_pubuse.dat.gz
# Download https://www2.census.gov/programs-surveys/cps/datasets/2019/march/asec2019_pubuse.dat.gz
# 
# See https://cengel.github.io/gearup2016/SULdataAccess.html for R packages to
# get CPS data. Main packages are impusr and icpsrdata
#
# NOTE: if using ipums data or icpsrdata cite
#   Sarah Flood, Miriam King, Renae Rodgers, Steven Ruggles and J. Robert
#   Warren. Integrated Public Use Microdata Series, Current Population Survey:
#   Version 7.0 [dataset]. Minneapolis, MN: IPUMS, 2020.
#   https://doi.org/10.18128/D030.V7.0
################################################################################
# Variable names and labels available at cps_00004.cbk
# vignette("ipums-cps", package = "ipumsr")
# NOTE: To load data, you must download both the extract's data and the DDI
################################################################################

################################################################################
# Reading
################################################################################
# Directory
if (! dir.exists(paste(INPDIR, 'cps2year/', sep = ''))) {
  new.dir <- paste(INPDIR, 'cps2year/', sep = '')
  dir.create(new.dir); cat('Created directory:', new.dir, '\n')
}

# Dependencies
require(ipumsr)

# Read file
file.out <- paste(INPDIR, 'cps2year/cps04.RData', sep = '')
if (! file.exists(file.out)) {
  cat('* Reading raw .fwf files...\n')
  ddi.file <- paste(INPDIR, 'cps2year/cps_00004.xml', sep = '')
  ddi.read <- ipumsr::read_ipums_ddi(ddi.file)
  raw <- ipumsr::read_ipums_micro(ddi.read)
  save(raw, file = file.out); cat('Done.\nSaved file:', file.out, '\n')
} else {
  cat('* Reading raw .RData... ')
  raw <- get(load(file = file.out)); cat('Done.\nLoaded file:', file.out, '\n')
}; rm(file.out)

################################################################################
# Formatting input + discarding of useless/repetitive variables
################################################################################
data <- as.data.frame(raw)
colnames(data) <- tolower(colnames(data))

# Remove all weights
data[, colnames(data)[grepl('repwt', colnames(data))]] <- NULL

# Choose valid samples
data <- data[which(data[, 'gq'] == 1), ]  # Only households

# Remove more irrelevant/repetitive/id variables
rm.vars.in <- c('asecflag', 'asecwth', 'mish', 'gq', 'hhintype', 'metarea',
  'statecensus', 'ownershp', 'faminc', 'durunem2', 'whyunemp', 'whyabsnt',
  'whyptlwk', 'wnftlook', 'asecoverp', 'wksunem2', 'fullpart', 'wkswork2',
  'wantjob', 'whyptly', 'payifabs', 'wnlwnilf', 'whynwly', 'actnlfly',
  'srcsurv1', 'srcsurv2', 'srcdisa1', 'srcdisa2', 'srcearn', 'srceduc',
  'srcunem', 'srcwelfr', 'srcwkcom', 'vetqa', 'whyssi1', 'whyssi2', 'gotvdisa',
  'gotveduc', 'gotvothe', 'gotvpens', 'gotvsurv', 'offpov', 'offpovuniv',
  'offreason', 'poverty', 'spmpov', 'whymove', 'caidly', #'grpwho1', 'dpwho1',
  'trccovly', 'champvaly', 'inhcovly', 'schiply', 'hiurule', 'vet1', 'vet2',
  'vet3', 'vet4', 'gotwic', 'kidcneed', 'eligorg', 'momrule', 'mom2rule',
  'poprule', 'pop2rule', 'sprule')  # In descriptive file
rm.vars.out <- c('county', 'metfips', 'individcc', 'heatval', 'stampval',
  'marbasecidh', 'hrhhid', 'hrhhid2', 'hseq', 'asecwt', 'momloc', 'momloc2',
  'poploc', 'poploc2', 'sploc', 'aspouse', 'pecohab', 'ahrsworkt', 'durunemp',
  'marbasecidp', 'earnwt', 'lnkfw1ywt', 'famid', 'asecfwt', 'mthwelfr',
  'fedretir', 'fedtax', 'depstat', 'statetax', 'offtotval', 'offcutoff',
  'cutoff', 'hiuid', 'hiufpginc', 'hiufpginc', 'hourwage', 'earnweek',
  'hiufpgbase', 'srcunemp', 'hhincome')  # Not in descriptive file
data[, unique(c(rm.vars.in, rm.vars.out))] <- NULL

################################################################################
# Adaptation of variables
################################################################################
# IDs
names(data)[names(data) == 'serial'] <- 'householdid'
data[, 'cpsid'] <- as.character(data[, 'cpsid'])
data[, 'cpsidp'] <- as.character(data[, 'cpsidp'])

# Oversample (this might be relevant, we don't know)
data[, 'oversample'] <- as.integer(data[, 'asecoverh'])
data[, 'asecoverh'] <- NULL

# Household tenure
hhtenure <- rep(NA, nrow(data))
hhtenure[data[, 'hhtenure'] == 1] <- 'own'
hhtenure[data[, 'hhtenure'] == 2] <- 'rent'
hhtenure[data[, 'hhtenure'] == 3] <- 'occupy'
data[, 'hhtenure'] <- as.factor(hhtenure)

# Region
region <- rep(NA, nrow(data))
region[data[, 'region'] == 11 | data[, 'region'] == 12] <- 'northeast'
region[data[, 'region'] == 21 | data[, 'region'] == 22] <- 'midwest'
region[data[, 'region'] == 31 | data[, 'region'] == 32] <- 'south'
region[data[, 'region'] == 41 | data[, 'region'] == 42] <- 'west'
region[data[, 'region'] == 33] <- 'south'
data[, 'region'] <- as.factor(region)

# State
state <- rep(NA, nrow(data))
state[data[, 'statefip'] ==  1] <- 'AL'
state[data[, 'statefip'] ==  2] <- 'AK'
state[data[, 'statefip'] ==  4] <- 'AZ'
state[data[, 'statefip'] ==  5] <- 'AR'
state[data[, 'statefip'] ==  6] <- 'CA'
state[data[, 'statefip'] ==  8] <- 'CO'
state[data[, 'statefip'] ==  9] <- 'CT'
state[data[, 'statefip'] == 10] <- 'DE'
state[data[, 'statefip'] == 11] <- 'DC'
state[data[, 'statefip'] == 12] <- 'FL'
state[data[, 'statefip'] == 13] <- 'GA'
state[data[, 'statefip'] == 15] <- 'HI'
state[data[, 'statefip'] == 16] <- 'ID'
state[data[, 'statefip'] == 17] <- 'IL'
state[data[, 'statefip'] == 18] <- 'IN'
state[data[, 'statefip'] == 19] <- 'IA'
state[data[, 'statefip'] == 20] <- 'KS'
state[data[, 'statefip'] == 21] <- 'KY'
state[data[, 'statefip'] == 22] <- 'LA'
state[data[, 'statefip'] == 23] <- 'ME'
state[data[, 'statefip'] == 24] <- 'MD'
state[data[, 'statefip'] == 25] <- 'MA'
state[data[, 'statefip'] == 26] <- 'MI'
state[data[, 'statefip'] == 27] <- 'MN'
state[data[, 'statefip'] == 28] <- 'MS'
state[data[, 'statefip'] == 29] <- 'MO'
state[data[, 'statefip'] == 30] <- 'MT'
state[data[, 'statefip'] == 31] <- 'NE'
state[data[, 'statefip'] == 32] <- 'NV'
state[data[, 'statefip'] == 33] <- 'NH'
state[data[, 'statefip'] == 34] <- 'NJ'
state[data[, 'statefip'] == 35] <- 'NM'
state[data[, 'statefip'] == 36] <- 'NY'
state[data[, 'statefip'] == 37] <- 'NC'
state[data[, 'statefip'] == 38] <- 'ND'
state[data[, 'statefip'] == 39] <- 'OH'
state[data[, 'statefip'] == 40] <- 'OK'
state[data[, 'statefip'] == 41] <- 'OR'
state[data[, 'statefip'] == 42] <- 'PA'
state[data[, 'statefip'] == 44] <- 'RI'
state[data[, 'statefip'] == 45] <- 'SC'
state[data[, 'statefip'] == 46] <- 'SD'
state[data[, 'statefip'] == 47] <- 'TN'
state[data[, 'statefip'] == 48] <- 'TX'
state[data[, 'statefip'] == 49] <- 'UT'
state[data[, 'statefip'] == 50] <- 'VT'
state[data[, 'statefip'] == 51] <- 'VA'
state[data[, 'statefip'] == 53] <- 'WA'
state[data[, 'statefip'] == 54] <- 'WV'
state[data[, 'statefip'] == 55] <- 'WI'
state[data[, 'statefip'] == 56] <- 'WY'
data[, 'state'] <- as.factor(state)
data[, 'statefip'] <- NULL

# Living in Metro area
metro <- rep(NA, nrow(data))
metro[data[, 'metro'] == 0] <- 'NA'
metro[data[, 'metro'] == 1] <- 'no'
metro[data[, 'metro'] == 2] <- 'in central city'
metro[data[, 'metro'] == 3] <- 'out central city'
metro[data[, 'metro'] == 4] <- 'unknown central city'
data[, 'metro'] <- metro
data <- data[which(data[, 'metro'] != 'NA'), ]
data[, 'metro'] <- as.factor(data[, 'metro'])

# Core-based statistical area size
areasize <- rep(NA, nrow(data))
areasize[data[, 'cbsasz'] == 0] <- '0 (<100k)'
areasize[data[, 'cbsasz'] == 1] <- '1 (100k-250k)'
areasize[data[, 'cbsasz'] == 2] <- '2 (250k-500k)'
areasize[data[, 'cbsasz'] == 3] <- '3 (500k-1M)'
areasize[data[, 'cbsasz'] == 4] <- '4 (1M-2.5M)'
areasize[data[, 'cbsasz'] == 5] <- '5 (2.5M-5M)'
areasize[data[, 'cbsasz'] == 6] <- '6 (>5)'
data[, 'areasize'] <- as.factor(areasize)
data[, 'cbsasz'] <- NULL

# Public help in housing/commodities
data[, 'pubhous'] <- ifelse(data[, 'pubhous'] == 2, 1, 0)
data[, 'rentsub'] <- ifelse(data[, 'rentsub'] == 2, 1, 0)
data[, 'heatsub'] <- ifelse(data[, 'heatsub'] == 2, 1, 0)

# Public help in food
data[, 'lunchsub'] <- ifelse(data[, 'lunchsub'] == 1, 1, 0)
data[which(data[, 'atelunch'] == 99), 'atelunch'] <- 0
data[which(data[, 'frelunch'] > 9), 'frelunch'] <- 0

# Units in structure
buildtype <- character(nrow(data))
buildtype[data[, 'unitsstr'] %in% c(1, 11)] <- 'single'
buildtype[data[, 'unitsstr'] %in% c(5, 6, 7, 12)] <- 'multiple'
data[, 'buildtype'] <- factor(buildtype)
data[, 'unitsstr'] <- NULL

# Access to telephone
data[, 'phone'] <- ifelse(data[, 'phone'] == 1, 0, 1)

# Relation to head of household
isrelhh <- relhh <- character(nrow(data))
relhh[data[, 'relate'] == 101] <- 'head'
relhh[data[, 'relate'] %in% c(201, 202, 203)] <- 'spouse'
relhh[data[, 'relate'] %in% c(301, 303)] <- 'child'
relhh[data[, 'relate'] == 501] <- 'parent'
relhh[data[, 'relate'] == 701] <- 'sibling'
relhh[data[, 'relate'] == 901] <- 'grandchild'
relhh[data[, 'relate'] == 1001] <- 'otherrel'
relhh[data[, 'relate'] %in% c(1113, 1114, 1116, 1117)] <- 'partner'
relhh[data[, 'relate'] == 1115] <- 'roomate'
relhh[data[, 'relate'] == 1241] <- 'lodger'
relhh[data[, 'relate'] == 1242] <- 'fosterchild'
relhh[data[, 'relate'] == 1260] <- 'othernonrel'
isrelhh[data[, 'relate'] <= 1001] <- 'related'
isrelhh[data[, 'relate'] >= 1113] <- 'unrelated'
data[, 'relhh'] <- factor(relhh)
data[, 'isrelhh'] <- factor(isrelhh)
data[, 'relate'] <- NULL

# Sex
data[, 'female'] <- ifelse(data[, 'sex'] == 2, 1, 0)
data[, 'sex'] <- NULL

# Race
race <- character(nrow(data))
race[data[, 'race'] == 100] <- 'White'
race[data[, 'race'] == 200] <- 'Black'
race[data[, 'race'] == 300] <- 'Native American'
race[data[, 'race'] == 651] <- 'Asian'
race[data[, 'race'] == 652] <- 'Hawaiian/Pacific Islander'
data[, 'race'] <- race
data <- data[which(data[, 'race'] != ''), ]
data[, 'race'] <- factor(data[, 'race'])

# Marital status
marital <- character(nrow(data))
marital[data[, 'marst'] ==1] <- 'Married'
marital[data[, 'marst'] ==2] <- 'Married'
marital[data[, 'marst'] ==3] <- 'Separated'
marital[data[, 'marst'] ==4] <- 'Divorced'
marital[data[, 'marst'] ==5] <- 'Widowed'
marital[data[, 'marst'] ==6] <- 'NeverMarried'
data[, 'marital'] <- marital
data[, 'marst'] <- NULL
data <- data[which(data[, 'marital'] != ''), ]
data[, 'marital'] <- factor(data[, 'marital'])

# Population status
data <- data[which(data[, 'popstat'] != 3), ]  # No children obs.
data[, 'popstat'] <- NULL

# Veteran
data[, 'veteran'] <- ifelse(data[, 'vetstat'] == 2, 1, 0)
data[, 'vetstat'] <- NULL

# Family unit
data <- data[which(data[, 'famunit'] <= 3), ]  # Suppress very weird obs.

# Age of children
data[which(data[, 'eldch'] == 99), 'eldch'] <- -1
data[which(data[, 'yngch'] == 99), 'yngch'] <- -1

# Family type
ftype <- rep(NA, nrow(data))
ftype[data[, 'ftype'] == 1] <- 'primary'
ftype[data[, 'ftype'] == 2] <- 'nonfamily'
ftype[data[, 'ftype'] == 3] <- 'related subfamily'
ftype[data[, 'ftype'] == 4] <- 'unrelated subfamily'
ftype[data[, 'ftype'] == 5] <- 'secondary'
data[, 'ftype'] <- as.factor(ftype)

# Reference kind of family
fkind <- character(nrow(data))
fkind[data[, 'famkind'] == 1] <- 'couple'
fkind[data[, 'famkind'] == 2] <- 'maleref'
fkind[data[, 'famkind'] == 3] <- 'femaleref'
data[, 'fkind'] <- as.factor(fkind)
data[, 'famkind'] <- NULL

# Family relationship
famrel <- rep(NA, nrow(data))
famrel[data[, 'famrel'] == 0] <- 'not member'
famrel[data[, 'famrel'] == 1] <- 'reference'
famrel[data[, 'famrel'] == 2] <- 'spouse'
famrel[data[, 'famrel'] == 3] <- 'child'
famrel[data[, 'famrel'] == 4] <- 'other'
data[, 'famrel'] <- as.factor(famrel)

# Birthplace
bpl <- rep(NA, nrow(data))
bpl[data[, 'bpl'] == 9900] <- 'US mainland'
bpl[data[, 'bpl'] >= 10000 & data[, 'bpl'] <= 12090] <- 'US overseas'
bpl[data[, 'bpl'] >= 15000 & data[, 'bpl'] <= 19900] <- 'North America'
bpl[data[, 'bpl'] >= 20000 & data[, 'bpl'] <= 21090] <- 'Central America'
bpl[data[, 'bpl'] >= 25000 & data[, 'bpl'] <= 26091] <- 'Caribbean'
bpl[data[, 'bpl'] >= 30005 & data[, 'bpl'] <= 31000] <- 'South America'
bpl[data[, 'bpl'] >= 40000 & data[, 'bpl'] <= 49900] <- 'Europe'
bpl[data[, 'bpl'] >= 50000 & data[, 'bpl'] <= 59900] <- 'Asia'
bpl[data[, 'bpl'] >= 60010 & data[, 'bpl'] <= 60099] <- 'Africa'
bpl[data[, 'bpl'] >= 70010 & data[, 'bpl'] <= 72000] <- 'Oceania'
bpl[data[, 'bpl'] == 96000] <- 'NA'
data[, 'bpl'] <- bpl
data <- data[which(bpl != 'NA'), ]
data[, 'bpl'] <- as.factor(data[, 'bpl'] )

# Time in the US
timeus <- rep(NA, nrow(data))
timeus[data[, 'yrimmig'] == 0] <- '0 (native)'
timeus[data[, 'year'] - data[, 'yrimmig'] < 71] <- '1 (>30yr)'
timeus[data[, 'year'] - data[, 'yrimmig'] < 31] <- '2 (16-30yr)'
timeus[data[, 'year'] - data[, 'yrimmig'] < 16] <- '3 (5-15yr)'
timeus[data[, 'year'] - data[, 'yrimmig'] <  5] <- '4 (<5yr)'
data[, 'timeus'] <- as.factor(timeus)
data[, 'yrimmig'] <- NULL

# Citizenship
citizen <- rep(NA, nrow(data))
citizen[data[, 'citizen'] == 1] <- 'Born US'
citizen[data[, 'citizen'] == 2] <- 'Born US'
citizen[data[, 'citizen'] == 3] <- 'Born abroad, US parents'
citizen[data[, 'citizen'] == 4] <- 'Citizen'
citizen[data[, 'citizen'] == 5] <- 'Not citizen'
data[, 'citizen'] <- as.factor(citizen)

# Mother's birthplace
mbpl <- rep(NA, nrow(data))
mbpl[data[, 'mbpl'] == 9900] <- 'US mainland'
mbpl[data[, 'mbpl'] >= 10000 & data[, 'mbpl'] <= 12090] <- 'US overseas'
mbpl[data[, 'mbpl'] >= 15000 & data[, 'mbpl'] <= 19900] <- 'North America'
mbpl[data[, 'mbpl'] >= 20000 & data[, 'mbpl'] <= 21090] <- 'Central America'
mbpl[data[, 'mbpl'] >= 25000 & data[, 'mbpl'] <= 26091] <- 'Caribbean'
mbpl[data[, 'mbpl'] >= 30005 & data[, 'mbpl'] <= 31000] <- 'South America'
mbpl[data[, 'mbpl'] >= 40000 & data[, 'mbpl'] <= 49900] <- 'Europe'
mbpl[data[, 'mbpl'] >= 50000 & data[, 'mbpl'] <= 59900] <- 'Asia'
mbpl[data[, 'mbpl'] >= 60010 & data[, 'mbpl'] <= 60099] <- 'Africa'
mbpl[data[, 'mbpl'] >= 70010 & data[, 'mbpl'] <= 72000] <- 'Oceania'
mbpl[data[, 'mbpl'] == 96000] <- 'NA'
data[, 'mbpl'] <- mbpl
data <- data[mbpl != 'NA', ]
data[, 'mbpl'] <- as.factor(data[, 'mbpl'])

# Father's birthplace
fbpl <- rep(NA, nrow(data))
fbpl[data[, 'fbpl'] == 9900] <- 'US mainland'
fbpl[data[, 'fbpl'] >= 10000 & data[, 'fbpl'] <= 12090] <- 'US overseas'
fbpl[data[, 'fbpl'] >= 15000 & data[, 'fbpl'] <= 19900] <- 'North America'
fbpl[data[, 'fbpl'] >= 20000 & data[, 'fbpl'] <= 21090] <- 'Central America'
fbpl[data[, 'fbpl'] >= 25000 & data[, 'fbpl'] <= 26091] <- 'Caribbean'
fbpl[data[, 'fbpl'] >= 30005 & data[, 'fbpl'] <= 31000] <- 'South America'
fbpl[data[, 'fbpl'] >= 40000 & data[, 'fbpl'] <= 49900] <- 'Europe'
fbpl[data[, 'fbpl'] >= 50000 & data[, 'fbpl'] <= 59900] <- 'Asia'
fbpl[data[, 'fbpl'] >= 60010 & data[, 'fbpl'] <= 60099] <- 'Africa'
fbpl[data[, 'fbpl'] >= 70010 & data[, 'fbpl'] <= 72000] <- 'Oceania'
fbpl[data[, 'fbpl'] == 96000] <- 'NA'
data[, 'fbpl'] <- fbpl
data <- data[which(fbpl != 'NA'), ]
data[, 'fbpl'] <- as.factor(data[, 'fbpl'])

# Foreign birthplace
nativity <- character(nrow(data))
nativity[data[, 'nativity'] == 0] <- 'NA'
nativity[data[, 'nativity'] == 1] <- 'Both US'
nativity[data[, 'nativity'] == 2] <- 'One US'
nativity[data[, 'nativity'] == 3] <- 'One US'
nativity[data[, 'nativity'] == 4] <- 'Both foreign'
nativity[data[, 'nativity'] == 5] <- 'Foreign born'
data[, 'nativity'] <- nativity
data <- data[nativity != 'NA', ]
data[, 'nativity'] <- as.factor(nativity)

# Hispanic
data[, 'hispanic'] <- ifelse(data[, 'hispan'] != 0, 1, 0)
data[, 'hispan'] <- NULL

# Employment status
employment <- character(nrow(data))
employment[data[, 'empstat'] == 01] <- 'Armed Forces'
employment[data[, 'empstat'] == 10] <- 'At work'
employment[data[, 'empstat'] == 12] <- 'Has job, not at work last week'
employment[data[, 'empstat'] == 20] <- 'Unemployed'
employment[data[, 'empstat'] == 21] <- 'Unemployed'
employment[data[, 'empstat'] == 22] <- 'Unemployed'
employment[data[, 'empstat'] == 30] <- 'NILF'
employment[data[, 'empstat'] == 31] <- 'NILF - housework'
employment[data[, 'empstat'] == 32] <- 'NILF'
employment[data[, 'empstat'] == 33] <- 'NILF - school'
employment[data[, 'empstat'] == 34] <- 'NILF'
employment[data[, 'empstat'] == 35] <- 'NILF - unpaid'
employment[data[, 'empstat'] == 36] <- 'NILF - retired'
data[, 'employment'] <- employment
data <- data[which(data[, 'employment'] != ''), ]
data[, 'employment'] <- as.factor(data[, 'employment'])
data[, 'empstat'] <- NULL

# Labor force
data <- data[which(data[, 'labforce'] == 2), ]  # Only people IN labor force
data[, 'labforce'] <- NULL
data[, 'lfproxy'] <- ifelse(data[, 'lfproxy'] == 2, 1, 0)

# Occupation
occ <- character(nrow(data))
occ[data[, 'occ2010'] < 500] <- 'management'
occ[data[, 'occ2010'] >=  500 & data[, 'occ2010'] <  800] <- 'business operations'
occ[data[, 'occ2010'] >=  800 & data[, 'occ2010'] < 1000] <- 'finance'
occ[data[, 'occ2010'] >= 1000 & data[, 'occ2010'] < 1300] <- 'computer/maths'
occ[data[, 'occ2010'] >= 1300 & data[, 'occ2010'] < 1550] <- 'architect/engineer'
occ[data[, 'occ2010'] >= 1550 & data[, 'occ2010'] < 1600] <- 'technician'
occ[data[, 'occ2010'] >= 1600 & data[, 'occ2010'] < 2000] <- 'science'
occ[data[, 'occ2010'] >= 2000 & data[, 'occ2010'] < 2100] <- 'social service'
occ[data[, 'occ2010'] >= 2100 & data[, 'occ2010'] < 2200] <- 'legal'
occ[data[, 'occ2010'] >= 2200 & data[, 'occ2010'] < 2600] <- 'education'
occ[data[, 'occ2010'] >= 2600 & data[, 'occ2010'] < 3000] <- 'arts/sports/media'
occ[data[, 'occ2010'] >= 3000 & data[, 'occ2010'] < 3700] <- 'health'
occ[data[, 'occ2010'] >= 3700 & data[, 'occ2010'] < 4000] <- 'protective'
occ[data[, 'occ2010'] >= 4000 & data[, 'occ2010'] < 4200] <- 'food'
occ[data[, 'occ2010'] >= 4200 & data[, 'occ2010'] < 4300] <- 'maintenance'
occ[data[, 'occ2010'] >= 4300 & data[, 'occ2010'] < 4700] <- 'personal care'
occ[data[, 'occ2010'] >= 4700 & data[, 'occ2010'] < 5000] <- 'sales'
occ[data[, 'occ2010'] >= 5000 & data[, 'occ2010'] < 6000] <- 'office'
occ[data[, 'occ2010'] >= 6000 & data[, 'occ2010'] < 6200] <- 'farming'
occ[data[, 'occ2010'] >= 6200 & data[, 'occ2010'] < 6800] <- 'construction'
occ[data[, 'occ2010'] >= 6800 & data[, 'occ2010'] < 7000] <- 'extraction'
occ[data[, 'occ2010'] >= 7000 & data[, 'occ2010'] < 7700] <- 'installation'
occ[data[, 'occ2010'] >= 7700 & data[, 'occ2010'] < 9000] <- 'production'
occ[data[, 'occ2010'] >= 9000 & data[, 'occ2010'] < 9800] <- 'transportation'
occ[data[, 'occ2010'] >= 9800 & data[, 'occ2010'] < 9900] <- 'military'
occ[data[, 'occ2010'] == 9920] <- 'unemployed/never worked'
#occ[data[, 'occ2010'] == 9999] <- 'unknown'
#data <- data[occ != '', ]
data[, 'choccly'] <- ifelse(data[, 'occ'] == data[, 'occly'], 0, 1)
data[, 'occ'] <- as.factor(occ)
data[, 'occly'] <- NULL
data[, 'occ2010'] <- NULL
data[, 'occ1990'] <- NULL
data[, 'occ1950'] <- NULL
data[, 'occ10ly'] <- NULL
data[, 'occ90ly'] <- NULL
data[, 'occ50ly'] <- NULL

# Industry (codified as sector)
ind <- character(nrow(data))
ind[data[, 'ind'] <= 770] <- 'primary'
ind[data[, 'ind'] >= 1070 & data[, 'ind'] <= 3990] <- 'secondary'
ind[data[, 'ind'] >= 4070 & data[, 'ind'] <= 9590] <- 'terciary'
ind[data[, 'ind'] == 1190] <- 'terciary'
ind[data[, 'ind'] == 9890] <- 'military'
data[, 'chindly'] <- ifelse(data[, 'ind'] == data[, 'indly'], 0, 1)
data[, 'ind'] <- ind
data <- data[which(ind != 'military'), ]
data[, 'ind'] <- factor(data[, 'ind'] )
data[, 'ind1990'] <- NULL
data[, 'ind1950'] <- NULL
data[, 'indly'] <- NULL
data[, 'ind90ly'] <- NULL
data[, 'ind50ly'] <- NULL

# Class of worker
classworker <- character(nrow(data))
classworker[data[, 'classwkr'] == 0] <- 'NIU'
classworker[data[, 'classwkr'] %in% c(10, 13, 14)] <- 'Self-employed'
classworker[data[, 'classwkr'] %in% c(20, 21, 22, 23, 24)] <- 'Wage/salary'
classworker[data[, 'classwkr'] %in% c(25, 26, 27, 28, 29)] <- 'Government employee'
classworker[data[, 'classwkr'] == 29] <- 'Unpaid family worker'
classwly <- character(nrow(data))
classwly[data[, 'classwly'] == 0] <- 'NIU'
classwly[data[, 'classwly'] %in% c(10, 13, 14)] <- 'Self-employed'
classwly[data[, 'classwly'] %in% c(20, 21, 22, 23, 24)] <- 'Wage/salary'
classwly[data[, 'classwly'] %in% c(25, 26, 27, 28, 29)] <- 'Government employee'
classwly[data[, 'classwly'] == 29] <- 'Unpaid family worker'
data[, 'classworker'] <- factor(classworker)
data[, 'chclasswkr'] <- as.numeric(classworker != classwly)
data[, 'classwkr'] <- NULL
data[, 'classwly'] <- NULL
data <- data[which(data[, 'classworker'] != 'NIU'), ]

# Hours worked
data <- data[which(data[, 'uhrsworkt'] > 0), ]
data <- data[which(data[, 'uhrsworkt'] < 997), ]
data[, 'hoursworked'] <- data[, 'uhrsworkt']
data[, 'pcthrsmainjob'] <- data[, 'uhrswork1'] / data[, 'uhrsworkt']
data[, 'uhrsworkt'] <- NULL
data[, 'uhrswork1'] <- NULL
data[, 'absent'] <- NULL  # Nothing relevant

# Weeks unemployed last year
data[, 'wksunemly'] <- ifelse(data[, 'wksunem1'] == 99, 0, data[, 'wksunem1'])
data[, 'wksunem1'] <- NULL

# Did they look for a job
data[, 'nwlookwkly'] <- ifelse(! data[, 'nwlookwk'] %in% c(0, 99), 1, 0)
data[, 'nwlookwk'] <- NULL

# Work status
wkstat <- character(nrow(data))
wkstat[data[, 'wkstat'] == 10] <- 'Full-time'
wkstat[data[, 'wkstat'] == 11] <- 'Full-time'
wkstat[data[, 'wkstat'] == 12] <- 'Full-time'
wkstat[data[, 'wkstat'] == 13] <- 'Full-time'
wkstat[data[, 'wkstat'] == 14] <- 'Part-time'
wkstat[data[, 'wkstat'] == 15] <- 'Part-time'
wkstat[data[, 'wkstat'] == 20] <- 'Part-time'
wkstat[data[, 'wkstat'] == 21] <- 'Part-time'
wkstat[data[, 'wkstat'] == 22] <- 'Part-time'
wkstat[data[, 'wkstat'] == 40] <- 'Part-time'
wkstat[data[, 'wkstat'] == 41] <- 'Part-time'
wkstat[data[, 'wkstat'] == 42] <- 'Part-time'
wkstat[data[, 'wkstat'] == 50] <- 'Unemployed'
wkstat[data[, 'wkstat'] == 60] <- 'Unemployed'
wkstat[data[, 'wkstat'] == 99] <- 'NIU'
data[, 'wkstat'] <- factor(wkstat)
# data <- data[which(data[, 'wkstat'] != 'Unemployed'), ]
# data <- data[which(data[, 'wkstat'] != 'NIU'), ]

# Usually works full time
data[, 'usftptlw'] <- ifelse(data[, 'usftptlw'] == 2, 1, 0)

# Worked last year
data[, 'workly'] <- ifelse(data[, 'workly'] == 2, 1, 0)

# Education
edu <- character(nrow(data))
edu[data[, 'educ'] < 73] <- 'No HSD'
edu[data[, 'educ'] == 73] <- 'HSD'  # High school diploma
edu[data[, 'educ'] > 73 & data[, 'educ']< 110] <- 'SC'  # Some college (up to 3 years of college)
edu[data[, 'educ'] >= 110 & data[, 'educ'] <= 122] <- 'CG' # College graduate or >3 years of college
edu[data[, 'educ'] %in% 123:125] <- 'Adv'  # MSc / prof. school / PhD
data[, 'edu'] <- edu
data <- data[which(data[, 'edu'] != ''), ]
data[, 'edu'] <- factor(data[, 'edu'])
data[, 'educ'] <- NULL
data[, 'educ99'] <- NULL

# Attending school or college
schlcoll <- character(nrow(data))
schlcoll[data[, 'schlcoll'] == 0] <- 'NIU'
schlcoll[data[, 'schlcoll'] %in% 1:2] <- 'High school'
schlcoll[data[, 'schlcoll'] %in% 3:4] <- 'College or university'
schlcoll[data[, 'schlcoll'] == 5] <- 'Does not attend school'
data[, 'schlcoll'] <- factor(schlcoll)

# Cognitive difficulties
data[, 'diffany'] <- ifelse(data[, 'diffany'] == 2, 1, 0)
data[, 'diffmob'] <- ifelse(data[, 'diffmob'] == 2, 1, 0)
data[, 'diffeye'] <- ifelse(data[, 'diffeye'] == 2, 1, 0)
data[, 'diffrem'] <- ifelse(data[, 'diffrem'] == 2, 1, 0)
data[, 'diffhear'] <- ifelse(data[, 'diffhear'] == 2, 1, 0)
data[, 'diffphys'] <- ifelse(data[, 'diffphys'] == 2, 1, 0)
data[, 'diffcare'] <- ifelse(data[, 'diffcare'] == 2, 1, 0)

# Pension plan
data[, 'pension'] <- ifelse(data[, 'pension'] == 3, 1, 0)

# Firm size
firmsize <- character(nrow(data))
firmsize[data[, 'firmsize'] == 0] <- '0'
firmsize[data[, 'firmsize'] >  0 & data[, 'firmsize'] <= 3] <- '1 (1-24)'
firmsize[data[, 'firmsize'] >= 5 & data[, 'firmsize'] <= 6] <- '2 (25-99)'
firmsize[data[, 'firmsize'] == 7] <- '3 (100-499)'
firmsize[data[, 'firmsize'] == 8] <- '4 (500-999)'
firmsize[data[, 'firmsize'] == 9] <- '5 (>1000)'
data[, 'firmsize'] <- firmsize
data <- data[which(data[, 'firmsize'] != '0'), ]
data[, 'firmsize'] <- factor(data[, 'firmsize'])

# Receives social security income
data[, 'ssocinc'] <- ifelse(data[, 'whyss1'] != 0, 1, 0)
data[, 'whyss1'] <- NULL
data[, 'whyss2'] <- NULL

# Tax filer status
taxstatus <- character(nrow(data))
taxstatus[which(data[, 'filestat'] %in% 1:3)] <- 'joint'
taxstatus[which(data[, 'filestat'] == 4)] <- 'headhh'
taxstatus[which(data[, 'filestat'] == 5)] <- 'single'
taxstatus[which(data[, 'filestat'] == 6)] <- 'nonfiler'
data[, 'taxstatus'] <- factor(taxstatus)
data[, 'filestat'] <- NULL

# Mortgage
data[, 'mortgage'] <- ifelse(data[, 'spmmort'] == 1, 1, 0)
data[, 'spmmort'] <- NULL

# State of residence (1 year ago)
movedstate <- character(nrow(data))
movedstate[data[, 'migsta1'] == 0] <- 'NIU'
movedstate[(data[, 'migsta1'] >= 1) & (data[, 'migsta1'] <= 60)] <- 'From US'
movedstate[data[, 'migsta1'] == 91] <- 'From out US'
movedstate[data[, 'migsta1'] == 99] <- 'No'
data[, 'movedstate'] <- factor(movedstate)
data[, 'migsta1'] <- NULL

# Place of residence (1 year ago)
movedhouse <- character(nrow(data))
movedhouse[data[, 'migrate1'] == 0] <- 'NIU'
movedhouse[data[, 'migrate1'] == 1] <- 'Same house'
movedhouse[data[, 'migrate1'] == 3] <- 'Moved within county'
movedhouse[data[, 'migrate1'] == 4] <- 'Moved within state'
movedhouse[data[, 'migrate1'] == 5] <- 'Moved between states'
movedhouse[data[, 'migrate1'] == 6] <- 'Abroad'
data[, 'movedhouse'] <- factor(movedhouse)
data[, 'migrate1'] <- NULL

# Work disability
data[, 'disabwrk'] <- ifelse(data[, 'disabwrk'] == 2, 1, 0)
data[, 'quitsick'] <- ifelse(data[, 'quitsick'] == 2, 1, 0)

# Health status
health <- character(nrow(data))
health[data[, 'health'] == 1] <- 'excellent'
health[data[, 'health'] == 2] <- 'very good'
health[data[, 'health'] == 3] <- 'good'
health[data[, 'health'] == 4] <- 'fair'
health[data[, 'health'] == 5] <- 'poor'
data[, 'health'] <- factor(health)

# Employer pays group health plan
emphplan <- character(nrow(data))
emphplan[data[, 'paidgh'] %in% c(0, 10)] <- 'no'
emphplan[data[, 'paidgh'] == 21] <- 'partly'
emphplan[data[, 'paidgh'] == 22] <- 'yes'
data[, 'emphplan'] <- factor(emphplan)
data[, 'paidgh'] <- NULL

# Health Insurance
data[, 'medicare'] <- ifelse(data[, 'himcaidly'] == 2, 1, 0)
data[, 'medicaid'] <- ifelse(data[, 'himcarely'] == 2, 1, 0)
data[, 'milhealth'] <- ifelse(data[, 'hichamp'] == 2, 1, 0)
data[, 'phinsur'] <- ifelse(data[, 'phinsur'] == 2, 1, 0)
data[, 'phiown'] <- ifelse(data[, 'phiown'] == 2, 1, 0)
data[, 'out'] <- ifelse(data[, 'out'] == 2, 1, 0)
data[, 'himcaidly'] <- NULL
data[, 'himcarely'] <- NULL
data[, 'hichamp'] <- NULL

# Health insurance of family members
data[, 'grpdeply'] <- ifelse(data[, 'grpdeply'] == 2, 1, 0)
data[, 'grpownly'] <- ifelse(data[, 'grpownly'] == 2, 1, 0)
data[, 'grpoutly'] <- ifelse(data[, 'grpoutly'] == 2, 1, 0)

# Type of coverage
grptyply <- character(nrow(data))
grptyply[data[, 'grptyply'] == 0] <- 'none'
grptyply[data[, 'grptyply'] == 1] <- 'family'
grptyply[data[, 'grptyply'] == 2] <- 'self'
grptyply[data[, 'grptyply'] == 3] <- 'plusone'
data[, 'grptyply'] <- factor(grptyply)

# Purchase insurance
data[, 'dpdeply'] <- ifelse(data[, 'dpdeply'] == 2, 1, 0)
data[, 'dpownly'] <- ifelse(data[, 'dpownly'] == 2, 1, 0)
data[, 'dpoutly'] <- ifelse(data[, 'dpoutly'] == 2, 1, 0)

# Type of coverage
dptyply <- character(nrow(data))
dptyply[data[, 'dptyply'] == 0] <- 'none'
dptyply[data[, 'dptyply'] == 1] <- 'family'
dptyply[data[, 'dptyply'] == 2] <- 'self'
dptyply[data[, 'dptyply'] == 3] <- 'plusone'
data[, 'dptyply'] <- factor(dptyply)

# # Veteran status
# data[, 'veteran'] <- ifelse(data[, 'vetlast'] <= 1, 0, 1)
data[, 'vetlast'] <- NULL

# Paid by the hour
data[, 'paidhour'] <- ifelse(data[, 'paidhour'] == 2, 0, 1)

# Union membership/coverage
data[, 'union'] <- ifelse(data[, 'union'] >= 2, 1, 0)

# Number of employers last year
data[, 'numempsly'] <- data[, 'numemps']
data[, 'numemps'] <- NULL

# Different income columns
inc.cols <- c('inctot', 'incwage', 'incbus', 'incfarm', 'incss', 'incwelfr',
  'incretir', 'incssi', 'incint', 'incunemp', 'incwkcom', 'incvet', 'incsurv',
  'incdisab', 'incdivid', 'incrent', 'inceduc', 'incchild', 'incasist',
  'incother', 'inclongj', 'oincbus', 'oincfarm', 'oincwage', 'incsurv1',
  'incsurv2', 'incdisa1', 'incdisa2', 'ftotval')
for (icol in inc.cols) {
  invalid <- c(99999998, 99999999)
  data[, icol] <- as.numeric(data[, icol])
  data[which(data[, icol] %in% invalid), icol] <- 0
  #data <- data[which(! data[, icol] %in% invalid), ]
}

# Renaming of some of them
data[, 'famincome'] <- data[, 'ftotval']  # Total family income
data[, 'ftotval'] <- NULL
data[, 'personalincome'] <- data[, 'inctot'] # Total personal income
data[, 'inctot'] <- NULL
data[, 'incomewage'] <- data[, 'incwage']
data[, 'incwage'] <- NULL

# Other numeric variables
data[, 'ctccrd'] <- as.numeric(data[, 'ctccrd'] + data[, 'actccrd'])
data[, 'actccrd'] <- NULL
data[, 'adjginc'] <- as.numeric(data[, 'adjginc'])
data[, 'eitcred'] <- as.numeric(data[, 'eitcred'])
data[, 'fedtaxac'] <- as.numeric(data[, 'fedtaxac'])
data[, 'fica'] <- as.numeric(data[, 'fica'])
data[, 'margtax'] <- as.numeric(data[, 'margtax'] / 10)
data[, 'stataxac'] <- as.numeric(data[, 'stataxac'])
data[, 'taxinc'] <- as.numeric(data[, 'taxinc'])
data[, 'foodstamp'] <- as.numeric(data[, 'foodstamp'])
data[, 'schllunch'] <- as.numeric(data[, 'schllunch'])

# Suppress some variables from Poverty Supplementary Questionnaire (SPM)
spm.in <- c('spmlunch', 'spmmedxpns', 'spmchxpns', 'spmcapxpns', 'spmchsup',
  'spmnadults', 'spmnchild', 'spmnpers', 'spmwkxpns')
spm.out <- c('spmcaphous', 'spmeqscale', 'spmwt', 'spmsttax', 'spmfedtaxac',
  'spmeitc', 'spmfica', 'spmfedtaxbc', 'spmwic', 'spmheat', 'spmsnap',
  'spmftotval', 'spmtotres', 'spmgeoadj', 'spmthresh', 'spmfamunit')
data[, spm.out] <- NULL
for (incol in spm.in) {
  data[, incol] <- as.numeric(data[, incol])
}

################################################################################
# Format final output
################################################################################
# Select valid cases
data <- data[which((data[, 'age'] >= 16) & (data[, 'age'] <= 65)), ]
data <- data[which(data[, 'employment'] == 'At work'), ]
data <- data[which(data[, 'classworker'] != 'NIU'), ]
data <- data[which(data[, 'nativity'] != 'NA'), ]

# Selected variables in order
var.order <- c('householdid', 'cpsid', 'cpsidp', 'oversample', 'numprec',
  'hhrespln', 'pernum', 'lineno', 'month', 'year', 'region', 'state',
  'metro', 'areasize', 'age', 'female', 'race', 'marital', 'edu', 'schlcoll',
  'hispanic', 'bpl', 'timeus', 'citizen', 'mbpl', 'fbpl', 'nativity',
  'movedstate', 'movedhouse', 'classworker', 'chclasswkr', 'employment',
  'lfproxy', 'occ', 'choccly', 'ind', 'chindly', 'firmsize', 'hoursworked',
  'wkstat', 'pcthrsmainjob', 'wkswork1', 'usftptlw', 'workly', 'uhrsworkly',
  'wksunemly', 'nwlookwkly', 'strechlk', 'ptweeks', 'paidhour', 'union',
  'numempsly', 'ssocinc', 'pension', 'emphplan', 'veteran', 'medicare',
  'medicaid', 'milhealth', 'hiunpers', 'phinsur', 'phiown', 'out', 'grpdeply',
  'grpownly', 'grpoutly', 'grptyply', 'grpwho1', 'dpdeply', 'dpownly',
  'dpoutly', 'dptyply', 'dpwho1', 'health', 'diffany', 'diffmob', 'diffeye',
  'diffrem', 'diffhear', 'diffphys', 'diffcare', 'disabwrk', 'quitsick',
  'taxstatus', 'ctccrd', 'adjginc', 'eitcred', 'fedtaxac', 'fica', 'margtax',
  'stataxac', 'taxinc', 'hhtenure', 'buildtype', 'mortgage', 'phone', 'nfams',
  'ncouples', 'nmothers', 'nfathers', 'famsize', 'nchild', 'nchlt5', 'eldch',
  'yngch', 'famunit', 'nsibs', 'ftype', 'fkind', 'famrel', 'relhh', 'isrelhh',
  'pubhous', 'rentsub', 'heatsub', 'lunchsub', 'atelunch', 'frelunch',
  'stampno', 'stampmo', 'foodstmp', 'foodstamp', 'schllunch', 'spmlunch',
  'spmmedxpns', 'spmchxpns', 'spmcapxpns', 'spmchsup', 'spmnadults',
  'spmnchild', 'spmnpers', 'spmwkxpns', 'famincome','personalincome',
  'incomewage', 'incbus', 'incfarm', 'incss', 'incwelfr', 'incretir', 'incssi',
  'incint', 'incunemp', 'incwkcom', 'incvet', 'incsurv', 'incdisab', 'incdivid',
  'incrent', 'inceduc', 'incchild', 'incasist', 'incother', 'inclongj',
  'oincbus', 'oincfarm', 'oincwage', 'incsurv1', 'incsurv2', 'incdisa1',
  'incdisa2', 'cpi99')

# Check variable list
chck1 <- colnames(data)[! colnames(data) %in% var.order]
chck2 <- var.order[! var.order %in% colnames(data)]
if (length(chck1 > 1) | length(chck2 > 1)) {
  stop('number of variables does not check out.')
}

# Save data
file.out <- paste(INPDIR, 'cps2year/salary_2010to2019.RData', sep = '')
salary <- as.data.frame(data[, var.order])
save(salary, file = file.out); cat('Saved file:', file.out, '\n')
rm(file.out, data, raw); gc()
# END OF SCRIPT
