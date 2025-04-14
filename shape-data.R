# Shaping data for parametric g-formula and the sdr lmtp ####
# November 3, 2020

library(here)
library(boxr)
library(survival)
box_auth()

# Get raw data and helper functions ####
exposure.lag <- 20
race.imputed.melt <- NULL
source("~/headRs/00-my-theme.R")
library(data.table)
# # paste(ls(), collapse = ",\n") %>% clipr::write_clip()
source(here::here("breaks.R"))

# Get data ####
source(here::here("get-data.R"))
# Clear out unnecessary gunk
rm(
	additional_outcomes,
	dta,
	exposure,
	get.coef,
	get.cohort2,
	get.exposure,
	get.facet_tikz,
	get.ggtab,
	get.hwse.ggtab,
	get.hwse2.coxph,
	get.hwse3.coxph,
	get.jobhist,
	get.tikz,
	icd_codes.function,
	jobhist,
	jobhist_py,
	jobhist_py.cast,
	ltab_age,
	ltab_calendar,
	seer,
	self_injury.function,
	spec_icd_codes
)

# Which outcome ?
outcomes.which <- grep("hodgkin", incidence.key$description, ignore.case = T)
incidence.key[outcomes.which,]

# Remove people who experienced event in first year
dat <- cohort_analytic[!studyno %in% c(103598, 132387)]
# dat <- copy(cohort_analytic)

# Lag time off
dat[,`:=`(
	off = shift(off, exposure.lag, 0),
	cum_off = shift(cum_off, exposure.lag, 0)
)]
dat <- dat[yin >= as.Date("1985-01-01") - exposure.lag * 365.25]
# Get data analytic ####
get.coxph(
	cohort_name = "dat",
	run_model = F,
	time_scale = "age",
	start.year = -Inf,
	employment_status.lag = exposure.lag)
dat <- copy(nhl.dat)

# Getting histology codes from cancer data
cancer <- box_read('512727972584')
setDT(cancer)

nhl_morphology <- merge(
	dat[status == 1, .(studyno, yoi)],
	unique(cancer[morph %in% c(
		9590:9597, 9670:9671, 9673, 9675, 9678:9680, 9684, 9687:9691,
		9695, 9698:9702, 9705, 9708:9709, 9712, 9714:9719, 9724:9729,
		9735, 9737:9738, 9811:9818, 9823, 9827, 9837, 9590:9597,
		9670:9671, 9673, 9675, 9678:9680, 9684, 9687, 9688, 9689:9691,
		9695, 9698:9702, 9705, 9708:9709, 9712, 9714:9719, 9724:9729,
		9735, 9737, 9738, 9811:9818, 9823, 9827, 9837
	),.(studyno, yoi = ddiag, morph)]),
	by = c("studyno", "yoi"))

dat <- merge(
	dat, nhl_morphology[,.(studyno, morph = as.integer(morph))],
	by = "studyno",
	all.x = T)

# TV covariates: calendar year, duration of employment, time spent off, cumulative MWF exposure up to ___ year, employment status, censoring/death
# Baseline covariates: year of hire, age at hire, sex, race, plant
# Outcome: NHL

# Index rows
setorder(dat, studyno, year)
dat[immortal == 0 & year >= 1985,`:=`(
	I = 1:.N,
	N = .N), by = .(studyno)]

# Some covariate engineering
dat[, `:=`(
	yob.gm = date.to.gm(yob),
	yin.gm = date.to.gm(yin),
	jobloss.date.gm = date.to.gm(jobloss.date),
	yod.gm = date.to.gm(yod)
)]
dat[,`:=`(
	Plant = which.max(table(Plant))
), by = .(studyno)]
dat[,`:=`(
	`age_at_hire` = yin.gm - yob.gm
)]
dat[I == N & status == 0 & year == floor(yod.gm),
		`Other event` := 1]
dat[is.na(`Other event`),
		`Other event` := 0]
# Lag employment status and duration of employment ####
dat[, `Employment status` := ifelse(
	year > floor(jobloss.date.gm) + exposure.lag, "Left work", "At work")]
dat[, employment.years := apply(
	data.frame(jobloss.date.gm, year + 1), 1, min) - yin.gm, studyno]
dat[, employment.years := shift(
	employment.years, exposure.lag, 0), studyno]

covariates_to_lag <- c(
	"Employment status", "employment.years"
)
dat[,covariates_to_lag := lapply(covariates_to_lag, function(x) {
	shift(get(x), exposure.lag, fill = get(x)[1])
}), studyno]

dat[I == N & status == 0 & (I < 1994 + exposure.lag - 1985 + 1), Censored := 1]
dat[is.na(Censored), Censored := 0]

mwf.breaks <- get.mwf.breaks(dat)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save to Box                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box_write(dat,
					"dat.rds",
					dir_id = 125518026124,
					description = paste0("Output of get.coxph(), lag of ", exposure.lag))

dat <- dat[yin >= as.Date("1938-01-01") &
					 	year <= exposure.lag + 1994 &
					 	yin.gm < 1982 &
					 	jobloss.date.gm - yin.gm >= 3,]
sum(dat$status)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read from Box                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- box_read(737540256124)

# # Exposure over time ####
# exposure.ggplot <- dat[, .(
# 	Soluble   = median(  soluble[soluble > 0]),
# 	Straight  = median( straight[straight > 0]),
# 	Synthetic = median(synthetic[synthetic > 0])
# ), .(year = year - exposure.lag)][year >= 1941] %>% melt(
# 	id.vars = "year",
# 	variable.name = "mwf") %>%
# 	ggplot(aes(x = year, y = value)) +
# 	geom_line()
#
# tikzDevice::tikz(here::here("resources", "exposure.tex"), height = 2.5, width = 6, standAlone = T)
# exposure.ggplot +
# 	facet_wrap(. ~ mwf) +
# 	labs(x = "Calendar year",
# 			 y = "Median exposure (mg/m$^3$) among exposed") +
# 	coord_cartesian(xlim = c(1938, 2000)) +
# 	mytheme
# dev.off()
# lualatex("^exposure\\.tex", here::here("resources"))
#
# tikzDevice::tikz(here::here("resources", "exposure_long.tex"), height = 5.5, width = 3, standAlone = T)
# exposure.ggplot +
# 	facet_wrap(. ~ mwf, nrow = 3) +
# 	labs(x = "Calendar year",
# 			 y = "Median exposure (mg/m$^3$) among exposed") +
# 	coord_cartesian(xlim = c(1917, 2005)) +
# 	mytheme
# dev.off()
# lualatex("^exposure_long\\.tex", here::here("resources"))

# # Population characteristics ####
# get.tab1 <- function(
		# 		df = dat[year >= 1985 & year <= 1994 + exposure.lag],
# 		table.engine = "pander",
# 		use_finrace = T,
# 		incidence = F,
# 		py = "py",
# 		mathmode = F,
# 		nrow_as_fu = T) {
#
# 	setorder(df, studyno, year)
# 	# Individual-level summary
# 	tab1.sum <- df[,.(
# 		'Race' = NA,
# 		'\\hspace{10pt}White' = if (use_finrace) {
# 			ifelse(finrace[1] == 1, 1, 0)} else {
# 				as.numeric(race[1] == "White")
# 			},
# 		'\\hspace{10pt}Black' = if (use_finrace) {
# 			ifelse(finrace[1] == 2, 1, 0)} else {
# 				as.numeric(race[1] == "Black")
# 			},
# 		'\\hspace{10pt}Unknown' = if (use_finrace) {
# 			ifelse(finrace[1] %in% c(0, 9), 1, 0)} else {
# 				as.numeric(race[1] == "Unknown")
# 			},
# 		'Sex' = NA,
# 		'\\hspace{10pt}Male' = ifelse(sex[1] == 'M', 1, 0),
# 		'\\hspace{10pt}Female' = ifelse(sex[1] == 'F', 1, 0),
# 		"Plant$^1$" = NA,
# 		"\\hspace{10pt}Plant 1" = as.numeric(Plant[1]) == 1,
# 		"\\hspace{10pt}Plant 2" = as.numeric(Plant[1]) == 2,
# 		"\\hspace{10pt}Plant 3" = as.numeric(Plant[1]) == 3,
# 		"Ever exposed to MWFs" = NA,
# 		# as.numeric(cum_soluble[.N] + cum_straight[.N] + cum_synthetic[.N] > 0)
# 		"\\hspace{10pt}Straight"  = as.numeric(cum_straight[.N] > 0),
# 		"\\hspace{10pt}Soluble"   = as.numeric(cum_soluble[.N] > 0),
# 		"\\hspace{10pt}Synthetic" = as.numeric(cum_synthetic[.N] > 0),
# 		# "Employment status in 1995" = NA,
# 		# "\\hspace{10pt}Left work" = as.numeric(jobloss.date[1] < as.Date("1995-01-01")),
# 		# "\\hspace{10pt}Still at work" = as.numeric(jobloss.date[1] >= as.Date("1995-01-01")),
# 		"Deceased by end of follow-up" = {
# 			if (!incidence) {
# 				as.numeric(!is.na(yod[1]) & yod[1] <= as.Date(paste0(year[.N], "-12-31")))
# 			} else {
# 				as.numeric(!is.na(ddiag_first[1]) & ddiag_first[1] <= as.Date(paste0(year[.N], "-12-31")))
# 			}
# 		},
# 		# "\\hspace{10pt}ALD" = as.numeric(max(`Alcohol-related Liver Disease`) > 0),
# 		# "\\hspace{10pt}Suicide" = as.numeric(max(Suicide) > 0),
# 		# "\\hspace{10pt}Overdose" = as.numeric(max(Overdose) > 0),
#
#
# 		# Years since follow-up
# 		"\\hline Years of follow-up" = if (nrow_as_fu) {.N} else {
# 			as.numeric(difftime(min(yod[1], yoc[1], as.Date(paste0(max(year), "-12-31")), na.rm = T),
# 													as.Date(paste0(min(year[immortal == 0]), "-01-01")), units = "days"))/365
# 		},
# 		'Year of birth' = yob.gm[1],
# 		'Year of hire' = yin.gm[1],
# 		'Age at hire (years)' = yin.gm[1] - yob.gm[1],
# 		'Year of leaving work$^2$' = {
# 			if (year(jobloss.date[1]) == 1995) {NaN} else {jobloss.date.gm[1]}},
# 		'Age at leaving work (years)$^2$' = {
# 			if (year(jobloss.date[1]) == 1995) {NaN} else {
# 				jobloss.date.gm[1] - yob.gm[1]
# 			}},
# 		"Years at work$^2$" = {
# 			if (year(jobloss.date[1]) == 1995) {
# 				NaN} else {as.numeric(difftime(jobloss.date[1], yin[1], units = "days"))/365}},
# 		'Year of death among deceased' = {
# 			if (!incidence) {
# 				as.numeric(year(yod[1]))} else {
# 					as.numeric(year(ddiag_first[1]))
# 				}},
# 		'Age at death (years) among deceased' = {
# 			if (!incidence) {
# 				as.numeric(difftime(yod[1], yob[1], units = "days"))/365} else {
# 					as.numeric(difftime(ddiag_first[1], yob[1], units = "days"))/365
# 				}},
# 		# "Maximum cumulative exposure$^3$ (mg/m$^3\\cdot$y)}" =
# 		# 	NA,
# 		# 	max(cum_soluble[cum_soluble != 0]) +
# 		# 	  max(cum_straight[cum_straight != 0]) +
# 		#   	max(cum_synthetic[cum_synthetic != 0])
#
# 		'Cumulative time off (years)' = {
# 			max(cum_off)
# 		},
#
# 		"Cumulative exposure$^3$ to MWFs (mg/m$^3\\cdot$y)" = NA,
#
# 		"\\hspace{10pt}Straight "  = as.numeric(ifelse(
# 			sum(cum_straight != 0) > 0,
# 			max(cum_straight[cum_straight != 0], na.rm = T),
# 			NA)),
#
# 		"\\hspace{10pt}Soluble "   = as.numeric(ifelse(
# 			sum(cum_soluble != 0) > 0,
# 			max(cum_soluble[cum_soluble != 0], na.rm = T),
# 			NA)),
#
# 		"\\hspace{10pt}Synthetic " = as.numeric(ifelse(
# 			sum(cum_synthetic != 0) > 0,
# 			max(cum_synthetic[cum_synthetic != 0], na.rm = T),
# 			NA))
#
# 	),
# 	by = .(studyno)][,-'studyno']
#
# 	# Correct names if looking at incidence outcome
# 	if (incidence) {
# 		names(tab1.sum) <- gsub(
# 			"Deceased by end of follow-up",
# 			"Diagnosed with cancer by end of follow-up",
# 			names(tab1.sum))
# 		names(tab1.sum) <- gsub(
# 			"Year of death among deceased",
# 			"Year of first cancer diagnosis",
# 			names(tab1.sum))
# 		names(tab1.sum) <- gsub(
# 			"Age at death \\(years\\) among deceased",
# 			"Age at first cancer diagnosis (years)",
# 			names(tab1.sum))
# 	}
#
# 	# Population-level summary
# 	table.break <- "\\hline Years of follow-"
# 	tab1 <- rbind(
# 		t(apply(tab1.sum[,1:(grep(table.break, names(tab1.sum)) - 1)], 2, function(x) {
# 			return(
# 				c(mean(x, na.rm = T) * sum(!is.na(x)),
# 					mean(x, na.rm = T) * 100,
# 					NA
# 					# NA,
# 					# NA
# 				)
# 			)
# 		})),
# 		t(apply(tab1.sum[,grep(table.break, names(tab1.sum)):ncol(tab1.sum)], 2,
# 						function(x) {
# 							return(
# 								c(
# 									# mean(x, na.rm = T),
# 									# sd(x, na.rm = T)
# 									# min(x, na.rm = T),
# 									median(as.numeric(x), na.rm = T),
# 									quantile(as.numeric(x), 0.25, na.rm = T),
# 									quantile(as.numeric(x), 0.75, na.rm = T)
# 									# max(x, na.rm = T)
# 								)
# 							)
# 						}))
# 	)
#
# 	tab1[!is.finite(as.matrix(tab1))] <- NA
#
# 	# Table names
# 	colnames(tab1) <- c('center', 'spread', " ")#, 'Minimum', 'Median', 'Maximum')
#
# 	# Digits
# 	tab1.digits <- matrix(2, ncol = 3,#5
# 												nrow = nrow(tab1))
# 	tab1.digits[grep("Age", rownames(tab1), ignore.case = T), ] <- 1
# 	tab1.digits[grep("years ", rownames(tab1), ignore.case = T), ] <- 1
# 	tab1.digits[grep("year ", rownames(tab1), ignore.case = T), ] <- 0
# 	tab1.digits[1:(grep(table.break, rownames(tab1), ignore.case = T) - 1), 1] <- 0
# 	tab1.digits[1:(grep(table.break, rownames(tab1), ignore.case = T) - 1), -1] <- 0
#
# 	# more digits if small %s
# 	tab1.digits[1:nrow(tab1.digits) < grep(table.break, rownames(tab1), ignore.case = T) &
# 								(tab1[,2] < 2 & tab1[,2] != 0), -1] <- 2
#
# 	tab1.digits[1:nrow(tab1.digits) < grep(table.break, rownames(tab1), ignore.case = T) &
# 								(tab1[,2] < 1 & tab1[,2] != 0), -1] <- 2
#
# 	tab1.digits[1:nrow(tab1.digits) < grep(table.break, rownames(tab1), ignore.case = T) &
# 								(tab1[,2] > 99 & tab1[,2] != 100), -1] <- 1
#
# 	# more digits if numbers
# 	tab1.digits[1:nrow(tab1.digits) >= grep(table.break, rownames(tab1), ignore.case = T) &
# 								(tab1[,2] < .1 & tab1[,2] != 0), -1] <- 2
#
# 	tab1.digits[1:nrow(tab1.digits) >= grep(table.break, rownames(tab1), ignore.case = T) &
# 								(tab1[,2] == 0), -1] <- 0
#
# 	tab1.which_year <- matrix(F, ncol = 3, nrow = nrow(tab1))
# 	tab1.which_year[grep("year ", rownames(tab1), ignore.case = T), ] <- T
#
# 	tab1 <- matrix(
# 		sapply(1:length(tab1), function(i) {
# 			if (is.na(as.vector(tab1)[i])) {NA} else {
# 				formatC(as.vector(tab1)[i], as.vector(tab1.digits)[i], format = "f",
# 								big.mark = if (!i %in% which(as.vector(tab1.which_year))) {
# 									if (mathmode) {"\\\\,"} else {","}} else {""})
# 			}}),
# 		ncol = ncol(tab1),
# 		nrow = nrow(tab1),
# 		dimnames = dimnames(tab1)
# 	)
#
# 	# Pretty percents
# 	tab1[1:(grep("Years of fo", rownames(tab1)) - 1), 2] <- sapply(
# 		tab1[1:(grep("Years of fo", rownames(tab1)) - 1), 2],
# 		function (x) {
# 			if (!is.na(x)) {
# 				paste0(x, ifelse(table.engine == "xtable", '\\%', "%"))
# 			} else {NA}
# 		})
#
# 	# Math mode
# 	if (mathmode) {
# 		tab1 <- matrix(
# 			sapply(1:length(tab1), function(i) {
# 				if (!is.na(as.vector(tab1)[i])) {
# 					paste0('$', as.vector(tab1)[i], '$')
# 				} else {NA}
# 			}),
# 			ncol = ncol(tab1),
# 			nrow = nrow(tab1),
# 			dimnames = dimnames(tab1))
# 	}
#
# 	tab1 <- rbind(
# 		matrix(c(
# 			paste0(
# 				if (mathmode) {"$"} else {""},
# 				prettyNum(n_distinct(df$studyno),
# 									ifelse(mathmode, '\\\\,', ",")),
# 				if (mathmode) {"$"} else {""}),
# 			paste0(
# 				if (mathmode) {"$"} else {""},
# 				formatC(ifelse(nrow_as_fu, nrow(df), sum(df[,py])),
# 								format = "f", digits = 0,
# 								big.mark = if (mathmode) {"\\\\,"} else {","}),
# 				if (mathmode) {"$"} else {""}),
# 			NA
# 			# , NA, NA
# 		),
# 		nrow = 1,
# 		dimnames = list(c("Study population size (person-years)"), NULL)),
# 		tab1)
#
# 	# Make column indicating stat type
# 	tab1 <- cbind(tab1, spread.which = c(
# 		rep(paste0(if (mathmode) {"$"} else {""}, "n", if (mathmode) {"$"} else {""}, ifelse(table.engine == "xtable", "\\%", "%")),
# 				grep("\\hline", rownames(tab1)) - 1),
# 		rep("median, Q1, Q3", nrow(tab1) - grep("\\hline", rownames(tab1)) + 1))
# 	)
#
# 	# Clean up rownames
# 	rownames(tab1)[duplicated(rownames(tab1))] <- paste0(rownames(tab1)[duplicated(rownames(tab1))], " ")
# 	tab1 <- as.data.frame(tab1, make.names = F)
#
# 	print(matrix(as.vector(sapply(tab1[,1:3], function(x) gsub("\\$|\\\\", "", x))),
# 							 ncol = 3,
# 							 dimnames = list(rownames(tab1), colnames(tab1)[1:3])))
#
# 	return(tab1)
# }
# tab1 <- get.tab1(dat[year >= 1985 & year <= 1994 + exposure.lag], mathmode = F, nrow_as_fu = T)
# saveRDS(tab1, here::here("resources", "tab1.rds"))
# nhl.tab1 <- get.tab1(dat[
# 	studyno %in% dat[status == 1 & year <= 1994 + exposure.lag, studyno] &
# 		year >= 1985 & year <= 1994 + exposure.lag],
# 	mathmode = F, nrow_as_fu = T)
# saveRDS(nhl.tab1, here::here("resources", "nhl.tab1.rds"))
#
# clean.tab1 <- function(tab1) {
# 	tab1[grepl("^n", tab1[,4]) & !is.na(tab1[,1]), 2] <- paste0(
# 		"(", tab1[grepl("^n", tab1[,4]) & !is.na(tab1[,1]), 2], ")"
# 	)
#
# 	tab1[grepl("median", tab1[,4]) & !is.na(tab1[,1]), 2] <- paste0(
# 		"(", tab1[grepl("median", tab1[,4]) & !is.na(tab1[,1]), 2], ", ",
# 		tab1[grepl("median", tab1[,4]) & !is.na(tab1[,1]), 3],
# 		")"
# 	)
# 	return(tab1[,1:2])
# }
# tab1 <- readRDS(here::here("resources", "tab1.rds"))
# nhl.tab1 <- readRDS(here::here("resources", "nhl.tab1.rds"))
#
# clipr::write_clip(
# 	kable(cbind(clean.tab1(tab1), NA, clean.tab1(nhl.tab1)))
# )
# clipr::write_clip(
# 	print(xtable(cbind(
# 		rownames(clean.tab1(tab1)),
# 		mutate(clean.tab1(tab1), spread = sanitize(spread)), NA,
# 		mutate(clean.tab1(nhl.tab1), spread = sanitize(spread))
# 	)))
# )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Refresh categorical coding                     ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat[finrace == 9, Race := "Unknown"]

setorder(dat, studyno, year)
dat[,I := NA]
dat[immortal == 0 & year >= 1985, I := 1:.N, studyno]

mwf.breaks[,cum_soluble := cum_soluble5]

soluble.cutpoints <- quantile(dat[status == 1 & soluble > 0, soluble], seq(0, 1, length.out = 10))
soluble.cutpoints[c(length(soluble.cutpoints))] <- max(dat$soluble, na.rm = T)
soluble.cutpoints <- unique(c(-Inf, 0, soluble.cutpoints))

cum_soluble.cutpoints <- quantile(dat[status == 1 & cum_soluble > 0.05, cum_soluble], seq(0, 1, length.out = 10))
cum_soluble.cutpoints[c(length(cum_soluble.cutpoints))] <- max(dat$cum_soluble, na.rm = T)
cum_soluble.cutpoints <- unique(c(-Inf, 0.05, cum_soluble.cutpoints))

cum_straight.cutpoints <- quantile(dat[status == 1 & cum_straight > 0, cum_straight], seq(0, 1, length.out = 10))
cum_straight.cutpoints[c(length(cum_straight.cutpoints))] <- max(dat$cum_straight, na.rm = T)
cum_straight.cutpoints <- unique(c(-Inf, 0, cum_straight.cutpoints))

cum_synthetic.cutpoints <- quantile(dat[status == 1 & cum_synthetic > 0, cum_synthetic], seq(0, 1, length.out = 10))
cum_synthetic.cutpoints[c(length(cum_synthetic.cutpoints))] <- max(dat$cum_synthetic, na.rm = T)
cum_synthetic.cutpoints <- unique(c(-Inf, 0, cum_synthetic.cutpoints))

cum_off.cutpoints <- quantile(dat[status == 1 & cum_off > 0, cum_off], seq(0, 1, length.out = 3))
cum_off.cutpoints[c(length(cum_off.cutpoints))] <- max(dat$cum_off, na.rm = T)
cum_off.cutpoints <- unique(c(-Inf, 0, cum_off.cutpoints))

age_hire.cutpoints <- quantile(dat[status == 1, yin.gm - yob.gm], seq(0, 1, length.out = 3))
age_hire.cutpoints[c(length(age_hire.cutpoints))] <- max(dat[, yin.gm - yob.gm], na.rm = T)
age_hire.cutpoints <- unique(c(-Inf, 0, age_hire.cutpoints))

# years_since_hire.cutpoints <- quantile(dat[status == 1, employment.years], seq(0, 1, length.out = 20))
years_since_hire.cutpoints <- seq(max(dat$employment.years), max(dat$employment.years), 5)
years_since_hire.cutpoints[c(length(years_since_hire.cutpoints))] <- max(dat$employment.years)
years_since_hire.cutpoints <- unique(c(years_since_hire.cutpoints))

dat[,`:=`(
	`Soluble` = cut(soluble, soluble.cutpoints, include.lowest = T, dig.lab = 2),
	`Cumulative soluble` = cut(cum_soluble, cum_soluble.cutpoints, include.lowest = T, dig.lab = 2),
	# `Cumulative straight` = cut(cum_straight, cum_straight.cutpoints, include.lowest = T, dig.lab = 2),
	# `Cumulative synthetic` = cut(cum_synthetic, cum_synthetic.cutpoints, include.lowest = T, dig.lab = 2),
	`Cumulative time off` =  cut(cum_off, cum_off.cutpoints, include.lowest = T, dig.lab = 2),
	`Age at hire` =  cut(yin.gm - yob.gm, age_hire.cutpoints, include.lowest = T, dig.lab = 2),
	`Duration of employment` =  cut(employment.years, years_since_hire.cutpoints, include.lowest = T, dig.lab = 2)
)]
dat[, `Cumulative soluble shifted` := shift(`Cumulative soluble`, fill = levels(`Cumulative soluble`)[1])]

table(dat$`Soluble`, useNA = 'always')
table(dat$`Cumulative soluble`, useNA = 'always')
table(dat$`Cumulative straight`, useNA = 'always')
table(dat$`Cumulative synthetic`, useNA = 'always')
table(dat$`Cumulative time off`, useNA = 'always')
table(dat$`Age at hire`, useNA = 'always')
table(dat$`Duration of employment`, useNA = 'always')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Followed intervention indicators               ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setorder(dat, studyno, year)
dat[, years_since_start := 1:.N, studyno]

followed.names <- paste0("followed_", c("200", "100", "50", "25", "05", "025"))

for (i in with(dat, min(years_since_start):max(years_since_start))) {

	dat.tmp <- dat[years_since_start == i, c("studyno", "soluble"), with = F]

	# followed in that year?
	dat.tmp[, `:=`(
		followed_200 = as.numeric(soluble < 2),
		followed_100  = as.numeric(soluble < 1),
		followed_50  = as.numeric(soluble < 0.5),
		followed_25  = as.numeric(soluble < 0.25),
		followed_05  = as.numeric(soluble < 0.05),
		followed_025 = as.numeric(soluble < 0.025)),
		studyno
	]

	setorder(dat.tmp, studyno)
	dat[years_since_start == i, (followed.names) := as.list(
		dat.tmp[, followed.names, with = F])]
}

dat[,(followed.names) := lapply(
	followed.names, function(followed) {
		as.numeric(cumsum(get(followed)) == years_since_start)
	}),
	studyno]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save full dataset to Box                       ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box_write(dat,
					file_name = "dat.analytic.rds",
					dir_id = 125518026124)
# file id: 1179471733155

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reduced dataset                                ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat[status == 1 & year <= 1994 + exposure.lag, I] %>% table

times <- seq(0, max(dat[status == 1, I], na.rm = T) + 1, by = 4)

data.frame(
	time = times[-1],
	cases = sapply(2:length(times), function (i = 2) {
		dat[year <= 1994 + exposure.lag, sum(
			status[I <= times[i] & I > times[i - 1]],
			na.rm = T)]})
)

reduce.data <- function(
		df = data.table::copy(dat),
		times) {

	df <- rbindlist(lapply(length(times):2, function(i = length(times)) {
		# Any rows in this interval? If so, pick the last row
		d <- df[immortal == 0 & I > times[i - 1] & I <= times[i],]
		d[,`:=`(I.which = max(I)), by = .(studyno)]
		d[, (followed.names) := lapply(followed.names, function(followed){
			all(get(followed) == 1)
		}), studyno]
		# d[,.(studyno, year, start, end, py)]
		d[, `:=`(
			end = as.Date(apply(data.frame(
				as.Date(paste0(year, "-12-31")),
				yod,
				yoi,
				yoc
			), 1, min, na.rm = T)),
			start = as.Date(paste0(year, "-01-01"))
		)]
		d[,`:=`(py = as.numeric(end + 1 - start)/365), studyno]
		d[,py := sum(py), studyno]
		return(d[I == I.which])
	}))
	setorder(df, studyno, year)
	df[immortal == 0, I := 1:.N, by = .(studyno)]
	return(df)
}

dat.reduced <- reduce.data(dat, times)
dat.reduced[,`Employment status` := factor(`Employment status`, levels = c("At work", "Left work"))]
dat.reduced[, Observed := as.numeric(factor(Censored, 1:0)) - 1]
dat.reduced$Observed %>% table(useNA = "ifany")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save reduced dataset to Box                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box_write(dat.reduced,
					file_name = "dat.reduced.rds",
					dir_id = 125518026124)
