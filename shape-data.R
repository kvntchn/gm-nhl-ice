# Shaping data for parametric g-formula and the sdr lmtp ####
# Kevin Chen
# November 3, 2020

library(here)
library(boxr)
library(survival)
box_auth()

# Get raw data and helper functions ####
exposure.lag <- 10
race.imputed.melt <- NULL
source("~/headRs/00-my-theme.R")
library(data.table)
# # paste(ls(), collapse = ",\n") %>% clipr::write_clip()
source(here::here("../gm-cancer-inc/breaks.R"))
source("blip.r")

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
# dat <- dat[yin <= as.Date("1983-01-01") & year >= 1985 & year <= 1994 + exposure.lag]
# Get data analytic ####
get.coxph(
	cohort_name = "dat",
	run_model = F,
	time_scale = "age",
	start.year = -Inf,
	employment_status.lag = exposure.lag)
dat <- copy(nhl.dat)

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
	`Age at hire` = yob.gm - yin.gm
)]
dat[I == N & status == 0 & year == floor(yod.gm),
		`Other event` := 1]
dat[is.na(`Other event`),
		`Other event` := 0]
# Lag employment status and duration of employment ####
dat[, `Employment status` := ifelse(year > floor(jobloss.date.gm) + exposure.lag, "Left work", "At work")]
dat[, employment.years := apply(data.frame(jobloss.date.gm, year + 1), 1, min) - yin.gm, studyno]
dat[, employment.years := shift(employment.years, exposure.lag, 0), studyno]
covariate.breaks <- get.covariate.breaks(dat)
dat[, `Duration of employment` := get.cut(employment.years, covariate.breaks)]

dat[I == N & status == 0 & (I < 1994 + exposure.lag - 1985 + 1), Censored := 1]
dat[is.na(Censored), Censored := 0]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save to Box                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box_write(dat,
					"dat.rds",
					dir_id = 125518026124)

dat <- box_read(737540256124)
mwf.breaks <- get.mwf.breaks(dat)
dat <- dat[yout >= as.Date("1941-01-01") & year <= exposure.lag + 1994 & yin.gm < 1982,]

# Exposure over time ####
exposure.ggplot <- dat[, .(
	Soluble   = median(  soluble[soluble > 0]),
	Straight  = median( straight[straight > 0]),
	Synthetic = median(synthetic[synthetic > 0])
), .(year = year - exposure.lag)] %>% melt(
	id.vars = "year",
	variable.name = "mwf") %>%
	ggplot(aes(x = year, y = value)) +
	geom_line()
tikzDevice::tikz(here::here("resources", "exposure.tex"), height = 2.5, width = 6, standAlone = T)
 exposure.ggplot +
	facet_wrap(. ~ mwf) +
	labs(x = "Calendar year",
			 y = "Median exposure (mg/m$^3$) among exposed") +
	coord_cartesian(xlim = c(1917, 2005)) +
	mytheme
dev.off()
lualatex("^exposure\\.tex", here::here("resources"))

tikzDevice::tikz(here::here("resources", "exposure_long.tex"), height = 5.5, width = 3, standAlone = T)
 exposure.ggplot +
	facet_wrap(. ~ mwf, nrow = 3) +
	labs(x = "Calendar year",
			 y = "Median exposure (mg/m$^3$) among exposed") +
	coord_cartesian(xlim = c(1917, 2005)) +
	mytheme
dev.off()
lualatex("^exposure_long\\.tex", here::here("resources"))

# Population characteristics ####
get.tab1 <- function(
	df = dat[year >= 1985 & year <= 1994 + exposure.lag],
	table.engine = "pander",
	use_finrace = T,
	incidence = F,
	py = "py",
	mathmode = F,
	nrow_as_fu = T) {

	setorder(df, studyno, year)
	# Individual-level summary
	tab1.sum <- df[,.(
		'Race' = NA,
		'\\hspace{10pt}White' = if (use_finrace) {
			ifelse(finrace[1] == 1, 1, 0)} else {
				as.numeric(race[1] == "White")
			},
		'\\hspace{10pt}Black' = if (use_finrace) {
			ifelse(finrace[1] == 2, 1, 0)} else {
				as.numeric(race[1] == "Black")
			},
		'\\hspace{10pt}Unknown' = if (use_finrace) {
			ifelse(finrace[1] %in% c(0, 9), 1, 0)} else {
				as.numeric(race[1] == "Unknown")
			},
		'Sex' = NA,
		'\\hspace{10pt}Male' = ifelse(sex[1] == 'M', 1, 0),
		'\\hspace{10pt}Female' = ifelse(sex[1] == 'F', 1, 0),
		"Plant$^1$" = NA,
		"\\hspace{10pt}Plant 1" = as.numeric(Plant[1]) == 1,
		"\\hspace{10pt}Plant 2" = as.numeric(Plant[1]) == 2,
		"\\hspace{10pt}Plant 3" = as.numeric(Plant[1]) == 3,
		"Ever exposed to MWFs" = NA,
		# as.numeric(cum_soluble[.N] + cum_straight[.N] + cum_synthetic[.N] > 0)
		"\\hspace{10pt}Straight"  = as.numeric(cum_straight[.N] > 0),
		"\\hspace{10pt}Soluble"   = as.numeric(cum_soluble[.N] > 0),
		"\\hspace{10pt}Synthetic" = as.numeric(cum_synthetic[.N] > 0),
		# "Employment status in 1995" = NA,
		# "\\hspace{10pt}Left work" = as.numeric(jobloss.date[1] < as.Date("1995-01-01")),
		# "\\hspace{10pt}Still at work" = as.numeric(jobloss.date[1] >= as.Date("1995-01-01")),
		"Deceased by end of follow-up" = {
			if (!incidence) {
				as.numeric(!is.na(yod[1]) & yod[1] <= as.Date(paste0(year[.N], "-12-31")))
			} else {
				as.numeric(!is.na(ddiag_first[1]) & ddiag_first[1] <= as.Date(paste0(year[.N], "-12-31")))
			}
		},
		# "\\hspace{10pt}ALD" = as.numeric(max(`Alcohol-related Liver Disease`) > 0),
		# "\\hspace{10pt}Suicide" = as.numeric(max(Suicide) > 0),
		# "\\hspace{10pt}Overdose" = as.numeric(max(Overdose) > 0),


		# Years since follow-up
		"\\hline Years of follow-up" = if (nrow_as_fu) {.N} else {
			as.numeric(difftime(min(yod[1], yoc[1], as.Date(paste0(max(year), "-12-31")), na.rm = T),
				as.Date(paste0(min(year[immortal == 0]), "-01-01")), units = "days"))/365
		},
		'Year of birth' = yob.gm[1],
		'Year of hire' = yin.gm[1],
		'Age at hire (years)' = yin.gm[1] - yob.gm[1],
		'Year of leaving work$^2$' = {
			if (year(jobloss.date[1]) == 1995) {NaN} else {jobloss.date.gm[1]}},
		'Age at leaving work (years)$^2$' = {
			if (year(jobloss.date[1]) == 1995) {NaN} else {
				jobloss.date.gm[1] - yob.gm[1]
			}},
		"Years at work$^2$" = {
			if (year(jobloss.date[1]) == 1995) {
				NaN} else {as.numeric(difftime(jobloss.date[1], yin[1], units = "days"))/365}},
		'Year of death among deceased' = {
			if (!incidence) {
				as.numeric(year(yod[1]))} else {
					as.numeric(year(ddiag_first[1]))
				}},
		'Age at death (years) among deceased' = {
			if (!incidence) {
				as.numeric(difftime(yod[1], yob[1], units = "days"))/365} else {
					as.numeric(difftime(ddiag_first[1], yob[1], units = "days"))/365
				}},
		# "Maximum cumulative exposure$^3$ (mg/m$^3\\cdot$y)}" =
		# 	NA,
		# 	max(cum_soluble[cum_soluble != 0]) +
		# 	  max(cum_straight[cum_straight != 0]) +
		#   	max(cum_synthetic[cum_synthetic != 0])

		'Cumulative time off (years)' = {
			max(cum_off)
			},

		"Cumulative exposure$^3$ to MWFs (mg/m$^3\\cdot$y)" = NA,

		"\\hspace{10pt}Straight "  = as.numeric(ifelse(
			sum(cum_straight != 0) > 0,
			max(cum_straight[cum_straight != 0], na.rm = T),
			NA)),

		"\\hspace{10pt}Soluble "   = as.numeric(ifelse(
			sum(cum_soluble != 0) > 0,
			max(cum_soluble[cum_soluble != 0], na.rm = T),
			NA)),

		"\\hspace{10pt}Synthetic " = as.numeric(ifelse(
			sum(cum_synthetic != 0) > 0,
			max(cum_synthetic[cum_synthetic != 0], na.rm = T),
			NA))

	),
	by = .(studyno)][,-'studyno']

	# Correct names if looking at incidence outcome
	if (incidence) {
		names(tab1.sum) <- gsub(
			"Deceased by end of follow-up",
			"Diagnosed with cancer by end of follow-up",
			names(tab1.sum))
		names(tab1.sum) <- gsub(
			"Year of death among deceased",
			"Year of first cancer diagnosis",
			names(tab1.sum))
		names(tab1.sum) <- gsub(
			"Age at death \\(years\\) among deceased",
			"Age at first cancer diagnosis (years)",
			names(tab1.sum))
	}

	# Population-level summary
	table.break <- "\\hline Years of follow-"
	tab1 <- rbind(
		t(apply(tab1.sum[,1:(grep(table.break, names(tab1.sum)) - 1)], 2, function(x) {
			return(
				c(mean(x, na.rm = T) * sum(!is.na(x)),
					mean(x, na.rm = T) * 100,
					NA
					# NA,
					# NA
				)
			)
		})),
		t(apply(tab1.sum[,grep(table.break, names(tab1.sum)):ncol(tab1.sum)], 2,
						function(x) {
							return(
								c(
									# mean(x, na.rm = T),
									# sd(x, na.rm = T)
									# min(x, na.rm = T),
									median(as.numeric(x), na.rm = T),
									quantile(as.numeric(x), 0.25, na.rm = T),
									quantile(as.numeric(x), 0.75, na.rm = T)
									# max(x, na.rm = T)
								)
							)
						}))
	)

	tab1[!is.finite(as.matrix(tab1))] <- NA

	# Table names
	colnames(tab1) <- c('center', 'spread', " ")#, 'Minimum', 'Median', 'Maximum')

	# Digits
	tab1.digits <- matrix(2, ncol = 3,#5
												nrow = nrow(tab1))
	tab1.digits[grep("Age", rownames(tab1), ignore.case = T), ] <- 1
	tab1.digits[grep("years ", rownames(tab1), ignore.case = T), ] <- 1
	tab1.digits[grep("year ", rownames(tab1), ignore.case = T), ] <- 0
	tab1.digits[1:(grep(table.break, rownames(tab1), ignore.case = T) - 1), 1] <- 0
	tab1.digits[1:(grep(table.break, rownames(tab1), ignore.case = T) - 1), -1] <- 0

	# more digits if small %s
	tab1.digits[1:nrow(tab1.digits) < grep(table.break, rownames(tab1), ignore.case = T) &
								(tab1[,2] < 2 & tab1[,2] != 0), -1] <- 2

	tab1.digits[1:nrow(tab1.digits) < grep(table.break, rownames(tab1), ignore.case = T) &
								(tab1[,2] < 1 & tab1[,2] != 0), -1] <- 2

	tab1.digits[1:nrow(tab1.digits) < grep(table.break, rownames(tab1), ignore.case = T) &
								(tab1[,2] > 99 & tab1[,2] != 100), -1] <- 1

	# more digits if numbers
	tab1.digits[1:nrow(tab1.digits) >= grep(table.break, rownames(tab1), ignore.case = T) &
								(tab1[,2] < .1 & tab1[,2] != 0), -1] <- 2

	tab1.digits[1:nrow(tab1.digits) >= grep(table.break, rownames(tab1), ignore.case = T) &
								(tab1[,2] == 0), -1] <- 0

	tab1.which_year <- matrix(F, ncol = 3, nrow = nrow(tab1))
	tab1.which_year[grep("year ", rownames(tab1), ignore.case = T), ] <- T

	tab1 <- matrix(
		sapply(1:length(tab1), function(i) {
			if (is.na(as.vector(tab1)[i])) {NA} else {
			formatC(as.vector(tab1)[i], as.vector(tab1.digits)[i], format = "f",
							big.mark = if (!i %in% which(as.vector(tab1.which_year))) {
								if (mathmode) {"\\\\,"} else {","}} else {""})
		}}),
		ncol = ncol(tab1),
		nrow = nrow(tab1),
		dimnames = dimnames(tab1)
	)

	# Pretty percents
	tab1[1:(grep("Years of fo", rownames(tab1)) - 1), 2] <- sapply(
		tab1[1:(grep("Years of fo", rownames(tab1)) - 1), 2],
		function (x) {
			if (!is.na(x)) {
				paste0(x, ifelse(table.engine == "xtable", '\\%', "%"))
			} else {NA}
		})

	# Math mode
	if (mathmode) {
	tab1 <- matrix(
		sapply(1:length(tab1), function(i) {
			if (!is.na(as.vector(tab1)[i])) {
				paste0('$', as.vector(tab1)[i], '$')
			} else {NA}
		}),
		ncol = ncol(tab1),
		nrow = nrow(tab1),
		dimnames = dimnames(tab1))
	}

	tab1 <- rbind(
		matrix(c(
			paste0(
				if (mathmode) {"$"} else {""},
				prettyNum(n_distinct(df$studyno),
									ifelse(mathmode, '\\\\,', ",")),
				if (mathmode) {"$"} else {""}),
			paste0(
				if (mathmode) {"$"} else {""},
				formatC(ifelse(nrow_as_fu, nrow(df), sum(df[,py])),
								format = "f", digits = 0,
								big.mark = if (mathmode) {"\\\\,"} else {","}),
				if (mathmode) {"$"} else {""}),
			NA
			# , NA, NA
		),
		nrow = 1,
		dimnames = list(c("Study population size (person-years)"), NULL)),
		tab1)

	# Make column indicating stat type
	tab1 <- cbind(tab1, spread.which = c(
		rep(paste0(if (mathmode) {"$"} else {""}, "n", if (mathmode) {"$"} else {""}, ifelse(table.engine == "xtable", "\\%", "%")),
				grep("\\hline", rownames(tab1)) - 1),
		rep("median, Q1, Q3", nrow(tab1) - grep("\\hline", rownames(tab1)) + 1))
	)

	# Clean up rownames
	rownames(tab1)[duplicated(rownames(tab1))] <- paste0(rownames(tab1)[duplicated(rownames(tab1))], " ")
	tab1 <- as.data.frame(tab1, make.names = F)

	print(matrix(as.vector(sapply(tab1[,1:3], function(x) gsub("\\$|\\\\", "", x))),
				 ncol = 3,
				 dimnames = list(rownames(tab1), colnames(tab1)[1:3])))

	return(tab1)
}
tab1 <- get.tab1(dat[year >= 1985 & year <= 1994 + exposure.lag], mathmode = F, nrow_as_fu = T)
saveRDS(tab1, here::here("resources", "tab1.rds"))
nhl.tab1 <- get.tab1(dat[
	studyno %in% dat[status == 1 & year <= 1994 + exposure.lag, studyno] &
	year >= 1985 & year <= 1994 + exposure.lag],
	mathmode = F, nrow_as_fu = T)
saveRDS(nhl.tab1, here::here("resources", "nhl.tab1.rds"))

clean.tab1 <- function(tab1) {
	tab1[grepl("^n", tab1[,4]) & !is.na(tab1[,1]), 2] <- paste0(
		"(", tab1[grepl("^n", tab1[,4]) & !is.na(tab1[,1]), 2], ")"
	)

	tab1[grepl("median", tab1[,4]) & !is.na(tab1[,1]), 2] <- paste0(
		"(", tab1[grepl("median", tab1[,4]) & !is.na(tab1[,1]), 2], ", ",
		tab1[grepl("median", tab1[,4]) & !is.na(tab1[,1]), 3],
		")"
	)
	return(tab1[,1:2])
}
tab1 <- readRDS(here::here("resources", "tab1.rds"))
nhl.tab1 <- readRDS(here::here("resources", "nhl.tab1.rds"))

clipr::write_clip(
	kable(cbind(clean.tab1(tab1), NA, clean.tab1(nhl.tab1)))
	)
clipr::write_clip(
	print(xtable(cbind(
		rownames(clean.tab1(tab1)),
		mutate(clean.tab1(tab1), spread = sanitize(spread)), NA,
		mutate(clean.tab1(nhl.tab1), spread = sanitize(spread))
	)))
	)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unknown race as separate category              ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat[finrace == 9, Race := "Unknown"]

setorder(dat, studyno, year)
dat[,I := NA]
dat[immortal == 0 & year >= 1985, I := 1:.N, studyno]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Post-shift-intervention data                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat[,`Cumulative soluble` := `Cumulative soluble 5`]
mwf.breaks[,cum_soluble := cum_soluble5]

# dat.shift <- get.shift_intervention(dta = dat[,.(studyno, immortal, I, year, straight, status, yod, yoi, yoc)], "soluble", delta = 0.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Scale down all three MWF so total under REL    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.total <- dat[,.(studyno, immortal, I, year, straight, soluble, synthetic, status, yoi, yoc,  yod)]
# Shift down if total above REL = 0.5 mg/m3
dat.total[soluble + straight + synthetic > 0.5, `:=`(
	straight = straight * 0.5 / (soluble + straight + synthetic),
	soluble = soluble * 0.5 / (soluble + straight + synthetic),
	synthetic = synthetic * 0.5 / (soluble + straight + synthetic)
	)]

dat.total[, `:=`(
	cum_straight = cumsum(straight),
	cum_soluble = cumsum(soluble),
	cum_synthetic = cumsum(synthetic)
	), by = .(studyno)]

dat.total[,`:=`(
			Straight = get.cut(straight, mwf.breaks, "straight", dig.lab = 3),
			`Cumulative straight` = get.cut(cum_straight, mwf.breaks, "cum_straight", dig.lab = 3),
			Soluble = get.cut(soluble, mwf.breaks, "soluble", dig.lab = 3),
			`Cumulative soluble` = get.cut(cum_soluble, mwf.breaks, "cum_soluble", dig.lab = 3),
			Synthetic = get.cut(soluble, mwf.breaks, "soluble", dig.lab = 3),
			`Cumulative synthetic` = get.cut(cum_synthetic, mwf.breaks, "cum_synthetic", dig.lab = 3)
				)]


dat.total2 <- dat[,.(studyno, immortal, I, year, straight, soluble, synthetic, status, yoi, yoc,  yod)]
# Shift down if total above REL = 0.5 mg/m3
dat.total2[soluble + straight + synthetic > 0.5/2, `:=`(
	straight = straight * 0.5/2 / (soluble + straight + synthetic),
	soluble = soluble * 0.5/2 / (soluble + straight + synthetic),
	synthetic = synthetic * 0.5/2 / (soluble + straight + synthetic)
	)]

dat.total2[, `:=`(
	cum_straight = cumsum(straight),
	cum_soluble = cumsum(soluble),
	cum_synthetic = cumsum(synthetic)
	), by = .(studyno)]

dat.total2[,`:=`(
			Straight = get.cut(straight, mwf.breaks, "straight", dig.lab = 3),
			`Cumulative straight` = get.cut(cum_straight, mwf.breaks, "cum_straight", dig.lab = 3),
			Soluble = get.cut(soluble, mwf.breaks, "soluble", dig.lab = 3),
			`Cumulative soluble` = get.cut(cum_soluble, mwf.breaks, "cum_soluble", dig.lab = 3),
			Synthetic = get.cut(soluble, mwf.breaks, "soluble", dig.lab = 3),
			`Cumulative synthetic` = get.cut(cum_synthetic, mwf.breaks, "cum_synthetic", dig.lab = 3)
				)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get post-intervention data (exposure limit)   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.limit_intervention <- function(
	df = data.table::copy(dat),
	mwf = 'straight',
	limit = 0.5
) {
	mwf <- grep(mwf, c("straight", "soluble", "synthetic"), value = T)

	df <- df[,c("studyno", "immortal", "I", "year",
							"straight", "soluble", "synthetic",
							"status", "yoi", "yoc", "yod"), with = F]
	df[get(mwf) > limit, (mwf) := list(limit)]
	df[, (paste0("cum_", mwf)) := cumsum(get(mwf)), by = .(studyno)]
	df[,(c(str_to_title(mwf), paste0("Cumulative ", mwf))) := list(
		get.cut(get(mwf), mwf.breaks, mwf, dig.lab = 3),
		get.cut(get(paste0('cum_', mwf)), mwf.breaks, paste0('cum_', mwf), dig.lab = 3)
				)]

	return(df)
}

get.limit_total_intervention <- function(
	df = data.table::copy(dat),
	mwf = 'straight',
	limit = 0.5
) {
	mwf <- grep(mwf, c("straight", "soluble", "synthetic"), value = T)

	df <- df[,.(studyno, immortal, I, year,
							straight, soluble, synthetic,
							status, yoi, yoc,  yod)]

	other_mwf <- c('straight', 'soluble', 'synthetic')
	other_mwf <- other_mwf[other_mwf != mwf]

	df[straight + soluble + synthetic > limit, (mwf) := lapply(
		list(get(mwf)), function(x) {
			x <- get(mwf)
			x <- limit - get(other_mwf[1]) - get(other_mwf[2])
			x[x < 0] <- 0
			return(x)
		})]

	df[straight + soluble + synthetic < limit, .(straight, studyno, year)]

	df[, (paste0("cum_", mwf)) := cumsum(get(mwf)), by = .(studyno)]
	df[,(c(str_to_title(mwf), paste0("Cumulative ", mwf))) := list(
		get.cut(get(mwf), mwf.breaks, mwf, dig.lab = 3),
		get.cut(get(paste0('cum_', mwf)), mwf.breaks, paste0('cum_', mwf), dig.lab = 3)
	)]

	return(df)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Straight to REL							                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat.str   <- get.limit_intervention(mwf = 'str', limit = 0.5)
dat.str2  <- get.limit_intervention(mwf = 'str', limit = 0.5/2)
dat.str10 <- get.limit_intervention(mwf = 'str', limit = 0.5/10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Straight so that total under REL	             ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.str_total  <- get.limit_total_intervention(mwf = "str", limit = 0.5)
dat.str_total2 <- get.limit_total_intervention(mwf = "str", limit = 0.5/2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Soluble to REL							                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.sol   <- get.limit_intervention(mwf = 'sol', limit = 0.5)
dat.sol2  <- get.limit_intervention(mwf = 'sol', limit = 0.5/2)
dat.sol10 <- get.limit_intervention(mwf = 'sol', limit = 0.5/10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Soluble so that total under REL	               ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.sol_total  <- get.limit_total_intervention(mwf = "sol", limit = 0.5)
dat.sol_total2 <- get.limit_total_intervention(mwf = "sol", limit = 0.5/2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Synthetic to REL							                 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.syn   <- get.limit_intervention(mwf = 'syn', limit = 0.5)
dat.syn2  <- get.limit_intervention(mwf = 'syn', limit = 0.5/2)
dat.syn10 <- get.limit_intervention(mwf = 'syn', limit = 0.5/10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Synthetic so that total under REL	             ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat.syn_total  <- get.limit_total_intervention(mwf = "syn", limit = 0.5)
dat.syn_total2 <- get.limit_total_intervention(mwf = "syn", limit = 0.5/2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compare exposure before/after intervention     ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (mwf in c("Straight", "Soluble", "Synthetic" )) {
	var.names <- c("studyno", "year", "status",
								 mwf, tolower(mwf), paste0("Cumulative ", tolower(mwf)),
								 paste0("cum_", tolower(mwf))
	)
	exposure.comparison <- rbindlist(list(
		dat[, var.names, with = F],
		get(paste0("dat.", tolower(substr(mwf, 1, 3))))[, var.names, with = F],
		get(paste0("dat.", tolower(substr(mwf, 1, 3)), "_total"))[, var.names, with = F],
		dat.total[, var.names, with = F]
	),
	idcol = "Intervention")
	exposure.comparison[,`:=`(Intervention = c(
		"No intervention",
		paste0(mwf, " to REL"),
		paste0(mwf, " so total under REL"),
		"Scale all")[Intervention]
		)]

	exposure.comparison[, year.max := max(year), studyno]

	print(table(
		exposure.comparison[
			status == 1,
			c(paste0("Cumulative ", tolower(mwf)), "Intervention"),
			with = F]
	))

	# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# # Plot post-shift-intervention data              ####
	# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# exposure.comparison[year == year.max] %>% ggplot(
	# 	aes(x = get(paste0("cum_", tolower(mwf))))
	# ) +
	# 	geom_density(aes(color = Intervention)) +
	# 	# geom_histogram(bins = 40) +
	# 	# facet_wrap(. ~ Intervention, ncol = 1) +
	# 	mytheme +
	# 	scale_x_continuous(
	# 		n.breaks = 6,
	# 		trans = "log",
	# 		labels = function(x) {sprintf("%.3f", x)}) +
	# 	labs(
	# 		y = "Density",
	# 		# y = "Count (person-time)",
	# 		x = "Cumulative exposure (mg/m$^3$)") + theme(
	# 			axis.title = element_text(margin = margin(10, 5, 10, 5)),
	# 			legend.title = element_blank(),
	# 			legend.position = "bottom",
	# 			legend.box.spacing = unit(2, "pt"),
	# 			legend.key.size = unit(10, "pt")
	# 		)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reduced dataset                                ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat[status == 1 & year <= 1994 + exposure.lag, I] %>% table

times <- c(0, 4, seq(
	8, max(dat$I, na.rm = T),
	2))
times[length(times)] <- max(times) + 1
times

data.frame(
	time = times[-1],
	cases = sapply(2:length(times), function (i = 2) {
		dat[year <= 1994 + exposure.lag, sum(
			status[I <= times[i] & I > times[i - 1]],
			na.rm = T)]}))

reduce.data <- function(
	df = data.table::copy(dat), times, observed.dat = data.table::copy(dat.reduced), merge_obs = T) {
	df <- rbindlist(lapply(length(times):2, function(i = length(times)) {
		# Any rows in this interval? If so, pick the last row
		d <- df[immortal == 0 & I > times[i - 1] & I <= times[i],]
		d[,I.which := max(I), by = .(studyno)]
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
		d[,`:=`(
			py = as.numeric(end + 1- start)/365
		),
		studyno]
		d[,py := sum(py), studyno]
		return(d[I == I.which])
	}))
	setorder(df, studyno, year)
	df[,.(studyno, I, year)]
	df[immortal == 0, I := 1:.N, by = .(studyno)]

	if (merge_obs) {
		df <- merge(
			df,
			observed.dat[,.(studyno, year, Observed)],
			by = c("studyno", "year")
		)
	}

	return(df)
}

dat.reduced <- reduce.data(dat, times, merge = F)
table(dat.reduced[, c("I", "year")])
table(dat.reduced[status == 1]$I)
sum(dat.reduced$status)

dat.reduced[,`Employment status` := factor(`Employment status`, levels = c("At work", "Left work"))]
dat.reduced[, Observed := as.numeric(factor(Censored, 1:0)) - 1]
dat.reduced$Observed %>% table(useNA = "ifany")

# Scale down all by constant factor
dat.total.reduced  <- reduce.data(dat.total, times)
dat.total.reduced2 <- reduce.data(dat.total2, times)

# Intervene on straight
dat.str.reduced   <- reduce.data(dat.str, times)
dat.str.reduced2  <- reduce.data(dat.str2, times)
dat.str.reduced10 <- reduce.data(dat.str10, times)

dat.str_total.reduced  <- reduce.data(dat.str_total, times)
dat.str_total.reduced2 <- reduce.data(dat.str_total2, times)

# Intervene on soluble
dat.sol.reduced   <- reduce.data(dat.sol, times)
dat.sol.reduced2  <- reduce.data(dat.sol2, times)
dat.sol.reduced10 <- reduce.data(dat.sol10, times)

dat.sol_total.reduced  <- reduce.data(dat.sol_total, times)
dat.sol_total.reduced2 <- reduce.data(dat.sol_total2, times)

# Intervene on synthetic
dat.syn.reduced   <- reduce.data(dat.syn, times)
dat.syn.reduced2  <- reduce.data(dat.syn2, times)
dat.syn.reduced10 <- reduce.data(dat.syn10, times)

dat.syn_total.reduced  <- reduce.data(dat.syn_total, times)
dat.syn_total.reduced2 <- reduce.data(dat.syn_total2, times)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save to Box                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box_write(dat.reduced,
					file_name = "dat.reduced.rds",
					dir_id = 125518026124)

# name        : dat.reduced.rds
# file id     : 928445783979

dat.intervention.names <- paste0("dat.",
			 paste0(rep(c("str", "sol", "syn"), each = 2),
			 			 rep(c("", "_total"), 3)),
			 ".reduced")

for (j in dat.intervention.names) {
	box_write(get(j),
						file_name = paste0(j, ".rds"),
						dir_id = 125518026124)
}

for (j in paste0(dat.intervention.names, "2")) {
	box_write(get(j),
						file_name = paste0(j, ".rds"),
						dir_id = 125518026124)
}

for (j in paste0(dat.intervention.names[
	- grep("total",
				 dat.intervention.names)], "10")) {
	box_write(get(j),
						file_name = paste0(j, ".rds"),
						dir_id = 125518026124)
}

# name          : dat.str.reduced.rds
# file id       : 928446426341
# name          : dat.str_total.reduced.rds
# file id       : 928447106294
# name          : dat.sol.reduced.rds
# file id       : 928445303042
# name          : dat.sol_total.reduced.rds
# file id       : 928447942153
# name          : dat.syn.reduced.rds
# file id       : 928451552327
# name          : dat.syn_total.reduced.rds
# file id       : 928447597463

# name          : dat.str.reduced2.rds
# file id       : 929952765078
# name          : dat.str_total.reduced2.rds
# file id       : 929957702377
# name          : dat.sol.reduced2.rds
# file id       : 929953044799
# name          : dat.sol_total.reduced2.rds
# file id       : 929954713019
# name          : dat.syn.reduced2.rds
# file id       : 929951401906
# name          : dat.syn_total.reduced2.rds
# file id       : 929954585841

# name          : dat.str.reduced10.rds
# file id       : 942770799114
# name          : dat.sol.reduced10.rds
# file id       : 930687554963
# name          : dat.syn.reduced10.rds
# file id       : 942772072287


box_write(dat.total.reduced,
						file_name = paste0("dat.total.reduced", ".rds"),
						dir_id = 125518026124)
# name        : dat.total.reduced.rds
# file id     : 928452001570


box_write(dat.total.reduced2,
						file_name = paste0("dat.total.reduced2", ".rds"),
						dir_id = 125518026124)
# name        : dat.total.reduced2.rds
# file id     : 929952423129


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Wide dat a for sdr lmtp                        ####
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# dat.reduced <- box_read(802974178895)
# dat.shift.reduced <- box_read(818286283716)
#
# # Relevant covariates
# covariates.names <- c(
# 	"studyno", "I", "year", "immortal",
# 	"status", "Observed",
# 	"age.year2", "year2",
# 	"Year", "Age",
# 	 "cum_soluble", "soluble",
# 	"Cumulative soluble", "Soluble",
# 	"straight", "synthetic",
# 	"Straight", "Synthetic",
# 	"cum_straight", "cum_synthetic",
# 	"Cumulative straight", "Cumulative synthetic",
# 	"employment.years", "cum_off",
# 	"Duration of employment", "Cumulative time off",
# 	"yin.gm",
# 	"Year of hire", "yout16",
# 	"Sex", "Race", "Plant",
# 	"Employment status")
#
# timevar.names <- c(
# 	"studyno", "I",
# 	"status", "Observed",
# 	"soluble", "cum_soluble",
# 	"Cumulative soluble", "Soluble",
# 	"year", "year2", "age.year2",
# 	"Year", "Age",
# 	"straight", "synthetic",
# 	"Straight", "Synthetic",
# 	"cum_straight", "cum_synthetic",
# 	"Cumulative straight", "Cumulative synthetic",
# 	"employment.years", "cum_off",
# 	"Duration of employment", "Cumulative time off",
# 	"Employment status")
#
# # Make data wide
# get.wide_data <- function(
# 	dat.long = copy(dat.reduced),
# 	covariates.names,
# 	timevar.names) {
#
# 	covariates.names <- covariates.names[covariates.names %in% names(dat.long)]
# 	timevar.names <- timevar.names[timevar.names %in% names(dat.long)]
#
# 	add.who <- dat.long[,.(I.max = max(I)), by = .(studyno)][I.max < max(dat.long$I), .(
# 		I = (I.max + 1):max(dat.long$I)
# 	), by = .(studyno)]
#
# 	dat.wide <- rbindlist(list(
# 		dat.long[,covariates.names, with = F],
# 		add.who
# 	), use.names = T, fill = T)
#
# 	# Impute in new rows
# 	if ("Observed" %in% covariates.names) {
# 	dat.wide[Observed == 0, status := NA]}
# 	if ("status" %in% covariates.names) {
# 	dat.wide[studyno %in% dat.wide[status == 1, studyno],
# 					 status := zoo::na.locf(status, na.rm = F),
# 					 .(studyno)]}
#
# 	if ("yout16" %in% covariates.names) {
# 		dat.wide[,`:=`(yout.gm = yout16 + 1900)]
# 		# # Make year of leaving work categorical
# 		# dat.wide[,`:=`(
# 		# 	Year_left_work = cut(yout.gm, unique(quantile(yout.gm, seq(0, 1, 1/6))), include.lowest = T)
# 		# )]
# 	}
#
# 	if ("Sex" %in% covariates.names) {
# 		# One-hot coding
# 		baseline_dat <- dat.wide[,.(
# 				Sex = as.numeric(Sex[1] == "Female"),
# 				Plant2 = as.numeric(which.max(table(Plant)) == 2),
# 				Plant3 = as.numeric(which.max(table(Plant)) == 3),
# 				yout = yout16[1],
# 				Race = as.numeric(Race[1]  == levels(Race)[1]),
# 				yin.gm = yin.gm[1],
# 				Year_of_hire = `Year of hire`[1]
# 			), by = .(studyno)]
#
# 		setDT(baseline_dat)
# 	} else {baseline_dat <- unique(dat.wide[, .(studyno)])}
#
# 	# cast wide
# 	timevar_dat <- dcast(dat.wide[, timevar.names, with = F],
# 											 studyno ~ I,
# 											 value.var = timevar.names[!grepl("studyno|I$", timevar.names)])
#
# 	if (length(timevar.names[!grepl("studyno|I$", timevar.names)]) == 1) {
# 		names(timevar_dat)[-1] <- paste0(timevar.names[!grepl("studyno|I$", timevar.names)], "_", names(timevar_dat)[-1])
# 	}
#
# 	# Merge
# 	dat.wide <- merge(baseline_dat, timevar_dat, by = "studyno", all = T)
#
# 	names(dat.wide) <- gsub(" ", "_", names(dat.wide))
#
# 	return(dat.wide)
# }
#
# dat.wide <- get.wide_data(copy(dat.reduced), covariates.names, timevar.names)
# # dat.shift.wide <- get.wide_data(
# # 	copy(dat.shift.reduced[,.(studyno, I, soluble, cum_soluble, `Cumulative soluble`)]),
# # 	covariates.names, timevar.names)
#
# # save to box
# box_write(dat.wide,
# 					"dat.reduced.wide.rds",
# 					dir_id = 125518026124)
#
# box_write(dat.shift.wide,
# 					"dat.shift.reduced.wide.rds",
# 					dir_id = 125518026124)