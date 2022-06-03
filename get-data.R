# Get data for cancer incidence ####
# Kevin Chen
# November 3, 2020

library(here)

library(survival)

if (!("exposure.lag" %in% ls())) {exposure.lag <- 10}

year.max <- 2015
employment_status.lag <- 0
yout.which <- "yout"

# Source 00-hello.R
if (!('cohort' %in% ls(envir = .GlobalEnv))) {
	source(here::here('../gm-wrangling/wrangling', '00-hello.R'))
}

# Source script for cancer incidence helper functions
source(here::here('../gm-cancer-inc', 'incidence.R'))

incidence.key <- data.table::fread(here::here("../gm-wrangling/cancer incidence", 'cancer-key.tsv'))
# Expand incidence key ####
incidence.key[, `:=`(
	var.name = paste0("canc_", code),
	date.name = paste0("ddiag_", code))]
incidence.key <- rbindlist(
	list(incidence.key,
			 data.table(
			 	code = c("first", "copd", "external"),
			 	description = c(
			 		"All cancers",
			 		"Chronic obstructive pulmonary disease",
			 		"All external causes"),
			 	var.name = c(
			 		"canc_first",
			 		"Chronic obstructive pulmonary disease",
			 		"All external causes"),
			 	date.name = c(
			 		"ddiag_first",
			 		"yod", "yod")
			 )))

# Get cohort analytic ####
if (!('cohort_analytic' %in% ls())) {
	outcome.type <- 'incidence'
	cohort_analytic <- get.cohort_analytic(
		outcome_type = outcome.type,
		exposure.lag = exposure.lag,
		deathage.max = NULL,
		year.max = year.max,
		hire.year.min = -Inf,
		use_seer = T,
	)
	setorder(cohort_analytic, studyno, year)
	cohort_analytic[, `:=`(yin.gm = date.to.gm(yin))]

	# Keep only people who appear in the exposure data
	cohort_analytic <- cohort_analytic[studyno %in% unique(exposure$studyno)]

	# PICK YOUT ####
	cohort_analytic[, jobloss.date := get(yout.which)]
	cohort_analytic[year > year(jobloss.date), employment := "Left work"]
	cohort_analytic[year <= year(jobloss.date), employment := "At work"]

	# Exposure after leaving work is 0
	cohort_analytic[year > (year(jobloss.date) + exposure.lag), `:=`(
		straight = 0,
		soluble = 0,
		synthetic = 0)]
	# NA fill
	cohort_analytic[year <= (year(jobloss.date) + exposure.lag), `:=`(
		straight = zoo::na.locf(straight),
		soluble = zoo::na.locf(soluble),
		synthetic = zoo::na.locf(synthetic)
		), by = .(studyno)]

	cohort_analytic[,`:=`(
		off = {
			off.gan[is.na(off.gan)] <- 0
			off.han[is.na(off.han)] <- 0
			off.san[is.na(off.san)] <- 0
			off <- off.gan + off.san + off.han
			off[off > 1] <- 1
			off}
	), by = .(studyno)]

	cohort_analytic[, `:=`(
		cum_straight = cumsum(straight),
		cum_soluble = cumsum(soluble),
		cum_synthetic = cumsum(synthetic),
		cum_off = cumsum(off)
	), by = .(studyno)]

	# Which columns ####
	col.names <- names(cohort_analytic[, c(
		"studyno",
		"age.year1",
		"age.year2",
		"year1",
		"year2",
		grep("canc\\_", names(cohort_analytic), value = T),
		paste0(c("cum_", ""), rep(c("straight", "soluble", "synthetic", "off"), each = 2)),
		"year",
		"yin.gm",
		"yin",
		"yrin",
		"yrin16",
		"race",
		"finrace",
		"plant",
		grep("ddiag", names(cohort_analytic), value = T),
		"yod",
		"yoc",
		"yob",
		"sex",
		"dateout.date",
		"employment_end.date",
		"employment_end.date.legacy",
		"yout", "yout16",
		"yout_recode",
		"jobloss.date",
		"All causes",
		"Chronic obstructive pulmonary disease",
		"All external causes",
		"nohist", "wh", "immortal", "right.censored",
		"possdiscr_new", "flag77", "oddend",
		"status15", "cancinccoh15_new"), with = F])

	# Drop unnecessary data ####
	cohort_analytic <- cohort_analytic[
		wh == 1 & nohist == 0 &
			# cancinccoh15_new == 1 &
			possdiscr_new == 0 &
			# immortal == 0 &
			right.censored == 0,
		col.names, with = F]
}

Sys.sleep(0)

# # Try out get.coxph() ####
# get.coxph(outcomes = 31)