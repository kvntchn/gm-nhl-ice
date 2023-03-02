# Parametric g-formula ####
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

# Relevant covariates ####
covariates.which <- c(
	"studyno", "I", "year", "immortal",
	"status", "Censored", "Other event", "All causes",
	"age.year2", "Age", "year2", "Year",
	"straight", "soluble", "synthetic",
	"Straight", "Soluble", "Synthetic",
	"Cumulative straight", "Cumulative soluble", "Cumulative synthetic",
	"cum_straight", "cum_soluble", "cum_synthetic",
	"Employment status",
	"Duration of employment", "employment.years",
	"Cumulative time off", "cum_off",
	"Year of hire", "yin.gm", "Sex", "Race", "Plant")

#####################################################
# Iterated conditional expectation g-computation ####
#####################################################

# Hazard-extended ICE implementation of g-formula ####
# Wen et al. (2020)
ice_gcomp <- function(
		a = 0,
		dta = copy(dat.reduced[,covariates.which, with = F]),
		dta.a1 = NULL,
		# dta.a1 = copy(dat.a),
		times = sort(unique(dat.reduced[status == 1]$I), decreasing = T),
		mwf.name = "soluble",
		quiet = F,
		total_effect = F
) {

	if (total_effect) {
		dta[,`Other event` := max(`Other event`), studyno]
	}

	if ("immortal" %in% names(dta.a1)) {
		dta.a1 <- dta.a1[,-c("immortal", "I", "I.which"), with = F]
	}

	if (a != 0) {
		# nrow(dta.a1) == nrow(dta)
		dta.a1 <- merge(
			dta.a1[,c("studyno", "year", grep(
				mwf.name, names(dta.a1), ignore.case = T, value = T)),
				with = F],
			dta[,covariates.which[!as.logical(
				apply(sapply(mwf.name, grepl, x = covariates.which, ignore.case = T), 1,
							sum)
			)], with = F],
			by = c("studyno", "year"))

	} else {dta.a1 <- copy(dta)}

	setorder(dta, studyno, year)
	dta[, I := 1:.N, by = studyno]

	setorder(dta.a1, studyno, year)
	dta.a1[, I := 1:.N, by = studyno]

	# Pick covariates
	timevar_covariates.which <- c(
		"Age",
		# "age.year2",
		"Cumulative straight", "Cumulative soluble", "Cumulative synthetic",
		# "Duration of employment",
		"Employment status",
		# "Year",
		"Cumulative time off"
	)
	baseline_covariates.which <- c(
		"Year of hire", "Race", "Plant", "Sex")

	# Check support
	support <- sapply(c(timevar_covariates.which, baseline_covariates.which), function(x) {
		length(table(dta[status == 1, x, with = F]))
	})
	no_support <- F

	if (sum(support <= 1)) {
		no_support <- T
		message("Check support")
	}

	# 1. Pooled logistic regression for discrete hazard ####
	# Total number of outcomes
	J <- max(times)
	if (!quiet) {
		cat(paste0("Iteration ", 0, " of ", J - 1, "...\n"))
		cat(paste0("\t Getting discrete hazards\n"))
	}
	require(fastglm, quietly = T)
	discrete.formula <- as.formula(paste(
		"status",
		"~",
		"Age +",
		"`Cumulative straight` +",
		"`Cumulative soluble` +",
		"`Cumulative synthetic` +",
		"`Employment status` +",
		"`Year of hire` +",
		"`Cumulative time off` +",
		"Sex +",
		"Race + Plant",
		NULL
	))
	h <- fastglm(
		model.matrix(discrete.formula, data = dta[immortal == 0 & Censored == 0]),
		dta[immortal == 0 & Censored == 0]$status,
		family = "binomial"
	)
	# summary(h)

	# Check positivity
	if (a != 0) {
		dta.g <- merge(
			dta.a1[immortal == 0 & I %in% (J - 1):J, c("studyno",
																								 paste0("Cumulative ", mwf.name)
			),
			with = F],
			dta[immortal == 0 & I %in% (J - 1):J & Censored == 0]
		)

		dta.g.tab <- dta.g[get(
			paste0("Cumulative ", mwf.name, ".x")
		) == get(
			paste0("Cumulative ", mwf.name, ".y")
		), .(
			n = n_distinct(studyno),
			overlap = sum(get(
				paste0("Cumulative ", mwf.name, ".x")
			) == get(
				paste0("Cumulative ", mwf.name, ".y")
			))
		),
		by = c(names(dta.g)[!grepl(
			paste0("studyno",
						 "|",
						 "h.pred",
						 "|",
						 paste0("\\.[xy]$")
			),
			names(dta.g)
		)])
		]

		# No overlap for
		message('\t No overlap for ',
						signif(mean(dta.g.tab[,overlap == 0]), 3) * 100,
						"% of covariate histories")
	}

	# 2. Get predicted values ####
	# get discrete hazard (including those who died of other causes)
	# dta.a1$status <- 0
	dta[immortal == 0,`:=`(h.discrete = predict(
		h, newdata = model.matrix(discrete.formula, data = dta.a1[immortal == 0]), type = "response"))]

	if (total_effect) {
		dta[`Other event` == 1, h.discrete := 0]
	}

	dta[immortal == 0,`:=`(h.pred = shift(h.discrete, -1, NA)), by = .(studyno)]

	# set iterator q
	q <- 1

	# Subset
	dta <- dta[I <= J - q]

	# Loop over times
	while (q < J & !no_support) {
		# q <- q - 1
		k <- J - q

		if (!quiet) {
			cat(paste0("Iteration ", q, " of ", J - 1, "...\n"))
			cat(paste0("\t Getting cumulative hazard over [", k, ", ", J, "]\n"))
		}

		# # Make new covariates, if necessary
		# covariate.breaks <- get.covariate.breaks(dta[immortal == 0 & I == k + 1])
		# dta[, `:=`(
		# 	Age = get.cut(age.year2 / 365, covariate.breaks, "age"),
		# 	Year = get.cut(year, covariate.breaks),
		# 	`Duration of employment` = get.cut(employment.years, covariate.breaks),
		# 	`Year of hire` = get.cut(yin.gm, covariate.breaks, "yin"),
		# 	`Cumulative time off` = get.cut(cum_off, covariate.breaks, "cum_off"),
		# 	`Cumulative straight` = get.cut(cum_straight, mwf.breaks),
		# 	`Cumulative soluble` = get.cut(cum_soluble, mwf.breaks),
		# 	`Cumulative synthetic` = get.cut(cum_synthetic, mwf.breaks)
		# )]

		support <- sapply(c(timevar_covariates.which, baseline_covariates.which), function(x) {
			length(table(dta[status == 1 & I == k, x, with = F]))
		})

		if (sum(support <= 1)) {
			no_support <- T
			cat(paste0("Inadequate support at k = ", k))
		}

		# 3. Regress hazard from previous step ####
		# Cast time varying-covariates
		dta.wide <- dcast(
			dta[studyno %in% dta[immortal == 0 & I == k & status == 0 & Censored == 0, studyno]],
			studyno ~ I,
			value.var = timevar_covariates.which)

		# Merge baseline covariates
		dta.wide <- merge(
			dta[I == k, .(h.pred, `Year of hire`, Sex, Race, Plant), studyno],
			dta.wide,
			by = "studyno",
			all.x = F, all.y = T)

		h <- fastglm(
			model.matrix(h.pred ~ . - studyno, data = dta.wide),
			dta.wide$h.pred,
			family = "quasibinomial")

		na_coef <- is.na(coef(h))
		if (!quiet) {
			if (sum(na_coef) > 0) {
				message("Covariates dropped: ", paste0(names(coef(h))[na_coef], collapse = ", "))
			}
		}

		# 4. Use predicted hazard to obtain that to be regressed next ####
		# Get post-intervention data
		dta.a <- dta.a1[studyno %in% dta[immortal == 0 & I == k, studyno]]
		dta.a <- merge(
			dta[,.(`Year of hire` = `Year of hire`[1], Sex = Sex[1],
						 Race = Race[1], Plant = Plant[1]), studyno],
			dcast(dta.a, studyno ~ I,
						value.var = c(timevar_covariates.which, "Censored", "status")),
			by = "studyno",
			all.x = F, all.y = T)

		# Check positivity ####
		if (a != 0 & k > 1) {
			if (!quiet) {
				cat(paste0("\t Checking positivity over (", k, ", ", J, "]\n"))
			}
			dta.g <- merge(
				dta.a[get(paste0("Censored_", k - 1)) == 0 & get(paste0("status_", k - 1)) == 0,
							c("studyno", paste0("Cumulative ", mwf.name, "_", c(k - 1, k))),
							with = F],
				dta.wide
			)

			dta.g.tab <- dta.g[
				get(paste0("Cumulative ", mwf.name, "_", k - 1, ".x")) == get(
					paste0("Cumulative ", mwf.name, "_", k - 1, ".y")),
				.(n = n_distinct(studyno),
					overlap = sum(
						get(paste0("Cumulative ", mwf.name, "_", k, ".x")) == get(
							paste0("Cumulative ", mwf.name, "_", k, ".y")))
				),
				by = c(names(dta.g)[
					!grepl(paste0(
						"studyno", "h.pred",
						paste0("_", k, "$"), paste0("_", k, "\\.[xy]$"), paste0("_", k - 1, "\\.[xy]$"),
						sep = "|"),
						names(dta.g)
					)])
			]

			# No overlap for
			message('\t No overlap for ',
							signif(mean(dta.g.tab[,overlap == 0]), 3) * 100,
							"% of covariate histories")
			message('\t No overlap for ',
							sum(dta.g.tab[overlap == 0, n]),
							" people")
		}

		if (a != 0 & k == 1) {
			if (!quiet) {
				cat(paste0("\t Checking positivity over (", k - 1, ", ", J, "]\n"))
			}
			dta.g <- merge(
				dta.a[,c("studyno", paste0("Cumulative ", mwf.name, "_", k)), with = F],
				dta.wide
			)

			dta.g.tab <- dta.g[, .(
				n = n_distinct(studyno),
				overlap = sum(
					get(paste0("Cumulative ", mwf.name, "_", k, ".x")) == get(
						paste0("Cumulative ", mwf.name, "_", k, ".y")))
			),
			by = c(names(dta.g)[!grepl(
				paste0("studyno",
							 "|",
							 "h.pred",
							 "|",
							 paste0("_", k, "$"),
							 "|",
							 paste0("_", k, "\\.[xy]$")
				),
				names(dta.g)
			)])
			]

			# No overlap for
			message('\t No overlap for ',
							signif(mean(dta.g.tab[,overlap == 0]), 3) * 100,
							"% of covariate histories")
			message('\t No overlap for ',
							sum(dta.g.tab[overlap == 0,n]),
							" people")
		}

		# get predicted values
		if ("h.predict" %in% names(dta.a)) {
			dta.a <- dta.a[,-c("h.predict"), with = F]
		}

		h.predict <- predict(
			h,
			model.matrix(
				~ . - studyno,
				data = dta.a[, names(dta.wide)[names(dta.wide) != "h.pred"], with = F]),
			"response")

		# Get iterated expectation
		if ("h.predict" %in% names(dta)) {
			dta <- dta[,-c("h.predict"), with = F]
		}
		dta <- merge(
			dta,
			dta.a[,.(studyno, h.predict = h.predict)],
			by = "studyno", all.x = T)

		dta[immortal == 0 & I == max(k - 1, 1),
				`:=`(h.pred = h.predict * (1 - h.discrete) + h.discrete)]

		# Subset
		dta <- dta[I <= max(k - 1, 1)]

		# Update iterator
		q <- q + 1

	} # End loop over times

	if (!no_support) {return(dta[,.(studyno, h.pred)])
	} else {return(data.table(studyno = NA, h.pred = NA))}
}

get.bs <- function(dat = copy(dat.reduced),
									 dat.a1 = list(),
									 a1.names = NULL,
									 B,
									 mwf.name = list("soluble"),
									 run_intervention = T,
									 run_natural = T
) {
	require(lubridate, quietly = T)
	require(foreach, quietly = T)
	require(doParallel, quietly = T)
	registerDoParallel(detectCores() - 1)
	# start time
	start <- Sys.time()
	message("BS started at ", start, "\n")
	# # Progress bar
	# pb <- txtProgressBar(min = 0, max = B, style = 3)
	bs <- foreach (b = 1:B) %dopar% {
		# who.mcmc <- data.table(
		# 	studyno = sample(
		# 		unique(dat$studyno),
		# 		n_distinct(dat$studyno), replace = T))
		# who.mcmc$id <- 1:nrow(who.mcmc)

		# Get mcmc data (this may take a long time)
		data.mcmc <- merge(
			who.mcmc[[b]],
			dat[studyno %in% unique(who.mcmc[[b]]$studyno), covariates.which, with = F],
			by = "studyno",
			all.x = T, allow.cartesian = T)
		Sys.sleep(0)
		data.mcmc[, studyno := id]
		# data.mcmc <- data.mcmc[,-"id", with = F]

		if (run_intervention) {
			dat.a1 <- lapply(dat.a1, function(x) {
				x <- merge(
					who.mcmc[[b]],
					x[studyno %in% unique(who.mcmc[[b]]$studyno),],
					all.x = T, allow.cartesian = T)
				Sys.sleep(0)
				x[, studyno := id]
				return(x[,-"id", with = F])
			})
		}

		out <- cbind(
			b = b,
			if (run_natural) {mean(ice_gcomp(dta = data.mcmc, dta.a1 = NULL, quiet = T)$h.pred)},
			do.call(cbind, lapply(1:length(dat.a1), function(i) {
				return(mean(ice_gcomp(
					a = 1,
					dta = data.mcmc,
					dta.a1 = dat.a1[[i]],
					mwf.name = mwf.name[[i]], quiet = T)$h.pred))
			}))
		)

		if (run_natural) {colnames(out)[2] <- "natural"}
		colnames(out)[-(if (run_natural) {1:2} else {1})] <- a1.names

		return(as.data.table(out))

		# # Set progressbar
		# setTxtProgressBar(pb, b)
	}
	# bs
	time.unit <- "minutes"
	since_start <- lubridate::time_length(difftime(Sys.time(), start), time.unit)
	if (since_start > 90) {
		since_start <- since_start/60
		time.unit <- "hours"
	}

	cat(paste0("\n", round(since_start, 2), " ", time.unit, " since get.bs() was called.\n"))
	return(rbindlist(bs))
}