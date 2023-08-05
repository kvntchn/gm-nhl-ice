# Parametric g-formula ####
# Kevin Chen
# November 3, 2020

library(here)
library(boxr); box_auth()
library(survival)

# Get raw data and helper functions ####
race.imputed.melt <- NULL
source("~/headRs/00-my-theme.R")
library(data.table)
# # paste(ls(), collapse = ",\n") %>% clipr::write_clip()
source(here::here("breaks.R"))

#####################################################
# Iterated conditional expectation g-computation ####
#####################################################
# Hazard-extended ICE implementation of g-formula ####
# Wen et al. (2020)
ice_gcomp <- function(
		a = "50",
		dta = copy(dat.reduced[!is.na(I), covariates.which, with = F]),
		times = sort(unique(dat.reduced[status == 1]$I), decreasing = T),
		quiet = F
) {

	# Indicator for natural course (everone followed)
	dta$followed_Any <- 1

	# Index follow-up periods
	setorder(dta, studyno, year)
	dta[, I := 1:.N, by = studyno]

	# Pick covariates
	timevar_covariates.which <- c(
		"Employment status",
		"Soluble",
		NULL
	)

	# Static covariates for discrete hazard and cumulative hazard
	static_covariates.which <- c(
		# "Employment status",
		"Cumulative straight",
		"Cumulative synthetic",
		"Baseline cumulative soluble",
		"Age",
		"Race", "Plant", "Sex",
		"Cumulative time off",
		"Duration of employment",
		NULL
	)

	# Intervention node name
	intervene.which <- "Cumulative soluble"

	# History to condition on for intervention
	dynamic.which <- c(
		"Employment status",
		"Cumulative straight",
		"Cumulative synthetic",
		"Cumulative time off",
		"Age",
		"Duration of employment",
		"Race",
		"Plant",
		"Sex",
		NULL
	)

	# 1. Regression for discrete hazard ####
	# Total number of outcomes
	J <- max(times)
	require(fastglm, quietly = T)

	# Formula for discrete hazard regression
	discrete.formula <- as.formula(paste(
		"status",
		"~",
		"`Employment status` +",
		"`Cumulative soluble` +",
		"`Cumulative straight` +",
		"`Cumulative synthetic` +",
		"`Cumulative time off` +",
		"Age +",
		"`Duration of employment` +",
		"Race +",
		"Plant +",
		"Sex",
		NULL
	))

	h <- fastglm(
		model.matrix(discrete.formula, data = dta)[with(dta, immortal == 0 & Censored == 0),],
		dta[immortal == 0 & Censored == 0]$status,
		family = "binomial"
	)
	# sum(h$fitted.values)
	# summary(h)

	# 2. Get predicted values ####
	# get discrete hazard (including those who died of other causes)
	dta[immortal == 0, `:=`(
		h.discrete = predict(
			h,
			newdata = model.matrix(discrete.formula, data = dta)[with(dta, immortal == 0),],
			type = "response"))]

	# dta[dta[,.(always_followed = as.numeric(all(get(paste0("followed_", a)) == 1))),
	# 		eval(dynamic.which)],
	# 		on = dynamic.which])

	if (a != "Any") {
		n_p <- length(levels(unlist(dta[, intervene.which, with = F])))
		dta.probs <- lapply(1:J, function(k = 1) {

			# Get exposures among followers
			dta.probs_followed <- dta[
				get(paste0("followed_", a)) == 1 & I == k, .(
					Exposure = get(intervene.which)),
				eval(dynamic.which)]

			# What is the highest level of rule-abiding exposure
			dta.probs_barely_followed <- dta.probs_followed[,.(
				Exposure = levels(Exposure)[max(as.numeric(Exposure))]
			),
			eval(dynamic.which)]

			# Get non-followers
			dta.probs_not_followed <- dta[
				get(paste0("followed_", a)) == 0 & I == k,
				c(dynamic.which, intervene.which), with = F]

			# Get lowest exposure among non-followers
			dta.probs_not_followed[,`:=`(
				Exposure = levels(get(intervene.which))[
					min(as.numeric(get(intervene.which)))]
			), eval(dynamic.which)]

			# Add rule-abiding exposure if positivity allows
			dta.probs_not_followed <- dta.probs_barely_followed[
				dta.probs_not_followed,
					on = dynamic.which
				]
			dta.probs_not_followed[is.na(Exposure), Exposure := i.Exposure]

			# Combine distribution to get the intervention distribution
			dta.probs <- rbindlist(list(
				dta.probs_followed, dta.probs_not_followed[,names(dta.probs_followed), with = F]
			))

			# Get counts within exposure levels
			dta.probs <- lapply(
				levels(dta.probs_followed$Exposure), function(x) {
					probs.tmp <- dta.probs_followed[, .(sum(Exposure == x)),
																 eval(dynamic.which)]
					names(probs.tmp)[length(dynamic.which) + 1] <- x
					return(probs.tmp)
				})

			# Reduce
			dta.probs <- as.data.frame(Reduce(function(...) {
				merge(..., all = T, by = dynamic.which)}, dta.probs))

			# Compute probabilities
			dta.probs[,-(1:length(dynamic.which))] <- dta.probs[,-(1:length(dynamic.which))] /
				matrix(rep(
					rowSums(dta.probs[,-(1:length(dynamic.which))]),
					ncol(dta.probs) - length(dynamic.which)), nrow = nrow(dta.probs))
			names(dta.probs)[-(1:length(dynamic.which))] <- paste0("p", 1:n_p)
			names(dta.probs)[1:length(dynamic.which)] <- dynamic.which
			dta.probs$I <- k
			return(dta.probs)
		})

		dta.probs <- rbindlist(dta.probs)
		# dta[, paste0(dynamic.which, "_og") := get(dynamic.which)]
		# for (k in 1:J) {
		# 	dta[I == k, (dynamic.which) :=  approxfun(
		# 		dta.probs[I == k, 1], rule = 2, method = "constant")(get(dynamic.which))]
		# }

		dta.join <- dta.probs[dta, on = c("I", paste0(dynamic.which))]
		dta.join[is.na(p1), paste0("p", 1:n_p) := lapply(
			1:n_p, function(k) {as.numeric(
				get(intervene.which) == levels(get(intervene.which))[k]
			)}
		)]

		# dta[, (dynamic.which) := get(paste0(dynamic.which, "_og"))]
		for (j in 1:length(levels(dta[[intervene.which]]))) {
			dta.join[, paste0(intervene.which) := levels(
				dta[[paste0(intervene.which)]])[j]]
			dta.join[, paste0('h.discrete', j) := list(
				predict(
					h,
					newdata = model.matrix(discrete.formula, data = dta.join),
					type = "response"))]
		}

		h.discrete <- rowSums(
			dta.join[,paste0("h.discrete", 1:n_p), with = F] *
				dta.join[,paste0("p", 1:n_p), with = F])
		# mean(h.discrete)/mean(dta.join$h.discrete)
		dta.join$h.discrete <- h.discrete
		dta.join <- dta.join[,c(
			"studyno", "I", paste0("h.discrete", c("", 1:n_p))), with = F]
		setorder(dta.join, studyno, I)

		if (all.equal(dta[,.(studyno, I)], dta.join[,.(studyno, I)])) {
			dta[, paste0("h.discrete", c("", 1:n_p)) := as.list(
				dta.join[, paste0("h.discrete", c("", 1:n_p)), with = F])]
		} else {stop("check code")}
	}

	dta[immortal == 0,`:=`(
		h.pred = shift(h.discrete, -1, NA)),
		studyno]

	# set iterator q
	q <- 1

	# Subset
	dta <- dta[I <= J - q]

	intervene.which <- "Soluble"

	# Loop over times ####
	total_py <- sum(dta$py)
	intervened_py <- 0
	while (q < J) {
		# q <- q - 1
		k <- J - q
		# sum(dta[I == k, h.pred], na.rm = T)
		# sum(dat.reduced[I >= k, status])

		if (intervene.which %in% timevar_covariates.which) {
			intervene.which.k <- paste0(intervene.which, "_", k)
		} else {
			intervene.which.k <- intervene.which
		}

		if (k == 1) {
			intervene.which.k <- "Baseline cumulative soluble"
			timevar_covariates.which <- timevar_covariates.which[
				!timevar_covariates.which == "Soluble"]
		}

		dynamic.which <- c(
			if (k != 1) {"Soluble"} else {NULL},
			if (k != 1) {"Baseline cumulative soluble"} else {NULL},
			if (k == 1) {"Duration of employment"} else {NULL},
			"Employment status",
			# if (k != 1) {"Baseline cumulative straight"} else {NULL},
			"Cumulative straight",
			# if (k != 1) {"Baseline cumulative synthetic"} else {NULL},
			"Cumulative synthetic",
			"Cumulative time off",
			"Age",
			"Race",
			"Plant",
			"Sex",
			NULL
		)

		dynamic.which.static <- dynamic.which[!dynamic.which %in% timevar_covariates.which]
		dynamic.which.timevar <- dynamic.which[dynamic.which %in% timevar_covariates.which]

		dynamic.which <- c(
			dynamic.which.static,
			sapply(dynamic.which.timevar, paste0, "_", 1:k)
		)

		names(dynamic.which) <- NULL

		dynamic.which <- dynamic.which[!dynamic.which %in% intervene.which.k]

		if (!quiet) {
			# cat(paste0("Iteration ", q, " of ", J - 1, "...\n"))
			cat(paste0("\t Getting cumulative hazard over [", k, ", ", J, "]\n"))
		}

		# 3. Regress hazard from previous step ####
		# Cast time varying-covariates
		if (length(timevar_covariates.which) == 0) {
			dta.wide <- dta[immortal == 0 & ((I == k - 1 & Censored == 0 & status == 0) | k == 1),
											.(studyno = unique(studyno))]
		} else {
			dta.wide <- dcast(
				dta[studyno %in% dta[immortal == 0 & ((I == k - 1 & Censored == 0 & status == 0) | k == 1), studyno]],
				studyno ~ I,
				value.var = timevar_covariates.which)
			if (length(timevar_covariates.which) == 1) {
				names(dta.wide)[-1] <- paste0(timevar_covariates.which,  "_", names(dta.wide)[-1])
			}}

		# Followers
		dta.followed <- dta[I == k & studyno %in% dta.wide$studyno, .(
			followed = as.numeric(all(get(paste0("followed_", a)) == 1))),
			studyno]

		# Merge baseline covariates
		dta.wide <- merge(
			dta[I == k, unique(c(
				"studyno", "py",
				"h.pred",
				static_covariates.which, dynamic.which.static)), with = F][
					dta.followed,
					on = "studyno"],
			dta.wide,
			by = "studyno",
			all.x = F, all.y = T)

		h <- fastglm(
			model.matrix(
				h.pred ~ .,
				data = dta.wide[
					!studyno %in% dta[I == k & (Censored == 1 | status == 1), studyno],
					c("h.pred", static_covariates.which,
						sapply(timevar_covariates.which, function(x) {
							grep(paste0("^", x, "_"), names(dta.wide), value = T)})),
					with = F]),
			dta.wide[!studyno %in% dta[I == k & (Censored == 1 | status == 1), studyno], h.pred],
			family = "quasibinomial")
		# sum(h$fitted.values)
		# summary(h)

		# 4. Use predicted hazard to obtain that to be regressed next ####
		dta.wide[is.na(h.pred), h.pred := 0]
		dta.wide$h.predict <- predict(h, newdata = model.matrix(
			h.pred ~ .,
			data = dta.wide[,c("h.pred", static_covariates.which,
												 sapply(timevar_covariates.which, function (x) {
												 	grep(paste0("^", x, "_"),
												 			 names(dta.wide), value = T)})), with = F]), type = "response")
		setorder(dta.wide, studyno)
		dta[immortal == 0 & I == max(k - 1, 1) & studyno %in% dta.wide[,studyno],
				h.pred := dta.wide$h.predict * (1 - h.discrete) + h.discrete]

		if (a != "Any") {
			n_p <- length(levels(unlist(dta.wide[, intervene.which.k, with = F])))

			# Get exposures among followers
			dta.probs_followed <- dta.wide[
				followed == 1, .(Exposure = get(intervene.which.k)),
				eval(dynamic.which)]

			# What is the highest level of rule-abiding exposure
			dta.probs_barely_followed <- dta.probs_followed[,.(
				Exposure = levels(Exposure)[max(as.numeric(Exposure))]
			),
			eval(dynamic.which)]

			# Get non-followers
			dta.probs_not_followed <- dta.wide[
				followed == 0 ,
				c(dynamic.which, intervene.which.k), with = F]

			# Get lowest exposure among non-followers
			dta.probs_not_followed[,`:=`(
				Exposure = levels(get(intervene.which.k))[
					min(as.numeric(get(intervene.which.k)))]
			), eval(dynamic.which)]

			# Add rule-abiding exposure if positivity allows
			dta.probs_not_followed <- dta.probs_barely_followed[
					dta.probs_not_followed,
					on = dynamic.which
				]
			dta.probs_not_followed[is.na(Exposure), Exposure := i.Exposure]

			# Combine distribution to get the intervention distribution
			dta.probs <- rbindlist(list(
				dta.probs_followed, dta.probs_not_followed[,names(dta.probs_followed), with = F]
			))

			dta.probs <- lapply(
				levels(dta.probs$Exposure), function(x) {
					probs.tmp <- dta.probs[, .(sum(Exposure == x)), eval(dynamic.which)]
					names(probs.tmp)[length(dynamic.which) + 1] <- x
					return(probs.tmp)
				})
			dta.probs <- as.data.frame(Reduce(function(...) {
				merge(..., all = T, by = dynamic.which)}, dta.probs))
			dta.probs[,-(1:length(dynamic.which))] <- dta.probs[,-(1:length(dynamic.which))] /
				matrix(rep(
					rowSums(dta.probs[,-(1:length(dynamic.which))]),
					n_p), nrow = nrow(dta.probs))
			names(dta.probs)[-(1:length(dynamic.which))] <- paste0("p", 1:n_p)
			names(dta.probs)[1:length(dynamic.which)] <- paste0(dynamic.which)
			dta.probs <- as.data.table(dta.probs)

			# if (k == 1) {
			# 	dta.wide[, `Duration of employment` :=  approxfun(
			# 		dta.probs[, `Duration of employment`], rule = 2, method = "constant")(
			# 			`Duration of employment`)]
			# }

			dta.join <- dta.probs[dta.wide[,-'h.predict'], on = c(dynamic.which)]
			intervened_py <- sum(dta.join[!is.na(p1), py])
			if (!quiet) {
				# cat(paste0("Iteration ", q, " of ", J - 1, "...\n"))
				cat(paste0("\t\t No intervention on ",
									 round(nrow(dta.join[is.na(p1)])/nrow(dta.join) * 100, 2),
									 "\n"))
			}
			# dta.join$`Duration of employment`
			dta.join[is.na(p1), paste0("p", 1:n_p) := lapply(
				1:n_p, function(k) {as.numeric(
					get(intervene.which.k) == levels(get(intervene.which.k))[k]
				)}
			)]

			# dta.wide[, (dynamic.which) := get(paste0(dynamic.which, "_og"))]

			for (j in 1:length(levels(dta.wide[[intervene.which.k]]))) {
				dta.join[, (intervene.which.k) := levels(
					dta.wide[[intervene.which.k]])[j]]
				dta.join[, paste0("h.predict", j) := predict(
					h, newdata = model.matrix(
						h.pred ~ .,
						data = dta.join[
							,
							c("h.pred", static_covariates.which,
								sapply(timevar_covariates.which, function(x) {
									grep(paste0("^", x, "_"),
											 names(dta.wide), value = T)
								})
							),
							with = F]), type = "response")]
			}

			dta.join <- dta.join[
				dta[immortal == 0 & I == max(k - 1, 1) & studyno %in% dta.wide[,studyno], c(
					"studyno",
					paste0("h.discrete", 1:n_p)
				), with = F],
				on = "studyno"]
			dta.join$h.pred <- rowSums(
				(
					dta.join[,paste0("h.predict", 1:n_p), with = F] * (
						1 - dta.join[,paste0("h.discrete", 1:n_p), with = F]
					) + dta.join[,paste0("h.discrete", 1:n_p), with = F]
				) *
					dta.join[,paste0("p", 1:n_p), with = F])
			dta.join <- dta.join[,.(studyno, h.pred)]
			setorder(dta.join, studyno)

			if (all.equal(dta.join$studyno,
										dta[immortal == 0 & I == max(k - 1, 1) & studyno %in% dta.join[,studyno], studyno])) {
				dta[immortal == 0 & I == max(k - 1, 1) & studyno %in% dta.join[,studyno],
						h.pred := dta.join$h.pred]
			} else {stop("check code")}

		}

		# Subset
		dta <- dta[I <= max(k - 1, 1)]

		# Update iterator
		q <- q + 1

	} # End loop over times

	return(dta[,.(studyno, h.pred, total_py = total_py, intervened_py = intervened_py)])
}

get.bs <- function(dat = copy(dat.reduced),
									 interventions = c("Any", "5", "25", "05"),
									 B = 10
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

		out <- cbind(
			b = b,
			do.call(cbind, lapply(1:length(interventions), function(i) {
				return(mean(ice_gcomp(
					dta = copy(data.mcmc),
					a = interventions[i],
					quiet = T)$h.pred))
			}))
		)

		colnames(out)[-1] <- unlist(interventions)

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
