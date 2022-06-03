# Blip function for stochastic intervention ####
# Kevin
# April 21, 2021

get.blip_intervention <- function(dta, mwf.intervention = mwf.which, a = 1, delta = 0.6) {

	if (a == 1) {
		for (mwf.which in mwf.intervention) {

			dta.a1 <- data.table::copy(dta)

			mwf <- unlist(dta.a1[, mwf.which, with = F])
			intervention.which <- which(mwf > min(mwf) + delta)

			mwf[intervention.which] <- mwf[intervention.which] - delta

			# Change continuous
			dta.a1[, (mwf.which) := list(mwf)]
			dta.a1[, (paste0("cum_", mwf.which)) := list(cumsum(get(mwf.which))), by = .(studyno)]

			dta.a1[,(c(str_to_title( mwf.which),
								 paste0("Cumulative ", mwf.which))) := list(
								 	get.cut(get(mwf.which), mwf.breaks, mwf.which),
								 	get.cut(get(paste0("cum_", mwf.which)), mwf.breaks, paste0("cum_", mwf.which))
								 )]

		}
		return(dta.a1)
	} else {return(dta)}
}

get.static_intervention <- function(dta, mwf.intervention = mwf.which, a = 0, delta = NULL) {

	a <- a + 1

	for (mwf.which in mwf.intervention) {

		dta.a1 <- data.table::copy(dta)
		mwf <- mwf.breaks[a:(a + 1), mwf.which, with = F]
		if (a == 1) {
			mwf <- max(mwf)
		} else {
			if (a >= length(table(dta.a1[,str_to_title(mwf.which), with = F]))) {
				mwf <- min(mwf)
			} else {
				mwf <- median(mwf)
			}
		}

		dta.a1[,(mwf.which) := list(mwf)]
		dta.a1[, (paste0("cum_", mwf.which)) := list(cumsum(get(mwf.which))), by = .(studyno)]
		# dta.a1[,.(studyno, soluble, Soluble, cum_soluble, `Cumulative soluble`)]
		dta.a1[,(c(str_to_title( mwf.which),
							 paste0("Cumulative ", mwf.which))) := list(
							 	get.cut(get(mwf.which), mwf.breaks, mwf.which),
							 	get.cut(get(paste0("cum_", mwf.which)), mwf.breaks, paste0("cum_", mwf.which))
							 )]

		dta[, A := as.numeric(as.numeric(get(str_to_title(mwf.which))) == a + 1)]

	}

	return(dta.a1)
}

get.shift_intervention <- function(
	dta = copy(dat.reduced), mwf.intervention = mwf.which, a = 1, delta = 0.1, dig.lab = 3) {

	if (a == 1) {
		for (mwf.which in mwf.intervention) {
			# mwf.which <- "soluble"

			mwf <- unlist(dta[, mwf.which, with = F])
			# mwf.cat <- unlist(dta[, str_to_title(mwf.which), with = F])
			intervention.which <- which(mwf > 0.05)

			new_mwf <- as.numeric(mwf[intervention.which]) - delta
			new_mwf[new_mwf < 0.05] <- 0.05
			mwf[intervention.which] <- new_mwf
				# unlist(mwf.breaks[is.finite(get(mwf.which)), mwf.which, with = F])[new_mwf]

			# Change continuous
			dta[, (mwf.which) := list(mwf)]
			dta[, (paste0("cum_", mwf.which)) := list(cumsum(get(mwf.which))), by = .(studyno)]

			dta[,(c(str_to_title( mwf.which),
								 paste0("Cumulative ", mwf.which))) := list(
								 	get.cut(get(mwf.which), mwf.breaks, mwf.which, dig.lab = dig.lab),
								 	get.cut(get(paste0("cum_", mwf.which)), mwf.breaks, paste0("cum_", mwf.which), dig.lab = dig.lab)
								 )]

		}
		return(dta.a1)
	} else {return(data.table::copy(dta))}
}