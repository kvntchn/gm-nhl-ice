# Define quantiles ####
get.covariate.breaks <- function(dat, cutpoints = c(100, 80, 60, 40) - 20) {
	covariate.breaks <- apply(dat[status == 1, .(
		year,
		yin = yin.gm,
		employment.years,
		age = age.year2 / 365,
		cum_off)], 2, function(x) {
			if (length(x) > cutpoints[1]) {
				breaks <- quantile(x, seq(0, 1, 1 / 6))
			} else if (length(x) > cutpoints[2]) {
				breaks <- quantile(x, seq(0, 1, 1 / 5))
			} else if (length(x) > cutpoints[3]) {
				breaks <- quantile(x, seq(0, 1, 1 / 4))
			} else if (length(x) > cutpoints[4]) {
				breaks <- quantile(x, seq(0, 1, 1 / 3))
			} else {
				breaks <- quantile(x, seq(0, 1, 1 / 2))
			}
			breaks[-length(breaks)] <- floor(breaks[-length(breaks)])
			breaks[length(breaks)] <- ceiling(breaks[length(breaks)])
			breaks[c(1, length(breaks))] <- c(-Inf, Inf)
			breaks
		})

	covariate.breaks <- as.data.table(covariate.breaks)

	return(covariate.breaks)
}

get.mwf.breaks <- function(dat) {
	mwf.breaks <- cbind(apply(dat[status == 1, .(
		straight,
		soluble,
		synthetic,
		off,
		cum_straight,
		# cum_soluble,
		cum_synthetic,
		cum_off)], 2, function(x) {
			x <- x[x > 0]
			if (length(x) > 40) {
				if (length(x) > 60) {
					probs <- seq(0, 1, 1 / 3)
				} else {
					probs <- seq(0, 1, 1 / 2)
				}
				breaks <- quantile(x, probs)
				breaks[c(1, length(probs))] <- c(0, Inf)
				breaks <- c(-Inf, breaks)
			} else {
				breaks <- c(-Inf, 0, Inf)
			}
			if (length(breaks) < 5) {
				breaks <- c(breaks, rep(NA, 5 - length(breaks)))
			}
			names(breaks) <- NULL
			breaks
		}),
		apply(dat[status == 1, .(cum_soluble)], 2, function(x) {
			x <- x[x > 0]
			if (length(x) > 40) {
				if (length(x) > 60) {
					probs <- seq(0, 1, 1 / 3)
				} else {
					probs <- seq(0, 1, 1 / 2)
				}
				breaks <- quantile(x[x > 0], probs)
				breaks[c(1, length(probs))] <- c(0, Inf)
				breaks <- c(-Inf, breaks)
			} else {
				breaks <- c(-Inf, 0, Inf)
			}
			if (length(breaks) < 5) {
				breaks <- c(breaks, rep(NA, 5 - length(breaks)))
			}
			breaks
		}),
		apply(dat[status == 1, .(cum_soluble5 = cum_soluble,
														 soluble5 = soluble)], 2, function(x) {
														 	x <- x[x > 0.05]
														 	if (length(x) > 40) {
														 		if (length(x) > 60) {
														 			probs <- seq(0, 1, 1 / 3)
														 		} else {
														 			probs <- seq(0, 1, 1 / 2)
														 		}
														 		breaks <- quantile(x[x > 0], probs)
														 		breaks[c(1, length(probs))] <- c(0.05, Inf)
														 		breaks <- c(-Inf, breaks)
														 	} else {
														 		breaks <- c(-Inf, 0.05, Inf)
														 	}
														 	if (length(breaks) < 5) {
														 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
														 	}
														 	breaks
														 }),
		apply(dat[status == 1, .(cum_soluble127 = cum_soluble,
														 soluble127 = soluble)], 2, function(x) {
														 	x <- x[x > 1.2742]
														 	if (length(x) > 40) {
														 		if (length(x) > 60) {
														 			probs <- seq(0, 1, 1 / 3)
														 		} else {
														 			probs <- seq(0, 1, 1 / 2)
														 		}
														 		breaks <- quantile(x[x > 0], probs)
														 		breaks[c(1, length(probs))] <- c(1.2742, Inf)
														 		breaks <- c(-Inf, breaks)
														 	} else {
														 		breaks <- c(-Inf, 1.2742, Inf)
														 	}
														 	if (length(breaks) < 5) {
														 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
														 	}
														 	breaks
														 }),
		apply(dat[status == 1, .(cum_soluble11 = cum_soluble,
														 soluble11 = soluble)], 2, function(x) {
														 	x <- x[x > 0.11]
														 	if (length(x) > 40) {
														 		if (length(x) > 60) {
														 			probs <- seq(0, 1, 1 / 3)
														 		} else {
														 			probs <- seq(0, 1, 1 / 2)
														 		}
														 		breaks <- quantile(x[x > 0], probs)
														 		breaks[c(1, length(probs))] <- c(0.11, Inf)
														 		breaks <- c(-Inf, breaks)
														 	} else {
														 		breaks <- c(-Inf, 0.11, Inf)
														 	}
														 	if (length(breaks) < 5) {
														 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
														 	}
														 	breaks
														 }),
		apply(dat[status == 1, .(cum_synthetic01 = cum_synthetic,
														 synthetic01 = synthetic)], 2, function(x) {
														 	x <- x[x > 0.0015]
														 	if (length(x) > 40) {
														 		if (length(x) > 60) {
														 			probs <- seq(0, 1, 1 / 3)
														 		} else {
														 			probs <- seq(0, 1, 1 / 2)
														 		}
														 		breaks <- quantile(x[x > 0], probs)
														 		breaks[c(1, length(probs))] <- c(0.0015, Inf)
														 		breaks <- c(-Inf, breaks)
														 	} else {
														 		breaks <- c(-Inf, 0.0015, Inf)
														 	}
														 	if (length(breaks) < 5) {
														 		breaks <- c(breaks, rep(NA, 5 - length(breaks)))
														 	}
														 	breaks
														 })
	)

	mwf.breaks <- as.data.table(mwf.breaks)

	return(mwf.breaks)}

get.cut <- function(x, breaks, y = NULL, include.lowest = T, dig.lab = 4) {
	if (is.null(y)) {y <- deparse(substitute(x))}
	cut(
		x,
		unlist(unique(na.exclude(breaks[, y, with = F]))),
		include.lowest = include.lowest,
		dig.lab = dig.lab
	)
}