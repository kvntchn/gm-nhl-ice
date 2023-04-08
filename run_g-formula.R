# Run parametric g-formula ####
# Kevin Chen
# November 3, 2020

library(here)
library(boxr); box_auth()

# # Shape data
# source(here::here("shape-data.R"))
# dat <- box_read(1179471733155)
# dat <- dat[yin.gm >= 1965]
dat.reduced <- box_read(928445783979)
N <- length(unique(dat.reduced[!is.na(I)]$studyno))
message(sum(dat.reduced[!is.na(I)]$status))
source(here("g-formula.R"))

followed.names <- c("followed_200", "followed_100", "followed_50",
										"followed_25",  "followed_05",  "followed_025")
covariates.which <- c(
	followed.names,
	"studyno", "I", "year", "immortal",
	"status", "Censored", "Other event", "All causes",
	"age.year2", "Age", "Age at hire",
	"year2", "Year",
	"straight", "soluble", "synthetic",
	"Straight", "Soluble", "Synthetic",
	"Cumulative straight", "Cumulative synthetic",
	"Cumulative soluble", "Cumulative soluble shifted",
	"Baseline cumulative soluble",
	"Baseline cumulative straight",
	"Baseline cumulative synthetic",
	"cum_straight", "cum_soluble", "cum_synthetic",
	"Employment status",
	"Duration of employment", "employment.years",
	"Cumulative time off", "cum_off",
	"Year of hire", "yin.gm", "Sex", "Race", "Plant", 'py')

dat.reduced[,`:=`(
	`Baseline cumulative soluble` = `Cumulative soluble`[I == 1],
	`Baseline cumulative straight` = `Cumulative straight`[I == 1],
	`Baseline cumulative synthetic` = `Cumulative synthetic`[I == 1]
), studyno]

#####################################################
# Point estimates                                ####
#####################################################

# Baseline
baseline <- ice_gcomp(a = "Any")
message(round(mean(baseline$h.pred) * N))

# Intervention
reduction.sol200 <- ice_gcomp(a = "200")
message(round(mean(reduction.sol200$h.pred) * N))

reduction.sol100 <- ice_gcomp(a = "100")
message(round(mean(reduction.sol100$h.pred) * N))

reduction.sol50 <- ice_gcomp(a = "50")
message(round(mean(reduction.sol50$h.pred) * N))

reduction.sol25 <- ice_gcomp(a = "25")
message(round(mean(reduction.sol25$h.pred) * N))

reduction.sol05 <- ice_gcomp(a = "05")
message(round(mean(reduction.sol05$h.pred) * N))

# Save results
estimates <- list(
	baseline = baseline$h.pred,
	reduction.sol200 = reduction.sol200$h.pred,
	reduction.sol100 = reduction.sol100$h.pred,
	reduction.sol50 = reduction.sol50$h.pred,
	reduction.sol25 = reduction.sol25$h.pred,
	reduction.sol05 = reduction.sol05$h.pred,
	observed = sum(dat.reduced$status)/n_distinct(dat.reduced$studyno)
)

save(estimates, file = here::here("resources", paste0("estimates.rdata")))

load(here::here("resources", paste0("estimates.rdata")))

sapply(estimates, mean) * N

#####################################################
# MCMC variance estimate                         ####
#####################################################
B <- 1e3

# # Get MCMC people
# get_mcmc_list <- function(B, id = unique(dat.reduced$studyno)) {
# 	pb <- txtProgressBar(min = 0, max = 1e3, style = 3)
# 	who.mcmc.list <- lapply(1:B, function(b) {
# 		who.mcmc <- data.table(studyno = sample(id, replace = T))
# 		who.mcmc$id <- 1:nrow(who.mcmc)
# 		setTxtProgressBar(pb, b)
# 		return(who.mcmc)
# 	})
# 	return(who.mcmc.list)
# }
#
# set.seed(232)
# who.mcmc.list <- get_mcmc_list(B)
# saveRDS(who.mcmc.list, here::here("resources/private", "who.mcmc.list.rds"))
# who.mcmc.list2 <- get_mcmc_list(B)
# saveRDS(who.mcmc.list2, here::here("resources/private", "who.mcmc.list2.rds"))

# who.mcmc <- readRDS(here::here("resources/private", "who.mcmc.list.rds"))
# # system.time({
# estimates.bs <- get.bs(
# 	dat = data.table::copy(dat.reduced),
# 	B = B,
# 	interventions = list("Any", "200", "100", "50", "25", "05"))
# # })
# saveRDS(estimates.bs, here::here("resources", paste0("estimates.bs.rds")))

# Examine BS sampling distribution ####
estimates.bs <- readRDS(here::here("resources", paste0("estimates.bs.rds")))


do.call(cbind, lapply(estimates.bs[,-(1:2)], function(x) {
	x/estimates.bs[,2]
})) %>% melt(measure.vars = names(.)) -> rr.bs

rr.point_estimates <- data.frame(
	variable = paste0(c(
		"200", "100", "50", "25", "05"
	), ".Any"),
	estimate = sapply(estimates[2:6], mean) / mean(estimates[[1]])
)

rr <- merge(rr.bs,
						rr.point_estimates,
						by = "variable")

rr[,.(estimate = estimate[1],
			lower = 2 * estimate[1] - quantile(value, 1 - 0.025),
			upper = 2 * estimate[1] - quantile(value, 0.025)
),
variable][c(3, 2, 5, 4, 1),]

do.call(cbind, lapply(estimates.bs[,-(1:2)], function(x) {
	x - estimates.bs[,2]
})) %>% melt(measure.vars = names(.)) -> rd.bs

rd.point_estimates <- data.frame(
	variable = paste0(c(
		"200", "100", "50", "25", "05"
	), ".Any"),
	estimate = sapply(estimates[2:6], mean) - mean(estimates[[1]])
)

rd <- merge(rd.bs,
						rd.point_estimates,
						by = "variable")

rd[,.(estimate = estimate[1] * N,
			lower = (2 * estimate[1] - quantile(value, 1 - 0.025)) * N,
			upper = (2 * estimate[1] - quantile(value, 0.025)) * N
			# wald.lower = (estimate[1] + qnorm(0.025) * sd(value)) * N,
			# wald.upper = (estimate[1] + qnorm(1 - 0.025) * sd(value)) * N
),
variable][c(3, 2, 5, 4, 1),]

# # Plot sampling distributions ####
# rr.bs %>% ggplot(
# 	aes(x = log(value))) +
# 	geom_histogram(bins = 50) +
# 	geom_vline(xintercept = 0, color = "salmon") +
# 	facet_wrap(. ~ variable) +
# 	mytheme
#
# rd.bs %>% ggplot(
# 	aes(x = value)) +
# 	geom_histogram(bins = 50) +
# 	geom_vline(xintercept = 0, color = "salmon") +
# 	facet_wrap(. ~ variable) +
# 	mytheme
