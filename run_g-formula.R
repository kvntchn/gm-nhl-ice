# Run parametric g-formula ####
# Kevin Chen
# November 3, 2020

library(here)
source(here("g-formula.R"))

# # Shape data
# source(here::here("shape-data.R"))

dat.reduced <- box_read(928445783979)

# name          : dat.sol.reduced.rds
dat.sol.reduced <- box_read(928445303042)
# name          : dat.sol.reduced2.rds
dat.sol.reduced2 <- box_read(929953044799)
# name          : dat.sol.reduced10.rds
dat.sol.reduced10 <- box_read(930687554963)

#####################################################
# Post-intervention distribution of exposure     ####
#####################################################
# # Proportion of person-years intervened upon
# interventions <- data.frame(
# 	dt.name = c(
# 		"dat.sol.reduced", "dat.sol.reduced2", "dat.sol.reduced10"
# 	)
# )
# interventions$mwf <- sapply(
# 	substr(interventions$dt.name, 5, 7),
# 	grep, x = c("straight", "soluble", "synthetic"),
# 	value = T)
#
# interventions$percent_py_intervened <- apply(interventions, 1, function(x) {
# 	dat <- get(x[1])
# 	mwf <- x[2]
# 	dat[,(paste0("Cumulative ", mwf, " (new)")):=list(get(paste0("Cumulative ", mwf)))]
# 	merge(dat.reduced,
# 				dat[,c("studyno", "I", paste0("Cumulative ", mwf, " (new)")), with = F],
# 				by = c("studyno", "I"))[,.(
# 					V1 = sum(py[get(paste0("Cumulative ", mwf)) != get(paste0("Cumulative ", mwf, " (new)"))])/sum(py)
# 				)]$V1
# })
#
# saveRDS(interventions, here::here("resources", "percent-py-intervened.rds"))

exposure.comparison <- rbindlist(c(
	list(dat.reduced[, .(
		studyno, I, status, year, py, `Cumulative soluble`, cum_soluble, soluble)]),
	lapply(list(dat.sol.reduced, dat.sol.reduced2, dat.sol.reduced10
	), function(x) {
		x[,.(
			studyno, I, status, year, py, `Cumulative soluble`, cum_soluble, soluble)]
	})),
	idcol = "Intervention"
)

exposure.comparison[Intervention == 1,
										Type := "None"]
exposure.comparison[Intervention %in% 2:4,
										Type := "Static"]
exposure.comparison[Intervention %in% 5:7,
										Type := "Dynamic"]

exposure.comparison <- rbindlist(list(
	exposure.comparison[Type == "None",.(
		Intervention, studyno, I, status, year, py,
		`Cumulative soluble`, cum_soluble, soluble, Type = "Static"
	)],
	exposure.comparison
))

exposure.comparison[Type == "None", Type := "Dynamic"]

exposure.comparison[,Type := factor(Type, c("Static", "Dynamic"), c("A) Static intervention", "B) Dynamic intervention"))]

exposure.comparison[,`:=`(Intervention = factor(
	Intervention, labels = c(
		"No limit",
		rep(c("0.5", "0.25", "0.05"), 1)))
)]

exposure.comparison[,`:=`(year.max = max(year), I.max = max(I)), studyno]

dcast(exposure.comparison[,.(py = sum(py)), .(`Cumulative soluble`, Intervention)],
			`Cumulative soluble` ~ Intervention,
			value.var = "py")
table(exposure.comparison[status == 1,.(`Cumulative soluble`, Intervention)])

tikzDevice::tikz(here::here("resources", "shift-cum-soluble.tex"), height = 2.75, width = 6, standAlone = T)
exposure.comparison[I == I.max] %>% ggplot(
	aes(x = cum_soluble)
	#, color = '#102451'
) +
	geom_density(aes(lty = Intervention)
							 # , color = '#102451'
	) +
	# geom_histogram(aes(fill = Intervention), bins = 40, alpha = 0.4, position = 'identity') +
	# facet_wrap(. ~ Intervention, ncol = 1) +
	scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash", "longdash", "dotted", "dotdash")) +
	mytheme +
	scale_x_continuous(
		n.breaks = 6,
		trans = "log",
		labels = function(x) {sprintf("%.2f", x)}) +
	geom_rug(data = exposure.comparison[status == 1, .(cum_soluble)],
					 alpha = 0.5, size = 0.05, length = grid::unit(3.5, "pt")
					 # color = '#102451'
	) +
	# facet_wrap(. ~ Type) +
	labs(
		y = "Density",
		x = "Cumulative exposure (mg/m$^3\\cdot\\,$years)",
		lty = "Hypothetical exposure limit (mg/m$^3\\cdot\\,$years): ") + theme(
			axis.title = element_text(margin = margin(10, 5, 10, 5)),
			# legend.title = element_blank(),
			legend.position = "bottom",
			legend.box.spacing = unit(2, "pt"),
			legend.key.size = unit(10, "pt"),
			strip.text = element_text(hjust = 0)
		) +
	theme(panel.grid = element_blank()) +
	coord_cartesian(ylim = c(0, 0.5),
									xlim = c(1e-4, 805))
dev.off()
lualatex("shift-cum-soluble\\.tex", here::here('resources'))

tikzDevice::tikz(here::here("resources", "shift-avg-soluble.tex"), height = 2.75, width = 6, standAlone = T)
exposure.comparison %>% ggplot(
	aes(x = soluble)
	# color = '#102451'
) +
	geom_density(aes(lty = Intervention),
							 bw = 0.02
	) +
	# geom_histogram(aes(fill = Intervention), bins = 40, alpha = 0.4, position = 'identity') +
	# facet_wrap(. ~ Intervention, ncol = 1) +
	scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash", "longdash", "dotted", "dotdash")) +
	mytheme +
	scale_x_continuous(
		n.breaks = 4,
		trans = "sqrt",
		breaks = c(0.01, 0.2, 0.5, 1, 2, 3),
		labels = function(x) {sprintf("%.2f", x)}) +
	facet_wrap(. ~ Type) +
	labs(
		y = "Density",
		x = "Average annual exposure (mg/m$^3$)",
		lty = "Hypothetical exposure limit (mg/m$^3$): ") + theme(
			axis.title = element_text(margin = margin(10, 5, 10, 5)),
			# legend.title = element_blank(),
			legend.position = "bottom",
			legend.box.spacing = unit(2, "pt"),
			legend.key.size = unit(10, "pt")
		) +
	theme(panel.grid = element_blank(),
				strip.text = element_text(hjust = 0)) +
	coord_cartesian(ylim = c(0, 3.5), xlim = c(-0.1, 3.5))
dev.off()
lualatex("shift-avg-soluble\\.tex", here::here('resources'))

# #####################################################
# # Point estimates                                ####
# #####################################################
#
# Baseline
baseline <- ice_gcomp()

# Intervention
reduction.sol <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol.reduced), mwf.name = "soluble")
reduction.sol2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol.reduced2), mwf.name = "soluble")
reduction.sol10 <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol.reduced10), mwf.name = "soluble")
observed <- sum(dat.reduced$status)/n_distinct(dat.reduced$studyno)

# Save results
estimates <- list(
	baseline = baseline,
	reduction.sol = reduction.sol,
	reduction.sol2 = reduction.sol2,
	reduction.sol10 = reduction.sol10,
	observed = sum(dat.reduced$status)/n_distinct(dat.reduced$studyno)
)

save(
	baseline,
	reduction.sol,
	reduction.sol2,
	reduction.sol10,
	observed, file = here::here("resources", paste0("estimates.rdata")))

h.baseline <- mean(baseline$h.pred)
h.reduction.sol <- mean(reduction.sol$h.pred)
h.reduction.sol2 <- mean(reduction.sol2$h.pred)
h.reduction.sol10 <- mean(reduction.sol10$h.pred)

save(
	h.baseline,
	h.reduction.sol,
	h.reduction.sol2,
	h.reduction.sol10,
	file = here::here("resources", paste0("h.estimates.rdata")))

load(here::here("resources", paste0("estimates.rdata")))
load(here::here("resources", paste0("h.estimates.rdata")))

# #####################################################
# # MCMC variance estimate                         ####
# #####################################################
B <- 1e3
#
# # # Get MCMC people
# # set.seed(232)
# # pb <- txtProgressBar(min = 0, max = 1e3, style = 3)
# # who.mcmc.list <- lapply(1:1e3, function(b) {
# # 	who.mcmc <- data.table(studyno = sample(unique(dat.reduced$studyno), replace = T))
# # 	who.mcmc$id <- 1:nrow(who.mcmc)
# # 	setTxtProgressBar(pb, b)
# # 	return(who.mcmc)
# # })
# # saveRDS(who.mcmc.list, here::here("resources/private", "who.mcmc.list.rds"))
# # saveRDS(who.mcmc.list2, here::here("resources/private", "who.mcmc.list2.rds"))
# #
# who.mcmc <- readRDS(here::here("resources/private", "who.mcmc.list.rds"))
#
# # system.time({
# estimates.bs <- get.bs(
# 	dat = data.table::copy(dat.reduced), B = B,
# 	run_natural = F,
# 	dat.a1 = list(
# 		dat.str.reduced, dat.str_total.reduced,
# 		dat.str.reduced2, dat.str_total.reduced2,
# 		dat.str.reduced10,
# 		dat.str_total.reduced10,
# 		dat.sol.reduced, dat.sol_total.reduced,
# 		dat.sol.reduced2, dat.sol_total.reduced2,
# 		dat.sol.reduced10,
# 		dat.sol_total.reduced10,
# 		dat.syn.reduced, dat.syn_total.reduced,
# 		dat.syn.reduced2, dat.syn_total.reduced2,
# 		dat.syn.reduced10,
# 		dat.syn_total.reduced10
# 	),
# 	a1.names = c(
# 		"straight", "total by straight",
# 		"straight2", "total by straight2",
# 		"straight10",
# 		"total by straight10",
# 		"soluble", "total by soluble",
# 		"soluble2", "total by soluble2",
# 		"soluble10",
# 		"total by soluble10",
# 		"synthetic", "total by synthetic",
# 		"synthetic2", "total by synthetic2",
# 		"synthetic10",
# 		"total by synthetic10"),
# 	mwf.name = list(
# 		"straight", "straight",
# 		"straight", "straight",
# 		"straight",
# 		"straight",
# 		"soluble", "soluble",
# 		"soluble", "soluble",
# 		"soluble",
# 		"soluble",
# 		"synthetic", "synthetic",
# 		"synthetic", "synthetic",
#     "synthetic",
# 		"synthetic"
# ))
# # })
#
# saveRDS(estimates.bs, here::here("resources", paste0("estimates.bs.rds")))

# Examine BS sampling distribution ####
estimates.bs <- readRDS(here::here("resources", paste0("estimates.bs.rds")))
do.call(cbind, lapply(estimates.bs[,-(1:2)], function(x) {
	x/estimates.bs[,2]
})) %>% melt(measure.vars = names(.)) -> rr.bs

rr.point_estimates <- data.frame(
	variable = paste0(c(
		"soluble", "soluble2", "soluble10"
	), ".natural"),
	estimate = c(
		h.reduction.sol, h.reduction.sol2, h.reduction.sol10
	)/h.baseline
)

rr <- merge(rr.bs,
						rr.point_estimates,
						by = "variable")

rr[,.(estimate = estimate[1],
			lower = quantile(value, 0.025) - mean(value) + estimate[1],
			upper = quantile(value, 1 - 0.025) - mean(value) + estimate[1]
			# wald.lower = exp(log(estimate[1]) + qnorm(0.025) * sd(log(value))),
			# wald.upper = exp(log(estimate[1]) - qnorm(0.025) * sd(log(value)))
),
variable]

do.call(cbind, lapply(estimates.bs[,-(1:2)], function(x) {
	x - estimates.bs[,2]
})) %>% melt(measure.vars = names(.)) -> rd.bs

rd.point_estimates <- data.frame(
	variable = paste0(c(
		"soluble", "soluble2", "soluble10"), ".natural"),
	estimate = c(
		h.reduction.sol, h.reduction.sol2, h.reduction.sol10) - h.baseline
)

rd <- merge(rd.bs,
						rd.point_estimates,
						by = "variable")

rd[,.(estimate = estimate[1] * 34738,
			lower = (quantile(value, 0.025) - mean(value) + estimate[1]) * 34738,
			upper = (quantile(value, 1 - 0.025) - mean(value) + estimate[1]) * 34738
),
variable]

# Plot sampling distributions ####
rr.bs %>% ggplot(
	aes(x = log(value))) +
	geom_histogram(bins = 50) +
	geom_vline(xintercept = 0, color = "salmon") +
	facet_wrap(. ~ variable) +
	mytheme

rd.bs %>% ggplot(
	aes(x = value)) +
	geom_histogram(bins = 50) +
	geom_vline(xintercept = 0, color = "salmon") +
	facet_wrap(. ~ variable) +
	mytheme
