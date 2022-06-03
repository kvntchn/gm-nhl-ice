# Run parametric g-formula ####
# Kevin Chen
# November 3, 2020

library(here)
source(here("g-formula.R"))

# # Shape data
# source(here::here("shape-data.R"))

dat.reduced <- box_read(928445783979)

# name          : dat.total.reduced.rds
dat.total.reduced <- box_read(928452001570)

# name          : dat.str.reduced.rds
dat.str.reduced <- box_read(928446426341)
# name          : dat.str_total.reduced.rds
dat.str_total.reduced <- box_read(928447106294)

# name          : dat.sol.reduced.rds
dat.sol.reduced <- box_read(928445303042)
# name          : dat.sol_total.reduced.rds
dat.sol_total.reduced <- box_read(928447942153)

# name          : dat.syn.reduced.rds
dat.syn.reduced <- box_read(928451552327)
# name          : dat.syn_total.reduced.rds
dat.syn_total.reduced <- box_read(928447597463)

# name          : dat.total.reduced2.rds
dat.total.reduced2 <- box_read(929952423129)

# name          : dat.str.reduced2.rds
dat.str.reduced2 <- box_read(929952765078)
# name          : dat.sol.reduced2.rds
dat.sol.reduced2 <- box_read(929953044799)
# name          : dat.syn.reduced2.rds
dat.syn.reduced2 <- box_read(929951401906)

# name          : dat.str_total.reduced2.rds
dat.str_total.reduced2 <- box_read(929957702377)
# name          : dat.sol_total.reduced2.rds
dat.sol_total.reduced2 <- box_read(929954713019)
# name          : dat.syn_total.reduced2.rds
dat.syn_total.reduced2 <- box_read(929954585841)

# name          : dat.str.reduced10.rds
dat.str.reduced10 <- box_read(942770799114)
# name          : dat.sol.reduced10.rds
dat.sol.reduced10 <- box_read(930687554963)
# name          : dat.syn.reduced10.rds
dat.syn.reduced10 <- box_read(942772072287)


#####################################################
# Post-intervention distribution of exposure     ####
#####################################################
# # Proportion of person-years intervened upon
# interventions <- data.frame(
# 	dt.name = c(
# 		"dat.sol.reduced", "dat.sol.reduced2", "dat.sol.reduced10",
# 		"dat.sol_total.reduced", "dat.sol_total.reduced2",
# 		"dat.str.reduced", "dat.str.reduced2", "dat.str.reduced10",
# 		"dat.str_total.reduced", "dat.str_total.reduced2",
# 		"dat.syn.reduced", "dat.syn.reduced2", "dat.syn.reduced10",
# 		"dat.syn_total.reduced", "dat.syn_total.reduced2"
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
#
exposure.comparison <- rbindlist(c(
	list(dat.reduced[, .(
		studyno, I, status, year, py, `Cumulative soluble`, cum_soluble)]),
	lapply(list(dat.sol.reduced, dat.sol.reduced2, dat.sol.reduced10), function(x) {
		x[,.(
		studyno, I, status, year, py, `Cumulative soluble`, cum_soluble)]
	})),
	idcol = "Intervention"
	)

exposure.comparison[,`:=`(Intervention = c(
	"No limit", "0.5", "0.25", "0.05")[Intervention])]

exposure.comparison[,`:=`(year.max = max(year), I.max = max(I)), studyno]

dcast(exposure.comparison[,.(py = sum(py)), .(`Cumulative soluble`, Intervention)],
			`Cumulative soluble` ~ Intervention)
table(exposure.comparison[status == 1,.(`Cumulative soluble`, Intervention)])

tikzDevice::tikz(here::here("resources", "shift-soluble.tex"), height = 2.75, width = 4, standAlone = T)
exposure.comparison[I == I.max] %>% ggplot(
	aes(x = cum_soluble)
	# color = '#102451'
) +
	geom_density(aes(lty = Intervention)
							 # , color = '#102451'
							 ) +
	# geom_histogram(aes(fill = Intervention), bins = 40, alpha = 0.4, position = 'identity') +
	# facet_wrap(. ~ Intervention, ncol = 1) +
	scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotted")) +
	mytheme +
	scale_x_continuous(
		n.breaks = 6,
		trans = "log",
		labels = function(x) {sprintf("%.3f", x)}) +
	geom_rug(data = exposure.comparison[status == 1, .(cum_soluble)],
					 alpha = 0.5, size = 0.05, length = grid::unit(3.5, "pt")
					 # color = '#102451'
					 ) +
	labs(
		y = "Density",
		x = "Cumulative exposure (mg/m$^3\\cdot\\,$years)",
		color = "Exposure limit") + theme(
			axis.title = element_text(margin = margin(10, 5, 10, 5)),
			legend.title = element_blank(),
			legend.position = "bottom",
			legend.box.spacing = unit(2, "pt"),
			legend.key.size = unit(10, "pt")
		) +
	theme(panel.grid = element_blank()) +
	coord_cartesian(ylim = c(0, 0.5))
dev.off()
lualatex("shift-soluble\\.tex", here::here('resources'))

dynamic_exposure.comparison <- rbindlist(c(
	list(dat.reduced[, .(
		studyno, I, status, year, py, `Cumulative soluble`, cum_soluble)]),
	lapply(list(dat.sol_total.reduced, dat.sol_total.reduced2), function(x) {
		x[,.(
		studyno, I, status, year, py, `Cumulative soluble`, cum_soluble)]
	})),
	idcol = "Intervention"
	)

dynamic_exposure.comparison[,`:=`(Intervention = c(
	"No limit", "0.5", "0.25")[Intervention])]

dynamic_exposure.comparison[,`:=`(year.max = max(year), I.max = max(I)), studyno]

dcast(dynamic_exposure.comparison[,.(py = sum(py)), .(`Cumulative soluble`, Intervention)],
			`Cumulative soluble` ~ Intervention)
table(dynamic_exposure.comparison[status == 1,.(`Cumulative soluble`, Intervention)])

tikzDevice::tikz(here::here("resources", "shift-soluble-dynamic.tex"), height = 2.75, width = 4, standAlone = T)
dynamic_exposure.comparison[I == I.max] %>% ggplot(
	aes(x = cum_soluble)
	# color = '#102451'
) +
	geom_density(aes(lty = Intervention)
							 # , color = '#102451'
							 ) +
	# geom_histogram(aes(fill = Intervention), bins = 40, alpha = 0.4, position = 'identity') +
	# facet_wrap(. ~ Intervention, ncol = 1) +
	scale_linetype_manual(values = c("solid", "longdash", "dotted")) +
	mytheme +
	scale_x_continuous(
		n.breaks = 6,
		trans = "log",
		labels = function(x) {sprintf("%.3f", x)}) +
	geom_rug(data = exposure.comparison[status == 1, .(cum_soluble)],
					 alpha = 0.5, size = 0.05, length = grid::unit(6.5, "pt")
					 # color = '#102451'
					 ) +
	labs(
		y = "Density",
		x = "Cumulative exposure (mg/m$^3\\cdot\\,$years)",
		color = "Exposure limit") + theme(
			axis.title = element_text(margin = margin(10, 5, 10, 5)),
			legend.title = element_blank(),
			legend.position = "bottom",
			legend.box.spacing = unit(2, "pt"),
			legend.key.size = unit(10, "pt")
		) +
	theme(panel.grid = element_blank()) +
	coord_cartesian(ylim = c(0, 0.5))
dev.off()
lualatex("shift-soluble-dynamic\\.tex", here::here('resources'))

#
# #####################################################
# # Point estimates                                ####
# #####################################################
#
# Baseline
baseline <- ice_gcomp()

# Intervention
reduction.str <- ice_gcomp(a = 1, dta.a1 = copy(dat.str.reduced), mwf.name = "straight")
reduction.str_total <- ice_gcomp(a = 1, dta.a1 = copy(dat.str_total.reduced), mwf.name = "straight")

reduction.str2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.str.reduced2), mwf.name = "straight")
reduction.str_total2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.str_total.reduced2), mwf.name = "straight")

reduction.str10 <- ice_gcomp(a = 1, dta.a1 = copy(dat.str.reduced10), mwf.name = "straight")

reduction.sol <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol.reduced), mwf.name = "soluble")
reduction.sol_total <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol_total.reduced), mwf.name = "soluble")

reduction.sol2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol.reduced2), mwf.name = "soluble")
reduction.sol_total2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol_total.reduced2), mwf.name = "soluble")

reduction.sol10 <- ice_gcomp(a = 1, dta.a1 = copy(dat.sol.reduced10), mwf.name = "soluble")

reduction.syn <- ice_gcomp(a = 1, dta.a1 = copy(dat.syn.reduced), mwf.name = "synthetic")
reduction.syn_total <- ice_gcomp(a = 1, dta.a1 = copy(dat.syn_total.reduced), mwf.name = "synthetic")

reduction.syn2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.syn.reduced2), mwf.name = "synthetic")
reduction.syn_total2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.syn_total.reduced2), mwf.name = "synthetic")

reduction.syn10 <- ice_gcomp(a = 1, dta.a1 = copy(dat.syn.reduced10), mwf.name = "synthetic")

# reduction.total <- ice_gcomp(a = 1, dta.a1 = copy(dat.total.reduced), mwf.name = c("straight", "soluble", "synthetic"))
# reduction.total2 <- ice_gcomp(a = 1, dta.a1 = copy(dat.total.reduced2), mwf.name = c("straight", "soluble", "synthetic"))

observed <- sum(dat.reduced$status)/n_distinct(dat.reduced$studyno)
#
# # Save results
# estimates <- list(
# 	baseline = baseline,
# 	reduction.str = reduction.str,
# 	reduction.str2 = reduction.str2,
# 	reduction.str_total = reduction.str_total,
# 	reduction.str_total2 = reduction.str_total2,
# 	reduction.str10 = reduction.str10,
# 	reduction.sol = reduction.sol,
# 	reduction.sol_total = reduction.sol_total,
# 	reduction.sol2 = reduction.sol2,
# 	reduction.sol_total2 = reduction.sol_total2,
# 	reduction.sol10 = reduction.sol10,
# 	reduction.syn = reduction.syn,
# 	reduction.syn_total = reduction.syn_total,
# 	reduction.syn2 = reduction.syn2,
# 	reduction.syn_total2 = reduction.syn_total2,
# 	reduction.syn10 = reduction.syn10,
# 	reduction.total = reduction.total,
# 	reduction.total2 = reduction.total2,
# 	observed = sum(dat.reduced$status)/n_distinct(dat.reduced$studyno)
# )
#
# save(
# 	baseline,
# 	reduction.str,
# 	reduction.str_total,
# 	reduction.str2,
# 	reduction.str_total2,
# 	reduction.str10,
# 	reduction.sol,
# 	reduction.sol_total,
# 	reduction.sol2,
# 	reduction.sol_total2,
# 	reduction.sol10,
# 	reduction.syn,
# 	reduction.syn_total,
# 	reduction.syn2,
# 	reduction.syn_total2,
# 	reduction.syn10,
# 	reduction.total,
# 	reduction.total2,
# 	observed, file = here::here("resources", paste0("estimates.rdata")))
#
# h.baseline <- mean(baseline$h.pred)
# h.reduction.str <- mean(reduction.str$h.pred)
# h.reduction.str_total <- mean(reduction.str_total$h.pred)
# h.reduction.sol <- mean(reduction.sol$h.pred)
# h.reduction.sol_total <- mean(reduction.sol_total$h.pred)
# h.reduction.syn <- mean(reduction.syn$h.pred)
# h.reduction.syn_total <- mean(reduction.syn_total$h.pred)
# h.reduction.total <- mean(reduction.total$h.pred)
# h.reduction.str2 <- mean(reduction.str2$h.pred)
# h.reduction.str_total2 <- mean(reduction.str_total2$h.pred)
# h.reduction.sol2 <- mean(reduction.sol2$h.pred)
# h.reduction.sol_total2 <- mean(reduction.sol_total2$h.pred)
# h.reduction.syn2 <- mean(reduction.syn2$h.pred)
# h.reduction.syn_total2 <- mean(reduction.syn_total2$h.pred)
# h.reduction.total2 <- mean(reduction.total2$h.pred)
# h.reduction.str10 <- mean(reduction.str10$h.pred)
# h.reduction.sol10 <- mean(reduction.sol10$h.pred)
# h.reduction.syn10 <- mean(reduction.syn10$h.pred)
#
# save(
# 	h.baseline,
# 	h.reduction.str,
# 	h.reduction.str_total,
# 	h.reduction.str2,
# 	h.reduction.str_total2,
# 	h.reduction.str10,
# 	h.reduction.sol,
# 	h.reduction.sol_total,
# 	h.reduction.sol2,
# 	h.reduction.sol_total2,
# 	h.reduction.sol10,
# 	h.reduction.syn,
# 	h.reduction.syn_total,
# 	h.reduction.syn2,
# 	h.reduction.syn_total2,
# 	h.reduction.syn10,
# 	h.reduction.total,
# 	h.reduction.total2,
# 	file = here::here("resources", paste0("h.estimates.rdata")))

load(here::here("resources", paste0("estimates.rdata")))
load(here::here("resources", paste0("h.estimates.rdata")))


# Print results
cat("Observed cumulative incidence:", observed,
		"\nIntervention cumulative incidence:",
		"\nNatural course cumulative incidence:\t\t", h.baseline,
		"\n\tReduce straight to REL:\t\t\t", h.reduction.str,
		"\n\tReduce straight for total under REL:\t", h.reduction.str_total,
		"\n\tReduce straight to REL/2:\t\t", h.reduction.str2,
		"\n\tReduce straight for total under REL/2:\t", h.reduction.str_total2,
		"\n\tReduce straight to REL/10:\t\t", h.reduction.str10/2,
		"\n\tReduce soluble to REL:\t\t\t", h.reduction.sol,
		"\n\tReduce soluble for total under REL:\t", h.reduction.sol_total,
		"\n\tReduce soluble to REL/2:\t\t", h.reduction.sol/2,
		"\n\tReduce soluble for total under REL/2:\t", h.reduction.sol_total/2,
		"\n\tReduce soluble to REL/10:\t\t", h.reduction.sol10/2,
		"\n\tReduce synthetic to REL:\t\t", h.reduction.syn,
		"\n\tReduce synthetic for total under REL:\t", h.reduction.syn_total,
		"\n\tReduce synthetic to REL/2:\t\t", h.reduction.syn/2,
		"\n\tReduce synthetic for total under REL/2:\t", h.reduction.syn_total/2,
		"\n\tReduce synthetic to REL/10:\t\t", h.reduction.syn10/2,
		"\n\tScale all down to REL:\t\t\t", h.reduction.total,
		"\n\tScale all down to REL/2:\t\t", h.reduction.total2,
		"\n",
		"Risk ratio:",
		"\n\tReduce straight to REL:\t\t\t", h.reduction.str/h.baseline,
		"\n\tReduce straight to REL/2:\t\t", h.reduction.str2/h.baseline,
		"\n\tReduce straight for total under REL:\t", h.reduction.str_total/h.baseline,
		"\n\tReduce straight for total under REL/2:\t", h.reduction.str_total2/h.baseline,
		"\n\tReduce straight to REL/10:\t\t", h.reduction.str10/h.baseline,
		"\n\tReduce soluble to REL:\t\t\t", h.reduction.sol/h.baseline,
		"\n\tReduce soluble to REL/2:\t\t", h.reduction.sol2/h.baseline,
		"\n\tReduce soluble to REL/10:\t\t", h.reduction.sol10/h.baseline,
		"\n\tReduce soluble for total under REL:\t", h.reduction.sol_total/h.baseline,
		"\n\tReduce soluble for total under REL/2:\t", h.reduction.sol_total2/h.baseline,
		"\n\tReduce synthetic to REL:\t\t", h.reduction.syn/h.baseline,
		"\n\tReduce synthetic to REL/2:\t\t", h.reduction.syn2/h.baseline,
		"\n\tReduce synthetic for total under REL:\t", h.reduction.syn_total/h.baseline,
		"\n\tReduce synthetic for total under REL/2:\t", h.reduction.syn_total2/h.baseline,
		"\n\tReduce synthetic to REL/10:\t\t", h.reduction.syn10/h.baseline,
		"\n\tScale all down to REL:\t\t\t", h.reduction.total/h.baseline,
		"\n\tScale all down to REL/2:\t\t", h.reduction.total2/h.baseline,
		"\n",
		"Risk difference:",
		"\n\tReduce straight to REL:\t\t\t", h.reduction.str - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.str - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce straight to REL2:\t\t", h.reduction.str2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.str2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce straight to REL/10:\t\t", h.reduction.str10 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.str10 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce straight for total under REL:\t", h.reduction.str_total - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.str_total - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce straight for total under REL/2:\t", h.reduction.str_total2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.str_total2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce soluble to REL:\t\t\t", h.reduction.sol - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.sol - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce soluble to REL/2:\t\t", h.reduction.sol2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.sol2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce soluble to REL/10:\t\t", h.reduction.sol10 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.sol10 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce soluble for total under REL:\t", h.reduction.sol_total - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.sol_total - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce soluble for total under REL/2:\t", h.reduction.sol_total2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.sol_total2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce synthetic to REL:\t\t", h.reduction.syn - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.syn - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce synthetic to REL/2:\t\t", h.reduction.syn2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.syn2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce synthetic to REL/10:\t\t", h.reduction.syn10 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.syn10 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce synthetic for total under REL:\t", h.reduction.syn_total - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.syn_total - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tReduce synthetic for total under REL/2:\t", h.reduction.syn_total2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.syn_total2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tScale all down to REL:\t\t\t", h.reduction.total - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.total - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		"\n\tScale all down to REL/2:\t\t", h.reduction.total2 - h.baseline,
		"\n\t\t\t\t\t\t", (h.reduction.total2 - h.baseline) * n_distinct(dat.reduced$studyno), "\tcases",
		NULL)

# #####################################################
# # MCMC variance estimate                         ####
# #####################################################
# B <- 1e3
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
# 	dat.a1 = list(
# 		dat.str.reduced, dat.str_total.reduced,
# 		dat.str.reduced2, dat.str_total.reduced2,
#     dat.str.reduced10,
# 		dat.sol.reduced, dat.sol_total.reduced,
# 		dat.sol.reduced2, dat.sol_total.reduced2,
#     dat.sol.reduced10,
# 		dat.syn.reduced, dat.syn_total.reduced,
# 		dat.syn.reduced2, dat.syn_total.reduced2,
#     dat.syn.reduced10,
# 		dat.total.reduced, dat.total.reduced2),
# 	a1.names = c(
# 		"straight", "total by straight",
# 		"straight2", "total by straight2",
#     "straight10",
# 		"soluble", "total by soluble",
# 		"soluble2", "total by soluble2",
#     "soluble10",
# 		"synthetic", "total by synthetic",
# 		"synthetic2", "total by synthetic2",
#     "synthetic10",
# 		"scale all", "scale all2"),
# 	mwf.name = list(
# 		"straight", "straight",
# 		"straight", "straight",
#     "straight",
# 		"soluble", "soluble",
# 		"soluble", "soluble",
# 		"soluble",
# 		"synthetic", "synthetic",
# 		"synthetic", "synthetic",
#     "synthetic",
# 		c("straight", "soluble", "synthetic"),
# 		c("straight", "soluble", "synthetic")))
# # })
#
# saveRDS(estimates.bs, here::here("resources", paste0("estimates.bs.rds")))
# estimates.bs <- get.bs(
# 	dat = data.table::copy(dat.reduced), B = 1e3,
# 	dat.a1 = list(
#     dat.str.reduced10,
# 		dat.syn.reduced10),
# 	a1.names = c(
#     "straight10",
# 		"synthetic10"),
# 	mwf.name = list(
# 		"straight",
# 		"synthetic"))
# # })
#
# saveRDS(estimates.bs, here::here("resources", paste0("estimates_10th-of-REL.bs.rds")))
#
# # # Get leave-one-out estimates ####
# # system.time({
# # estimates.loo <- get.loo(
# # 	dat = data.table::copy(dat.reduced),
# # 	leave_out = unique(dat.reduced$studyno)[1:10],
# # 	dat.a1 = list(
# # 		dat.str.reduced, dat.str_total.reduced,
# # 		dat.str.reduced2, dat.str_total.reduced2,
# # 		dat.sol.reduced, dat.sol_total.reduced,
# # 		dat.sol.reduced2, dat.sol_total.reduced2,
# # 		dat.syn.reduced, dat.syn_total.reduced,
# # 		dat.syn.reduced2, dat.syn_total.reduced2,
# # 		dat.total.reduced, dat.total.reduced2),
# # 	a1.names = c(
# # 		"straight", "total by straight",
# # 		"straight2", "total by straight2",
# # 		"soluble", "total by soluble",
# # 		"soluble2", "total by soluble2",
# # 		"synthetic", "total by synthetic",
# # 		"synthetic2", "total by synthetic2",
# # 		"scale all", "scale all2"),
# # 	mwf.name = list(
# # 		"straight", "straight",
# # 		"straight", "straight",
# # 		"soluble", "soluble",
# # 		"soluble", "soluble",
# # 		"synthetic", "synthetic",
# # 		"synthetic", "synthetic",
# # 		c("straight", "soluble", "synthetic"),
# # 		c("straight", "soluble", "synthetic")))
# # })
# #
# # saveRDS(estimates.loo, here::here("resources", paste0("estimates.loo.rds")))
# #
# # estimates.loo <- readRDS(here::here("resources", paste0("estimates.loo.rds")))

# Examine BS sampling distribution ####
estimates.bs <- readRDS(here::here("resources", paste0("estimates.bs.rds")))
do.call(cbind, lapply(estimates.bs[,-(1:2)], function(x) {
	x/estimates.bs[,2]
})) %>% melt(measure.vars = names(.)) -> rr.bs

rr.point_estimates <- data.frame(
	variable = paste0(c(
		"straight", "total by straight", "straight2", "total by straight2",
		"straight10",
		"soluble", "total by soluble", "soluble2", "total by soluble2",
		"soluble10",
		"synthetic", "total by synthetic", "synthetic2", "total by synthetic2",
		"synthetic10",
		"scale all", "scale all2"), ".natural"),
	estimate = c(
		h.reduction.str, h.reduction.str_total, h.reduction.str2, h.reduction.str_total2,
		h.reduction.str10,
		h.reduction.sol, h.reduction.sol_total, h.reduction.sol2, h.reduction.sol_total2,
		h.reduction.sol10,
		h.reduction.syn, h.reduction.syn_total, h.reduction.syn2, h.reduction.syn_total2,
		h.reduction.syn10,
		h.reduction.total, h.reduction.total)/h.baseline
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
		"straight", "total by straight", "straight2", "total by straight2",
		"straight10",
		"soluble", "total by soluble", "soluble2", "total by soluble2",
		"soluble10",
		"synthetic", "total by synthetic", "synthetic2", "total by synthetic2",
		"synthetic10",
		"scale all", "scale all2"), ".natural"),
	estimate = c(
		h.reduction.str, h.reduction.str_total, h.reduction.str2, h.reduction.str_total2,
		h.reduction.str10,
		h.reduction.sol, h.reduction.sol_total, h.reduction.sol2, h.reduction.sol_total2,
		h.reduction.sol10,
		h.reduction.syn, h.reduction.syn_total, h.reduction.syn2, h.reduction.syn_total2,
		h.reduction.syn10,
		h.reduction.total, h.reduction.total) - h.baseline
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
