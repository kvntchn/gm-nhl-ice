# Figures for illustrating intervention ####
# Kevin Chen
# May 25, 2023

library(here)
library(boxr); box_auth()

source("~/HeadRs/00-my-theme.R")
options(
	tikzLualatexPackages = c(
		# "\\usepackage[utf8]{inputenc}",
		"\\usepackage{amssymb}",
		"\\usepackage{tikz}\n",
		"\\RequirePackage[T1]{fontenc}\n",
		"\\RequirePackage{lmodern}",
		"\\usepackage[active,tightpage,psfixbb]{preview}\n",
		"\\PreviewEnvironment{pgfpicture}\n",
		"\\setlength\\PreviewBorder{0pt}\n",
		# "\\input{\\string~/HeadRs/common_supplement.tex}\n",
		# "\\input{\\string~/HeadRs/stathead.sty}\n",
		NULL
	)
)

# # Shape data
# source(here::here("shape-data.R"))
# dat <- box_read(1179471733155)
# dat <- dat[yin.gm >= 1965]
dat.reduced <- box_read(928445783979)
N <- length(unique(dat.reduced[!is.na(I)]$studyno))
message(sum(dat.reduced[!is.na(I)]$status))
source(here("g-formula.R"))
dat.reduced[,.(year1 = min(year), year2 = max(year)), I]

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

k <- 2
dta <- dat.reduced[I <= k]
timevar_covariates.which <- c(
	"Employment status",
	"soluble", "Soluble",
	NULL
)
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
intervene.which.k <- c(paste0("soluble_", k), paste0("Soluble_", k))
dynamic.which.static <- dynamic.which[!dynamic.which %in% timevar_covariates.which]
dynamic.which.timevar <- dynamic.which[dynamic.which %in% timevar_covariates.which]
dynamic.which <- c(
	dynamic.which.static,
	sapply(dynamic.which.timevar, paste0, "_", 1:k)
)
names(dynamic.which) <- NULL
dynamic.which <- dynamic.which[!dynamic.which %in% intervene.which.k]
dta.wide <- dcast(
	dta[studyno %in% dta[
		immortal == 0 & (I == 1 & Censored == 0 & status == 0), studyno]],
	studyno ~ I,
	value.var = timevar_covariates.which)
dta.wide <- merge(
	dta[I == k, unique(c(
		"studyno", "py",
		static_covariates.which,
		dynamic.which.static)), with = F],
	dta.wide,
	by = "studyno",
	all.x = F, all.y = T)

a <- 0.25
post_intervention <- dta.wide[,.(soluble = {
	if (sum(get(paste0('soluble_', k)) <= a) > 0) {
		max(get(paste0('soluble_', k))[get(paste0('soluble_', k)) <= a])
	} else {min(get(paste0('soluble_', k)))}
}), eval(dynamic.which)]

soluble_distribution <- post_intervention[
	dta.wide[,c('studyno', dynamic.which, intervene.which.k), with = F],
	on = eval(dynamic.which)
][,.(Observed = get(paste0('soluble_', k)), `Post-intervention` = soluble,
		 studyno)]

soluble_distribution[is.na(`Post-intervention`), `Post-intervention` := Observed]
soluble_distribution[Observed <= a, `Post-intervention` := Observed]



soluble_combos <- dta.wide[, c('studyno', dynamic.which), with = F][
	soluble_distribution, on = 'studyno'
]

soluble_combo_distributions <- rbindlist(list(
	soluble_combos[soluble_combos[,.N, dynamic.which][order(N, decreasing = T)][1,], on = dynamic.which],
	soluble_combos[soluble_combos[,.N, dynamic.which][order(N, decreasing = T)][8,], on = dynamic.which],
	soluble_combos[soluble_combos[,.N, dynamic.which][order(N, decreasing = T)][46,], on = dynamic.which]
), idcol = T)

soluble_combo_distributions[,.id := factor(.id, labels = c(
	"Target exposure limit is supported",
	"Supportable exposure limit below target",
	"Supportable exposure limit above target"
))]

mytheme <- theme_bw() +
	theme(legend.position = 'bottom',
				legend.box = 'vertical',
				legend.margin = margin(t = 5, b = 5),
				panel.grid = element_blank(),
				panel.border = element_rect(color = "black", linewidth = 1 / .pt),
				strip.background = element_rect(
					color = 'black', fill="lightgrey", linewidth = 1 / .pt),
				axis.ticks = element_line(linewidth = 1 / .pt),
				legend.key.size = unit(8, 'pt'),
				legend.spacing.x = unit(2, 'pt'),
				legend.spacing.y = unit(2.5, 'pt'),
				legend.title = element_text(size = 8, colour = 'black'),
				legend.text  = element_text(size = 8, colour = 'black'),
				axis.text    = element_text(size = 8, colour = 'black'),
				axis.title.y = element_text(size = 8, margin = margin(r = 6)),
				axis.title.x = element_text(size = 8, margin = margin(t = 6)),
				plot.margin = margin(t = 6, b = 6, l = 6, r = 6)
	)

tikz('reports/SER (2023)/resources/intervention-density.tex',
		 width = 3.75, height = 2.7,
		 documentDeclaration = "\\documentclass{beamer}",
		 bareBones = T, standAlone = F)
melt(soluble_distribution[Observed > 0, -'studyno'],
		 measure.vars = names(soluble_distribution[,-'studyno'])) -> ggdat
ggdat[grep("Post", variable), variable := "Post-intervention\\ \\ \\ \\ "]
ggdat[,variable := factor(variable, levels = rev(unique(variable)))]
print(
	ggplot(ggdat,
				 aes(x = value, col = variable, fill = variable,
				 		alpha = variable, linewidth = variable)) +
		# geom_density(fill = NA, bw = 0.1, linewidth = 8 / .pt) +
		geom_histogram(position = 'identity',
									 aes(y = after_stat(density)),
									 breaks = exp(seq(log(1e-3), log(10), 0.5))
		) +
		geom_vline(aes(lty = '0.25 mg/m\\textsuperscript{3}', xintercept = a),
							 linewidth = 4 / .pt) +
		scale_linetype_manual(values = 2) +
		labs(x = "Exposure to soluble MWF (mg/m\\textsuperscript{3})",
				 y = "Density",
				 col = "Distribution:\\ \\ \\ ",
				 fill = "Distribution:\\ \\ \\ ",
				 alpha = "Distribution:\\ \\ \\ ",
				 linewidth = "Distribution:\\ \\ \\ ",
				 lty = "Target exposure limit:") +
		scale_color_manual(values = c("#A3CFEB", "#ba8123")) +
		scale_fill_manual(values = c("#A3CFEB", "#ffffff")) +
		scale_alpha_manual(values = 1:0) +
		scale_linewidth_manual(values = c(0, 8 / .pt)) +
		scale_x_continuous(trans = 'log',
											 breaks = c(1e-3, 0.05, 1, 10),
											 labels = c(1e-3, 0.05, 1, 10)) +
		guides(linewidth = guide_legend(override.aes = list(
			linewidth = 4))) +
		mytheme
)
dev.off()
# lualatex(directory = 'reports/private/SER (2023)/resources')

for (i in 1:3) {
	fig.name <- paste0('reports/paper/resources/intervention-density-witnin-combos-',
										 letters[i], '.tex')
	tikz(fig.name,
			 width = 3.5, height = 2.5,
			 # width = 3.75, height = 1.85 * 3 + 0.5,
			 # documentDeclaration = "\\documentclass{beamer}",
			 bareBones = F, standAlone = T)
	melt(soluble_combo_distributions[Observed > 0, .(Observed, `Post-intervention`, .id)],
			 measure.vars = c('Observed', 'Post-intervention')) -> ggdat
	ggdat[grep("Post", variable), variable := "Post-intervention"]
	ggdat[,variable := factor(variable, levels = rev(unique(variable)))]
	print(cowplot::plot_grid(
		ggplot(ggdat[.id == levels(.id)[i],],
					 aes(x = value, col = variable, fill = variable,
					 		alpha = variable, linewidth = variable)) +
			geom_histogram(position = 'identity',
										 aes(y = after_stat(density)),
										 breaks = seq(0, 0.65, 0.025)) +
			geom_vline(aes(lty = '0.25 mg/m$^3$', xintercept = a),
								 linewidth = 1 / .pt) +
			scale_linetype_manual(values = 5) +
			labs(x = "Exposure to soluble MWF (mg/m$^3$)",
					 y = "Density",
					 col = "Distribution:",
					 fill = "Distribution:",
					 alpha = "Distribution:",
					 linewidth = "Distribution:",
					 lty = "Target exposure limit:") +
			scale_color_manual(values = c("lightgrey", "black")) +
			scale_fill_manual(values = c("lightgrey", "#ffffff")) +
			scale_alpha_manual(values = 1:0) +
			scale_linewidth_manual(values = c(0, 1 / .pt)) +
			# facet_wrap(. ~ .id, ncol = 3) +
			mytheme + theme(
				legend.position = c(ifelse(i != 3, 1, 0.01), 0.8),
				legend.justification = ifelse(i != 3, "right", "left"),
				legend.margin = margin(t = 0, r = 0, b = 2.5, l = 0),
				legend.box.margin = margin(t = 0, r = 2, b = 0, l = 2),
				plot.margin = margin(t = 15, l = 12, 6, 6)
			),
		labels = paste0(letters[i], ")"), label_fontface = "plain", label_size = 10,
		hjust = -0.25, vjust = 1.5
	))
	dev.off()
	}
lualatex(directory = 'reports/paper/resources')


# tikz('reports/private/SER (2023)/resources/intervention-density-combined.tex',
# 		 width = (3.5 - 0.5) * 4 + 0.5, height = 2.5,
# 		 documentDeclaration = "\\documentclass{beamer}",
# 		 bareBones = T, standAlone = F)
# print(
# 	ggplot(rbindlist(list(
# 		cbind('.id' = "Marginal distribution",
# 					melt(soluble_distribution[Observed > 0, .(Observed, `Post-intervention`)],
# 				 measure.vars = c('Observed', 'Post-intervention'))),
# 		melt(soluble_combo_distributions[Observed > 0, .(Observed, `Post-intervention`, .id)],
# 							 measure.vars = c('Observed', 'Post-intervention'))
# 	)),
# 	aes(x = value, col = variable)) +
# 		# geom_density(fill = NA, bw = 0.1, linewidth = 8 / .pt) +
# 		facet_wrap(. ~ .id, ncol = 4) +
# 		geom_histogram(linewidth = 8 / .pt, position = 'identity', alpha = 0.2,
# 									 aes(y = after_stat(density), fill = variable),
# 									 breaks = seq(0, 0.65, 0.025)) +
# 		geom_vline(aes(lty = '0.25 mg/m\\textsuperscript{3}', xintercept = a)) +
# 		scale_linetype_manual(values = 2) +
# 		# geom_vline(lty = 2, xintercept = a) +
# 		labs(x = "Exposure to soluble MWF (mg/m\\textsuperscript{3})",
# 				 y = "Density",
# 				 col = "Distribution:\\ \\ \\ ",
# 				 fill = "Distribution:\\ \\ \\ ",
# 				 lty = "Target exposure limit:") +
# 		# scale_x_continuous(trans = 'log',
# 		# 									 breaks = c(1e-3, 0.05, 1, 10),
# 		# 									 labels = c(1e-3, 0.05, 1, 10)) +
# 		theme_bw() +
# 		scale_color_manual(values = c("#ba8123", "#003c61")) +
# 		mytheme)
# )
# dev.off()