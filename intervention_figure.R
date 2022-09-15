# Making figure for explaining interventions ####
# Kevin Chen
# June 28, 2022

library(here)
source("~/headrs/00-my-theme.R")
library(data.table)

eg.dat <- data.table(
	Exposure = c(0.16, 0.1, 0.02),
	`Other exposure` = c(0.04, 0.1, 0.02)
)

eg.dat[, Static := {
	x <- Exposure
	x[x > 0.05] <- 0.05
	x
}]

eg.dat[Exposure + `Other exposure` > 0.05, Dynamic := {
	x <- Exposure
	x <- 0.05 - `Other exposure`
	x[x < 0] <- 0
	x
}]

eg.dat[Exposure + `Other exposure` < 0.05, Dynamic := Exposure]

eg.melt <- data.table(
	mid = c(eg.dat$Exposure, eg.dat$Static,
					eg.dat$Exposure, eg.dat$Dynamic),
	upper = c(eg.dat$`Other exposure`,
						eg.dat$`Other exposure`,
						eg.dat$`Other exposure` + eg.dat$Exposure,
						eg.dat$`Other exposure` + eg.dat$Dynamic),
	exposure = c(rep("Observed exposure", 3), rep("Post-intervention exposure", 3),
							 rep("Observed exposure", 3), rep("Post-intervention exposure", 3)),
	intervention = c(rep("Static", 6), rep("Dynamic", 6))
)

eg.melt[,upper := mid + rep(eg.dat$`Other exposure`, 4)]

eg.melt[,intervention := factor(
	intervention,
	c("Static", "Dynamic"),
	c("Static intervention", "Dynamic intervention"))]

eg.melt[,`:=`(
	xmin = rep(1:3 - 0.4, 4),
	xmax = rep(1:3 + 0.4, 4)
)]

# Simplified layout?
eg.melt <- eg.melt[c(1:6, 10:12)]
eg.melt$type <- rep(c(
	"No intervention",
	"Static intervention",
	"Dynamic intervention"), each = 3)
eg.melt[,type:= factor(type, c(
	"No intervention",
	"Static intervention",
	"Dynamic intervention"),
	c(
	"A) No intervention",
	"B) Static intervention",
	"C) Dynamic intervention"))]

eg.melt %>% ggplot(aes(
	xmin = xmin,
	xmax = xmax,
)) +
	geom_rect(
		aes(ymin = 0, ymax = mid, fill = "Exposure to soluble MWF"),
		color = "black") +
	geom_rect(
		aes(ymin = mid, ymax = upper, fill = "Other MWF  exposure"),
		color = "black") +
	geom_rect(
		ymin = 0, ymax = 0, xmin = 0, xmax = 0,
		fill = "white", color = "black",
		aes(lty = "Hypothetical exposure limit")) +
	geom_hline(aes(yintercept = 0.05), lty = 2) +
	scale_y_continuous(breaks = seq(0, 0.2, 0.05),
										 labels = c(0, "", "", "", "")) +
	scale_x_continuous(breaks = 1:3,
										 labels = paste0("p-yr ", 1:3)) +
	facet_wrap(. ~ type, nrow = 1) +
	scale_fill_manual(values = c("white", "lightgrey")) +
	scale_linetype_manual(values = 2) +
	coord_cartesian(ylim = c(0, 0.2)) +
	mytheme +
	labs(y = "Average annual exposure (mg/m$^3$)") +
	theme(
		axis.line.y = element_line(),
		axis.line.x = element_blank(),
		# axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		plot.background = element_blank(),
		panel.border = element_blank(),
		strip.background = element_rect(color = "white", size = 3/.pt),
		strip.text = element_text(hjust = 0),
		legend.position = "bottom",
		legend.title = element_blank(),
		legend.key.size = unit(8, "pt"),
		# legend.direction = "vertical",
		# legend.box = "vertical",
		# legend.spacing = unit(5, "pt"),
		# # legend.box.margin = margin(),
		# legend.margin = margin(t = -10),
		# legend.box.just = "left"
	) -> fig

print(fig)

tikzDevice::tikz(here::here("resources", "interventions.tex"),
								 height = 3, width = 6.5, standAlone = T)
print(fig)
dev.off()
lualatex("interventions\\.tex", here::here('resources'))
